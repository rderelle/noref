//! Exact k-mer occupancy counting: the number of genomes containing each canonical
//! k-mer, kept only for the k-mers found in at least `min_needed` of the `n` genomes.
//!
//! Pass 1 adds the k-mers of each genome to [`SHARDS`] hash tables, one table at a time
//! under its lock, so that several threads can add genomes at the same time. A table
//! counts a genome once per k-mer (it remembers the last genome counted) and only holds
//! the k-mers that can still reach `min_needed` genomes:
//!
//! - a k-mer found in at least `min_needed` genomes is found in at least one of any
//!   `n - min_needed + 1` of them (pigeonhole principle), so only the first
//!   `n - min_needed + 1` genomes added to a table create entries; later genomes only
//!   increase the counts of the k-mers already there;
//! - once `merged` genomes are in, a k-mer counted `c` times can reach at most
//!   `c + n - merged` genomes: k-mers that can no longer reach `min_needed` are dropped
//!   along the way, and all those below it once every genome is in.
//!
//! Counts are therefore exact and do not depend on the order in which genomes are added
//! (nor on the number of threads). Pass 2 looks up the remaining (shared) k-mers in
//! [`KmerCounts`]; the others read as 0, which gives the same anchors as their actual
//! counts since anchors only compare occupancies of shared k-mers with each other or
//! with lower ones (see `find_segment_anchor_pairs`).
//!
//! Tables store a hash of each k-mer in place of the k-mer: [`kmer_hash`] is a
//! bijection, so this loses nothing, and the hash also selects the table and the slot.

use hashbrown::HashMap;
use std::sync::Mutex;

/// log2 of the number of pass 1 tables: with many tables, threads rarely wait for a
/// lock, and growing or pruning a table is quick.
const SHARD_BITS: u32 = 10;
const SHARDS: usize = 1 << SHARD_BITS;
/// Smallest pass 1 table (slots).
const MIN_CAPACITY: usize = 16;
/// K-mers looked up in one call by pass 2 (32 KiB of hashes).
pub const BATCH: usize = 1 << 12;
/// Number of k-mers between prefetching the memory of a k-mer and using it.
const AHEAD: usize = 32;

/// SplitMix64 finalizer: a bijection of `u64` with good avalanche.
#[inline]
fn mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

/// Hash of a canonical k-mer code, stored in the tables in place of the k-mer (a
/// bijection, so distinct k-mers have distinct hashes). It is 0 only for `u64::MAX`,
/// which is never a canonical k-mer (for `k = 32` its reverse complement is 0, and
/// shorter k-mers are smaller), so 0 marks empty slots.
#[inline]
pub(crate) fn kmer_hash(key: u64) -> u64 {
    mix64(!key)
}

const EMPTY: u64 = 0;

/// Pass 1 table of a hash: its top bits (the slot uses the low bits).
#[inline]
fn shard_of(hash: u64) -> usize {
    (hash >> (u64::BITS - SHARD_BITS)) as usize
}

/// Ask the CPU to start loading the cache line of `target`. This is only a hint: it
/// never changes results, and does nothing on other architectures.
#[inline(always)]
fn prefetch<T>(target: &T) {
    let p: *const u8 = (target as *const T).cast();
    // SAFETY (both blocks): a prefetch is not a memory access for the program (it
    // cannot fault), and `p` comes from a valid reference anyway.
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
        _mm_prefetch(p.cast(), _MM_HINT_T0);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        std::arch::asm!(
            "prfm pldl1keep, [{p}]",
            p = in(reg) p,
            options(nostack, preserves_flags, readonly)
        );
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = p;
}

/// One k-mer of a pass 1 table (16 bytes): its hash ([`EMPTY`] for an empty slot, whose
/// count is 0), the number of genomes containing it, and the last genome counted.
#[derive(Clone, Copy, Debug, Default)]
struct Slot {
    hash: u64,
    count: u16,
    last: u16,
}

/// Pass 1 table, to which `merged` genomes have been added: open addressing with linear
/// probing, a power-of-two number of slots, at most three-quarters full.
#[derive(Debug)]
struct Table {
    slots: Vec<Slot>,
    len: usize,
    merged: usize,
}

impl Table {
    fn new() -> Self {
        Self {
            slots: vec![Slot::default(); MIN_CAPACITY],
            len: 0,
            merged: 0,
        }
    }

    /// First slot probed for a hash.
    #[inline]
    fn home(&self, hash: u64) -> usize {
        hash as usize & (self.slots.len() - 1)
    }

    /// Slot holding `hash`, or the empty slot where it would be inserted.
    #[inline]
    fn find(&self, hash: u64) -> usize {
        let mask = self.slots.len() - 1;
        let mut i = self.home(hash);
        loop {
            let slot_hash = self.slots[i].hash;
            if slot_hash == hash || slot_hash == EMPTY {
                return i;
            }
            i = (i + 1) & mask;
        }
    }

    /// Count of a k-mer hash (0 when absent).
    #[cfg(test)]
    fn count(&self, hash: u64) -> u16 {
        self.slots[self.find(hash)].count
    }

    /// The same table with `capacity` slots and only the k-mers counted at least
    /// `min_count` times.
    fn rebuilt(&self, min_count: u16, capacity: usize) -> Self {
        let mut table = Self {
            slots: vec![Slot::default(); capacity],
            len: 0,
            merged: self.merged,
        };
        for slot in &self.slots {
            if slot.hash != EMPTY && slot.count >= min_count {
                let i = table.find(slot.hash);
                table.slots[i] = *slot;
                table.len += 1;
            }
        }
        table
    }

    /// Count genome `isolate` once for each of its k-mers in this table (`hashes`,
    /// duplicates included).
    fn add(&mut self, isolate: u16, hashes: &[u64], n: usize, min_needed: usize) {
        // A k-mer first seen now can reach at most `n - merged` genomes.
        let can_create = n - self.merged >= min_needed;
        for (j, &hash) in hashes.iter().enumerate() {
            if let Some(&next) = hashes.get(j + AHEAD) {
                prefetch(&self.slots[self.home(next)]);
            }
            let i = self.find(hash);
            let slot = &mut self.slots[i];
            if slot.hash == hash {
                if slot.last != isolate {
                    slot.count += 1;
                    slot.last = isolate;
                }
            } else if can_create {
                *slot = Slot {
                    hash,
                    count: 1,
                    last: isolate,
                };
                self.len += 1;
                if 4 * self.len > 3 * self.slots.len() {
                    *self = self.rebuilt(1, 2 * self.slots.len());
                }
            }
        }
        self.merged += 1;

        // Drop the k-mers that can no longer reach `min_needed`, about 16 times over the
        // run and always after the last genome (no k-mer is created any more then), and
        // leave the table at most half full.
        let remaining = n - self.merged;
        let step = (n / 16).max(1);
        if remaining < min_needed && (self.merged.is_multiple_of(step) || remaining == 0) {
            let min_count = (min_needed - remaining) as u16;
            let kept = self.slots.iter().filter(|s| s.count >= min_count).count();
            let capacity = (2 * kept).next_power_of_two().max(MIN_CAPACITY);
            *self = self.rebuilt(min_count, capacity);
        }
    }
}

/// The canonical k-mers of one genome as hashes, grouped by pass 1 table (duplicates
/// included). Reused from genome to genome.
pub(crate) struct KmerBuckets(Vec<Vec<u64>>);

impl KmerBuckets {
    pub(crate) fn new() -> Self {
        Self(vec![Vec::new(); SHARDS])
    }

    pub(crate) fn clear(&mut self) {
        self.0.iter_mut().for_each(Vec::clear);
    }

    #[inline]
    pub(crate) fn push(&mut self, key: u64) {
        let hash = kmer_hash(key);
        self.0[shard_of(hash)].push(hash);
    }
}

/// Pass 1: exact counts of the k-mers of `n` genomes, added by several threads (see
/// the module documentation).
pub(crate) struct KmerCounter {
    tables: Vec<Mutex<Table>>,
    n: usize,
    min_needed: usize,
}

impl KmerCounter {
    pub(crate) fn new(n: usize, min_needed: usize) -> Self {
        assert!(n <= u16::MAX as usize && (1..=n).contains(&min_needed));
        Self {
            tables: (0..SHARDS).map(|_| Mutex::new(Table::new())).collect(),
            n,
            min_needed,
        }
    }

    /// Add the k-mers of genome `isolate` (each genome exactly once).
    pub(crate) fn add_genome(&self, isolate: u16, kmers: &KmerBuckets) {
        // Tables held by other threads are done last, to wait for them less often.
        let mut busy = Vec::new();
        for (s, table) in self.tables.iter().enumerate() {
            match table.try_lock() {
                Ok(mut table) => table.add(isolate, &kmers.0[s], self.n, self.min_needed),
                Err(_) => busy.push(s),
            }
        }
        for s in busy {
            let mut table = self.tables[s].lock().expect("k-mer table mutex poisoned");
            table.add(isolate, &kmers.0[s], self.n, self.min_needed);
        }
    }

    /// The counts of the shared k-mers, once every genome has been added.
    pub(crate) fn finish(self) -> KmerCounts {
        let tables = self
            .tables
            .into_iter()
            .map(|table| table.into_inner().expect("k-mer table mutex poisoned"));
        KmerCounts::new(tables.collect())
    }
}

/// Hash bits stored in a [`KmerCounts`] slot, below the count.
const LOW_BITS: u32 = 48;
const LOW_MASK: u64 = (1 << LOW_BITS) - 1;
/// Fewest [`KmerCounts`] buckets: with more than `2^(64 - LOW_BITS)` buckets, the hashes
/// of a bucket are less than `2^LOW_BITS` apart, so their low bits tell them apart.
const MIN_BUCKETS: usize = (1 << (u64::BITS - LOW_BITS)) + 1;
/// Slots per bucket (one cache line), and average number of k-mers per bucket.
const BUCKET_SLOTS: usize = 8;
const KMERS_PER_BUCKET: usize = 3;

#[repr(C, align(64))]
#[derive(Clone, Copy, Debug, Default)]
struct Bucket([u64; BUCKET_SLOTS]);

/// Number of genomes containing each k-mer found in at least `min_needed` genomes.
///
/// Pass 2 looks up every k-mer of every genome, so a lookup reads one cache line and
/// does not branch on its content (mispredicted branches stall the processor as much
/// as memory does). Hashes are spread over buckets of [`BUCKET_SLOTS`] slots, each
/// holding the low [`LOW_BITS`] bits of a hash and its count above them (0 when empty;
/// a count is never 0). The k-mers of full buckets (rare) are kept in `overflow`.
pub(crate) struct KmerCounts {
    buckets: Vec<Bucket>,
    overflow: HashMap<u64, u16>,
}

impl KmerCounts {
    /// Lookup table of the k-mers of the pass 1 tables (dropped one by one).
    fn new(tables: Vec<Table>) -> Self {
        let len: usize = tables.iter().map(|t| t.len).sum();
        let mut counts = Self {
            buckets: vec![Bucket::default(); (len / KMERS_PER_BUCKET).max(MIN_BUCKETS)],
            overflow: HashMap::new(),
        };
        for table in tables {
            for slot in table.slots.iter().filter(|s| s.hash != EMPTY) {
                let bucket = counts.bucket_of(slot.hash);
                let entry = (slot.hash & LOW_MASK) | (u64::from(slot.count) << LOW_BITS);
                match counts.buckets[bucket].0.iter_mut().find(|e| **e == 0) {
                    Some(free) => *free = entry,
                    None => _ = counts.overflow.insert(slot.hash, slot.count),
                }
            }
        }
        counts
    }

    /// Bucket of a hash (its top bits, scaled to the number of buckets).
    #[inline]
    fn bucket_of(&self, hash: u64) -> usize {
        ((u128::from(hash) * self.buckets.len() as u128) >> u64::BITS) as usize
    }

    /// Count of a k-mer hash (0 when absent).
    #[inline]
    fn count(&self, hash: u64) -> u16 {
        let bucket = &self.buckets[self.bucket_of(hash)].0;
        let low = hash & LOW_MASK;
        // Empty slots (0) add nothing, even when `low` is 0.
        let count = bucket.iter().fold(0, |count, &entry| {
            count
                | if entry & LOW_MASK == low {
                    entry >> LOW_BITS
                } else {
                    0
                }
        });
        if count == 0 && bucket[BUCKET_SLOTS - 1] != 0 {
            return self.overflow.get(&hash).copied().unwrap_or(0);
        }
        count as u16
    }

    /// Write the number of genomes containing the k-mer of `hashes[i]` (see
    /// [`kmer_hash`]) into `out[i]`: 0 for the k-mers found in fewer than `min_needed`
    /// genomes.
    pub(crate) fn occupancy_batch(&self, hashes: &[u64], out: &mut [u16]) {
        for (i, (&hash, out)) in hashes.iter().zip(out).enumerate() {
            if let Some(&next) = hashes.get(i + AHEAD) {
                prefetch(&self.buckets[self.bucket_of(next)]);
            }
            *out = self.count(hash);
        }
    }

    /// Number of k-mers found in at least `min_needed` genomes.
    #[cfg(test)]
    fn len(&self) -> usize {
        let in_buckets = self.buckets.iter().flat_map(|b| b.0).filter(|&e| e != 0);
        in_buckets.count() + self.overflow.len()
    }

    /// Exact counts of the k-mers of `genomes` (one list of canonical codes per genome,
    /// duplicates allowed), single-threaded.
    #[cfg(test)]
    pub(crate) fn from_genomes(genomes: &[Vec<u64>], min_needed: usize) -> Self {
        let counter = KmerCounter::new(genomes.len(), min_needed);
        let mut kmers = KmerBuckets::new();
        for (isolate, genome) in genomes.iter().enumerate() {
            kmers.clear();
            genome.iter().for_each(|&key| kmers.push(key));
            counter.add_genome(isolate as u16, &kmers);
        }
        counter.finish()
    }

    #[cfg(test)]
    fn occupancy(&self, key: u64) -> u16 {
        let mut out = [0];
        self.occupancy_batch(&[kmer_hash(key)], &mut out);
        out[0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hashbrown::HashSet;

    /// `n` genomes drawn from a pool of k-mers in 8 frequency classes (rare to core),
    /// some of them twice.
    fn random_genomes(n: usize, pool: u64, seed: u64) -> Vec<Vec<u64>> {
        let mut state = seed;
        let mut next = move || {
            state = state.wrapping_add(0x9e3779b97f4a7c15);
            mix64(state)
        };
        (0..n)
            .map(|_| {
                let mut genome = Vec::new();
                for i in 0..pool {
                    let key = i * 2654435761 % (1 << 34);
                    if next() % 8 <= i % 8 {
                        genome.push(key);
                        if next() % 10 == 0 {
                            genome.push(key);
                        }
                    }
                }
                genome
            })
            .collect()
    }

    fn exact_counts(genomes: &[Vec<u64>]) -> HashMap<u64, u16> {
        let mut counts = HashMap::new();
        for genome in genomes {
            for &key in genome.iter().collect::<HashSet<_>>() {
                *counts.entry(key).or_default() += 1;
            }
        }
        counts
    }

    /// Inverse of [`kmer_hash`]: the k-mer code with this hash.
    fn kmer_with_hash(hash: u64) -> u64 {
        let inverse = |c: u64| {
            (0..6).fold(c, |x, _| {
                x.wrapping_mul(2u64.wrapping_sub(c.wrapping_mul(x)))
            })
        };
        let unshift = |y: u64, s: u32| (0..64 / s).fold(y, |x, _| y ^ (x >> s));
        let x = unshift(hash, 31).wrapping_mul(inverse(0x94d049bb133111eb));
        let x = unshift(x, 27).wrapping_mul(inverse(0xbf58476d1ce4e5b9));
        !unshift(x, 30)
    }

    #[test]
    fn hash_is_a_bijection_and_empty_only_for_a_non_canonical_code() {
        for key in (0..1000u64).map(|i| mix64(i) >> (i % 64)) {
            assert_eq!(kmer_with_hash(kmer_hash(key)), key);
        }
        assert_eq!(kmer_with_hash(EMPTY), u64::MAX);
    }

    #[test]
    fn full_buckets_and_hashes_without_low_bits() {
        // 30 k-mers in bucket 5 (10 of them shared), one whose stored hash bits are 0,
        // and many others; k-mers not shared are in 2 of the 4 genomes.
        let bucket_start = |b: u128| (b << 64).div_ceil(MIN_BUCKETS as u128) as u64;
        let mut hashes: Vec<u64> = (0..30).map(|i| bucket_start(5) + i * 1_000_003).collect();
        hashes.push(7 << LOW_BITS);
        hashes.extend((1..3000u64).map(mix64));
        let keys: Vec<u64> = hashes.iter().map(|&h| kmer_with_hash(h)).collect();
        let rare: Vec<u64> = keys.iter().step_by(3).copied().collect();
        let genomes = vec![keys.clone(), keys.clone(), rare.clone(), rare];
        let counts = KmerCounts::from_genomes(&genomes, 3);
        assert_eq!(counts.buckets.len(), MIN_BUCKETS);
        assert_eq!(counts.bucket_of(hashes[29]), 5);
        assert_eq!(counts.overflow.len(), 2);
        assert_eq!(counts.len(), keys.len().div_ceil(3));
        for (i, &key) in keys.iter().enumerate() {
            assert_eq!(counts.occupancy(key), if i % 3 == 0 { 4 } else { 0 });
        }
        // Absent k-mers: in the full bucket, and with the stored bits of a present one.
        assert_eq!(counts.occupancy(kmer_with_hash(bucket_start(5) + 1)), 0);
        assert_eq!(counts.occupancy(kmer_with_hash(8 << LOW_BITS)), 0);
    }

    #[test]
    fn counts_are_exact_for_kmers_reaching_the_threshold_only() {
        for (n, min_needed) in [(1, 1), (5, 3), (12, 7), (12, 12), (40, 21)] {
            let genomes = random_genomes(n, 20_000, n as u64);
            let counts = KmerCounts::from_genomes(&genomes, min_needed);
            let mut shared = 0;
            for (&key, &count) in &exact_counts(&genomes) {
                let want = count * u16::from(count as usize >= min_needed);
                assert_eq!(counts.occupancy(key), want, "n={n} m={min_needed}");
                shared += usize::from(want > 0);
            }
            assert_eq!(counts.len(), shared);
            // Absent k-mers read as 0.
            assert_eq!(counts.occupancy(1 << 40), 0);
        }
    }

    #[test]
    fn counts_do_not_depend_on_genome_order_or_threads() {
        let (n, min_needed) = (16, 10);
        let genomes = random_genomes(n, 30_000, 7);
        let queries: Vec<u64> = genomes.iter().flatten().map(|&k| kmer_hash(k)).collect();
        let lookup = |counts: &KmerCounts| {
            let mut out = vec![0; queries.len()];
            counts.occupancy_batch(&queries, &mut out);
            out
        };
        let expected = lookup(&KmerCounts::from_genomes(&genomes, min_needed));
        assert!(expected.iter().filter(|&&c| c > 0).count() > queries.len() / 3);
        for threads in [2, 3] {
            let counter = KmerCounter::new(n, min_needed);
            std::thread::scope(|scope| {
                for t in 0..threads {
                    let (counter, genomes) = (&counter, &genomes);
                    scope.spawn(move || {
                        let mut kmers = KmerBuckets::new();
                        // Genomes in reverse order, interleaved between threads.
                        for isolate in (0..n).rev().skip(t).step_by(threads) {
                            kmers.clear();
                            genomes[isolate].iter().for_each(|&key| kmers.push(key));
                            counter.add_genome(isolate as u16, &kmers);
                        }
                    });
                }
            });
            assert_eq!(lookup(&counter.finish()), expected, "{threads} threads");
        }
    }

    #[test]
    fn tables_drop_kmers_that_cannot_reach_the_threshold() {
        // 4 genomes, k-mers need 3 of them.
        let (core, early, late) = (kmer_hash(1), kmer_hash(2), kmer_hash(3));
        let mut table = Table::new();
        table.add(0, &[core, early, core], 4, 3);
        table.add(1, &[core], 4, 3);
        assert_eq!((table.count(core), table.count(early)), (2, 1));
        // 1 genome left: `early` can reach 2 at most, and `late` comes too late.
        table.add(2, &[core, late], 4, 3);
        assert_eq!(table.len, 1);
        table.add(3, &[core, late, late], 4, 3);
        assert_eq!((table.len, table.count(core), table.count(late)), (1, 4, 0));
    }
}
