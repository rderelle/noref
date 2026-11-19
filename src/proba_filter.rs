//! Probabilistic k-mer occupancy counting using:
//! - a per-genome Bloom filter, so each canonical k-mer is counted once per genome;
//! - a global Count-Min Sketch (CMS) estimating how many genomes contain each k-mer.
//!
//! Both structures are cache-blocked: all probes of a key fall in one 64-bit word
//! (Bloom) or one 64-byte cache line (CMS). Keys are processed in batches, and the
//! word or line of a key is prefetched [`AHEAD`] keys before it is used: these tables
//! are far larger than the CPU caches, so most of the time would otherwise be spent
//! waiting for memory, one key at a time. Threads count into the CMS region by region
//! (see [`SharedSketch`]).
//!
//! False positive trade-offs:
//! - Bloom false positives cause undercounting (k-mers may look already seen),
//!   which is conservative.
//! - CMS collisions overestimate counts, which is permissive.

use std::sync::Mutex;

/// Per-genome Bloom filter size (16 MiB): keeps false positives (which slightly
/// undercount occupancies) well below 1% for bacterial-sized genomes.
pub const BLOOM_BITS: usize = 1 << 27;
/// Number of 64-byte CMS lines (8 Mi lines = 512 MiB of `u16` counters).
pub const CMS_LINES: usize = 1 << 23;

/// Counters per cache line, and probes (one per 8-counter quarter) per key.
const LINE: usize = 32;
const PROBES: usize = 4;

const SEEDS: [u64; 4] = [
    0x9e3779b97f4a7c15,
    0xbf58476d1ce4e5b9,
    0x94d049bb133111eb,
    0x3c79ac492ba7b653,
];

fn mix64(mut x: u64) -> u64 {
    // SplitMix64 finalizer: cheap avalanche hash for already-encoded k-mer keys.
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

/// Keys handed to the filters in one call (32 KiB of keys).
pub const BATCH: usize = 1 << 12;

/// Number of keys between prefetching the memory of a key and using it.
const AHEAD: usize = 32;

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

pub struct BloomFilter {
    words: Vec<u64>,
    mask: u64,
}

impl BloomFilter {
    pub fn new(bits: usize) -> Self {
        let words = (bits / 64).max(1).next_power_of_two();
        Self {
            words: vec![0; words],
            mask: words as u64 - 1,
        }
    }

    /// Word and two bits of a key.
    #[inline]
    fn probe(&self, key: u64) -> (usize, u64) {
        let h = mix64(key ^ SEEDS[0]);
        let bits = (1u64 << (h >> 58)) | (1u64 << ((h >> 52) & 63));
        ((h & self.mask) as usize, bits)
    }

    /// Set the key's two bits; returns true when both were already set.
    #[inline]
    pub fn test_and_set(&mut self, key: u64) -> bool {
        let (word, bits) = self.probe(key);
        let word = &mut self.words[word];
        let seen = *word & bits == bits;
        *word |= bits;
        seen
    }

    /// Keep, in order, the keys that were not seen before (same decisions as
    /// [`Self::test_and_set`] on each key in turn).
    pub fn retain_unseen(&mut self, keys: &mut Vec<u64>) {
        let mut kept = 0;
        for i in 0..keys.len() {
            if let Some(&next) = keys.get(i + AHEAD) {
                prefetch(&self.words[self.probe(next).0]);
            }
            let key = keys[i];
            if !self.test_and_set(key) {
                keys[kept] = key;
                kept += 1;
            }
        }
        keys.truncate(kept);
    }
}

/// CMS hash of a key: its low bits select the line, its top 12 bits the counters.
#[inline]
fn cms_hash(key: u64) -> u64 {
    mix64(key ^ SEEDS[1])
}

/// Four distinct counter slots (one per 8-counter quarter) of a CMS hash.
#[inline]
fn slots(h: u64) -> [usize; PROBES] {
    let s = (h >> 52) as usize;
    [
        s & 7,
        8 + ((s >> 3) & 7),
        16 + ((s >> 6) & 7),
        24 + ((s >> 9) & 7),
    ]
}

#[repr(C, align(64))]
#[derive(Clone)]
struct Line([u16; LINE]);

/// Count each hash once in `lines` (a power-of-two number of lines, selected by
/// the low bits of the hash).
fn increment_lines(lines: &mut [Line], hashes: &[u64]) {
    let mask = lines.len() as u64 - 1;
    for (i, &h) in hashes.iter().enumerate() {
        if let Some(&next) = hashes.get(i + AHEAD) {
            prefetch(&lines[(next & mask) as usize]);
        }
        let line = &mut lines[(h & mask) as usize].0;
        for s in slots(h) {
            line[s] = line[s].saturating_add(1);
        }
    }
}

pub struct CountMinSketch {
    lines: Vec<Line>,
    mask: u64,
}

impl CountMinSketch {
    pub fn new(lines: usize) -> Self {
        let n = lines.max(1).next_power_of_two();
        Self {
            lines: vec![Line([0; LINE]); n],
            mask: n as u64 - 1,
        }
    }

    /// Count each key once (single writer).
    #[cfg(test)]
    pub fn increment_batch(&mut self, keys: &[u64]) {
        let hashes: Vec<u64> = keys.iter().map(|&key| cms_hash(key)).collect();
        increment_lines(&mut self.lines, &hashes);
    }

    /// Write the estimated number of genomes containing `keys[i]` into `out[i]`
    /// (minimum of the key's counters).
    pub fn estimate_batch(&self, keys: &[u64], out: &mut [u16]) {
        for (i, (&key, out)) in keys.iter().zip(out).enumerate() {
            if let Some(&next) = keys.get(i + AHEAD) {
                prefetch(&self.lines[(cms_hash(next) & self.mask) as usize]);
            }
            let h = cms_hash(key);
            let line = &self.lines[(h & self.mask) as usize].0;
            *out = slots(h).iter().map(|&s| line[s]).min().unwrap_or(0);
        }
    }

    /// Share the sketch between threads for counting.
    pub fn shared(&mut self) -> SharedSketch<'_> {
        let len = self.lines.len() / REGIONS.min(self.lines.len());
        SharedSketch {
            regions: self.lines.chunks_mut(len).map(Mutex::new).collect(),
            shift: len.trailing_zeros(),
        }
    }

    #[cfg(test)]
    fn estimate(&self, key: u64) -> u16 {
        let mut out = [0];
        self.estimate_batch(&[key], &mut out);
        out[0]
    }
}

/// Number of regions (blocks of consecutive lines) of a [`SharedSketch`].
const REGIONS: usize = 256;
/// Keys a [`SketchWriter`] buffers per region (2 MiB in all).
const REGION_BUFFER: usize = 1 << 10;

/// A [`CountMinSketch`] that several threads count into without atomic counters:
/// each region has its own lock, and threads add their keys region by region
/// through a [`SketchWriter`]. Increments commute, so the counts do not depend on
/// the order of the keys (nor on the number of threads); updating one region at a
/// time also makes memory accesses much more local.
pub struct SharedSketch<'a> {
    regions: Vec<Mutex<&'a mut [Line]>>,
    /// log2 of the number of lines per region.
    shift: u32,
}

impl<'a> SharedSketch<'a> {
    pub fn writer(&self) -> SketchWriter<'_, 'a> {
        SketchWriter {
            buffers: (0..self.regions.len())
                .map(|_| Vec::with_capacity(REGION_BUFFER))
                .collect(),
            sketch: self,
        }
    }
}

/// Buffers the keys of one thread per region of a [`SharedSketch`]; a full buffer
/// is counted under the lock of its region. Call [`Self::finish`] when done.
pub struct SketchWriter<'s, 'a> {
    sketch: &'s SharedSketch<'a>,
    /// CMS hashes of the buffered keys, per region.
    buffers: Vec<Vec<u64>>,
}

impl SketchWriter<'_, '_> {
    /// Count each key once.
    pub fn add(&mut self, keys: &[u64]) {
        let region_mask = self.buffers.len() as u64 - 1;
        for &key in keys {
            let h = cms_hash(key);
            let region = ((h >> self.sketch.shift) & region_mask) as usize;
            self.buffers[region].push(h);
            if self.buffers[region].len() == REGION_BUFFER {
                self.flush(region);
            }
        }
    }

    fn flush(&mut self, region: usize) {
        let mut lines = self.sketch.regions[region]
            .lock()
            .expect("sketch region mutex poisoned");
        increment_lines(&mut lines, &self.buffers[region]);
        self.buffers[region].clear();
    }

    /// Count the keys still buffered.
    pub fn finish(mut self) {
        for region in 0..self.buffers.len() {
            self.flush(region);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slots_use_one_counter_per_quarter() {
        for key in 0..10_000u64 {
            for (q, &s) in slots(cms_hash(key)).iter().enumerate() {
                assert!((8 * q..8 * q + 8).contains(&s));
            }
        }
    }

    #[test]
    fn bloom_reports_previously_inserted_keys() {
        let mut bloom = BloomFilter::new(1 << 16);
        let fresh = (0..500).filter(|&k| !bloom.test_and_set(mix64(k))).count();
        assert!(fresh > 490);
        assert!((0..500).all(|k| bloom.test_and_set(mix64(k))));
    }

    #[test]
    fn batched_bloom_makes_the_same_decisions_as_single_keys() {
        // Small filter and repeated keys: many keys are (or look) already seen.
        let keys: Vec<u64> = (0..5000u64).map(|i| mix64(i % 1777)).collect();
        let mut single = BloomFilter::new(1 << 12);
        let expected: Vec<u64> = keys
            .iter()
            .copied()
            .filter(|&k| !single.test_and_set(k))
            .collect();
        assert!(expected.len() < 1777);
        let mut batched = BloomFilter::new(1 << 12);
        let mut kept = Vec::new();
        for chunk in keys.chunks(700) {
            let mut chunk = chunk.to_vec();
            batched.retain_unseen(&mut chunk);
            kept.extend(chunk);
        }
        assert_eq!(kept, expected);
        assert_eq!(batched.words, single.words);
    }

    /// Count `keys` (then key 42 often enough to saturate) with `threads` writers.
    fn shared_counts(lines: usize, keys: &[u64], threads: usize) -> CountMinSketch {
        let mut cms = CountMinSketch::new(lines);
        let shared = cms.shared();
        std::thread::scope(|scope| {
            for t in 0..threads {
                let shared = &shared;
                scope.spawn(move || {
                    let mut writer = shared.writer();
                    for chunk in keys.chunks(97).skip(t).step_by(threads) {
                        writer.add(chunk);
                    }
                    for _ in 0..65536 / threads + 1 {
                        writer.add(&[42]);
                    }
                    writer.finish();
                });
            }
        });
        drop(shared);
        cms
    }

    #[test]
    fn shared_and_single_writer_counts_match_including_collisions_and_saturation() {
        let keys: Vec<u64> = (0..5000).map(|i| mix64(i % 311)).collect();
        let queries: Vec<u64> = (0..1200).map(mix64).chain([42]).collect();
        for lines in [1, 4, 1024, 1 << 12] {
            let mut plain = CountMinSketch::new(lines);
            plain.increment_batch(&keys);
            plain.increment_batch(&[42; 65536]);
            let mut expected = vec![0; queries.len()];
            plain.estimate_batch(&queries, &mut expected);
            for (&q, &est) in queries.iter().zip(&expected) {
                assert_eq!(est, plain.estimate(q));
            }
            assert_eq!(plain.estimate(42), u16::MAX);
            for threads in [1, 3] {
                let shared = shared_counts(lines, &keys, threads);
                let mut found = vec![0; queries.len()];
                shared.estimate_batch(&queries, &mut found);
                assert_eq!(found, expected);
            }
        }
    }

    #[test]
    fn estimates_never_undercount() {
        let mut cms = CountMinSketch::new(64);
        for round in 0..7 {
            let keys: Vec<u64> = (0..5000u64).filter(|k| k % 7 > round).collect();
            cms.increment_batch(&keys);
        }
        for key in 0..5000u64 {
            assert!(cms.estimate(key) >= (key % 7) as u16);
        }
    }
}
