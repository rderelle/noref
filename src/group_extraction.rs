//! Raw variant group extraction from genome assemblies.
//!
//! Pass 1 counts, for every canonical `k`-mer, the number of genomes containing it.
//! Pass 2 rescans each contig, finds shared `k`-mer anchors bracketing local drops in
//! occupancy, and stores the bases between adjacent anchors. These occurrences are
//! merged, genome by genome, into variant groups (one per anchor pair) by the thread
//! that extracted them, in [`SHARDS`] independent tables.
//!
//! Contig orientation in an assembly is arbitrary, so anchor pairs are detected on both
//! strands and stored in a canonical orientation: the same locus therefore yields the
//! same anchor pair and the same sequence whatever the strand it was assembled on.
//!
//! Terms used throughout the crate:
//! - *anchor*: a canonical `k`-mer shared by enough genomes, at the edge of a local
//!   drop in occupancy;
//! - *middle*: the bases between two adjacent anchors;
//! - *occurrence*: an anchor pair and its middle in one genome ([`LocalOccurrence`]
//!   before merging, [`BlockOccurrence`] once merged into its variant group);
//! - *variant group*: the occurrences of all genomes for one anchor pair, which becomes a
//!   set of alignment columns once selected and filtered;
//! - *block*: a middle stored in a [`BlockArena`], packed with 2 bits per base.
use anyhow::{bail, Result};
use hashbrown::HashMap;
use rayon::{prelude::*, ThreadPoolBuilder};
use seq_io::fasta::{Reader as FastaReader, Record};
use std::io::{self, IsTerminal, Write};
use std::sync::Mutex;

use crate::dna::{encode, normalize, revcomp_code, Roller, INVALID};
use crate::io::{open_fasta, IsolateInput};
use crate::proba_filter::{BloomFilter, CountMinSketch, BATCH, BLOOM_BITS, CMS_LINES};
use crate::run_log::{say, RunLog};

const START_CLUSTER_SPAN: u32 = 10;

/// Number of independent group tables (a power of two). A genome is merged into one
/// table at a time, under the lock of that table, so that several threads can merge
/// their genomes at the same time.
const SHARDS: usize = 256;

/// Progress bar over the genome files, drawn only when the standard error is a
/// terminal (so that redirected output and job logs are not filled with it).
struct Progress {
    current: Mutex<usize>,
    total: usize,
    shown: bool,
}

impl Progress {
    fn new(total: usize) -> Self {
        let progress = Self {
            current: Mutex::new(0),
            total,
            shown: io::stderr().is_terminal(),
        };
        progress.draw(0);
        progress
    }

    fn increment(&self) {
        let mut current = self.current.lock().expect("progress mutex poisoned");
        *current += 1;
        self.draw(*current);
    }

    fn draw(&self, current: usize) {
        if !self.shown {
            return;
        }
        const WIDTH: usize = 45;
        let filled = (current * WIDTH).checked_div(self.total).unwrap_or(WIDTH);
        eprint!(
            "\r  \x1b[32m[{}{}]\x1b[0m {}/{} genome files",
            "=".repeat(filled),
            " ".repeat(WIDTH - filled),
            current,
            self.total
        );
        let _ = io::stderr().flush();
    }

    fn finish(&self) {
        if self.shown {
            eprintln!();
        }
    }
}

/// Sequences of up to this many bases are held by their [`Block`] itself.
const INLINE_BASES: usize = 24;

/// Handle to one sequence of a [`BlockArena`], in 8 bytes: length (16 bits) and,
/// above it, the base offset of the sequence in the arena or, for sequences of at
/// most [`INLINE_BASES`] bases (the most common), the bases themselves. Short
/// sequences thus take no arena space, and are equal if their blocks are.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Block(u64);

impl Block {
    #[inline]
    pub(crate) fn len(self) -> usize {
        (self.0 & 0xffff) as usize
    }

    #[inline]
    fn is_inline(self) -> bool {
        self.len() <= INLINE_BASES
    }

    /// Offset of a sequence stored in the arena.
    #[inline]
    fn start(self) -> usize {
        debug_assert!(!self.is_inline());
        (self.0 >> 16) as usize
    }

    /// The same sequence, `bases` further in the arena (see [`BlockArena::append`]).
    #[inline]
    fn shifted(self, bases: usize) -> Self {
        if self.is_inline() {
            self
        } else {
            Self(self.0 + ((bases as u64) << 16))
        }
    }
}

/// Append-only storage of ACGT-only sequences, packed with 2 bits per base.
#[derive(Clone, Default, Debug)]
pub(crate) struct BlockArena {
    words: Vec<u64>,
    /// Number of stored bases.
    len: usize,
}

impl BlockArena {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Append base codes (0..=3) as one sequence.
    fn push_codes(&mut self, codes: impl Iterator<Item = u8>) -> Result<Block> {
        let start = self.len;
        // Bases are packed in a word first, appended every 32 bases.
        let (mut bits, mut n) = (0u64, 0);
        for code in codes {
            debug_assert!(code < 4);
            bits |= (code as u64) << (2 * n);
            n += 1;
            if n == 32 {
                self.push_chunk(bits, n);
                (bits, n) = (0, 0);
            }
        }
        if self.len == start && n <= INLINE_BASES {
            return Ok(Block((bits << 16) | n as u64));
        }
        self.push_chunk(bits, n);
        let len = self.len - start;
        if len > u16::MAX as usize {
            bail!("block length {len} exceeds u16::MAX");
        }
        Ok(Block(((start as u64) << 16) | len as u64))
    }

    /// Store an ACGT-only sequence.
    pub(crate) fn push(&mut self, seq: &[u8]) -> Result<Block> {
        self.push_codes(seq.iter().map(|&b| encode(b)))
    }

    /// Store the reverse complement of an ACGT-only sequence.
    fn push_revcomp(&mut self, seq: &[u8]) -> Result<Block> {
        self.push_codes(seq.iter().rev().map(|&b| 3 - encode(b)))
    }

    /// Append the content of another arena, starting on a word boundary. Returns the
    /// number of bases by which the blocks of `other` must be [`Block::shifted`].
    fn append(&mut self, other: &BlockArena) -> usize {
        let shift = 32 * self.words.len();
        self.words.extend_from_slice(&other.words);
        self.len = shift + other.len;
        shift
    }

    /// Append `n <= 32` bases packed in `bits` as by [`Self::chunk`] (higher bits zero).
    #[inline]
    fn push_chunk(&mut self, bits: u64, n: usize) {
        let shift = 2 * (self.len % 32);
        if shift == 0 {
            if n > 0 {
                self.words.push(bits);
            }
        } else {
            *self.words.last_mut().unwrap() |= bits << shift;
            if shift + 2 * n > 64 {
                self.words.push(bits >> (64 - shift));
            }
        }
        self.len += n;
    }

    /// Copy a sequence from another arena.
    fn copy_from(&mut self, other: &BlockArena, block: Block) -> Result<Block> {
        if block.is_inline() {
            return Ok(block);
        }
        let start = self.len;
        for i in (0..block.len()).step_by(32) {
            let n = (block.len() - i).min(32);
            self.push_chunk(other.chunk(block.start() + i, n), n);
        }
        Ok(Block(((start as u64) << 16) | block.len() as u64))
    }

    #[inline]
    fn code(&self, i: usize) -> u8 {
        ((self.words[i / 32] >> (2 * (i % 32))) & 3) as u8
    }

    /// Base codes (0..=3) of a sequence.
    pub(crate) fn codes(&self, block: Block) -> impl Iterator<Item = u8> + '_ {
        (0..block.len()).map(move |i| {
            if block.is_inline() {
                ((block.0 >> (16 + 2 * i)) & 3) as u8
            } else {
                self.code(block.start() + i)
            }
        })
    }

    /// Uppercase bases of a sequence.
    pub(crate) fn bases(&self, block: Block) -> impl Iterator<Item = u8> + '_ {
        self.codes(block).map(|code| b"ACGT"[code as usize])
    }

    /// Decoded copy of a sequence (tests only).
    #[cfg(test)]
    pub(crate) fn get(&self, block: Block) -> Vec<u8> {
        self.bases(block).collect()
    }

    /// `n <= 32` bases starting at base `start`, packed in one word.
    #[inline]
    fn chunk(&self, start: usize, n: usize) -> u64 {
        let (word, shift) = (start / 32, 2 * (start % 32));
        let mut bits = self.words[word] >> shift;
        if shift != 0 && word + 1 < self.words.len() {
            bits |= self.words[word + 1] << (64 - shift);
        }
        if n < 32 {
            bits &= (1u64 << (2 * n)) - 1;
        }
        bits
    }

    /// Whether `a` (in `self`) and `b` (in `other`) hold the same sequence.
    fn same(&self, a: Block, other: &BlockArena, b: Block) -> bool {
        if a.is_inline() || b.is_inline() {
            return a == b;
        }
        a.len() == b.len()
            && (0..a.len()).step_by(32).all(|i| {
                let n = (a.len() - i).min(32);
                self.chunk(a.start() + i, n) == other.chunk(b.start() + i, n)
            })
    }
}

/// Left/right forward-directed anchor codes (in canonical orientation).
pub(crate) type AnchorPair = (u64, u64);

/// One occurrence of a variant group in one genome, once merged (8 bytes).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BlockOccurrence {
    pub(crate) isolate: u16,
    /// Index of the sequence in [`RawGroup::sequences`].
    pub(crate) sequence_id: u16,
    /// Genome coordinate of the middle of the sequence (see [`genome_position`]).
    pub(crate) pos: u32,
}

/// 0-based coordinate, in the genome's contigs concatenated in file order, of the
/// middle base of the region between two anchors (forward strand, so independent of
/// the orientation in which the sequence is stored). It is at most `length / 2`
/// bases from any variable site of the sequence.
#[inline]
fn genome_position(segment_start: usize, left_end: usize, right_start: usize) -> u32 {
    (segment_start + left_end + (right_start - left_end) / 2) as u32
}

#[derive(Debug, Default)]
pub(crate) struct RawGroup {
    /// Unique sequences, in first-occurrence order.
    pub(crate) sequences: Vec<Block>,
    /// Every occurrence; the order within each isolate is deterministic.
    pub(crate) occurrences: Vec<BlockOccurrence>,
    /// Number of distinct genomes with at least one occurrence.
    n_isolates: u16,
    /// More than `u16::MAX` distinct sequences: the group is discarded.
    overflow: bool,
}

/// Table of the groups with these anchors (see [`SHARDS`]).
#[inline]
fn shard_of(key: AnchorPair) -> usize {
    let mixed = (key.0 ^ key.1.rotate_left(32)).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    (mixed >> (u64::BITS - SHARDS.trailing_zeros())) as usize
}

/// One occurrence of a variant group found in a genome, before merging across genomes.
#[derive(Clone, Copy, Debug)]
struct LocalOccurrence {
    key: AnchorPair,
    pos: u32,
    seq: Block,
}

/// A shared `k`-mer of a segment: a possible anchor (8 bytes).
#[derive(Clone, Copy, Debug)]
struct KmerHit {
    /// Position of the first base of the k-mer in the segment.
    start: u32,
    /// Estimated number of genomes containing the k-mer.
    occupancy: u16,
}

impl KmerHit {
    /// The same k-mer seen on the reverse strand of a segment of length `len`.
    #[inline]
    fn mirror(self, len: usize, k: usize) -> Self {
        Self {
            start: (len - k) as u32 - self.start,
            occupancy: self.occupancy,
        }
    }
}

/// Keep the best start (highest occupancy, then left-most) per local cluster.
fn select_best_starts_per_local_cluster(starts: &[KmerHit]) -> Vec<KmerHit> {
    let mut selected = Vec::new();
    let Some((&first, rest)) = starts.split_first() else {
        return selected;
    };
    let (mut cluster_start, mut best) = (first.start, first);
    for &candidate in rest {
        if candidate.start < cluster_start + START_CLUSTER_SPAN {
            if candidate.occupancy > best.occupancy {
                best = candidate;
            }
        } else {
            selected.push(best);
            cluster_start = candidate.start;
            best = candidate;
        }
    }
    selected.push(best);
    selected
}

/// Pair each selected start with its best right anchor (highest occupancy, then
/// left-most) located at most `length_middle` bases downstream.
fn pair_anchors(
    starts: &[KmerHit],
    ends: &[KmerHit],
    k: usize,
    length_middle: usize,
    mut emit: impl FnMut(KmerHit, KmerHit),
) {
    let mut first_end = 0;
    for left in select_best_starts_per_local_cluster(starts) {
        let left_end = left.start as usize + k;
        while first_end < ends.len() && ends[first_end].start as usize <= left_end {
            first_end += 1;
        }
        let mut best: Option<KmerHit> = None;
        for &right in &ends[first_end..] {
            if right.start as usize - left_end > length_middle {
                break;
            }
            if best.is_none_or(|b| right.occupancy > b.occupancy) {
                best = Some(right);
            }
        }
        if let Some(right) = best {
            emit(left, right);
        }
    }
}

/// Reusable per-thread buffers for segment scanning.
#[derive(Default)]
struct Scratch {
    /// Canonical k-mers of the current batch and their occupancies.
    keys: Vec<u64>,
    occupancies: Vec<u16>,
    starts: Vec<KmerHit>,
    ends: Vec<KmerHit>,
    mirrored_starts: Vec<KmerHit>,
    mirrored_ends: Vec<KmerHit>,
    /// Start positions of the left and right anchors of each anchor pair found.
    pairs: Vec<(u32, u32)>,
}

struct SegmentContext<'a> {
    k: usize,
    cms: &'a CountMinSketch,
    min_needed: u32,
    length_middle: usize,
}

/// Find the anchor pairs of an ACGT-only segment, on both strands, as the
/// forward-strand start positions of the `(left, right)` anchors, sorted and
/// deduplicated.
fn find_segment_anchor_pairs(seg: &[u8], ctx: &SegmentContext, s: &mut Scratch) {
    let (k, min_needed) = (ctx.k, ctx.min_needed);
    s.starts.clear();
    s.ends.clear();
    s.pairs.clear();

    // Occupancies of the k-mers of the segment are looked up in batches. Start
    // candidates are shared k-mers followed by a lower-occupancy k-mer; end
    // candidates are shared k-mers preceded by a lower-occupancy k-mer.
    let mut roller = Roller::new(k);
    for &b in &seg[..k - 1] {
        roller.push(encode(b));
    }
    let mut previous: Option<KmerHit> = None;
    for (batch, chunk) in seg[k - 1..].chunks(BATCH).enumerate() {
        s.keys.clear();
        s.keys.extend(chunk.iter().map(|&b| {
            roller.push(encode(b));
            roller.canonical()
        }));
        s.occupancies.resize(chunk.len(), 0);
        ctx.cms.estimate_batch(&s.keys, &mut s.occupancies);
        for (i, &occupancy) in s.occupancies.iter().enumerate() {
            let hit = KmerHit {
                start: (batch * BATCH + i) as u32,
                occupancy,
            };
            if let Some(prev) = previous {
                if u32::from(prev.occupancy) >= min_needed && prev.occupancy > hit.occupancy {
                    s.starts.push(prev);
                }
                if u32::from(hit.occupancy) >= min_needed && prev.occupancy < hit.occupancy {
                    s.ends.push(hit);
                }
            }
            previous = Some(hit);
        }
    }

    // Forward strand.
    let pairs = &mut s.pairs;
    pair_anchors(&s.starts, &s.ends, k, ctx.length_middle, |l, r| {
        pairs.push((l.start, r.start))
    });

    // Reverse strand: the occupancy profile is mirrored, so forward ends become
    // reverse starts and vice versa. Pairs are mapped back to forward coordinates,
    // which makes the set of anchor pairs independent of the contig orientation.
    let len = seg.len();
    s.mirrored_starts.clear();
    s.mirrored_starts
        .extend(s.ends.iter().rev().map(|h| h.mirror(len, k)));
    s.mirrored_ends.clear();
    s.mirrored_ends
        .extend(s.starts.iter().rev().map(|h| h.mirror(len, k)));
    pair_anchors(
        &s.mirrored_starts,
        &s.mirrored_ends,
        k,
        ctx.length_middle,
        |l, r| pairs.push((r.mirror(len, k).start, l.mirror(len, k).start)),
    );

    pairs.sort_unstable();
    pairs.dedup();
}

/// Forward and reverse-complement codes of an ACGT-only k-mer.
fn kmer_codes(kmer: &[u8]) -> (u64, u64) {
    let fwd = kmer
        .iter()
        .fold(0, |code, &b| (code << 2) | encode(b) as u64);
    (fwd, revcomp_code(fwd, kmer.len()))
}

/// Extract the occurrences of one ACGT-only segment starting at genome coordinate `offset`.
fn extract_segment(
    seg: &[u8],
    offset: usize,
    ctx: &SegmentContext,
    scratch: &mut Scratch,
    occurrences: &mut [Vec<LocalOccurrence>],
    arena: &mut BlockArena,
) -> Result<()> {
    find_segment_anchor_pairs(seg, ctx, scratch);
    for &(left, right) in &scratch.pairs {
        let (left_start, right_start) = (left as usize, right as usize);
        let left_end = left_start + ctx.k;
        // Anchors are stored as positions only: their codes are computed here.
        let (left_fwd, left_rc) = kmer_codes(&seg[left_start..left_end]);
        let (right_fwd, right_rc) = kmer_codes(&seg[right_start..right_start + ctx.k]);
        let fwd_key = (left_fwd, right_fwd);
        let rev_key = (right_rc, left_rc);
        // Store each occurrence in the orientation with the smallest anchor pair so both
        // strands of the same locus share one key. Self-complementary pairs are
        // skipped as their orientation cannot be resolved.
        let middle = &seg[left_end..right_start];
        let (key, seq) = if fwd_key < rev_key {
            (fwd_key, arena.push(middle)?)
        } else if rev_key < fwd_key {
            (rev_key, arena.push_revcomp(middle)?)
        } else {
            continue;
        };
        occurrences[shard_of(key)].push(LocalOccurrence {
            key,
            pos: genome_position(offset, left_end, right_start),
            seq,
        });
    }
    Ok(())
}

struct IsolateExtraction {
    arena: BlockArena,
    /// Occurrences of the genome for each group table, in genome order.
    occurrences: Vec<Vec<LocalOccurrence>>,
    /// Genome coordinate of the end of each contig, in file order.
    contig_ends: Vec<u32>,
}

/// Merges genome extractions one at a time into the groups of one table.
///
/// Once `merged` genomes are in, a group seen in `s` of them can reach at most
/// `s + (n - merged)` genomes. Groups that can no longer reach `min_needed` are
/// dropped and new groups are no longer created when too few genomes remain. These
/// groups would be rejected by the support filter anyway, so pruning them changes
/// nothing but memory and time.
struct Merger {
    groups: HashMap<AnchorPair, RawGroup>,
    arena: BlockArena,
    n: usize,
    min_needed: usize,
    merged: usize,
}

impl Merger {
    fn new(n: usize, min_needed: usize) -> Self {
        Self {
            groups: HashMap::new(),
            arena: BlockArena::new(),
            n,
            min_needed,
            merged: 0,
        }
    }

    /// Merge the occurrences of genome `isolate` that belong to this table (all of them in
    /// one call), their sequences being stored in `arena`.
    fn merge(
        &mut self,
        isolate: u16,
        occurrences: &[LocalOccurrence],
        arena: &BlockArena,
    ) -> Result<()> {
        // A group first seen now can reach at most `n - merged` genomes.
        let can_create = self.n - self.merged >= self.min_needed;
        for occurrence in occurrences {
            let group = if can_create {
                self.groups.entry(occurrence.key).or_default()
            } else if let Some(group) = self.groups.get_mut(&occurrence.key) {
                group
            } else {
                continue;
            };
            if group.overflow {
                continue;
            }
            // Occurrences of one genome are merged together, so a new genome for
            // this group is necessarily the current one.
            if group
                .occurrences
                .last()
                .is_none_or(|last| last.isolate != isolate)
            {
                group.n_isolates += 1;
            }
            let sequence_id = match group
                .sequences
                .iter()
                .position(|&b| self.arena.same(b, arena, occurrence.seq))
            {
                Some(id) => id,
                None => {
                    group
                        .sequences
                        .push(self.arena.copy_from(arena, occurrence.seq)?);
                    group.sequences.len() - 1
                }
            };
            let Ok(sequence_id) = u16::try_from(sequence_id) else {
                group.overflow = true;
                group.occurrences = Vec::new();
                group.sequences = Vec::new();
                continue;
            };
            group.occurrences.push(BlockOccurrence {
                isolate,
                sequence_id,
                pos: occurrence.pos,
            });
        }
        self.merged += 1;

        // Prune about 16 times over the run, and always after the last genome.
        let remaining = self.n - self.merged;
        let step = (self.n / 16).max(1);
        if remaining < self.min_needed && (self.merged.is_multiple_of(step) || remaining == 0) {
            let min_seen = self.min_needed - remaining;
            self.groups
                .retain(|_, g| !g.overflow && g.n_isolates as usize >= min_seen);
            self.groups.shrink_to_fit();
        }
        Ok(())
    }
}

/// Pass 1: report the distinct canonical k-mers of a genome, in batches of up to
/// [`BATCH`] keys.
fn count_isolate_kmers(
    input: &IsolateInput,
    k: usize,
    mut flush: impl FnMut(&[u64]),
) -> Result<()> {
    let mut reader = FastaReader::new(open_fasta(&input.path)?);
    let mut bloom = BloomFilter::new(BLOOM_BITS);
    let mut roller = Roller::new(k);
    // Canonical k-mers in genome order, then only those seen for the first time.
    let mut keys = Vec::with_capacity(BATCH);
    while let Some(rec) = reader.next() {
        let rec = rec?;
        roller.reset();
        for &b in rec.seq() {
            if b.is_ascii_whitespace() {
                continue;
            }
            match encode(b) {
                INVALID => roller.reset(),
                code => {
                    if roller.push(code) {
                        keys.push(roller.canonical());
                        if keys.len() == BATCH {
                            bloom.retain_unseen(&mut keys);
                            flush(&keys);
                            keys.clear();
                        }
                    }
                }
            }
        }
    }
    bloom.retain_unseen(&mut keys);
    flush(&keys);
    Ok(())
}

/// Pass 2: extract all occurrences of one genome.
fn extract_isolate_occurrences(
    input: &IsolateInput,
    ctx: &SegmentContext,
) -> Result<IsolateExtraction> {
    let mut reader = FastaReader::new(open_fasta(&input.path)?);
    let mut occurrences = vec![Vec::new(); SHARDS];
    let mut arena = BlockArena::new();
    let mut scratch = Scratch::default();
    let mut contig = Vec::new();
    let mut contig_ends = Vec::new();
    let min_segment_len = 2 * ctx.k + 1;
    // Contigs are placed end to end in one genome coordinate system; k-mers and
    // anchor pairs are still searched within each contig only.
    let mut contig_start = 0usize;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        contig.clear();
        contig.extend(
            rec.seq()
                .iter()
                .filter(|b| !b.is_ascii_whitespace())
                .map(|&b| normalize(b)),
        );
        let contig_end = contig_start + contig.len();
        if contig_end > u32::MAX as usize {
            bail!("genome {} exceeds 4.29 Gb", input.name);
        }
        // Ambiguous bases split contigs into independent ACGT-only segments.
        let mut offset = contig_start;
        for seg in contig.split(|&b| b == b'N') {
            if seg.len() >= min_segment_len {
                extract_segment(seg, offset, ctx, &mut scratch, &mut occurrences, &mut arena)?;
            }
            offset += seg.len() + 1;
        }
        contig_ends.push(contig_end as u32);
        contig_start = contig_end;
    }
    Ok(IsolateExtraction {
        arena,
        occurrences,
        contig_ends,
    })
}

/// Variant groups as extracted, handed to [`crate::group_sorting`].
pub(crate) struct RawExtractedGroups {
    pub(crate) isolate_names: Vec<String>,
    /// Blocks (middles) referenced by the occurrences of the groups.
    pub(crate) arena: BlockArena,
    /// Variant groups: occurrences grouped by canonical anchor pair (in no particular order).
    pub(crate) groups: Vec<(AnchorPair, RawGroup)>,
    /// Per genome: coordinate of the end of each contig in the concatenated genome.
    pub(crate) contig_ends: Vec<Vec<u32>>,
    pub(crate) k: usize,
    pub(crate) min_needed: usize,
    pub(crate) n_isolates: usize,
}

/// Extract the variant groups found in at least `min_needed` genomes.
pub(crate) fn extract_groups(
    inputs: &[IsolateInput],
    k: usize,
    min_needed: usize,
    length_middle: usize,
    num_threads: usize,
    log: &RunLog,
) -> Result<RawExtractedGroups> {
    let n = inputs.len();
    if n > u16::MAX as usize {
        bail!("too many genomes ({n}); maximum is {}", u16::MAX);
    }
    let min_needed = min_needed as u32;
    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new().num_threads(n_threads).build()?;

    // Pass 1: count genome occupancies. Threads add k-mers to the sketch region by
    // region (see `SharedSketch`), so the counts do not depend on the number of threads.
    say!(log, " . pass 1: estimate k-mer occupancies");
    let progress = Progress::new(n);
    let mut cms = CountMinSketch::new(CMS_LINES);
    let shared = cms.shared();
    let result = pool.install(|| {
        inputs.par_iter().try_for_each(|input| {
            let mut writer = shared.writer();
            count_isolate_kmers(input, k, |keys| writer.add(keys))?;
            writer.finish();
            progress.increment();
            Ok::<_, anyhow::Error>(())
        })
    });
    drop(shared);
    progress.finish();
    result?;

    // Pass 2: each thread extracts the occurrences of a genome, then merges them into
    // the group tables, one table at a time.
    say!(log, " . pass 2: extract sequences");
    let ctx = SegmentContext {
        k,
        cms: &cms,
        min_needed,
        length_middle,
    };
    let progress = Progress::new(n);
    let mergers: Vec<Mutex<Merger>> = (0..SHARDS)
        .map(|_| Mutex::new(Merger::new(n, min_needed as usize)))
        .collect();
    let contig_ends = pool.install(|| {
        inputs
            .par_iter()
            .enumerate()
            .map(|(isolate, input)| {
                let extraction = extract_isolate_occurrences(input, &ctx)?;
                for (merger, occurrences) in mergers.iter().zip(&extraction.occurrences) {
                    merger.lock().expect("merger mutex poisoned").merge(
                        isolate as u16,
                        occurrences,
                        &extraction.arena,
                    )?;
                }
                progress.increment();
                Ok(extraction.contig_ends)
            })
            .collect::<Result<Vec<_>>>()
    });
    progress.finish();
    let contig_ends = contig_ends?;
    drop(cms);

    // Gather the groups of all tables, and their sequences in one arena.
    let mut arena = BlockArena::new();
    let mut groups = Vec::new();
    for merger in mergers {
        let merger = merger.into_inner().expect("merger mutex poisoned");
        let shift = arena.append(&merger.arena);
        groups.extend(merger.groups.into_iter().map(|(key, mut group)| {
            for block in &mut group.sequences {
                *block = block.shifted(shift);
            }
            (key, group)
        }));
    }

    Ok(RawExtractedGroups {
        isolate_names: inputs.iter().map(|x| x.name.clone()).collect(),
        arena,
        groups,
        contig_ends,
        k,
        min_needed: min_needed as usize,
        n_isolates: n,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::revcomp;

    fn hit(start: u32, occupancy: u16) -> KmerHit {
        KmerHit { start, occupancy }
    }

    #[test]
    fn block_arena_packs_forward_and_reverse_complement_sequences() {
        let mut arena = BlockArena::new();
        let a = arena.push(b"ACGTA").unwrap();
        let b = arena.push_revcomp(b"AACGC").unwrap();
        assert_eq!(arena.get(a), b"ACGTA");
        assert_eq!(arena.get(b), b"GCGTT");
        assert_eq!((a.len(), b.len()), (5, 5));
        assert_eq!(arena.codes(b).collect::<Vec<_>>(), [2, 1, 2, 3, 3]);
    }

    #[test]
    fn block_arena_compares_across_word_boundaries() {
        // Sequences of various lengths at various bit offsets in two arenas.
        let seqs: Vec<Vec<u8>> = (0..40)
            .map(|i| {
                (0..(i * 7 % 75))
                    .map(|j| b"ACGT"[(i * j + j / 3) % 4])
                    .collect()
            })
            .collect();
        let (mut a, mut b) = (BlockArena::new(), BlockArena::new());
        b.push(&[b'C'; 37]).unwrap(); // shift every stored block of `b`
        let blocks_a: Vec<Block> = seqs.iter().map(|s| a.push(s).unwrap()).collect();
        let blocks_b: Vec<Block> = seqs.iter().map(|s| b.push(s).unwrap()).collect();
        for (i, &x) in blocks_a.iter().enumerate() {
            assert_eq!(a.get(x), seqs[i]);
            for (j, &y) in blocks_b.iter().enumerate() {
                assert_eq!(a.same(x, &b, y), seqs[i] == seqs[j], "{i} vs {j}");
            }
        }
        let mut c = BlockArena::new();
        for (i, &x) in blocks_a.iter().enumerate() {
            let copied = c.copy_from(&a, x).unwrap();
            assert!(c.same(copied, &a, x));
            assert_eq!(c.get(copied), seqs[i]);
        }
        assert!(a.push(&vec![b'A'; 65536]).is_err());
    }

    #[test]
    fn short_sequences_take_no_arena_space() {
        let mut arena = BlockArena::new();
        let short = arena.push(&[b'T'; INLINE_BASES]).unwrap();
        let empty = arena.push(b"").unwrap();
        assert_eq!(arena.len, 0);
        let long = arena.push(&[b'T'; INLINE_BASES + 1]).unwrap();
        assert_eq!(arena.len, INLINE_BASES + 1);
        assert_eq!(arena.get(short), [b'T'; INLINE_BASES]);
        assert_eq!((arena.get(empty), short.shifted(9)), (vec![], short));
        assert_eq!(arena.get(long), [b'T'; INLINE_BASES + 1]);
    }

    #[test]
    fn starts_are_clustered_and_rights_prefer_occupancy_then_position() {
        let starts = [hit(0, 5), hit(3, 6), hit(9, 6), hit(20, 5)];
        assert_eq!(
            select_best_starts_per_local_cluster(&starts)
                .iter()
                .map(|h| h.start)
                .collect::<Vec<_>>(),
            vec![3, 20]
        );
        let ends = [hit(8, 5), hit(10, 7), hit(12, 7), hit(40, 9)];
        let mut pairs = Vec::new();
        pair_anchors(&[hit(0, 5)], &ends, 3, 10, |l, r| {
            pairs.push((l.start, r.start))
        });
        // 10 beats 8 on occupancy and 12 on position; 40 is beyond length_middle.
        assert_eq!(pairs, vec![(0, 10)]);
    }

    fn mirror_roundtrip(len: usize, k: usize) {
        let h = hit(4, 3);
        let m = h.mirror(len, k);
        // The k-mer starts where it ended on the other strand.
        assert_eq!(m.start as usize, len - (4 + k));
        let back = m.mirror(len, k);
        assert_eq!((back.start, back.occupancy), (h.start, h.occupancy));
    }

    #[test]
    fn kmer_codes_are_forward_and_reverse_complement() {
        let kmer = b"ACCGTTAGCATTGCA";
        let (fwd, rc) = kmer_codes(kmer);
        assert_eq!(fwd, crate::dna::pack(kmer));
        assert_eq!(rc, crate::dna::pack(&revcomp(kmer)));
        assert_eq!(kmer_codes(&revcomp(kmer)), (rc, fwd));
    }

    #[test]
    fn mirroring_is_an_involution() {
        mirror_roundtrip(30, 5);
        mirror_roundtrip(9, 5);
    }

    /// CMS where every k-mer of `shared` sequences reaches `n` and others stay at 0.
    fn sketch(shared: &[&[u8]], k: usize, n: u32) -> CountMinSketch {
        let mut cms = CountMinSketch::new(1 << 12);
        let mut keys = hashbrown::HashSet::new();
        for seq in shared {
            let mut roller = Roller::new(k);
            for &b in *seq {
                if roller.push(encode(b)) {
                    keys.insert(roller.canonical());
                }
            }
        }
        let keys: Vec<u64> = keys.into_iter().collect();
        for _ in 0..n {
            cms.increment_batch(&keys);
        }
        cms
    }

    fn extract(seg: &[u8], ctx: &SegmentContext) -> (Vec<LocalOccurrence>, BlockArena) {
        let mut groups = vec![Vec::new(); SHARDS];
        let mut arena = BlockArena::new();
        extract_segment(
            seg,
            7,
            ctx,
            &mut Scratch::default(),
            &mut groups,
            &mut arena,
        )
        .unwrap();
        (groups.concat(), arena)
    }

    const LEFT: &[u8] = b"GATTACACCGTAGGCATTACGGA";
    const RIGHT: &[u8] = b"TTGACCGGTAACGTCAGGATCCA";

    #[test]
    fn snp_occurrence_is_identical_on_both_strands() {
        let k = 7;
        let reference = [LEFT, b"C", RIGHT].concat();
        let variant = [LEFT, b"T", RIGHT].concat();
        // Only the reference k-mers are shared, so the SNP region forms one occurrence.
        let cms = sketch(&[&reference], k, 2);
        let ctx = SegmentContext {
            k,
            cms: &cms,
            min_needed: 2,
            length_middle: 20,
        };
        let (fwd_groups, fwd_arena) = extract(&variant, &ctx);
        let (rev_groups, rev_arena) = extract(&revcomp(&variant), &ctx);
        assert_eq!((fwd_groups.len(), rev_groups.len()), (1, 1));
        let (a, b) = (fwd_groups[0], rev_groups[0]);
        assert_eq!(a.key, b.key);
        // The stored middle is the SNP alone, in the same canonical orientation.
        let stored = fwd_arena.get(a.seq);
        assert_eq!(stored, rev_arena.get(b.seq));
        assert!(stored == b"T" || stored == b"A");
        // Both point at the SNP on the forward strand (segment offset is 7).
        let snp = LEFT.len();
        assert_eq!(a.pos as usize, 7 + snp);
        assert_eq!(b.pos as usize, 7 + variant.len() - 1 - snp);
    }

    /// Merge the occurrences `(anchors, position, sequence)` of genome `isolate`.
    fn merge(merger: &mut Merger, isolate: u16, rows: &[(AnchorPair, u32, &[u8])]) {
        let mut arena = BlockArena::new();
        let occurrences: Vec<LocalOccurrence> = rows
            .iter()
            .map(|&(key, pos, seq)| LocalOccurrence {
                key,
                pos,
                seq: arena.push(seq).unwrap(),
            })
            .collect();
        merger.merge(isolate, &occurrences, &arena).unwrap();
    }

    #[test]
    fn merge_interns_sequences_and_keeps_locations() {
        let key = (3, 5);
        // Long enough to be stored in the arena (see `INLINE_BASES`).
        let (acg, act) = (
            &[b"ACG" as &[u8]; 10].concat(),
            &[b"ACT" as &[u8]; 10].concat(),
        );
        let mut merger = Merger::new(3, 1);
        merge(&mut merger, 0, &[(key, 10, acg), (key, 11, acg)]);
        merge(&mut merger, 1, &[(key, 20, act)]);
        merge(&mut merger, 2, &[(key, 30, acg)]);
        let group = &merger.groups[&key];
        assert_eq!(group.sequences.len(), 2);
        assert_eq!(group.n_isolates, 3);
        assert_eq!(merger.arena.len, 60);
        assert_eq!(
            group
                .occurrences
                .iter()
                .map(|o| (o.isolate, o.pos, o.sequence_id))
                .collect::<Vec<_>>(),
            vec![(0, 10, 0), (0, 11, 0), (1, 20, 1), (2, 30, 0)]
        );
    }

    #[test]
    fn appended_arenas_keep_their_sequences() {
        let (mut a, mut b) = (BlockArena::new(), BlockArena::new());
        // More than one word, so that the arenas hold a partial last word.
        let first_seq: Vec<u8> = b"ACGTACGTA".iter().cycle().take(40).copied().collect();
        let first = a.push(&first_seq).unwrap();
        let long: Vec<u8> = b"GATTACA".iter().cycle().take(45).copied().collect();
        let blocks = [b.push(b"TTG").unwrap(), b.push(&long).unwrap()];
        let shift = a.append(&b);
        assert_eq!(shift, 64);
        assert_eq!(a.get(first), first_seq);
        assert_eq!(a.get(blocks[0].shifted(shift)), b"TTG");
        assert_eq!(a.get(blocks[1].shifted(shift)), long);
        // The arena can still grow, and be appended to an empty one.
        let last = a.push(b"CCA").unwrap();
        assert_eq!(a.get(last), b"CCA");
        assert_eq!(a.get(blocks[1].shifted(shift)), long);
        let mut c = BlockArena::new();
        assert_eq!(c.append(&a), 0);
        assert_eq!(c.get(last), b"CCA");
    }

    #[test]
    fn anchor_pairs_are_spread_over_all_tables() {
        let mut seen = [false; SHARDS];
        for left in 0..200u64 {
            for right in 0..200u64 {
                seen[shard_of((left * 7919, right * 104_729))] = true;
            }
        }
        assert!(seen.iter().all(|&s| s));
    }

    #[test]
    fn positions_are_middle_of_region_in_concatenated_genome() {
        // Middle of [left_end, right_start), shifted by the segment's genome offset.
        assert_eq!(genome_position(0, 10, 11), 10);
        assert_eq!(genome_position(100, 10, 14), 112);
        assert_eq!(genome_position(100, 10, 15), 112);

        // SNP occurrence in the second contig of a genome, after an N-split segment.
        let k = 7;
        let reference = [LEFT, b"C", RIGHT].concat();
        let variant = [LEFT, b"T", RIGHT].concat();
        let cms = sketch(&[&reference], k, 2);
        let ctx = SegmentContext {
            k,
            cms: &cms,
            min_needed: 2,
            length_middle: 20,
        };
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("g.fna");
        let second = [b"ACGTANNAC".as_slice(), &variant].concat();
        let mut fasta = b">c0 first\nACGTACGT\nAC\n>c1\n".to_vec();
        fasta.extend_from_slice(&second[..20]);
        fasta.push(b'\n');
        fasta.extend_from_slice(&second[20..]);
        fasta.extend_from_slice(b"\n>empty\n>c3\nACGT\n");
        std::fs::write(&path, fasta).unwrap();
        let input = IsolateInput {
            name: "g".into(),
            path,
        };
        let e = extract_isolate_occurrences(&input, &ctx).unwrap();
        let (c0, c1) = (10, second.len() as u32);
        assert_eq!(e.contig_ends, vec![c0, c0 + c1, c0 + c1, c0 + c1 + 4]);
        let occurrences = e.occurrences.concat();
        assert_eq!(occurrences.len(), 1);
        // SNP index in the concatenated genome: contig 0, then "ACGTANNAC", then LEFT.
        assert_eq!(occurrences[0].pos as usize, 10 + 9 + LEFT.len());
    }

    #[test]
    fn merge_prunes_groups_that_cannot_reach_the_threshold() {
        // 4 genomes, groups need 3 of them.
        let (core, early, late) = ((1, 1), (2, 2), (3, 3));
        let seq: &[u8] = b"ACGT";
        let mut merger = Merger::new(4, 3);
        merge(&mut merger, 0, &[(core, 0, seq), (early, 0, seq)]);
        merge(&mut merger, 1, &[(core, 0, seq)]);
        // 2 genomes left: `early` (1 genome) can still reach 3; `late` is created.
        assert!(merger.groups.contains_key(&early));
        merge(&mut merger, 2, &[(core, 0, seq), (late, 0, seq)]);
        // 1 genome left: `early` can reach 2 at most and `late` 2 at most.
        assert!(merger.groups.contains_key(&core));
        assert!(!merger.groups.contains_key(&early));
        assert!(!merger.groups.contains_key(&late));
        // Too late to create new groups; core keeps growing.
        merge(&mut merger, 3, &[(core, 0, seq), (late, 0, seq)]);
        assert_eq!(merger.groups.len(), 1);
        assert_eq!(merger.groups[&core].n_isolates, 4);
    }

    #[test]
    fn counting_is_strand_independent_and_splits_at_ambiguous_bases() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("g.fna");
        let seq = b"ACGTTGCAAGGCTTACGATCGRTCGGATCCATGC";
        let rc = revcomp(seq);
        let count = |data: &[u8]| {
            std::fs::write(&path, data).unwrap();
            let input = IsolateInput {
                name: "g".into(),
                path: path.clone(),
            };
            let mut kmers = Vec::new();
            count_isolate_kmers(&input, 5, |keys| kmers.extend_from_slice(keys)).unwrap();
            kmers.sort_unstable();
            kmers
        };
        let fwd = count(&[b">c1\n".as_slice(), seq, b"\n"].concat());
        let rev = count(&[b">c1\n".as_slice(), &rc, b"\n"].concat());
        assert_eq!(fwd, rev);
        // 20 + 12 valid bases on either side of R: 16 + 8 windows (some may repeat).
        assert!(fwd.len() <= 24 && fwd.len() >= 20);
        let wrapped = count(b">c1\nACGTTGCAAG\nGCTTACGATCG\r\nRTCGGATCCATGC\n");
        assert_eq!(wrapped, fwd);
    }
}
