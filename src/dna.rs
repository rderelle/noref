//! Nucleotide encoding and rolling canonical k-mers.
//!
//! Bases are packed with 2 bits (`A=0, C=1, G=2, T=3`), so the complement of a
//! code is `3 - code`. Every other byte (IUPAC ambiguity codes, `N`, gaps...) is
//! treated as missing data.

/// Code returned by [`encode`] for any non-ACGT byte.
pub(crate) const INVALID: u8 = 4;

/// Largest k-mer that fits in a `u64` with 2 bits per base.
pub(crate) const MAX_K: usize = 32;

static CODES: [u8; 256] = {
    let mut codes = [INVALID; 256];
    let bases = *b"ACGT";
    let mut i = 0;
    while i < 4 {
        codes[bases[i] as usize] = i as u8;
        codes[bases[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    codes
};

/// Encode one byte as a 2-bit base code, or [`INVALID`].
#[inline]
pub(crate) fn encode(b: u8) -> u8 {
    CODES[b as usize]
}

/// Uppercase ACGT, or `N` for anything else.
#[inline]
pub(crate) fn normalize(b: u8) -> u8 {
    match encode(b) {
        INVALID => b'N',
        c => b"ACGT"[c as usize],
    }
}

/// Complement of an uppercase base (`N` and gaps are left unchanged).
#[cfg(test)]
pub(crate) fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Rolling forward and reverse-complement k-mer codes over a stream of bases.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Roller {
    k: usize,
    mask: u64,
    rc_shift: u32,
    have: usize,
    /// Forward-strand code of the last `k` bases.
    pub(crate) fwd: u64,
    /// Reverse-complement code of the last `k` bases.
    pub(crate) rc: u64,
}

impl Roller {
    pub(crate) fn new(k: usize) -> Self {
        debug_assert!((1..=MAX_K).contains(&k));
        Self {
            k,
            mask: if k == MAX_K {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            },
            rc_shift: 2 * (k as u32 - 1),
            have: 0,
            fwd: 0,
            rc: 0,
        }
    }

    #[inline]
    pub(crate) fn reset(&mut self) {
        self.have = 0;
        self.fwd = 0;
        self.rc = 0;
    }

    /// Add one base code (0..=3); returns true once a complete k-mer is available.
    #[inline]
    pub(crate) fn push(&mut self, code: u8) -> bool {
        debug_assert!(code < 4);
        self.fwd = ((self.fwd << 2) | code as u64) & self.mask;
        self.rc = (self.rc >> 2) | (((3 - code) as u64) << self.rc_shift);
        if self.have < self.k {
            self.have += 1;
        }
        self.have == self.k
    }

    /// Strand-independent k-mer code.
    #[inline]
    pub(crate) fn canonical(&self) -> u64 {
        self.fwd.min(self.rc)
    }
}

/// Reverse-complement code of the code of a `k`-mer.
#[inline]
pub(crate) fn revcomp_code(code: u64, k: usize) -> u64 {
    // Complement all 32 bases, reverse the order of the 2-bit groups (bytes, then
    // nibbles, then bases), and drop the bases beyond the k-mer.
    let x = (!code).swap_bytes();
    let x = ((x >> 4) & 0x0f0f_0f0f_0f0f_0f0f) | ((x & 0x0f0f_0f0f_0f0f_0f0f) << 4);
    let x = ((x >> 2) & 0x3333_3333_3333_3333) | ((x & 0x3333_3333_3333_3333) << 2);
    x >> (64 - 2 * k)
}

#[cfg(test)]
pub(crate) fn pack(seq: &[u8]) -> u64 {
    seq.iter().fold(0, |v, &b| (v << 2) | encode(b) as u64)
}

#[cfg(test)]
pub(crate) fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoding_is_case_insensitive_and_ambiguity_is_invalid() {
        assert_eq!(b"ACGTacgt".map(encode), [0, 1, 2, 3, 0, 1, 2, 3]);
        for &b in b"NnRYKMSWBDHVX-*U." {
            assert_eq!(encode(b), INVALID);
            assert_eq!(normalize(b), b'N');
        }
        assert_eq!(b"acgtR".map(normalize), *b"ACGTN");
    }

    #[test]
    fn roller_matches_direct_forward_and_reverse_complement_codes() {
        let seq = b"ACGTTGCAAGGCTTACGATCGATCGGATCCATGCAAATTTGGGCCC";
        for k in [1, 5, 15, 31, 32] {
            let mut roller = Roller::new(k);
            for (i, &b) in seq.iter().enumerate() {
                let full = roller.push(encode(b));
                assert_eq!(full, i + 1 >= k);
                if full {
                    let kmer = &seq[i + 1 - k..=i];
                    assert_eq!(roller.fwd, pack(kmer));
                    assert_eq!(roller.rc, pack(&revcomp(kmer)));
                    assert_eq!(revcomp_code(roller.fwd, k), roller.rc);
                }
            }
            roller.reset();
            assert!(!roller.push(0) || k == 1);
        }
    }

    #[test]
    fn canonical_code_is_shared_by_both_strands() {
        let kmer = b"ACCGTTAGCATTGCA";
        let canonical = |s: &[u8]| {
            let mut roller = Roller::new(s.len());
            s.iter().for_each(|&b| {
                roller.push(encode(b));
            });
            roller.canonical()
        };
        assert_eq!(canonical(kmer), canonical(&revcomp(kmer)));
    }
}
