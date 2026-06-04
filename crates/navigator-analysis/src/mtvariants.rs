//! Derive mtDNA variants by comparing a sample sequence to a reference (rCRS, the revised
//! Cambridge Reference Sequence, NC_012920.1, 16,569 bp). Substitutions *and* indels are
//! derived via a banded global alignment (Needleman–Wunsch, unit edit costs); the band is
//! sized to the length difference plus slack, which is ample for mtDNA's few small indels.
//! Positions where either base is `N`/ambiguous are not called as substitutions.
//!
//! Pure: callers provide both sequences; the rCRS reference itself is supplied externally.

/// The kind of difference from the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MtVariantKind {
    Substitution,
    Insertion,
    Deletion,
}

/// A single mtDNA variant relative to the reference (1-based rCRS coordinates).
///
/// - substitution: `reference`/`alternate` are the single bases (e.g. A→G at 263).
/// - insertion: `reference` empty, `alternate` the inserted bases; `position` is the rCRS
///   base they follow (mtDNA `.1` convention, e.g. `315.1C`).
/// - deletion: `reference` the deleted bases, `alternate` empty; `position` is the first
///   deleted rCRS base (e.g. `8281-8289d`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MtVariant {
    pub position: i64,
    pub reference: String,
    pub alternate: String,
    pub kind: MtVariantKind,
}

impl MtVariant {
    /// Compact mtDNA notation, e.g. `263A>G`, `315.1C`, `8281-8289d`.
    pub fn notation(&self) -> String {
        match self.kind {
            MtVariantKind::Substitution => format!("{}{}>{}", self.position, self.reference, self.alternate),
            MtVariantKind::Insertion => format!("{}.1{}", self.position, self.alternate),
            MtVariantKind::Deletion => {
                let len = self.reference.len() as i64;
                if len <= 1 {
                    format!("{}d", self.position)
                } else {
                    format!("{}-{}d", self.position, self.position + len - 1)
                }
            }
        }
    }
}

fn is_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

const INF: i32 = i32::MAX / 2;

/// One alignment step, in forward (5'→3') order. Indices are 0-based into ref/sample.
enum Op {
    /// Aligned column (match or mismatch): ref[ri] vs sample[sj].
    Diag { ri: usize, sj: usize },
    /// Reference base deleted (gap in sample).
    Del { ri: usize },
    /// Sample base inserted (gap in reference); `after` ref bases precede it.
    Ins { sj: usize, after: usize },
}

/// The revised Cambridge Reference Sequence (rCRS, NC_012920.1, 16,569 bp), bundled as the
/// canonical mtDNA reference — used to self-generate a coordinate map onto a divergent
/// assembly's `chrM` (e.g. CHM13) so the rCRS-coordinate haplotree can be lifted onto it.
pub fn rcrs() -> &'static str {
    include_str!("../resources/rcrs.seq").trim()
}

/// Derive variants of `sample` relative to `reference` via banded global alignment.
pub fn derive(reference: &str, sample: &str) -> Vec<MtVariant> {
    let r = reference.as_bytes();
    let s = sample.as_bytes();
    if r.is_empty() || s.is_empty() {
        return Vec::new();
    }
    let ops = align_ops(r, s, 32);
    build_variants(r, s, &ops)
}

/// Aligned `(reference_index, sample_index)` pairs (0-based) from a banded global alignment —
/// i.e. the position map between two sequences. A wider band than [`derive`] tolerates the
/// extra indels seen aligning two *different* references (rCRS vs a CHM13-style `chrM`).
pub fn align_positions(reference: &str, sample: &str) -> Vec<(usize, usize)> {
    let r = reference.as_bytes();
    let s = sample.as_bytes();
    if r.is_empty() || s.is_empty() {
        return Vec::new();
    }
    align_ops(r, s, 256)
        .into_iter()
        .filter_map(|op| match op {
            Op::Diag { ri, sj } => Some((ri, sj)),
            _ => None,
        })
        .collect()
}

/// Banded Needleman–Wunsch (unit edit costs) of `r` vs `s`; returns the forward (5'→3') op
/// sequence. `slack` widens the band beyond the length difference to admit internal indels.
fn align_ops(r: &[u8], s: &[u8], slack: usize) -> Vec<Op> {
    let m = r.len();
    let n = s.len();

    let diff = n as isize - m as isize;
    let band = diff.unsigned_abs() + slack; // slack for internal indels
    let dmin = diff.min(0) - band as isize;
    let dmax = diff.max(0) + band as isize;
    let width = (dmax - dmin + 1) as usize;
    let col = |d: isize| (d - dmin) as usize;
    let in_band = |d: isize| d >= dmin && d <= dmax;

    // dp[i][col(d)] = min edit cost to align r[..i] with s[..i+d]; tb the chosen move.
    let mut dp = vec![vec![INF; width]; m + 1];
    let mut tb = vec![vec![0u8; width]; m + 1]; // 0=diag, 1=del(up), 2=ins(left)

    for d in dmin..=dmax {
        let j = d; // i = 0 → j = d
        if (0..=n as isize).contains(&j) {
            dp[0][col(d)] = j as i32; // all insertions
            tb[0][col(d)] = 2;
        }
    }
    dp[0][col(0)] = 0;

    for i in 1..=m {
        for d in dmin..=dmax {
            let j = i as isize + d;
            if j < 0 || j > n as isize {
                continue;
            }
            let j = j as usize;
            let mut best = INF;
            let mut mv = 0u8;
            if j >= 1 {
                let prev = dp[i - 1][col(d)]; // diag: (i-1, j-1), same d
                if prev < INF {
                    let cost = if r[i - 1].eq_ignore_ascii_case(&s[j - 1]) { 0 } else { 1 };
                    if prev + cost < best {
                        best = prev + cost;
                        mv = 0;
                    }
                }
            }
            if in_band(d + 1) {
                let prev = dp[i - 1][col(d + 1)]; // del: (i-1, j), d+1
                if prev < INF && prev + 1 < best {
                    best = prev + 1;
                    mv = 1;
                }
            }
            if j >= 1 && in_band(d - 1) {
                let prev = dp[i][col(d - 1)]; // ins: (i, j-1), d-1
                if prev < INF && prev + 1 < best {
                    best = prev + 1;
                    mv = 2;
                }
            }
            dp[i][col(d)] = best;
            tb[i][col(d)] = mv;
        }
    }

    // Traceback from (m, n).
    let mut i = m;
    let mut j = n;
    let mut ops: Vec<Op> = Vec::new();
    while i > 0 || j > 0 {
        let d = j as isize - i as isize;
        let mv = if i == 0 {
            2 // only insertions remain
        } else if j == 0 {
            1 // only deletions remain
        } else {
            tb[i][col(d)]
        };
        match mv {
            0 => {
                ops.push(Op::Diag { ri: i - 1, sj: j - 1 });
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(Op::Del { ri: i - 1 });
                i -= 1;
            }
            _ => {
                ops.push(Op::Ins { sj: j - 1, after: i });
                j -= 1;
            }
        }
    }
    ops.reverse();
    ops
}

fn up(b: u8) -> char {
    b.to_ascii_uppercase() as char
}

/// Right-align (3'-shift) a deletion of `len` reference bases starting at 0-based `start`,
/// the mtDNA/Phylotree convention: slide through a tandem/homopolymer run to the rightmost
/// equivalent position. Returns the new start and the (uppercased) deleted bases.
fn right_align_deletion(r: &[u8], mut start: usize, len: usize) -> (usize, String) {
    while start + len < r.len() && r[start].eq_ignore_ascii_case(&r[start + len]) {
        start += 1;
    }
    let deleted = r[start..start + len].iter().map(|&b| up(b)).collect();
    (start, deleted)
}

/// Right-align (3'-shift) an insertion of `ins` whose insertion point is before 0-based
/// reference index `after`. Sliding right rotates the inserted string. Returns the new
/// `after` (the rCRS position the insertion follows, 1-based) and the rotated bases.
fn right_align_insertion(r: &[u8], mut after: usize, ins: &str) -> (usize, String) {
    let mut bases: Vec<u8> = ins.bytes().map(|b| b.to_ascii_uppercase()).collect();
    while after < r.len() && r[after].to_ascii_uppercase() == bases[0] {
        bases.rotate_left(1);
        after += 1;
    }
    (after, String::from_utf8(bases).expect("ascii"))
}

fn build_variants(r: &[u8], s: &[u8], ops: &[Op]) -> Vec<MtVariant> {
    let mut variants = Vec::new();
    let mut k = 0;
    while k < ops.len() {
        match ops[k] {
            Op::Diag { ri, sj } => {
                let (rb, sb) = (r[ri], s[sj]);
                if !rb.eq_ignore_ascii_case(&sb) && is_base(rb) && is_base(sb) {
                    variants.push(MtVariant {
                        position: (ri + 1) as i64,
                        reference: up(rb).to_string(),
                        alternate: up(sb).to_string(),
                        kind: MtVariantKind::Substitution,
                    });
                }
                k += 1;
            }
            Op::Del { ri } => {
                let mut len = 0;
                while let Some(Op::Del { .. }) = ops.get(k) {
                    len += 1;
                    k += 1;
                }
                // Right-align through any tandem/homopolymer run (mtDNA convention).
                let (start, deleted) = right_align_deletion(r, ri, len);
                variants.push(MtVariant {
                    position: (start + 1) as i64,
                    reference: deleted,
                    alternate: String::new(),
                    kind: MtVariantKind::Deletion,
                });
            }
            Op::Ins { after, .. } => {
                let anchor = after;
                let mut inserted = String::new();
                while let Some(Op::Ins { sj, .. }) = ops.get(k) {
                    inserted.push(up(s[*sj]));
                    k += 1;
                }
                let (after, inserted) = right_align_insertion(r, anchor, &inserted);
                variants.push(MtVariant {
                    position: after as i64,
                    reference: String::new(),
                    alternate: inserted,
                    kind: MtVariantKind::Insertion,
                });
            }
        }
    }
    variants
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sub(pos: i64, r: &str, a: &str) -> MtVariant {
        MtVariant { position: pos, reference: r.into(), alternate: a.into(), kind: MtVariantKind::Substitution }
    }

    #[test]
    fn finds_substitutions_and_skips_ns() {
        //           1234567
        let refseq = "ACGTACG";
        let sample = "AGGTNCG"; // pos2 C>G; pos5 A vs N -> skipped
        let v = derive(refseq, sample);
        assert_eq!(v, vec![sub(2, "C", "G")]);
        assert_eq!(v[0].notation(), "2C>G");
    }

    #[test]
    fn identical_sequences_have_no_variants() {
        assert!(derive("ACGTACGT", "acgtacgt").is_empty()); // case-insensitive
    }

    #[test]
    fn detects_an_insertion() {
        // sample has an extra C inserted after ref position 4 (the 4th base).
        let v = derive("ACGTACGT", "ACGTCACGT");
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].kind, MtVariantKind::Insertion);
        assert_eq!(v[0].notation(), "4.1C");
    }

    #[test]
    fn detects_a_deletion() {
        // ref bases 5-6 (AC) are absent from the sample.
        let v = derive("ACGTACGT", "ACGTGT");
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].kind, MtVariantKind::Deletion);
        assert_eq!(v[0].reference, "AC");
        assert_eq!(v[0].notation(), "5-6d");
    }

    #[test]
    fn right_aligns_insertion_in_a_homopolymer() {
        // C-run at positions 2-4; the extra C is reported after the last C (315.1C-style).
        let v = derive("GCCCG", "GCCCCG");
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].kind, MtVariantKind::Insertion);
        assert_eq!(v[0].notation(), "4.1C"); // rightmost C of the run is position 4
    }

    #[test]
    fn right_aligns_deletion_in_a_homopolymer() {
        // deleting one C from the run is reported at the rightmost C.
        let v = derive("GCCCG", "GCCG");
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].kind, MtVariantKind::Deletion);
        assert_eq!(v[0].notation(), "4d");
    }

    #[test]
    fn substitution_plus_indel() {
        // pos2 C>T, and a single-base deletion of ref pos 6 (C).
        let v = derive("ACGTACGT", "ATGTAGT");
        let kinds: Vec<_> = v.iter().map(|x| x.kind).collect();
        assert!(kinds.contains(&MtVariantKind::Substitution));
        assert!(kinds.contains(&MtVariantKind::Deletion));
        assert!(v.iter().any(|x| x.notation() == "2C>T"));
        assert!(v.iter().any(|x| x.notation() == "6d"));
    }

    #[test]
    fn align_positions_maps_coordinates_across_an_indel() {
        // sample = ref with "TT" inserted after ref index 3 (0-based) → downstream shifts by 2.
        //   ref:    A C G T A C G T   (idx 0..7)
        //   sample: A C G T T T A C G T   (idx 0..9)
        let pairs = align_positions("ACGTACGT", "ACGTTTACGT");
        let map: std::collections::HashMap<usize, usize> = pairs.into_iter().collect();
        // Before the insertion: identity. (The lone T at ref idx 3 is ambiguous within the
        // inserted T-run, so it's not asserted.)
        assert_eq!(map[&0], 0);
        assert_eq!(map[&2], 2);
        // After the 2-base insertion: downstream coordinates shift by +2.
        assert_eq!(map[&4], 6); // 'A', uniquely placed right after the T-run
        assert_eq!(map[&5], 7);
        assert_eq!(map[&7], 9);
    }

    #[test]
    fn bundled_rcrs_is_the_full_mitogenome() {
        assert_eq!(rcrs().len(), 16_569);
        assert!(rcrs().starts_with("GATCACAGG")); // canonical rCRS 5' start
    }
}
