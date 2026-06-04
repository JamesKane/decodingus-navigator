//! Light local realignment around candidate indels (plan §4b mitigation).
//!
//! Ambiguous indels in homopolymers/repeats make BWA place the same insertion
//! differently across reads, smearing bases onto neighbouring positions — e.g. on
//! HG002 chrM a +1C in the 16295–16301 C-run makes ~47 reads put a spurious C on the
//! reference T at 16302, a false T>C SNP. GATK avoids this by local reassembly; here we
//! re-fit each read's bases over an active window back onto the reference with a
//! consistent gap model, so the homopolymer bases land in one place and the spurious
//! substitution disappears.
//!
//! The aligner is a **fitting alignment**: the read substring is fully consumed, with
//! free end gaps on the reference window (so reads starting/ending inside the window
//! still align). This module is pure and unit-tested; [`crate::caller`] drives it.

/// One aligned column between a read substring (query) and a reference window (target).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Query base aligned to a reference base (match or mismatch).
    Aligned,
    /// Query base with no reference base (insertion).
    Insertion,
    /// Reference base with no query base (deletion).
    Deletion,
}

// Linear-gap scoring. Single-base homopolymer indels dominate, so affine gaps add no
// resolving power here; a clear mismatch penalty drives the insertion choice.
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;
const GAP: i32 = -3;

#[derive(Clone, Copy)]
enum Move {
    Diag,
    Up,   // consume query only (insertion)
    Left, // consume target only (deletion)
    Stop,
}

/// Fit `query` fully into `target`, free end gaps on `target`. Returns the target
/// offset where the alignment starts and the column ops (query left-to-right).
pub fn fitting_align(query: &[u8], target: &[u8]) -> (usize, Vec<Op>) {
    let (n, m) = (query.len(), target.len());
    if n == 0 {
        return (0, Vec::new());
    }
    // score[i][j], i over query (0..=n), j over target (0..=m).
    let mut score = vec![vec![0i32; m + 1]; n + 1];
    let mut tb = vec![vec![Move::Stop; m + 1]; n + 1];

    // Query must be fully consumed: leading query bases against empty target cost gaps.
    for i in 1..=n {
        score[i][0] = GAP * i as i32;
        tb[i][0] = Move::Up;
    }
    // Row 0: free leading target gap (query can start anywhere in target) -> stays 0.

    for i in 1..=n {
        for j in 1..=m {
            let s = if query[i - 1].eq_ignore_ascii_case(&target[j - 1]) { MATCH } else { MISMATCH };
            let diag = score[i - 1][j - 1] + s;
            let up = score[i - 1][j] + GAP; // insertion (query base, no target)
            let left = score[i][j - 1] + GAP; // deletion (target base, no query)
            let mut best = diag;
            let mut mv = Move::Diag;
            if up > best {
                best = up;
                mv = Move::Up;
            }
            if left > best {
                best = left;
                mv = Move::Left;
            }
            score[i][j] = best;
            tb[i][j] = mv;
        }
    }

    // Free trailing target gap: best score across the last query row.
    let mut end_j = 0;
    let mut best = i32::MIN;
    for (j, &s) in score[n].iter().enumerate() {
        if s >= best {
            best = s;
            end_j = j;
        }
    }

    // Traceback to row 0 (free leading target gap).
    let mut ops = Vec::new();
    let (mut i, mut j) = (n, end_j);
    while i > 0 {
        match tb[i][j] {
            Move::Diag => {
                ops.push(Op::Aligned);
                i -= 1;
                j -= 1;
            }
            Move::Up => {
                ops.push(Op::Insertion);
                i -= 1;
            }
            Move::Left => {
                ops.push(Op::Deletion);
                j -= 1;
            }
            Move::Stop => break,
        }
    }
    ops.reverse();
    (j, ops) // j is now the target start offset
}

/// Project a realigned read onto reference positions. `query`/`quals` are the read's
/// bases/qualities over the window; `target_start` is the 0-based offset within the
/// window where the alignment begins; `window_start` is the window's reference index.
/// Returns `(ref_index, base, qual)` for each aligned (diagonal) column.
pub fn project(
    query: &[u8],
    quals: &[u8],
    window_start: usize,
    target_start: usize,
    ops: &[Op],
) -> Vec<(usize, u8, u8)> {
    let mut out = Vec::new();
    let mut qi = 0usize;
    let mut tj = target_start;
    for op in ops {
        match op {
            Op::Aligned => {
                out.push((window_start + tj, query[qi], quals.get(qi).copied().unwrap_or(0)));
                qi += 1;
                tj += 1;
            }
            Op::Insertion => qi += 1,
            Op::Deletion => tj += 1,
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aligned_string(query: &[u8], target: &[u8]) -> (usize, String) {
        let (start, ops) = fitting_align(query, target);
        let s: String = ops
            .iter()
            .map(|o| match o {
                Op::Aligned => 'M',
                Op::Insertion => 'I',
                Op::Deletion => 'D',
            })
            .collect();
        (start, s)
    }

    #[test]
    fn exact_match_is_all_aligned() {
        assert_eq!(aligned_string(b"ACGT", b"ACGT"), (0, "MMMM".into()));
    }

    #[test]
    fn query_fits_into_a_substring_of_target() {
        // query aligns to target[2..6]; free leading/trailing target gaps.
        assert_eq!(aligned_string(b"CGTA", b"AACGTACG"), (2, "MMMM".into()));
    }

    #[test]
    fn homopolymer_insertion_is_an_insertion_not_a_substitution() {
        // read has an extra C in a C-run then T; ref is CCC...T. The extra C must be an
        // insertion so the read's T aligns to the ref T (not a C smeared onto T).
        let (_start, ops) = fitting_align(b"CCCCCCCCT", b"CCCCCCCT");
        assert_eq!(ops.iter().filter(|o| **o == Op::Insertion).count(), 1);
        // last column aligns the trailing T.
        assert_eq!(*ops.last().unwrap(), Op::Aligned);

        // Projection puts the T on the last reference position, never a C.
        let quals = vec![40u8; 9];
        let proj = project(b"CCCCCCCCT", &quals, 100, _start, &ops);
        let last = proj.last().unwrap();
        assert_eq!(last.1, b'T');
        assert_eq!(last.0, 100 + 7); // ref window position of the T (8th target base)
    }

    #[test]
    fn read_ending_in_homopolymer_does_not_reach_the_trailing_base() {
        // read ends in the C-run (no T); after realignment it should not place any base
        // on the reference T position — the spurious-SNP fix.
        let ref_window = b"CCCCCCCT"; // ref C-run + T
        let read = b"CCCCCCCC"; // 8 C's, no T (read ended in the homopolymer)
        let (start, ops) = fitting_align(read, ref_window);
        let quals = vec![40u8; read.len()];
        let proj = project(read, &quals, 100, start, &ops);
        // The T position is window_start+7; the read must not cover it (it has no T).
        assert!(proj.iter().all(|(pos, _, _)| *pos != 100 + 7));
    }
}
