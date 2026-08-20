//! Light local realignment around candidate indels (plan §4b mitigation).
//!
//! An indel in a homopolymer or a repeat is ambiguous. So BWA puts the same insertion at a
//! different place in different reads. That spreads bases onto the positions beside it.
//!
//! Here is a real case. On HG002 chrM, a +1C in the C-run at 16295 to 16301 makes about 47 reads
//! put a false C onto the reference T at 16302. That reads as a T>C SNP, and it is not one.
//!
//! GATK avoids this with a local reassembly. This module instead puts the bases of each read, over
//! an active window, back onto the reference, with one consistent gap model. The homopolymer bases
//! then land in one place, and the false substitution goes away.
//!
//! The aligner does a **fit of the read into the window**. It consumes the whole part of the read
//! that it looks at, and the end gaps on the reference window are free. So a read that starts or
//! ends inside the window still aligns. This module is pure, and unit tests cover it.
//! [`crate::caller`] drives it.

/// One aligned column, between a part of a read, which is the query, and a reference window,
/// which is the target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Query base aligned to a reference base (match or mismatch).
    Aligned,
    /// Query base with no reference base (insertion).
    Insertion,
    /// Reference base with no query base (deletion).
    Deletion,
}

// The score uses a linear gap. An indel of one base in a homopolymer is the case that dominates
// here, so an affine gap separates nothing more. A clear mismatch penalty is what drives the
// choice of where the insertion goes.
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

    // The alignment must consume the whole query. So a query base at the start, against an empty
    // target, costs a gap.
    for i in 1..=n {
        score[i][0] = GAP * i as i32;
        tb[i][0] = Move::Up;
    }
    // Row 0. A target gap at the start is free, so the query can start anywhere in the target.
    // The row stays at 0.

    for i in 1..=n {
        for j in 1..=m {
            let s = if query[i - 1].eq_ignore_ascii_case(&target[j - 1]) {
                MATCH
            } else {
                MISMATCH
            };
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

    // A target gap at the end is free, so take the best score across the last query row.
    let mut end_j = 0;
    let mut best = i32::MIN;
    for (j, &s) in score[n].iter().enumerate() {
        if s >= best {
            best = s;
            end_j = j;
        }
    }

    // Trace back to row 0, because a target gap at the start is free.
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
        // The query aligns to target[2..6]. The target gaps at the start and at the end are
        // free.
        assert_eq!(aligned_string(b"CGTA", b"AACGTACG"), (2, "MMMM".into()));
    }

    #[test]
    fn homopolymer_insertion_is_an_insertion_not_a_substitution() {
        // read has an extra C in a C-run then T; ref is CCC...T. The extra C must be an
        // insertion so the read's T aligns to the ref T (not a C smeared onto T).
        let (_start, ops) = fitting_align(b"CCCCCCCCT", b"CCCCCCCT");
        assert_eq!(ops.iter().filter(|o| **o == Op::Insertion).count(), 1);
        // The last column aligns the T at the end.
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
        // The read ends inside the C-run, and it holds no T. After the realignment it must put no
        // base at all on the reference T position. That is the fix for the false SNP.
        let ref_window = b"CCCCCCCT"; // ref C-run + T
        let read = b"CCCCCCCC"; // 8 C's, no T (read ended in the homopolymer)
        let (start, ops) = fitting_align(read, ref_window);
        let quals = vec![40u8; read.len()];
        let proj = project(read, &quals, 100, start, &ops);
        // The T position is window_start+7; the read must not cover it (it has no T).
        assert!(proj.iter().all(|(pos, _, _)| *pos != 100 + 7));
    }
}
