//! Index batch size — the memory control for the whole module.
//!
//! minimap2 splits a reference into index *parts* of at most `batch_size` bases (the CLI's `-I`).
//! One part is resident at a time, so this single number decides peak RAM. Measured on CHM13v2
//! (3.1 Gbase) with the `sr` preset:
//!
//! | batch size | index build peak | mapping peak |
//! |-----------:|-----------------:|-------------:|
//! | whole genome (1 part) | 19.2 GiB | 10.25 GiB |
//! | 1 Gbase (4 parts) | 11.7 GiB | 5.4 GiB |
//! | 400 Mbase | 8.7 GiB | — |
//! | 200 Mbase | 7.5 GiB | — |
//!
//! Wall time was flat across all of them, so bounding memory here is close to free. Building one
//! monolithic index is the failure mode to avoid — it is what made an early estimate conclude the
//! module needed ~19 GB and could not run on a normal desktop.
//!
//! **Bigger is better, within budget.** A split index costs a little MAPQ fidelity: a read's
//! second-best hit can fall in another part and go uncounted, so MAPQ comes out slightly *too
//! high* at multi-mapping loci (measured at 7 of 5,045 records against a ~5-part split, with every
//! locus identical). So this picks the largest batch that fits, never the smallest that works.
//!
//! ## Why there is no RAM auto-detection here
//!
//! Reading physical memory needs `/proc/meminfo` on Linux, `sysctl` on macOS, and a Win32 call on
//! Windows. The workspace has no dependency that provides it and forbids shelling out, so a
//! detector here would mean either a new platform-specific dependency or a subprocess — both worse
//! than the alternative. Instead the budget is explicit: callers that know the machine's RAM pass
//! it to [`BatchSize::for_ram_gib`], and everyone else gets a default that fits a 16 GB desktop.

/// Bases per index part.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BatchSize(u64);

/// One gigabase.
const GBASE: u64 = 1_000_000_000;

/// The default: 1 Gbase, measured at 11.7 GiB to build and 5.4 GiB to map against CHM13. Chosen to
/// fit a 16 GB machine with room for the OS and the rest of the app.
const DEFAULT_BASES: u64 = GBASE;

impl Default for BatchSize {
    fn default() -> Self {
        // `NAVIGATOR_ALIGN_BATCH_MBASE` overrides, in megabases, for a machine this crate cannot
        // measure and a caller that hasn't been told either.
        let bases = std::env::var("NAVIGATOR_ALIGN_BATCH_MBASE")
            .ok()
            .and_then(|s| s.parse::<u64>().ok())
            .map(|mbase| mbase.saturating_mul(1_000_000))
            .unwrap_or(DEFAULT_BASES);
        Self(bases.max(1_000_000))
    }
}

impl BatchSize {
    /// An explicit batch size in bases.
    pub fn new(bases: u64) -> Self {
        Self(bases.max(1_000_000))
    }

    pub fn bases(self) -> u64 {
        self.0
    }

    /// The largest batch that fits a machine with `ram_gib` of physical memory, from the measured
    /// table in the module docs.
    ///
    /// The thresholds leave headroom deliberately: the numbers in that table are the mapper's peak
    /// alone, and a realignment job is also holding a sort buffer, the revert's scratch, and a
    /// desktop application. Below 8 GiB nothing here is comfortable, so the smallest step is
    /// offered rather than refusing outright — the preflight decides whether to proceed, not this.
    pub fn for_ram_gib(ram_gib: u64) -> Self {
        let bases = match ram_gib {
            0..=7 => 200_000_000,
            8..=15 => 400_000_000,
            16..=31 => GBASE,
            // Above 32 GiB a single part is affordable and costs no MAPQ fidelity at all, which is
            // the one thing splitting gives up.
            _ => u64::MAX / 2,
        };
        Self(bases)
    }

    /// Whether a reference of `total_bases` will be split into more than one part — i.e. whether
    /// the cross-part merge and its MAPQ caveat come into play at all.
    pub fn splits(self, total_bases: u64) -> bool {
        total_bases > self.0
    }

    /// An **upper bound** on how many parts `total_bases` will produce, for sizing a progress bar.
    ///
    /// Deliberately not exact. A part accumulates whole sequences until the running total
    /// *exceeds* the batch, so parts overshoot by up to one sequence and the real count comes in
    /// at or below this — measured, a 3 Mbase reference at a 1 Mbase batch yields 2 parts where
    /// this returns 3. A progress bar that finishes early is fine; one that runs past its own
    /// maximum is not.
    pub fn part_estimate(self, total_bases: u64) -> usize {
        if self.0 == 0 {
            return 1;
        }
        total_bases.div_ceil(self.0).max(1) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sizing table's whole purpose: a 16 GB desktop must land on a batch that fits it, and a
    /// small machine must land on a smaller one.
    #[test]
    fn ram_maps_to_the_largest_batch_that_fits() {
        assert_eq!(BatchSize::for_ram_gib(8).bases(), 400_000_000);
        assert_eq!(BatchSize::for_ram_gib(16).bases(), GBASE);
        assert!(
            BatchSize::for_ram_gib(64).bases() > BatchSize::for_ram_gib(16).bases(),
            "a large machine takes a single part, which costs no MAPQ fidelity"
        );
        assert!(BatchSize::for_ram_gib(4).bases() < BatchSize::for_ram_gib(8).bases());
    }

    /// Monotonicity: more memory must never select a smaller batch, or the table has a hole.
    #[test]
    fn more_ram_never_means_a_smaller_batch() {
        let mut previous = 0;
        for gib in [4u64, 8, 12, 16, 24, 32, 64, 128] {
            let bases = BatchSize::for_ram_gib(gib).bases();
            assert!(bases >= previous, "{gib} GiB regressed to {bases}");
            previous = bases;
        }
    }

    /// CHM13 is 3.1 Gbase; the default must split it (that is the point) into a handful of parts.
    #[test]
    fn the_default_splits_a_human_genome_into_a_few_parts() {
        let chm13 = 3_100_000_000u64;
        let default = BatchSize::default();
        assert!(default.splits(chm13));
        assert_eq!(default.part_estimate(chm13), 4);
    }

    #[test]
    fn a_reference_smaller_than_the_batch_is_one_part() {
        let b = BatchSize::new(GBASE);
        assert!(!b.splits(50_000_000));
        assert_eq!(b.part_estimate(50_000_000), 1);
    }

    /// A zero or absurdly small batch would produce an unbounded part count and thrash; the floor
    /// keeps a bad caller from turning the mapper into a no-op.
    #[test]
    fn the_batch_size_has_a_floor() {
        assert!(BatchSize::new(0).bases() >= 1_000_000);
    }
}
