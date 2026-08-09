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
//! ## Sizing itself
//!
//! [`BatchSize::for_this_machine`] reads the machine's physical memory and picks from the table
//! above. This is the path the app should use: the target user clicks "Realign" and gets a job
//! sized to their hardware, rather than being asked for a number in bases that nothing in their
//! experience equips them to choose. A wrong answer here is not a preference, it is either an
//! out-of-memory failure or a needlessly split index.
//!
//! It sizes from **total** memory, not currently-available memory, which is deliberate. The `.mmi`
//! is cached and reused for every later job against that build, so a machine that happens to be
//! busy at the moment of the first click would otherwise bake a more-split index — and its
//! permanent MAPQ cost — into the cache. Available memory is still reported by
//! [`detect_memory`], because deciding whether to start *right now* is a different question from
//! how to build the artifact, and belongs to the preflight.

/// Bases per index part.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BatchSize(u64);

/// One gigabase.
const GBASE: u64 = 1_000_000_000;

/// Bytes in a gibibyte.
const GIB: u64 = 1024 * 1024 * 1024;

/// What the machine has to work with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineMemory {
    /// Physical memory installed, in bytes.
    pub total: u64,
    /// Memory the OS reports as available right now, in bytes. Fluctuates; use it to decide
    /// whether to start a job, not to size a cached artifact.
    pub available: u64,
}

impl MachineMemory {
    pub fn total_gib(self) -> u64 {
        self.total / GIB
    }

    pub fn available_gib(self) -> u64 {
        self.available / GIB
    }
}

/// Read the machine's memory.
///
/// `None` if the platform will not say — sysinfo supports every desktop target Navigator ships to,
/// but reporting nothing is preferable to reporting a fabricated number that would then silently
/// size a multi-hour job.
pub fn detect_memory() -> Option<MachineMemory> {
    let mut system = sysinfo::System::new();
    system.refresh_memory();
    let total = system.total_memory();
    if total == 0 {
        return None;
    }
    Some(MachineMemory {
        total,
        available: system.available_memory(),
    })
}

/// The default: 1 Gbase, measured at 11.7 GiB to build and 5.4 GiB to map against CHM13. Chosen to
/// fit a 16 GB machine with room for the OS and the rest of the app.
const DEFAULT_BASES: u64 = GBASE;

/// "Do not split" — 8 Gbase, which is also minimap2's own default `batch_size`. Any human
/// reference fits in one part at this size, so it expresses the intent without a magic sentinel,
/// and it still renders as a real number wherever the choice is reported.
const UNSPLIT: u64 = 8 * GBASE;

/// `NAVIGATOR_ALIGN_BATCH_MBASE`, in megabases — the escape hatch for unusual hardware, and how
/// tests pin the value without depending on the machine they run on.
fn env_override() -> Option<u64> {
    std::env::var("NAVIGATOR_ALIGN_BATCH_MBASE")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .map(|mbase| mbase.saturating_mul(1_000_000).max(1_000_000))
}

/// The conservative default: honours the override, otherwise assumes a 16 GB desktop. Prefer
/// [`BatchSize::for_this_machine`], which actually looks.
impl Default for BatchSize {
    fn default() -> Self {
        Self(env_override().unwrap_or(DEFAULT_BASES).max(1_000_000))
    }
}

impl BatchSize {
    /// An explicit batch size in bases.
    pub fn new(bases: u64) -> Self {
        Self(bases.max(1_000_000))
    }

    /// The batch size for the machine this is running on — the button-click path.
    ///
    /// Precedence, highest first: the `NAVIGATOR_ALIGN_BATCH_MBASE` override, then detected
    /// physical memory, then the 16 GB-desktop default. The override comes first so a user on
    /// unusual hardware, or a test, can pin the value without having to defeat the detector.
    pub fn for_this_machine() -> Self {
        if let Some(bases) = env_override() {
            return Self(bases);
        }
        match detect_memory() {
            Some(memory) => Self::for_ram_gib(memory.total_gib()),
            None => Self(DEFAULT_BASES),
        }
    }

    /// Why [`BatchSize::for_this_machine`] chose what it did, for a log line or a UI tooltip.
    ///
    /// A realignment is a multi-hour job whose memory profile the user cannot see; when it is
    /// sized automatically, the sizing has to be inspectable rather than a mystery.
    pub fn explain() -> String {
        if let Some(bases) = env_override() {
            return format!(
                "index batch {} (from NAVIGATOR_ALIGN_BATCH_MBASE)",
                Self(bases).describe()
            );
        }
        match detect_memory() {
            Some(memory) => format!(
                "index batch {} (detected {} GiB RAM, {} GiB free)",
                Self::for_ram_gib(memory.total_gib()).describe(),
                memory.total_gib(),
                memory.available_gib(),
            ),
            None => format!(
                "index batch {} (could not detect RAM; using the default)",
                Self(DEFAULT_BASES).describe()
            ),
        }
    }

    /// The batch size as a person would say it.
    pub fn describe(self) -> String {
        if self.0 >= UNSPLIT {
            return "unsplit (one index part)".to_string();
        }
        if self.0 >= GBASE {
            let tenths = self.0 / 100_000_000;
            return format!("{}.{} Gbase", tenths / 10, tenths % 10);
        }
        format!("{} Mbase", self.0 / 1_000_000)
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
            // Above 32 GiB a single part is affordable, and a whole index costs no MAPQ fidelity
            // at all — the one thing splitting gives up. `UNSPLIT` rather than a saturating
            // sentinel: this number reaches a log line and a UI tooltip, and "9223372036854
            // Mbase" is not something to show a user who was promised a button.
            _ => UNSPLIT,
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

    /// Detection has to work on whatever machine this runs on — that is the entire point of taking
    /// the dependency. The assertions are about plausibility rather than a specific number, since
    /// the test cannot know the host.
    #[test]
    fn the_machine_reports_its_own_memory() {
        let memory = detect_memory().expect("every desktop target sysinfo supports reports memory");
        assert!(memory.total > 0);
        assert!(
            memory.total_gib() >= 1,
            "a machine running this test suite has at least 1 GiB"
        );
        assert!(
            memory.available <= memory.total,
            "available {} exceeded total {}",
            memory.available,
            memory.total
        );
    }

    /// The button-click path must always yield a usable batch, whatever the host, and must land on
    /// a value the sizing table actually produces rather than something improvised.
    #[test]
    fn sizing_for_this_machine_lands_on_a_table_value() {
        let chosen = BatchSize::for_this_machine();
        assert!(chosen.bases() >= 1_000_000);

        let from_table = detect_memory()
            .map(|m| BatchSize::for_ram_gib(m.total_gib()))
            .unwrap_or_default();
        assert_eq!(chosen, from_table, "detection and the table must agree");
    }

    /// The override is the escape hatch for hardware the table does not suit, so it has to beat
    /// detection rather than merely fill in for it.
    ///
    /// Serialized with the other env-reading test: `set_var` is process-global, and Rust runs tests
    /// in threads by default.
    #[test]
    fn the_env_override_beats_detection() {
        let _guard = ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: guarded by ENV_LOCK so no other test reads the environment concurrently.
        unsafe { std::env::set_var("NAVIGATOR_ALIGN_BATCH_MBASE", "250") };
        let chosen = BatchSize::for_this_machine();
        unsafe { std::env::remove_var("NAVIGATOR_ALIGN_BATCH_MBASE") };

        assert_eq!(chosen.bases(), 250_000_000);
    }

    #[test]
    fn the_explanation_names_the_reason_for_the_choice() {
        let _guard = ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: guarded by ENV_LOCK so no other test reads the environment concurrently.
        unsafe { std::env::set_var("NAVIGATOR_ALIGN_BATCH_MBASE", "400") };
        let explained = BatchSize::explain();
        unsafe { std::env::remove_var("NAVIGATOR_ALIGN_BATCH_MBASE") };
        assert!(explained.contains("400 Mbase"), "{explained}");
        assert!(explained.contains("NAVIGATOR_ALIGN_BATCH_MBASE"), "{explained}");

        // Without the override it reports what it detected, so a support log says why.
        let detected = BatchSize::explain();
        assert!(detected.contains("RAM") || detected.contains("default"), "{detected}");
    }

    /// This string reaches a log line and a UI tooltip, so no branch of the sizing table may render
    /// as a raw sentinel — the large-RAM case used to come out as "9223372036854 Mbase".
    #[test]
    fn every_table_choice_describes_itself_readably() {
        for gib in [4u64, 8, 16, 24, 32, 64, 128, 512] {
            let described = BatchSize::for_ram_gib(gib).describe();
            assert!(
                described.len() <= 32 && !described.contains("922337"),
                "{gib} GiB rendered as {described:?}"
            );
        }
        assert_eq!(BatchSize::for_ram_gib(8).describe(), "400 Mbase");
        assert_eq!(BatchSize::for_ram_gib(16).describe(), "1.0 Gbase");
        assert_eq!(BatchSize::for_ram_gib(128).describe(), "unsplit (one index part)");
    }

    /// `set_var` mutates process-global state, so the tests that touch it must not overlap.
    static ENV_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
}
