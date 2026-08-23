//! Index batch size: the memory control for the whole module.
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
//! Wall time was flat across all of them, so a limit on memory here is almost free. One
//! monolithic index is the failure mode to avoid. It is what made an early estimate say the module
//! needed ~19 GB, and could not run on a normal desktop.
//!
//! **Bigger is better, inside the budget.** A split index costs a little MAPQ fidelity. A read's
//! second-best hit can fall in another part, where the count misses it. So MAPQ comes out a little
//! *too high* at a locus with more than one hit. The measurement was 7 of 5,045 records against a
//! ~5-part split, and every locus was identical. So this code chooses the largest batch that fits,
//! and never the smallest one that works.
//!
//! ## How the code chooses the size
//!
//! [`BatchSize::for_this_machine`] reads the machine's physical memory and chooses from the table
//! above. This is the path the app must use. The target user clicks "Realign" and gets a job that
//! fits their hardware. Nobody asks them for a number in bases, because nothing in their
//! experience prepares them to choose one. A wrong answer here is not a preference. It is an
//! out-of-memory failure, or an index with more parts than it needs.
//!
//! It reads **total** memory, and not the memory that is free at that moment. That is deliberate.
//! The cache keeps the `.mmi`, and every later job against that build uses it again. A machine
//! that is busy at the moment of the first click would otherwise put an index with more parts into
//! the cache. That index carries a permanent MAPQ cost.
//!
//! [`detect_memory`] still reports the free memory. The question "can the job start *right now*"
//! belongs to the preflight. That is a different question from "how do we build the artifact".

/// Bases in each index part.
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
/// `None` if the platform will not say. sysinfo supports every desktop target that Navigator
/// ships to. But it is better to report nothing than to report an invented number. Such a number
/// would then set the size of a multi-hour job, with no warning.
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

/// "Do not split": 8 Gbase, which is also minimap2's own default `batch_size`. Any human
/// reference fits in one part at this size. So the value shows the intent, and it is not a special
/// sentinel. It also renders as a real number everywhere the code reports the choice.
const UNSPLIT: u64 = 8 * GBASE;

/// `NAVIGATOR_ALIGN_BATCH_MBASE`, in megabases. This is the override for unusual hardware. It is
/// also how tests pin the value, so that they do not depend on the machine they run on.
fn env_override() -> Option<u64> {
    std::env::var("NAVIGATOR_ALIGN_BATCH_MBASE")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .map(|mbase| mbase.saturating_mul(1_000_000).max(1_000_000))
}

/// The conservative default: it obeys the override, and if there is none it assumes a 16 GB
/// desktop. Prefer [`BatchSize::for_this_machine`], which reads the machine.
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

    /// The batch size for the machine this code runs on. This is the button-click path.
    ///
    /// Precedence, highest first: the `NAVIGATOR_ALIGN_BATCH_MBASE` override, then the physical
    /// memory that the code found, then the 16 GB-desktop default. The override comes first, so
    /// that a user on unusual hardware, or a test, can pin the value and not defeat the detector.
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
    /// A realignment is a multi-hour job, and the user can not see its memory profile. The code
    /// chooses the size without help, so the user must be able to see how it chose.
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
    /// The thresholds leave headroom on purpose. The numbers in that table are the peak of the
    /// mapper alone. A realignment job also holds a sort buffer, the scratch of the revert, and a
    /// desktop application. Below 8 GiB no value here is comfortable, so this code gives the
    /// smallest step and does not refuse. The preflight decides whether to go on, and this does
    /// not.
    pub fn for_ram_gib(ram_gib: u64) -> Self {
        let bases = match ram_gib {
            0..=7 => 200_000_000,
            8..=15 => 400_000_000,
            16..=31 => GBASE,
            // Above 32 GiB the machine can hold one part, and a whole index costs no MAPQ
            // fidelity at all. MAPQ fidelity is the one thing a split index gives up. Use
            // `UNSPLIT`, and not a sentinel at the top of the number range. This number reaches a
            // log line and a UI tooltip. "9223372036854 Mbase" is not a thing to show a user who
            // clicked one button.
            _ => UNSPLIT,
        };
        Self(bases)
    }

    /// True when a reference of `total_bases` needs more than one part. That is also when the
    /// cross-part merge and its MAPQ caveat apply at all.
    pub fn splits(self, total_bases: u64) -> bool {
        total_bases > self.0
    }

    /// An **upper bound** on how many parts `total_bases` makes. Use it to set the length of a
    /// progress bar.
    ///
    /// It is not exact, on purpose. A part collects whole sequences until the total *goes above*
    /// the batch. So a part overshoots by as much as one sequence, and the real count is this
    /// number or less. The measurement: a 3 Mbase reference at a 1 Mbase batch makes 2 parts,
    /// where this returns 3. A progress bar that ends early is acceptable. One that goes past its
    /// own maximum is not.
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

    /// The whole purpose of the table: a 16 GB desktop must land on a batch that fits it. A small
    /// machine must land on a smaller batch.
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
        // This reads the environment, so it takes ENV_LOCK too. The guard has no value unless
        // the readers hold it and the writers hold it.
        let _guard = ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());
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

    /// A batch of zero, or a very small batch, would make a part count with no limit, and the
    /// machine would thrash. The floor stops a bad caller who would make the mapper do nothing.
    #[test]
    fn the_batch_size_has_a_floor() {
        assert!(BatchSize::new(0).bases() >= 1_000_000);
    }

    /// Detection must work on any machine this runs on. That is the whole reason for the
    /// dependency. The assertions are about a plausible range, and not about a specific number,
    /// because the test can not know the host.
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

    /// The button-click path must always give a batch that works, on any host. It must also land
    /// on a value that the table makes, and not on an improvised one.
    #[test]
    fn sizing_for_this_machine_lands_on_a_table_value() {
        let _guard = ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let chosen = BatchSize::for_this_machine();
        assert!(chosen.bases() >= 1_000_000);

        let from_table = detect_memory()
            .map(|m| BatchSize::for_ram_gib(m.total_gib()))
            .unwrap_or_default();
        assert_eq!(chosen, from_table, "detection and the table must agree");
    }

    /// The override is for hardware that the table does not suit, so it must win over detection.
    /// It must not only fill in when detection gives nothing.
    ///
    /// This test runs in sequence with the other test that reads the environment: `set_var` is
    /// process-global, and Rust runs tests in threads by default.
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

    /// This string reaches a log line and a UI tooltip. So no branch of the table may render as a
    /// raw sentinel. The large-RAM case used to come out as "9223372036854 Mbase".
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
