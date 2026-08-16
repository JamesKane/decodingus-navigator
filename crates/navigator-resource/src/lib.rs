//! Watching what a long stage is doing to the machine.
//!
//! A realignment runs for hours on a machine its user is still trying to use, and the failure that
//! motivated this module was not the one anybody was watching for. On 2026-08-13 a WGS-scale sort
//! was killed six hours in, and the obvious suspect was memory. It was not: the compressor was
//! idle, no jetsam report was filed, and anonymous memory sat at 32 GB of 128. What the operating
//! system *did* complain about was writes — it filed a resource notice against the process for
//! dirtying 549 GB of file-backed memory at nearly 1.4x the system's sustained write-back limit.
//! WindowServer's main thread then missed a 40-second watchdog check-in, the watchdog killed it,
//! and every process in the login session went down with it, this job included.
//!
//! So the thing worth watching is not one number. Memory is cheap to sample and would have
//! exonerated itself immediately; the write rate is the number that was actually extreme, and
//! nothing was recording it. This samples both, on a cadence, and hands the caller a
//! [`ResourceSample`] to log.
//!
//! ## It reports; it does not intervene
//!
//! Nothing here aborts a job. A stage that is writing hard is doing its job — the sort *is* a
//! hundreds-of-GB write — and a watchdog that killed a six-hour run for going fast would be worse
//! than the problem. Bounding the damage belongs where the writes happen, on a byte cadence
//! ([`PacedFile`]); this exists so that the next time something goes wrong there is a record of
//! what the machine looked like, instead of an inference from a crash report.
//!
//! ## Why a crate of its own
//!
//! Because a byte counter only means anything if there is exactly one of it, and the pipeline's
//! writers are split across crates that do not depend on each other: `navigator-align` maps,
//! `navigator-analysis` reverts, sorts, marks, and compresses. This first shipped as a module in
//! the latter, which meant the mapping stage — the single longest in the job, and the one that
//! writes the ~60 GB `mapped.bam` — was neither paced nor counted. The run log read `0 MB/s`
//! straight through it, which is not a quiet stage but an unmeasured one.
//!
//! ## Portability
//!
//! Every probe here is `sysinfo`, which binds the platform APIs through pure-Rust crates on all
//! three desktop targets — the same reason `navigator-align` picked it for RAM detection. There is
//! deliberately no `fcntl`/`ioctl`/`/proc` in this crate: the guard has to hold on Windows, and a
//! guard that only arms itself on macOS would have been no guard at all for most users.

use std::fs::File;
use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

/// Bytes handed to disk by the pipeline's paced writers since the process started.
///
/// Self-accounted rather than read back from the OS: every platform exposes per-process I/O
/// counters differently (and Windows' are not in `sysinfo`'s default surface), whereas the writers
/// already know exactly how much they wrote. It costs one relaxed add per buffer.
///
/// One counter for the whole process, which is the reason this lives in a crate of its own: a
/// realignment's writes come from two crates that cannot see each other, and two counters would
/// have been two half-answers.
static BYTES_WRITTEN: AtomicU64 = AtomicU64::new(0);

/// Account `n` bytes written. Called from the write path; must stay this cheap.
pub fn record_bytes_written(n: u64) {
    BYTES_WRITTEN.fetch_add(n, Ordering::Relaxed);
}

/// Total bytes the post-processing writers have written this process.
pub fn bytes_written() -> u64 {
    BYTES_WRITTEN.load(Ordering::Relaxed)
}

/// Bytes written between forced flushes to disk. See [`PacedFile`].
const DEFAULT_SYNC_MB: u64 = 256;

/// How much may be left dirty before a paced stream pushes it to disk.
///
/// `NAVIGATOR_IO_SYNC_MB=0` turns the pacing off and restores the unbounded behaviour.
fn sync_interval() -> u64 {
    std::env::var("NAVIGATOR_IO_SYNC_MB")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(DEFAULT_SYNC_MB)
        * 1024
        * 1024
}

/// A file that will not let an unbounded amount of its output sit dirty in the page cache.
///
/// Without this, a stage writes as fast as it can into the page cache and leaves write-back to the
/// operating system, which sounds like the right division of labour and is, until the volume gets
/// far enough out of scale. One WGS realignment dirtied 549 GB of file-backed memory — enough that
/// macOS filed a disk-writes resource notice against the process for exceeding its sustained
/// write-back limit by 1.4x, and enough that WindowServer's main thread missed a 40-second watchdog
/// check-in and was killed, taking the login session and a six-hour job with it.
///
/// Flushing on a byte cadence caps how much can be outstanding at once. The write path pays for its
/// own I/O as it goes instead of handing the machine a debt to settle later, which is also why this
/// is not obviously a slowdown: the same bytes reach the same disk, in steadier instalments rather
/// than in storms.
///
/// It lives here rather than beside any one stage because every stage that writes tens of GB wants
/// it: the mapper's output and its 8.93 GB minimizer index, the revert's spill runs and FASTQ
/// (where the scratch peak actually is), the sort's runs and merged output, and the final CRAM. The
/// byte accounting that [`ResourceWatch`] reports comes from the same place, so a stream that is
/// paced is also a stream that is counted — and a writer nobody wrapped is a writer that shows up
/// in neither.
///
/// `sync_data` rather than `sync_all` — the contents must be durable, the metadata need not be, and
/// on a stream this size that is many thousands of inode updates. It is std's portable spelling:
/// `fdatasync` where there is one, `FlushFileBuffers` on Windows.
pub struct PacedFile {
    file: File,
    since_sync: u64,
    interval: u64,
}

impl PacedFile {
    pub fn new(file: File) -> Self {
        Self {
            file,
            since_sync: 0,
            interval: sync_interval(),
        }
    }

    /// Push everything still outstanding to disk.
    ///
    /// Called at the end of a stream whose end-of-file marker a later run may trust: a marker
    /// sitting in the page cache is a promise the disk has not made.
    pub fn sync(&self) -> io::Result<()> {
        self.file.sync_data()
    }
}

impl Write for PacedFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = self.file.write(buf)?;
        record_bytes_written(written as u64);

        if self.interval > 0 {
            self.since_sync += written as u64;
            if self.since_sync >= self.interval {
                self.file.sync_data()?;
                self.since_sync = 0;
            }
        }

        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}

/// The budget an external-sort stage holds in memory before spilling a run to disk.
///
/// Both of the pipeline's spill-to-disk stages — the revert's collator and the coordinate sort —
/// were sized by a constant, 256 MB and 512 MB, chosen when nobody had run a WGS through them.
/// The measured cost of that: a 30x WGS coordinate sort spilled **688 runs**, all of which the
/// merge then opens at once. It is bounded memory by design and it does work, but on a machine with
/// 128 GB of RAM it is a lot of fan-in bought for no reason, and the constant that produced it was
/// the same on a laptop that genuinely needed it.
///
/// So the number comes from the machine. `var` still wins when it is set — an explicit MB count is
/// the escape hatch for a run that has to be reproduced or squeezed — and the sizing is otherwise:
///
/// - **A quarter of installed RAM.** Total rather than free, because free fluctuates with whatever
///   the user happens to have open, and a stage whose run count depends on the browser is a stage
///   whose behaviour cannot be reproduced from a bug report.
/// - **Never below 512 MB.** That is the sort's existing default, so no machine sorts with less
///   than it does today.
/// - **Never above 8 GB.** Past that the returns are gone — 88 runs against 44 is nothing next to
///   688 against 88 — while the costs are not: the stable sort allocates half the buffer again as
///   scratch, and growing the record vector doubles its allocation while still holding the old one.
/// - **Never more than half of what is free right now.** The stable part of the rule assumes a
///   machine that is otherwise idle. When it is not, spilling an extra run is cheap and swapping
///   the buffer is not.
///
/// Memory the platform will not report reads as zero, and zero means *unknown*, which must not be
/// read as "no memory" — the same rule [`classify`] follows. An unknown machine gets the floor.
pub fn spill_budget(var: &str) -> u64 {
    if let Some(mb) = std::env::var(var).ok().and_then(|s| s.parse::<u64>().ok()) {
        return mb.max(1) * 1024 * 1024;
    }
    let mut system = sysinfo::System::new();
    system.refresh_memory();
    budget(system.total_memory(), system.available_memory())
}

/// What one heap allocation costs beyond the bytes asked for: the allocator's own bookkeeping plus
/// rounding up to a size class.
///
/// It lives beside [`spill_budget`] because the two are one contract. A budget is only as honest as
/// the tally that fills it, and a stage that counts only payload bytes will hold well over its
/// budget in real memory — which is fine against a hand-picked 512 MB constant chosen with a margin
/// nobody wrote down, and not fine against a fraction of the machine.
///
/// Sixteen bytes is the conventional figure for the allocators on the three desktop targets. This
/// is a budget estimate rather than an audit; the point is that a record with four small vectors
/// costs meaningfully more than the sum of their lengths.
pub const ALLOCATION_OVERHEAD: usize = 16;

/// The floor, and the answer for a machine that will not say how much memory it has.
const MIN_SPILL_BUDGET: u64 = 512 << 20;
/// The ceiling. See [`spill_budget`] for why bigger stops paying.
const MAX_SPILL_BUDGET: u64 = 8 << 30;

/// The sizing decision, split from the probe so it is testable on any machine — the same split
/// [`classify`] makes, and for the same reason.
fn budget(total_memory: u64, available_memory: u64) -> u64 {
    if total_memory == 0 {
        return MIN_SPILL_BUDGET;
    }
    let budget = (total_memory / 4).clamp(MIN_SPILL_BUDGET, MAX_SPILL_BUDGET);
    if available_memory == 0 {
        return budget;
    }
    // The floor holds even here: a machine this short of memory would have taken 512 MB under the
    // old constant anyway, so honouring the busy-machine guard past that point would be a
    // regression dressed as caution.
    budget.min((available_memory / 2).max(MIN_SPILL_BUDGET))
}

/// How hard the machine is being leaned on.
///
/// Bands, not a single threshold, because the interesting reading is the trend: a stage that
/// spends an hour at [`Pressure::Elevated`] is a different story from one that touches it once.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pressure {
    Normal,
    Elevated,
    Critical,
}

impl Pressure {
    pub fn label(self) -> &'static str {
        match self {
            Pressure::Normal => "normal",
            Pressure::Elevated => "elevated",
            Pressure::Critical => "critical",
        }
    }
}

/// Free memory below this fraction of total is [`Pressure::Elevated`].
const ELEVATED_FREE_FRACTION: f64 = 0.15;
/// Free memory below this fraction of total is [`Pressure::Critical`].
const CRITICAL_FREE_FRACTION: f64 = 0.07;
/// Swap growth (bytes, since the watch started) that counts as [`Pressure::Elevated`].
///
/// Growth, not absolute use: a desktop that has been up for a week has swap in use for reasons
/// that have nothing to do with this job, and blaming the job for it would cry wolf on every run.
const ELEVATED_SWAP_GROWTH: u64 = 2 << 30;
/// Swap growth that counts as [`Pressure::Critical`].
const CRITICAL_SWAP_GROWTH: u64 = 8 << 30;

/// One reading of the machine.
#[derive(Debug, Clone)]
pub struct ResourceSample {
    /// Time since the watch started.
    pub elapsed: Duration,
    /// Physical memory installed.
    pub total_memory: u64,
    /// Memory the OS says is available for allocation.
    pub available_memory: u64,
    /// Swap in use, and how much of it appeared since the watch started.
    pub used_swap: u64,
    pub swap_growth: u64,
    /// Bytes written by the pipeline's paced writers, in total and since the previous sample.
    pub bytes_written: u64,
    pub write_rate: f64,
    pub pressure: Pressure,
}

impl ResourceSample {
    /// A one-line rendering for a log.
    ///
    /// Deliberately compact: at a 30-second cadence a six-hour stage produces 720 of these, and
    /// they have to stay skimmable next to the stage timings they sit between.
    pub fn summary(&self) -> String {
        format!(
            "[{:>7.0}s] mem {:.1}/{:.1} GB free, swap +{:.1} GB, wrote {:.1} GB @ {:.0} MB/s — {}",
            self.elapsed.as_secs_f64(),
            gib(self.available_memory),
            gib(self.total_memory),
            gib(self.swap_growth),
            gib(self.bytes_written),
            self.write_rate / (1024.0 * 1024.0),
            self.pressure.label(),
        )
    }
}

fn gib(bytes: u64) -> f64 {
    bytes as f64 / (1024.0 * 1024.0 * 1024.0)
}

/// Classify a reading. Split from the probe so the decision is testable without a machine in a
/// particular state.
fn classify(total_memory: u64, available_memory: u64, swap_growth: u64) -> Pressure {
    // A platform that will not report its memory reports zero; that is "unknown", and unknown must
    // not read as "critical" — see `navigator_app::realign_job::has_room` for the same rule.
    let free_fraction = if total_memory == 0 {
        1.0
    } else {
        available_memory as f64 / total_memory as f64
    };

    if free_fraction < CRITICAL_FREE_FRACTION || swap_growth >= CRITICAL_SWAP_GROWTH {
        Pressure::Critical
    } else if free_fraction < ELEVATED_FREE_FRACTION || swap_growth >= ELEVATED_SWAP_GROWTH {
        Pressure::Elevated
    } else {
        Pressure::Normal
    }
}

/// How often to sample, when the caller does not say.
///
/// Thirty seconds is chosen against the thing being watched: the WindowServer watchdog fires at
/// 40 seconds, so a cadence slower than that could step over the entire window in which the
/// machine was in trouble.
pub const DEFAULT_INTERVAL: Duration = Duration::from_secs(30);

/// A background sampler, running for as long as it is held.
///
/// Dropping it stops the thread and joins it, so a stage cannot outlive its own watch or leave a
/// thread behind on the way out of an error.
pub struct ResourceWatch {
    stop: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
}

impl ResourceWatch {
    /// Start sampling every `interval`, passing each reading to `report`.
    ///
    /// `report` runs on the sampler thread and should do no real work — writing a line is the
    /// intended use.
    pub fn start(interval: Duration, mut report: impl FnMut(ResourceSample) + Send + 'static) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let flag = Arc::clone(&stop);

        let handle = thread::spawn(move || {
            let started = Instant::now();
            let mut system = sysinfo::System::new();
            system.refresh_memory();

            let baseline_swap = system.used_swap();
            let mut last_bytes = bytes_written();
            let mut last_at = started;

            // Poll the stop flag far more often than the sample cadence, so dropping the watch
            // returns promptly instead of blocking a job's teardown for up to `interval`.
            const TICK: Duration = Duration::from_millis(250);

            while !flag.load(Ordering::Relaxed) {
                let due = last_at + interval;
                while Instant::now() < due {
                    if flag.load(Ordering::Relaxed) {
                        return;
                    }
                    thread::sleep(TICK);
                }

                system.refresh_memory();
                let now = Instant::now();
                let total = bytes_written();
                let seconds = now.duration_since(last_at).as_secs_f64().max(f64::EPSILON);
                let used_swap = system.used_swap();
                let swap_growth = used_swap.saturating_sub(baseline_swap);
                let total_memory = system.total_memory();
                let available_memory = system.available_memory();

                report(ResourceSample {
                    elapsed: now.duration_since(started),
                    total_memory,
                    available_memory,
                    used_swap,
                    swap_growth,
                    bytes_written: total,
                    write_rate: total.saturating_sub(last_bytes) as f64 / seconds,
                    pressure: classify(total_memory, available_memory, swap_growth),
                });

                last_bytes = total;
                last_at = now;
            }
        });

        Self {
            stop,
            handle: Some(handle),
        }
    }
}

impl Drop for ResourceWatch {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plenty_of_memory_and_no_swap_is_normal() {
        assert_eq!(classify(128 << 30, 100 << 30, 0), Pressure::Normal);
    }

    #[test]
    fn a_nearly_full_machine_is_critical() {
        assert_eq!(classify(128 << 30, 4 << 30, 0), Pressure::Critical);
        assert_eq!(classify(128 << 30, 12 << 30, 0), Pressure::Elevated);
    }

    /// Swap that appeared *during* the job is the signal; swap that was already there is not, and
    /// the sample carries growth rather than absolute use for exactly that reason.
    #[test]
    fn swap_growth_raises_pressure_on_an_otherwise_idle_machine() {
        assert_eq!(classify(128 << 30, 100 << 30, 3 << 30), Pressure::Elevated);
        assert_eq!(classify(128 << 30, 100 << 30, 9 << 30), Pressure::Critical);
    }

    /// A platform that will not report memory must not read as a machine in trouble.
    #[test]
    fn unknown_memory_is_not_pressure() {
        assert_eq!(classify(0, 0, 0), Pressure::Normal);
    }

    /// The case the autosizing exists for: a big machine should stop spilling hundreds of runs.
    #[test]
    fn a_large_machine_gets_the_ceiling() {
        assert_eq!(budget(128 << 30, 100 << 30), MAX_SPILL_BUDGET);
    }

    #[test]
    fn an_ordinary_machine_gets_a_quarter_of_it() {
        assert_eq!(budget(16 << 30, 12 << 30), 4 << 30);
    }

    /// A machine whose memory is already spoken for gets a smaller buffer, because an extra spilled
    /// run costs a file and swapping the buffer costs the run.
    #[test]
    fn a_busy_machine_is_held_to_half_of_what_is_free() {
        assert_eq!(budget(64 << 30, 6 << 30), 3 << 30);
    }

    /// Never below what the sort used before any of this existed.
    #[test]
    fn a_small_or_busy_machine_never_goes_under_the_old_default() {
        assert_eq!(budget(2 << 30, 2 << 30), MIN_SPILL_BUDGET);
        assert_eq!(budget(64 << 30, 100 << 20), MIN_SPILL_BUDGET);
    }

    /// Unknown is not zero. A platform that will not report memory must not be sized as if it had
    /// none — and must not be sized as if it had plenty either.
    #[test]
    fn unknown_memory_gets_the_floor() {
        assert_eq!(budget(0, 0), MIN_SPILL_BUDGET);
        assert_eq!(budget(64 << 30, 0), MAX_SPILL_BUDGET.min(16 << 30));
    }

    /// The escape hatch has to win, or a run cannot be reproduced on a different machine.
    #[test]
    fn an_explicit_override_beats_the_machine() {
        let var = "NAVIGATOR_TEST_SPILL_MB_OVERRIDE";
        std::env::set_var(var, "7");
        assert_eq!(spill_budget(var), 7 * 1024 * 1024);
        std::env::remove_var(var);
    }

    /// A nonsense override must not produce a zero-byte budget, which would spill one run per
    /// record.
    #[test]
    fn a_zero_override_is_floored_at_one_megabyte() {
        let var = "NAVIGATOR_TEST_SPILL_MB_ZERO";
        std::env::set_var(var, "0");
        assert_eq!(spill_budget(var), 1024 * 1024);
        std::env::remove_var(var);
    }

    #[test]
    fn written_bytes_accumulate() {
        let before = bytes_written();
        record_bytes_written(1024);
        assert_eq!(bytes_written(), before + 1024);
    }

    /// The watch must stop when dropped, rather than leaving a sampler thread behind.
    #[test]
    fn dropping_the_watch_stops_it() {
        let watch = ResourceWatch::start(Duration::from_millis(50), |_| {});
        thread::sleep(Duration::from_millis(120));
        drop(watch);
    }
}
