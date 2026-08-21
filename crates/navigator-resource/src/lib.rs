//! A record of the load that a long stage puts on the machine.
//!
//! A realignment runs for hours. The user continues to use the machine during this time. The fault
//! that made this module necessary was not the fault that we monitored.
//!
//! On 2026-08-13, a WGS-scale sort stopped after six hours. We first thought that memory was the
//! cause. It was not. The compressor was idle. The system wrote no jetsam report. Anonymous memory
//! stayed at 32 GB of 128 GB.
//!
//! The operating system reported a different problem. The problem was the writes. The system filed
//! a resource notice against the process, because the process made 549 GB of file-backed memory
//! dirty. That rate is almost 1.4 times the sustained write-back limit of the system.
//!
//! The main thread of WindowServer then missed a 40-second watchdog check. The watchdog stopped
//! WindowServer. All processes in the login session stopped, and this job was one of them.
//!
//! So one number is not enough. Memory is easy to sample, and a memory sample shows immediately
//! that memory is not the cause. The write rate was the extreme value, but no code recorded it.
//! This module samples both values at an interval. It gives a [`ResourceSample`] to the caller, and
//! the caller writes the sample to the log.
//!
//! ## This module reports. It does not intervene.
//!
//! No code here stops a job. A stage that writes at a high rate does the correct work. The sort
//! must write hundreds of GB. A watchdog that stopped a six-hour run for high speed is worse than
//! the initial problem.
//!
//! To limit the damage, control the writes at the point where the code makes them. [`PacedFile`]
//! does this at a byte interval. This module gives you a record of the state of the machine. Before
//! this module, you could only make a deduction from a crash report.
//!
//! ## Why this is a crate
//!
//! A byte counter is correct only if there is one counter. The writers of the pipeline are in
//! crates that do not depend on each other. `navigator-align` maps the reads. `navigator-analysis`
//! reverts, sorts, marks, and compresses them.
//!
//! At first, this code was a module in `navigator-analysis`. So the mapping stage had no pace
//! control and no counter. That stage is the longest stage in the job, and it writes the 60 GB
//! `mapped.bam` file. The run log showed `0 MB/s` for the full stage. The stage was not quiet. No
//! code measured it.
//!
//! ## Portability
//!
//! Each probe uses `sysinfo`. That crate binds the platform APIs with pure-Rust code on all three
//! desktop targets. `navigator-align` uses `sysinfo` to find the quantity of RAM for the same
//! reason.
//!
//! This crate has no `fcntl` code, no `ioctl` code, and no `/proc` code. The guard must work on
//! Windows. A guard that starts only on macOS is not a guard for most users.

use std::fs::File;
use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

/// The quantity of bytes that the paced writers of the pipeline sent to disk. The count starts
/// when the process starts.
///
/// This crate counts the bytes. It does not read a count from the operating system. Each platform
/// gives the I/O counters for a process in a different form, and `sysinfo` does not show the
/// Windows counters by default. The writers already know the exact quantity. The count costs one
/// relaxed add for each buffer.
///
/// There is one counter for the full process. This is the reason for a separate crate. The writes
/// of a realignment come from two crates that can not see each other. Two counters would give two
/// incomplete answers.
static BYTES_WRITTEN: AtomicU64 = AtomicU64::new(0);

/// Add `n` bytes to the count. The write path calls this function, so it must stay small.
pub fn record_bytes_written(n: u64) {
    BYTES_WRITTEN.fetch_add(n, Ordering::Relaxed);
}

/// The total quantity of bytes that the post-process writers wrote in this process.
pub fn bytes_written() -> u64 {
    BYTES_WRITTEN.load(Ordering::Relaxed)
}

/// Bytes written between forced flushes to disk. See [`PacedFile`].
const DEFAULT_SYNC_MB: u64 = 256;

/// The quantity of data that can stay dirty before a paced stream sends it to disk.
///
/// `NAVIGATOR_IO_SYNC_MB=0` stops the pace control. The stream then has no limit.
fn sync_interval() -> u64 {
    std::env::var("NAVIGATOR_IO_SYNC_MB")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(DEFAULT_SYNC_MB)
        * 1024
        * 1024
}

/// A file that limits the quantity of its output that stays dirty in the page cache.
///
/// Without this control, a stage writes to the page cache at its maximum speed. The operating
/// system then does the write-back. This division of work is correct until the volume becomes very
/// large.
///
/// One WGS realignment made 549 GB of file-backed memory dirty. The process went past its
/// sustained write-back limit by 1.4 times, so macOS filed a disk-writes resource notice. The main
/// thread of WindowServer then missed a 40-second watchdog check and stopped. The login session
/// and a six-hour job stopped with it.
///
/// A flush at a byte interval limits the quantity of data that waits in the cache. The write path
/// does its own I/O during the run. It does not give the machine a debt to pay later. For this
/// reason the control is not a large delay: the same bytes go to the same disk at a more constant
/// rate.
///
/// This type is in this crate, not beside one stage, because each stage that writes tens of GB
/// needs it. These stages are the mapper with its 8.93 GB minimizer index, and the revert stage
/// with its spill runs and FASTQ. The sort with its runs and merged output and the last CRAM also
/// need it. The byte
/// count that [`ResourceWatch`] reports comes from the same place. So a paced stream is also a
/// counted stream, and a writer without this type is in neither total.
///
/// This type calls `sync_data`, not `sync_all`. The contents must be durable, but the metadata does
/// not need to be durable. On a stream of this size, the metadata is many thousands of inode
/// updates. `sync_data` is the portable name in std. It calls `fdatasync` where that call exists,
/// and `FlushFileBuffers` on Windows.
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

    /// Send all data that waits in the cache to disk.
    ///
    /// Call this at the end of a stream when a later run can trust the end-of-file marker of that
    /// stream. A marker in the page cache is a promise that the disk did not make.
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

/// The quantity of memory that an external-sort stage holds before it spills a run to disk.
///
/// The pipeline has two spill-to-disk stages: the collator of the revert stage and the coordinate
/// sort. A constant gave the size of each one, 256 MB and 512 MB. Nobody had run a WGS sample
/// through the two stages when they chose these values.
///
/// The measured result was bad. A 30x WGS coordinate sort spilled **688 runs**, and the merge
/// opens all of the runs at the same time. The memory has a limit by design, and the code works.
/// But on a machine with 128 GB of RAM, this fan-in is very large and gives no advantage. A laptop
/// that needs a small value uses the same constant.
///
/// So the machine gives the number. An explicit MB count in `var` has priority, because a run that
/// you must repeat or limit needs a fixed value. If `var` is empty, the rules are:
///
/// - **One quarter of the installed RAM.** Use the total RAM, not the free RAM. The free value
///   changes with the other applications of the user. The run count of a stage must not depend on
///   the browser. You can not repeat such behaviour from a bug report.
/// - **512 MB minimum.** This value is the current default of the sort. So no machine sorts with a
///   smaller value than it uses today.
/// - **8 GB maximum.** A larger value gives almost no advantage. The step from 688 runs to 88 runs
///   is large, but the step from 88 runs to 44 runs is small. The costs stay. The stable sort takes
///   half of the buffer again for scratch space. A larger record vector needs a second allocation,
///   and the first stays in memory at the same time.
/// - **Half of the free memory, maximum.** The rules above are for an idle machine. On a busy
///   machine, one more spilled run costs little, and a swap of the buffer costs much.
///
/// A platform that does not report the memory returns zero. Zero means *unknown*. Do not read zero
/// as "no memory". [`classify`] uses the same rule. An unknown machine gets the minimum value.
pub fn spill_budget(var: &str) -> u64 {
    if let Some(mb) = std::env::var(var).ok().and_then(|s| s.parse::<u64>().ok()) {
        return mb.max(1) * 1024 * 1024;
    }
    let mut system = sysinfo::System::new();
    system.refresh_memory();
    budget(system.total_memory(), system.available_memory())
}

/// The cost of one heap allocation above the quantity of bytes that the caller asked for. The cost
/// is the internal data of the allocator, and the increase to the next size class.
///
/// This constant is beside [`spill_budget`] because the two are one contract. A budget is correct
/// only if the total that fills it is correct. A stage that counts only the payload bytes holds
/// much more than its budget in real memory. This error was acceptable against a fixed 512 MB
/// constant with an unwritten margin. It is not acceptable against a fraction of the machine.
///
/// Sixteen bytes is the usual value for the allocators on the three desktop targets. This value is
/// an estimate, not an exact audit. It shows that a record with four small vectors costs much more
/// than the sum of the lengths of those vectors.
pub const ALLOCATION_OVERHEAD: usize = 16;

/// The floor, and the answer for a machine that will not say how much memory it has.
const MIN_SPILL_BUDGET: u64 = 512 << 20;
/// The maximum value. [`spill_budget`] gives the reason why a larger value has no advantage.
const MAX_SPILL_BUDGET: u64 = 8 << 30;

/// The size decision. It is separate from the probe, so a test can call it on any machine.
/// [`classify`] has the same separation for the same reason.
fn budget(total_memory: u64, available_memory: u64) -> u64 {
    if total_memory == 0 {
        return MIN_SPILL_BUDGET;
    }
    let budget = (total_memory / 4).clamp(MIN_SPILL_BUDGET, MAX_SPILL_BUDGET);
    if available_memory == 0 {
        return budget;
    }
    // The minimum value applies here also. With the old constant, a machine with this little
    // memory took 512 MB. To apply the busy-machine guard below that value gives a worse result,
    // and only looks careful.
    budget.min((available_memory / 2).max(MIN_SPILL_BUDGET))
}

/// The level of load on the machine.
///
/// There are bands, not one threshold, because the trend is the important measurement. A stage that
/// stays at [`Pressure::Elevated`] for one hour is not the same as a stage that reaches it once.
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
/// This value is the growth, not the absolute use. A desktop that ran for a week has swap in use
/// for other reasons. If the code reports that swap as a fault of this job, every run gives a false
/// alarm.
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
    /// One line of text for a log.
    ///
    /// The line is short by design. At an interval of 30 seconds, a six-hour stage makes 720 of
    /// these lines. A user must be able to read them quickly between the stage times.
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
    // A platform that does not report its memory returns zero. Zero means "unknown". Do not read
    // "unknown" as "critical". `navigator_app::realign_job::has_room` uses the same rule.
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
/// The value of 30 seconds comes from the event that we must see. The WindowServer watchdog starts
/// at 40 seconds. A longer interval can miss the full period in which the machine has a problem.
pub const DEFAULT_INTERVAL: Duration = Duration::from_secs(30);

/// A background sampler. It operates while the caller holds it.
///
/// When the caller drops it, the thread stops and joins. So a stage can not continue after its own
/// watch stops. An error can not leave a thread behind.
pub struct ResourceWatch {
    stop: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
}

impl ResourceWatch {
    /// Take a sample at each `interval` and give each sample to `report`.
    ///
    /// `report` runs on the sampler thread, so it must do very little work. The intended use is to
    /// write one line.
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

            // Read the stop flag much more often than the sample interval. A drop of the watch
            // then returns quickly. If not, the exit of a job waits for up to `interval`.
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

    /// New swap during the job is the signal. Swap that was present before the job is not a signal.
    /// For this reason the sample carries the growth, not the absolute use.
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

    /// The reason for the automatic size: a large machine must not spill hundreds of runs.
    #[test]
    fn a_large_machine_gets_the_ceiling() {
        assert_eq!(budget(128 << 30, 100 << 30), MAX_SPILL_BUDGET);
    }

    #[test]
    fn an_ordinary_machine_gets_a_quarter_of_it() {
        assert_eq!(budget(16 << 30, 12 << 30), 4 << 30);
    }

    /// A machine with little free memory gets a smaller buffer. One more spilled run costs one
    /// file. A swap of the buffer costs the full run.
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

    /// Unknown is not zero. For a platform that does not report the memory, do not calculate a
    /// size for a machine with no memory. Do not calculate a size for a large machine either.
    #[test]
    fn unknown_memory_gets_the_floor() {
        assert_eq!(budget(0, 0), MIN_SPILL_BUDGET);
        assert_eq!(budget(64 << 30, 0), MAX_SPILL_BUDGET.min(16 << 30));
    }

    /// The manual value must have priority. If not, you can not repeat a run on another machine.
    #[test]
    fn an_explicit_override_beats_the_machine() {
        let var = "NAVIGATOR_TEST_SPILL_MB_OVERRIDE";
        std::env::set_var(var, "7");
        assert_eq!(spill_budget(var), 7 * 1024 * 1024);
        std::env::remove_var(var);
    }

    /// An invalid manual value must not give a budget of zero bytes. Such a budget spills one run
    /// for each record.
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

    /// The watch must stop at the drop. It must not leave a sampler thread in the process.
    #[test]
    fn dropping_the_watch_stops_it() {
        let watch = ResourceWatch::start(Duration::from_millis(50), |_| {});
        thread::sleep(Duration::from_millis(120));
        drop(watch);
    }
}
