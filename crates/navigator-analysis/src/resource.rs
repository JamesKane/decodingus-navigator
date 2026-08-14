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
//! ([`crate::postprocess::bamio`]); this exists so that the next time something goes wrong there is
//! a record of what the machine looked like, instead of an inference from a crash report.
//!
//! ## Portability
//!
//! Every probe here is `sysinfo`, which binds the platform APIs through pure-Rust crates on all
//! three desktop targets — the same reason `navigator-align` picked it for RAM detection. There is
//! deliberately no `fcntl`/`ioctl`/`/proc` in this module: the guard has to hold on Windows, and a
//! guard that only arms itself on macOS would have been no guard at all for most users.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

/// Bytes handed to disk by the post-processing writers since the process started.
///
/// Self-accounted rather than read back from the OS: every platform exposes per-process I/O
/// counters differently (and Windows' are not in `sysinfo`'s default surface), whereas the writers
/// already know exactly how much they wrote. It costs one relaxed add per buffer.
static BYTES_WRITTEN: AtomicU64 = AtomicU64::new(0);

/// Account `n` bytes written. Called from the write path; must stay this cheap.
pub fn record_bytes_written(n: u64) {
    BYTES_WRITTEN.fetch_add(n, Ordering::Relaxed);
}

/// Total bytes the post-processing writers have written this process.
pub fn bytes_written() -> u64 {
    BYTES_WRITTEN.load(Ordering::Relaxed)
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
    /// Bytes written by the post-processing writers, in total and since the previous sample.
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
