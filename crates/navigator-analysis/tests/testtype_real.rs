//! A check of the test-type inference on a real BAM. It carries `#[ignore]`, and it needs a local
//! indexed BAM.
//!
//! Run it against one BAM:
//!
//! ```text
//! TESTTYPE_BAM=/path/to.bam [TESTTYPE_PLATFORM=ILLUMINA] [TESTTYPE_VENDOR=FamilyTreeDNA] \
//!   cargo test -p navigator-analysis --test testtype_real -- --ignored --nocapture
//! ```
//!
//! It prints the coverage profile from the BAI, and the test-type code that it inferred. A real
//! Big Y, Y Elite or mtFull BAM must come out with that type, and a WGS BAM must not.

use navigator_analysis::testtype::{coverage_profile_from_bai, infer_test_type};

#[test]
#[ignore]
fn print_inferred_test_type() {
    let Ok(bam) = std::env::var("TESTTYPE_BAM") else {
        eprintln!("set TESTTYPE_BAM to run");
        return;
    };
    let platform = std::env::var("TESTTYPE_PLATFORM").ok();
    let vendor = std::env::var("TESTTYPE_VENDOR").ok();

    // Show what the header probe takes out: the platform, and the hint about the vendor. An
    // environment variable can override each one, for a test.
    let probe = navigator_analysis::probe::probe_alignment(std::path::Path::new(&bam)).ok();
    let probe_platform = probe.as_ref().and_then(|p| p.platform.clone());
    let probe_vendor = probe.as_ref().and_then(|p| p.vendor_hint.clone());
    let big_y_code = probe.as_ref().and_then(|p| p.big_y_code.clone());
    println!("BAM: {bam}");
    println!("probe platform={probe_platform:?} vendor_hint={probe_vendor:?} big_y_code={big_y_code:?}");

    let platform = platform.or(probe_platform);
    let vendor = vendor.or(probe_vendor);
    let profile = coverage_profile_from_bai(std::path::Path::new(&bam), None);
    println!("profile: {profile:?}");
    let tt = infer_test_type(
        profile.as_ref(),
        platform.as_deref(),
        vendor.as_deref(),
        None,
        big_y_code.as_deref(),
    );
    println!("inferred test_type: {tt:?}");
    assert!(
        profile.is_some(),
        "expected a BAI-derived coverage profile for an indexed BAM"
    );
}
