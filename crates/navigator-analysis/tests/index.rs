//! Build a coordinate index, a `.bai` or a `.crai`, for a fixture that has none. Then confirm two
//! things: a region query on that file now works, and the index matches the one that the repo holds
//! beside the fixture.

use std::fs;
use std::path::{Path, PathBuf};

use navigator_analysis::index::{ensure_index, index_path_for};
use navigator_analysis::reader::has_region_index;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// A unique scratch dir under the system temp dir (this crate deliberately avoids a `tempfile` dep).
fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-index-test-{}-{tag}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// Copy the alignment into `dir`, and not the index that sits beside it. `ensure_index` must then
/// build one.
fn copy_without_index(name: &str, dir: &Path) -> PathBuf {
    let dst = dir.join(name);
    fs::copy(fixtures().join(name), &dst).unwrap();
    dst
}

#[test]
fn builds_bai_when_missing_and_matches_fixture() {
    let dir = scratch("bam");
    let bam = copy_without_index("coverage.bam", &dir);
    assert!(!has_region_index(&bam), "no index should exist in the scratch copy");

    let mut ticks = 0usize;
    let built = ensure_index(&bam, None, &mut |_done, total| {
        assert!(total.is_some(), "BAM progress should be determinate");
        ticks += 1;
    })
    .unwrap();

    let expected = index_path_for(&bam); // coverage.bam.bai
    assert_eq!(built.as_deref(), Some(expected.as_path()));
    assert!(has_region_index(&bam), "index should now be present");
    assert!(ticks >= 1, "at least the final 100% tick should fire");

    // Byte-identical to the fixture's checked-in .bai (same records, same sort order).
    let ours = fs::read(&expected).unwrap();
    let reference = fs::read(fixtures().join("coverage.bam.bai")).unwrap();
    assert_eq!(ours, reference, "built .bai should match the fixture index");

    // A second call is a no-op (index already present).
    let again = ensure_index(&bam, None, &mut |_, _| panic!("should not re-index")).unwrap();
    assert_eq!(again, None);

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn builds_crai_when_missing() {
    let dir = scratch("cram");
    let cram = copy_without_index("coverage.cram", &dir);
    assert!(!has_region_index(&cram));

    let built = ensure_index(&cram, None, &mut |_done, total| {
        assert!(total.is_none(), "CRAM progress should be indeterminate");
    })
    .unwrap();

    let expected = index_path_for(&cram); // coverage.cram.crai
    assert_eq!(built.as_deref(), Some(expected.as_path()));
    assert!(has_region_index(&cram), "crai should now be present");
    assert!(expected.exists());

    let _ = fs::remove_dir_all(&dir);
}
