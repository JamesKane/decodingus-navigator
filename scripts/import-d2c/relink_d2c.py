#!/usr/bin/env python3
"""Repoint workspace alignments whose files moved when the D2C collection was restructured.

The collection was reorganized from a kit-keyed layout
(`D2C/FTDNA/<kit>/<test>/CP086569.2/chrYM.cram`) to the UUID-keyed one the manifest describes
(`D2C/<subject>/<run>/CP086569.2/chrYM.cram`). Alignment rows still point at the old paths, so the
files read as missing: analysis silently degrades and `AlignmentArtifacts` treats every cached
artifact as un-stale-able ("the file is gone, so trust the cache").

`biosample_map.tsv` is the authority for where a kit's CRAM lives now, so this rewrites from the
manifest rather than by guessing at the directory tree.

    ./relink_d2c.py            # dry run (the default): classify everything, touch nothing
    ./relink_d2c.py --apply    # back the workspace up, then rewrite the safe subset

Only **relocations** are applied — a row whose reference build matches the file it is being
repointed at. A row whose build does *not* match is reported and left alone: repointing a GRCh38
alignment at a CHM13 CRAM would not fix a path, it would make the row lie about its reference.
Those need a real ingest of the CRAM (a new sequencing run), which `load_d2c.py` does properly.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

DEFAULT_ROOT = Path("/Volumes/Genomics/D2C")
DEFAULT_DB = Path.home() / ".decodingus" / "navigator-rs.db"
SERVER_PREFIX = "/mnt/md0/Repo"

# Build tokens that mean the same reference, so a relocation isn't rejected over spelling.
BUILD_ALIASES = {
    "chm13v2.0": "chm13",
    "chm13v2": "chm13",
    "chm13": "chm13",
    "hs1": "chm13",
    "t2t-chm13v2.0": "chm13",
    "grch38": "grch38",
    "hg38": "grch38",
    "grch37": "grch37",
    "hg19": "grch37",
    "b37": "grch37",
}


def canonical_build(b: str | None) -> str:
    return BUILD_ALIASES.get((b or "").strip().lower(), (b or "").strip().lower())


def build_of_path(p: str) -> str:
    """The reference a D2C artifact path implies. The collection encodes it as a directory:
    `CP086569.2` is the CHM13 chrY/chrM accession; `b38` is GRCh38."""
    parts = p.split("/")
    if "CP086569.2" in parts:
        return "chm13"
    if "b38" in parts:
        return "grch38"
    if "b37" in parts:
        return "grch37"
    return ""


def backup_db(db: Path) -> Path:
    """Consistent copy via SQLite's own backup API — a plain `cp` would tear the .db from its -wal."""
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dest = db.with_name(f"{db.name}.bak-{stamp}")
    need = db.stat().st_size
    free = shutil.disk_usage(db.parent).free
    if free < need * 1.1:
        sys.exit(f"error: need ~{need/2**30:.1f} GiB to back up {db.name}, only {free/2**30:.1f} GiB free")
    print(f"backing up {db} -> {dest} ({need/2**30:.1f} GiB) …", flush=True)
    started = time.time()
    subprocess.run(["sqlite3", str(db), f".backup '{dest}'"], check=True)
    subprocess.run(["sqlite3", str(dest), "PRAGMA quick_check;"], check=True, stdout=subprocess.DEVNULL)
    print(f"backup ok ({time.time()-started:.0f}s) — restore with:  cp '{dest}' '{db}'\n", flush=True)
    return dest


def manifest_index(manifest: Path) -> dict[tuple[str, str], dict]:
    """(lab, kit) → manifest row. The kit is how an old-layout path identifies itself."""
    idx: dict[tuple[str, str], dict] = {}
    with manifest.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            idx.setdefault((row["lab"], row["kit"]), row)
    return idx


KIT_PATTERNS = (
    re.compile(r"/D2C/FTDNA/([^/]+)/"),  # the old kit-keyed collection layout
    re.compile(r"/Downloads/FTDNA/([^/]+)/"),  # a local vendor download, since cleaned up
)


def kit_of(path: str) -> str | None:
    for pat in KIT_PATTERNS:
        m = pat.search(path)
        if m:
            return m.group(1)
    return None


def classify(db: Path, root: Path, index: dict) -> list[dict]:
    """Decide what to do with every alignment whose file is missing. Pure: reads only."""
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    rows = con.execute(
        """SELECT a.id, a.bam_path, a.reference_build, sr.biosample_guid, b.donor_identifier
             FROM alignment a
             JOIN sequence_run sr ON sr.id = a.sequence_run_id
             JOIN biosample b ON b.guid = sr.biosample_guid
            WHERE a.bam_path IS NOT NULL"""
    ).fetchall()

    out = []
    for aln_id, path, build, guid, donor in rows:
        if os.path.exists(path):
            continue
        item = {
            "alignment": aln_id,
            "subject": guid,
            "donor": donor,
            "old_path": path,
            "build": build,
            "new_path": "",
            "action": "",
            "why": "",
        }
        kit = kit_of(path)
        if not kit:
            item["action"], item["why"] = "unresolved", "no kit in the old path"
            out.append(item)
            continue
        row = index.get(("FTDNA", kit))
        if not row:
            item["action"], item["why"] = "unresolved", f"kit {kit} not in the manifest"
            out.append(item)
            continue
        cram = (row.get("cram") or "").strip()
        if not cram:
            item["action"], item["why"] = "unresolved", "manifest lists no CRAM for this kit"
            out.append(item)
            continue

        new = cram.replace(SERVER_PREFIX, str(root))
        item["new_path"] = new
        if not os.path.exists(new):
            item["action"], item["why"] = "unresolved", "manifest CRAM is not on disk either"
        elif canonical_build(build) != build_of_path(new):
            # Not a moved file: a different alignment of the same donor, against another reference.
            item["action"] = "reingest"
            item["why"] = f"build mismatch — row is {build}, target is {build_of_path(new) or '?'}"
        else:
            item["action"], item["why"] = "relink", "same build; the file moved"
        out.append(item)
    return out


def apply_relinks(db: Path, plan: list[dict]) -> int:
    relinks = [i for i in plan if i["action"] == "relink"]
    if not relinks:
        return 0
    con = sqlite3.connect(db)
    with con:  # one transaction — all or nothing
        for item in relinks:
            con.execute("UPDATE alignment SET bam_path = ? WHERE id = ?", (item["new_path"], item["alignment"]))
    con.close()
    return len(relinks)


def find_navigator(explicit: Path | None) -> str:
    """The navigator binary to drive.

    Picks the **most recently built** of release/debug rather than preferring release outright: a
    stale release build silently lacks migrations the workspace has already applied, and sqlx then
    refuses to open it ("migration N was previously applied but is missing in the resolved
    migrations"). Rejects a binary older than the newest migration for the same reason — better a
    clear "rebuild" than a run that fails on every subject.
    """
    if explicit:
        if not explicit.is_file():
            sys.exit(f"error: {explicit} not found")
        return str(explicit)
    repo = Path(__file__).resolve().parents[2]
    cands = [p for p in (repo / "target/release/navigator", repo / "target/debug/navigator") if p.is_file()]
    if not cands:
        sys.exit("error: navigator binary not found — build it (cargo build --release -p navigator-ui) or pass --navigator")
    chosen = max(cands, key=lambda p: p.stat().st_mtime)
    migrations = repo / "crates/navigator-store/migrations"
    newest = max((p.stat().st_mtime for p in migrations.glob("*.up.sql")), default=0)
    if chosen.stat().st_mtime < newest:
        sys.exit(
            f"error: {chosen} predates the newest migration in {migrations.name}/.\n"
            "       Rebuild (cargo build --release -p navigator-ui) or pass --navigator."
        )
    print(f"using {chosen}")
    return str(chosen)


def kit_for_subject(con: sqlite3.Connection, guid: str) -> str | None:
    row = con.execute(
        "SELECT external_id FROM external_id WHERE biosample_guid = ? AND source = 'FTDNA'", (guid,)
    ).fetchone()
    return row[0] if row else None


def apply_reingest(db: Path, navigator: str, plan: list[dict], project: str | None, limit: int | None) -> tuple[int, int]:
    """Ingest the CHM13 CRAM for each build-mismatched orphan, as a *new* sequencing run.

    The dead GRCh38 row is deliberately left in place: it records a file the workspace once had, and
    removing an alignment cascades into its derived analysis. Prune it separately once the new run
    is confirmed good.
    """
    targets = [i for i in plan if i["action"] == "reingest"]
    if limit:
        targets = targets[:limit]
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    ok = failed = 0

    def registered(item) -> bool:
        return (
            con.execute(
                """SELECT COUNT(*) FROM alignment a JOIN sequence_run sr ON sr.id = a.sequence_run_id
                    WHERE sr.biosample_guid = ? AND a.bam_path = ?""",
                (item["subject"], item["new_path"]),
            ).fetchone()[0]
            > 0
        )

    for n, item in enumerate(targets, 1):
        # A subject often owns *two* orphans — the dead vendor BAM and the moved CRAM — and the
        # relink pass above already repointed the latter at this exact file. Navigator would then
        # skip the import as an already-linked alignment and still exit 0, so checking first is the
        # difference between reporting "ingested" and reporting the truth.
        if registered(item):
            print(f"  [{n}/{len(targets)}] {item['donor']}: already linked (the relink pass covered it) — skipped")
            continue
        kit = kit_for_subject(con, item["subject"])
        args = [navigator, "ingest", "--db", str(db)]
        if kit:
            # Resolve by vendor id so the CRAM lands on the subject that already exists.
            args += ["--external-id", kit, "--id-source", "FTDNA", "--skip-unmatched"]
        else:
            args += ["--subject", item["donor"]]
        # A CRAM ships no .bai, so the coverage-shape detector would fall back to WGS.
        args += ["--test-type", "Big Y"]
        if project:
            args += ["--project", project]
        args.append(item["new_path"])

        print(f"  [{n}/{len(targets)}] {item['donor']} (kit {kit or '—'}) …", flush=True)
        res = subprocess.run(args, capture_output=True, text=True)
        if res.returncode != 0:
            failed += 1
            print(f"      FAILED: {(res.stderr or res.stdout).strip().splitlines()[-1:]}", flush=True)
            continue
        # Exit 0 is not proof of work: `--skip-unmatched` exits 0 when no subject matched, and an
        # already-linked alignment is skipped the same way. Confirm the row exists before claiming it.
        if registered(item):
            ok += 1
        else:
            failed += 1
            print("      reported success but no alignment appeared — check the subject match", flush=True)
    con.close()
    return ok, failed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--d2c-root", type=Path, default=DEFAULT_ROOT, help="collection root")
    ap.add_argument("--manifest", type=Path, help="default: <d2c-root>/_manifests/biosample_map.tsv")
    ap.add_argument("--db", type=Path, default=DEFAULT_DB, help="workspace database")
    ap.add_argument("--apply", action="store_true", help="write the relinks (default: dry run)")
    ap.add_argument(
        "--reingest",
        action="store_true",
        help="also ingest the CHM13 CRAM for build-mismatched orphans, as a new sequencing run",
    )
    ap.add_argument("--navigator", type=Path, help="navigator binary (default: the repo's build)")
    ap.add_argument("--project", default="D2C", help='project to file re-ingested runs under ("" to disable)')
    ap.add_argument("--limit", type=int, help="re-ingest at most N CRAMs (a bounded smoke test)")
    ap.add_argument("--no-backup", action="store_true", help="skip the pre-flight backup (NOT recommended)")
    ap.add_argument("--report", type=Path, default=Path("d2c_relink_report.tsv"), help="per-row decision")
    args = ap.parse_args()

    manifest = args.manifest or args.d2c_root / "_manifests" / "biosample_map.tsv"
    for p in (args.db, manifest):
        if not p.is_file():
            sys.exit(f"error: {p} not found")

    plan = classify(args.db, args.d2c_root, manifest_index(manifest))
    with args.report.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(plan[0].keys()) if plan else ["alignment"])
        w.writeheader()
        w.writerows(plan)

    counts = Counter(i["action"] for i in plan)
    print(f"missing alignment files: {len(plan)}")
    for action in ("relink", "reingest", "unresolved"):
        if counts[action]:
            print(f"  {counts[action]:5d}  {action}")
    for why, n in Counter(i["why"] for i in plan if i["action"] != "relink").most_common():
        print(f"          · {n:4d}  {why}")
    print(f"\nper-row decisions -> {args.report}")

    if counts["reingest"]:
        print(
            f"\n{counts['reingest']} row(s) point at a *different* alignment of the same donor (another\n"
            "reference), not a moved file. Rewriting those paths would make the row misreport its\n"
            "build, so they are left alone — ingest the CRAM properly with load_d2c.py instead."
        )

    if not args.apply:
        print("\ndry run — nothing written. Re-run with --apply to relink" + (", --reingest to import the CRAMs." if not args.reingest else "."))
        return
    if not counts["relink"] and not (args.reingest and counts["reingest"]):
        print("\nnothing to do.")
        return
    if not args.no_backup:
        backup_db(args.db)
    if counts["relink"]:
        print(f"relinked {apply_relinks(args.db, plan)} alignment(s).")
    if args.reingest and counts["reingest"]:
        print(f"\ningesting {counts['reingest']} CHM13 CRAM(s) as new sequencing runs:")
        ok, failed = apply_reingest(args.db, find_navigator(args.navigator), plan, args.project or None, args.limit)
        print(f"ingested {ok}, failed {failed}.")
        print(
            "the dead GRCh38 rows are left in place — they record files the workspace once had, and\n"
            "dropping an alignment cascades into its derived analysis. Prune them once you're happy."
        )


if __name__ == "__main__":
    main()
