#!/usr/bin/env python3
"""Re-import D2C Big Y VCFs so their variant sets carry per-call evidence.

Variant sets imported before migration 0042 kept only contig/position/ref/alt/rsID/genotype — QUAL,
FILTER, DP, GQ and AD were parsed and discarded. Nothing downstream could then tell a 40x hom-alt
call from a one-read artefact, which is exactly the judgement a private-Y engine has to make. The
importer now captures all of it, but only for *newly* imported sets: the stored rows can't be
back-filled, because the evidence was never written. The source VCFs are still on disk, so re-import
is the fix.

    ./refresh_variants_d2c.py              # dry run (the default)
    ./refresh_variants_d2c.py --apply      # back up, then refresh
    ./refresh_variants_d2c.py --limit 20   # bounded smoke test first

This is NOT a re-run of `load_d2c.py`. That loader would skip nearly every one of these subjects
("matched subject already has sequencing data"), and variant-set import is not content-idempotent —
anything that did re-run would *duplicate* its set rather than upgrade it. So this refreshes
deliberately, one subject at a time.

Order matters: **import first, delete second.** Ingesting produces a second set for the subject, and
only once the new evidence-bearing set is confirmed present is the old one removed. A crash mid-run
therefore leaves a duplicate — visible and easy to clean — never a subject with no Y variant data.
"""

from __future__ import annotations

import argparse
import csv
import os
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
BIG_Y_LABEL = "FTDNA Big Y (aengine)"
CALL_SCHEMA_EVIDENCE = 2


def backup_db(db: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dest = db.with_name(f"{db.name}.bak-{stamp}")
    need = db.stat().st_size
    if shutil.disk_usage(db.parent).free < need * 1.1:
        sys.exit(f"error: not enough space to back up {db.name}")
    print(f"backing up {db} -> {dest} ({need/2**30:.1f} GiB) …", flush=True)
    started = time.time()
    subprocess.run(["sqlite3", str(db), f".backup '{dest}'"], check=True)
    subprocess.run(["sqlite3", str(dest), "PRAGMA quick_check;"], check=True, stdout=subprocess.DEVNULL)
    print(f"backup ok ({time.time()-started:.0f}s) — restore with:  cp '{dest}' '{db}'\n", flush=True)
    return dest


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


def require_schema_column(db: Path) -> None:
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    cols = {r[1] for r in con.execute("PRAGMA table_info(variant_set)")}
    con.close()
    if "call_schema" not in cols:
        sys.exit(
            "error: this workspace predates migration 0042 (no variant_set.call_schema).\n"
            "       Open the app once (or run any navigator command) to apply migrations, then re-run."
        )


def plan(db: Path, root: Path, manifest: Path) -> list[dict]:
    """Every evidence-less Big Y set, paired with the VCF that can replace it. Reads only."""
    # The manifest is indexed two ways, because `load_d2c.py` identified subjects two ways: a kit it
    # already knew became an `external_id`, but most rows created a *new* subject named `[lab]-[kit]`
    # (the manifest's own `name` column), which lands in `donor_identifier` with no external id at
    # all. Keying on the kit alone finds only the minority that were matched.
    by_kit: dict[str, str] = {}
    by_name: dict[str, str] = {}
    with manifest.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["lab"] == "FTDNA" and row["y_tier"] == "vcf_aengine_native":
                artifact = row["y_artifact"] or ""
                by_kit.setdefault(row["kit"], artifact)
                by_name.setdefault(row["name"], artifact)

    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    rows = con.execute(
        """SELECT vs.id, vs.biosample_guid, b.donor_identifier, e.external_id
             FROM variant_set vs
             JOIN biosample b ON b.guid = vs.biosample_guid
        LEFT JOIN external_id e ON e.biosample_guid = vs.biosample_guid AND e.source = 'FTDNA'
            WHERE vs.source_label = ? AND vs.call_schema < ?
            ORDER BY vs.id""",
        (BIG_Y_LABEL, CALL_SCHEMA_EVIDENCE),
    ).fetchall()
    con.close()

    out = []
    for set_id, guid, donor, kit in rows:
        item = {
            "set": set_id,
            "subject": guid,
            "donor": donor,
            "kit": kit or "",
            "vcf": "",
            "action": "",
            "why": "",
        }
        # Prefer the vendor id when the subject carries one (that is how `load_d2c.py` addressed it);
        # otherwise fall back to the donor identifier the loader minted from the manifest's name.
        artifact = by_kit.get(kit or "") or by_name.get(donor or "")
        if not artifact:
            item["action"] = "skip"
            item["why"] = "no aengine VCF row in the manifest for this subject"
        else:
            vcf = artifact.replace(SERVER_PREFIX, str(root))
            item["vcf"] = vcf
            if not os.path.exists(vcf):
                item["action"], item["why"] = "skip", "source VCF not on disk"
            else:
                item["action"], item["why"] = "refresh", ""
        out.append(item)
    return out


def backfill_source_paths(db: Path, root: Path, manifest: Path, apply: bool) -> Counter:
    """Record `source_path` on Big Y sets that lack it, from the manifest's aengine VCF row.

    A set that already carries evidence needs no re-import to gain a path — the calls are identical
    either way — so this is the cheap equivalent of re-running the refresh purely to populate the
    column that lets the VCF be re-read for tree-position genotyping.
    """
    by_kit: dict[str, str] = {}
    by_name: dict[str, str] = {}
    with manifest.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["lab"] == "FTDNA" and row["y_tier"] == "vcf_aengine_native":
                artifact = row["y_artifact"] or ""
                by_kit.setdefault(row["kit"], artifact)
                by_name.setdefault(row["name"], artifact)

    con = sqlite3.connect(db)
    rows = con.execute(
        """SELECT vs.id, b.donor_identifier, e.external_id
             FROM variant_set vs
             JOIN biosample b ON b.guid = vs.biosample_guid
        LEFT JOIN external_id e ON e.biosample_guid = vs.biosample_guid AND e.source = 'FTDNA'
            WHERE vs.source_label = ? AND (vs.source_path IS NULL OR vs.source_path = '')""",
        (BIG_Y_LABEL,),
    ).fetchall()

    stat, updates = Counter(), []
    for set_id, donor, kit in rows:
        artifact = by_kit.get(kit or "") or by_name.get(donor or "")
        if not artifact:
            stat["no manifest row"] += 1
            continue
        path = artifact.replace(SERVER_PREFIX, str(root))
        if not os.path.exists(path):
            stat["VCF not on disk"] += 1
            continue
        updates.append((path, set_id))
        stat["resolved"] += 1
    if apply and updates:
        with con:
            con.executemany("UPDATE variant_set SET source_path = ? WHERE id = ?", updates)
    con.close()
    return stat


def sets_for(con: sqlite3.Connection, guid: str) -> dict[int, int]:
    """variant_set id -> call_schema, for this subject's Big Y sets."""
    return {
        r[0]: r[1]
        for r in con.execute(
            "SELECT id, call_schema FROM variant_set WHERE biosample_guid = ? AND source_label = ?",
            (guid, BIG_Y_LABEL),
        )
    }


def delete_set(con: sqlite3.Connection, set_id: int) -> None:
    """Children first, then the parent — the same order the store's `variant_set::delete` uses."""
    with con:
        con.execute("DELETE FROM variant_call WHERE variant_set_id = ?", (set_id,))
        con.execute("DELETE FROM variant_set WHERE id = ?", (set_id,))


def refresh(db: Path, navigator: str, items: list[dict]) -> Counter:
    stat = Counter()
    con = sqlite3.connect(db)
    for n, item in enumerate(items, 1):
        before = set(sets_for(con, item["subject"]))
        args = [navigator, "ingest", "--db", str(db)]
        if item["kit"]:
            args += ["--external-id", item["kit"], "--id-source", "FTDNA", "--skip-unmatched"]
        else:
            # The subject exists (it owns the set we're replacing), so an exact donor match finds it
            # rather than creating anything.
            args += ["--subject", item["donor"]]
        args.append(item["vcf"])
        res = subprocess.run(args, capture_output=True, text=True)
        if res.returncode != 0:
            stat["import-failed"] += 1
            print(f"  [{n}/{len(items)}] {item['donor']}: FAILED — {(res.stderr or res.stdout).strip()[-160:]}")
            continue

        after = sets_for(con, item["subject"])
        fresh = [sid for sid in set(after) - before if after[sid] >= CALL_SCHEMA_EVIDENCE]
        if not fresh:
            # The import reported success but produced nothing we can trust — leave the old set alone.
            stat["no-evidence-set"] += 1
            print(f"  [{n}/{len(items)}] {item['donor']}: imported but no evidence-bearing set appeared; old set kept")
            continue
        delete_set(con, item["set"])
        stat["refreshed"] += 1
        if n % 100 == 0 or n == len(items):
            print(f"  [{n}/{len(items)}] refreshed={stat['refreshed']} failed={stat['import-failed']}", flush=True)
    con.close()
    return stat


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--d2c-root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--manifest", type=Path, help="default: <d2c-root>/_manifests/biosample_map.tsv")
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--navigator", type=Path, help="navigator binary (default: the repo's build)")
    ap.add_argument("--apply", action="store_true", help="write (default: dry run)")
    ap.add_argument("--no-backup", action="store_true", help="skip the pre-flight backup (NOT recommended)")
    ap.add_argument("--limit", type=int, help="refresh at most N subjects (a bounded smoke test)")
    ap.add_argument(
        "--backfill-source-path",
        action="store_true",
        help="only record variant_set.source_path from the manifest, without re-importing. The calls "
        "themselves are unchanged by a re-import once they already carry evidence, so this reaches the "
        "same end state in seconds instead of re-parsing every VCF",
    )
    ap.add_argument("--report", type=Path, default=Path("d2c_refresh_report.tsv"))
    args = ap.parse_args()

    manifest = args.manifest or args.d2c_root / "_manifests" / "biosample_map.tsv"
    for p in (args.db, manifest):
        if not p.is_file():
            sys.exit(f"error: {p} not found")
    require_schema_column(args.db)

    if args.backfill_source_path:
        stat = backfill_source_paths(args.db, args.d2c_root, manifest, args.apply)
        print(f"Big Y sets missing source_path: {sum(stat.values())}")
        for k, v in stat.most_common():
            print(f"  {v:5d}  {k}")
        print("\ndry run — nothing written." if not args.apply else f"\nrecorded {stat['resolved']} source path(s).")
        return

    items = plan(args.db, args.d2c_root, manifest)
    with args.report.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(items[0].keys()) if items else ["set"])
        w.writeheader()
        w.writerows(items)

    counts = Counter(i["action"] for i in items)
    print(f"evidence-less Big Y variant sets: {len(items)}")
    print(f"  {counts['refresh']:5d}  refreshable")
    for why, n in Counter(i["why"] for i in items if i["action"] == "skip").most_common():
        print(f"  {n:5d}  skip — {why}")
    print(f"\nper-row plan -> {args.report}")

    todo = [i for i in items if i["action"] == "refresh"]
    if args.limit:
        todo = todo[: args.limit]
    if not args.apply:
        print(f"\ndry run — nothing written. --apply would refresh {len(todo)} subject(s).")
        return
    if not todo:
        print("\nnothing to refresh.")
        return
    if not args.no_backup:
        backup_db(args.db)
    print(f"refreshing {len(todo)} subject(s) — import first, delete the stale set only on success:")
    stat = refresh(args.db, find_navigator(args.navigator), todo)
    for k, v in stat.most_common():
        print(f"  {v:5d}  {k}")
    print(
        "\nY placement is derived from these calls, so re-run the subjects' Y analysis (or the\n"
        "project-wide analyze) to reconcile against the refreshed sets."
    )


if __name__ == "__main__":
    main()
