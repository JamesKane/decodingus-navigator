#!/usr/bin/env python3
"""Bulk-load the D2C sequencing collection into a Navigator workspace.

Drives `navigator ingest` once per sample, resolving each manifest row to the local files that
back it and to the workspace subject it belongs to.

    biosample_map.tsv row  ──resolve artifact──>  local dir/file
                           ──resolve subject───>  existing biosample, or a new [lab]-[kit]
                           ──navigator ingest──>  workspace

Subject resolution (the matching rule):
  1. An existing subject carrying this lab's kit as an `external_id` (source = the lab, e.g.
     FTDNA) is MATCHED and ingested into via `ingest --external-id` — never re-created.
  2. Failing that, an existing subject whose `donor_identifier` IS the raw kit id is matched.
  3. Otherwise a NEW subject is created, named `[lab]-[kit]` — taken from the manifest's own
     `name` column, which already renders that pattern and sanitizes the ~80 rows whose kit is a
     UUID or a stray file path (`Dante-0e74a433`, not `Dante-184b9e17-19f0-...`).

When several manifest rows land on the same subject, one wins by lab precedence
(FTDNA > FGC > YSEQ > Dante > Nebula > other), then by artifact richness. The rest are skipped.

Artifacts are resolved to the LOCAL tree: the manifest records server paths (`/mnt/md0/Repo/...`),
which are rewritten onto --d2c-root. The `ytree/flat` and `ytree/work/recover` Y artifacts have no
local counterpart, so those rows fall back to the per-run GATK output that backs them
(`<subject>/<run>/CP086569.2/gatk4/chrY.g.vcf.gz`).

Rows with a CRAM (or a recovered GVCF) are ingested as a sample DIRECTORY — `navigator ingest`
routes a directory through `add_sample_dir`, the sidecar fast path that places Y/mt from the GVCF
and fills sex / read-metrics / coverage from the text sidecars WITHOUT decoding the CRAM. Rows
whose only artifact is a VCF / STR / SNP csv are ingested as that single file.

Default is a DRY RUN. Pass --apply to write (and to back the database up first).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

# Lab precedence when several manifest rows resolve to the same subject: earlier wins.
LAB_PRECEDENCE = ["FTDNA", "FGC", "YSEQ", "Dante", "Nebula"]
SERVER_REPO_PREFIX = "/mnt/md0/Repo/"


def lab_rank(lab: str) -> int:
    """Precedence rank; every unlisted lab sorts last as "other", stably by name."""
    try:
        return LAB_PRECEDENCE.index(lab)
    except ValueError:
        return len(LAB_PRECEDENCE)


def artifact_rank(row: dict) -> int:
    """Richness of a row's data, for breaking a tie between two rows of the same lab.
    An alignment + GVCF beats a bare VCF, which beats an STR/SNP export, which beats nothing."""
    if row["cram"]:
        return 0
    tier = row["y_tier"]
    if tier in ("gvcf_haploid", "derived_rescue"):
        return 1
    if tier == "vcf_aengine_native":
        return 2
    if tier in ("str_only", "snp_panel"):
        return 3
    return 4


# ── artifact resolution ────────────────────────────────────────────────────────────────────────


def build_index(root: Path, cache: Path | None, refresh: bool) -> set[str]:
    """Every file under the D2C root, as an absolute-path set.

    One bulk walk, cached: the collection lives on an external volume where ~60k individual
    `stat` calls cost minutes, but a single traversal costs seconds.
    """
    if cache and cache.is_file() and not refresh:
        return {ln.rstrip("\n") for ln in cache.open() if ln.strip()}

    files: set[str] = set()
    for subject_dir in root.iterdir():
        # `_manifests`, `_repo_raw`, `_staging_review`, … are not sample directories.
        if not subject_dir.is_dir() or subject_dir.name.startswith(("_", ".")):
            continue
        for dirpath, dirnames, filenames in os.walk(subject_dir):
            dirnames[:] = [d for d in dirnames if not d.startswith(".")]
            for f in filenames:
                files.add(os.path.join(dirpath, f))

    if cache:
        cache.write_text("".join(f"{p}\n" for p in sorted(files)))
    return files


def to_local(server_path: str, root: Path) -> str | None:
    """Rewrite a manifest server path onto the local root. Only `/mnt/md0/Repo/...` has a
    counterpart; the `ytree/...` paths are server-side derivations and resolve elsewhere."""
    if server_path and server_path.startswith(SERVER_REPO_PREFIX):
        return str(root / server_path[len(SERVER_REPO_PREFIX) :])
    return None


def resolve_artifact(row: dict, root: Path, index: set[str], gvcfs: dict[str, list[str]]) -> tuple[str | None, str, str]:
    """Locate the local data backing one manifest row.

    Returns (path, mode, how): mode is "dir" (sidecar fast path) or "file"; path is None when the
    row has no data on disk.
    """
    cram = to_local(row["cram"], root)
    if cram and cram in index:
        # The CRAM's directory holds the whole staged sample: cram + coverage/stats/wgs_metrics
        # sidecars + gatk3/callable_status.bed + gatk4/chrY.g.vcf.gz.
        return os.path.dirname(cram), "dir", "cram"

    artifact = to_local(row["y_artifact"], root)
    if artifact and artifact in index:
        # An in-repo artifact: the b38 aengine variants.vcf.gz, or an STR / SNP csv.
        return artifact, "file", "repo-artifact"

    if row["y_artifact"]:
        # A `ytree/flat` or `ytree/work/recover` path: server-side derivations with no local copy.
        # The per-run GATK GVCF they were derived FROM is what we have, so ingest its directory.
        found = sorted(gvcfs.get(row["subject"], []))
        if found:
            return os.path.dirname(os.path.dirname(found[0])), "dir", "recovered-gvcf"
        return None, "", "unresolved"

    return None, "", "no-artifact"


# ── workspace state ────────────────────────────────────────────────────────────────────────────


def read_workspace(db: Path) -> tuple[dict[tuple[str, str], str], dict[str, str], set[str], set[str]]:
    """Read-only snapshot of what the workspace already knows: vendor ids, donor names, and the
    alignment paths already linked (so a re-run doesn't duplicate them)."""
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        ext = {(s.upper(), e): g for s, e, g in conn.execute("SELECT source, external_id, biosample_guid FROM external_id")}
        donors = {d: g for d, g in conn.execute("SELECT donor_identifier, guid FROM biosample")}
        linked = {p for (p,) in conn.execute("SELECT bam_path FROM alignment WHERE bam_path IS NOT NULL")}
        with_data = {
            g for (g,) in conn.execute(
                "SELECT DISTINCT sr.biosample_guid FROM sequence_run sr JOIN alignment a ON a.sequence_run_id = sr.id"
            )
        }
        return ext, donors, linked, with_data
    finally:
        conn.close()


def prune_empty_subject(db: Path, donor: str) -> bool:
    """Drop a subject this run created but failed to load anything into.

    `navigator ingest` creates the subject BEFORE it reads the file, so a corrupt or unrecognized
    input leaves an empty subject behind — and re-running would litter another one. Only a subject
    carrying no data of any kind is removed, so this can never take a real sample with it.
    """
    conn = sqlite3.connect(db)
    try:
        row = conn.execute(
            """SELECT guid FROM biosample b WHERE b.donor_identifier = ?
                 AND NOT EXISTS (SELECT 1 FROM sequence_run    WHERE biosample_guid = b.guid)
                 AND NOT EXISTS (SELECT 1 FROM variant_set     WHERE biosample_guid = b.guid)
                 AND NOT EXISTS (SELECT 1 FROM str_profile     WHERE biosample_guid = b.guid)
                 AND NOT EXISTS (SELECT 1 FROM chip_profile    WHERE biosample_guid = b.guid)
                 AND NOT EXISTS (SELECT 1 FROM haplogroup_call WHERE biosample_guid = b.guid)""",
            (donor,),
        ).fetchone()
        if not row:
            return False
        conn.execute("DELETE FROM biosample_project WHERE biosample_guid = ?", (row[0],))
        conn.execute("DELETE FROM biosample WHERE guid = ?", (row[0],))
        conn.commit()
        return True
    finally:
        conn.close()


def backup_db(db: Path) -> Path:
    """Consistent copy of the workspace via SQLite's own backup API — safe against the WAL, which
    a plain `cp` of the .db file would tear away from its -wal sidecar."""
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dest = db.with_name(f"{db.name}.bak-{stamp}")
    free = shutil.disk_usage(db.parent).free
    need = db.stat().st_size
    if free < need * 1.1:
        sys.exit(f"error: need ~{need/2**30:.1f} GiB to back up {db.name}, only {free/2**30:.1f} GiB free")

    print(f"backing up {db} -> {dest} ({need/2**30:.1f} GiB) …", flush=True)
    started = time.time()
    subprocess.run(["sqlite3", str(db), f".backup '{dest}'"], check=True)
    subprocess.run(["sqlite3", str(dest), "PRAGMA quick_check;"], check=True, stdout=subprocess.DEVNULL)
    print(f"backup ok ({time.time()-started:.0f}s) — restore with:  cp '{dest}' '{db}'\n", flush=True)
    return dest


# ── planning ───────────────────────────────────────────────────────────────────────────────────


def plan(rows, root, index, gvcfs, ext, donors, linked, with_data, done):
    """Decide, for every manifest row, what to do with it. Pure — touches nothing."""
    # The sample directories whose CRAM is already recorded as an alignment.
    linked_dirs = {os.path.dirname(p) for p in linked}

    resolved = []
    for row in rows:
        path, mode, how = resolve_artifact(row, root, index, gvcfs)
        lab, kit = row["lab"], row["kit"]

        # Match an existing subject: first by this lab's kit as a vendor id, then — FTDNA only — by a
        # subject whose donor_identifier IS the raw kit, which is how the GAP-loaded members are
        # named. The donor fallback is scoped to FTDNA deliberately: donor_identifier carries no
        # vendor namespace, so an unscoped match reads a bare "21418" as the same person whether it
        # is an FGC kit or an FTDNA kit. `target` is the dedup key — the subject, however named.
        guid = ext.get((lab.upper(), kit))
        if guid:
            target, match = guid, "ext"
        elif lab.upper() == "FTDNA" and kit in donors:
            guid, target, match = donors[kit], donors[kit], "donor"
        else:
            # New subject, named [lab]-[kit] — the manifest's `name` column already renders that
            # pattern and sanitizes the rows whose kit is a UUID or a stray path.
            target, match = row["name"], "new"

        resolved.append(dict(row=row, path=path, mode=mode, how=how, target=target, match=match, guid=guid))

    # One subject, one ingest: when several rows land on the same target, the highest-precedence
    # lab wins, then the richest artifact, then the subject uuid for a stable tie-break.
    by_target = defaultdict(list)
    for item in resolved:
        by_target[item["target"]].append(item)

    for group in by_target.values():
        if len(group) == 1:
            continue
        # A row with nothing on disk can never win: it would supersede — and so discard — the only
        # loadable copy of the sample. Among rows that DO have data, lab precedence decides.
        group.sort(
            key=lambda i: (i["path"] is None, lab_rank(i["row"]["lab"]), artifact_rank(i["row"]), i["row"]["subject"])
        )
        winner = group[0]
        for loser in group[1:]:
            loser["action"] = "skip-superseded"
            loser["why"] = (
                f'{loser["row"]["lab"]}/{loser["row"]["subject"][:8]} yields to '
                f'{winner["row"]["lab"]}/{winner["row"]["subject"][:8]}'
            )

    for item in resolved:
        if item.get("action"):
            continue
        row = item["row"]
        if item["path"] is None:
            item["action"] = "skip-no-data"
            item["why"] = "manifest lists no artifact" if item["how"] == "no-artifact" else "artifact not on disk"
        elif row["subject"] in done:
            item["action"] = "skip-done"
            item["why"] = "already ingested (state ledger)"
        elif item["match"] != "new" and item["guid"] in with_data:
            # A matched subject that already carries sequencing data: the files are present, so
            # there is nothing to add. (This is the "…added when the sequencing files are not
            # present" rule — we never stack a second dataset onto an already-sequenced subject.)
            item["action"] = "skip-present"
            item["why"] = "matched subject already has sequencing data"
        elif item["mode"] == "dir" and item["path"] in linked_dirs:
            item["action"] = "skip-present"
            item["why"] = "alignment already linked"
        elif row["y_tier"] == "snp_panel":
            # FTDNA's SNPs.csv is a COMMA-separated Y-SNP panel with the verdict in column 2.
            # `looks_like_ysnp_panel` (navigator-domain/src/filetype.rs) only splits on TAB and only
            # reads a verdict from column 3, so it never matches and `ingest` fails with "could not
            # recognize the data". Skip rather than create an empty subject on a guaranteed failure.
            # Extending that detector to accept the comma form would unlock these rows.
            item["action"] = "skip-unsupported"
            item["why"] = "comma-separated FTDNA SNPs.csv — navigator's Y-SNP detector is tab-only"
        else:
            item["action"] = "ingest"
            item["why"] = ""

    return resolved


def ingest_args(item, project: str | None) -> list[str]:
    """The `navigator ingest` invocation for one planned row."""
    row = item["row"]
    args = ["ingest"]
    if item["match"] == "ext":
        # Resolve by vendor id — the subject must already exist; `ingest` never creates one here.
        args += ["--external-id", row["kit"], "--id-source", row["lab"].upper(), "--skip-unmatched"]
    elif item["match"] == "donor":
        # An existing subject whose donor_identifier is the raw kit; --subject finds it by exact match.
        args += ["--subject", row["kit"]]
    else:
        # No such subject: --subject creates it, named [lab]-[kit].
        args += ["--subject", row["name"]]
    # A CRAM ships no .bai, so the coverage-shape detector can't tell Big Y from WGS and would fall
    # back to WGS. For FTDNA the layout names the test definitively, so force it.
    if row["lab"] == "FTDNA" and row["cram"]:
        args += ["--test-type", "Big Y"]
    if project:
        args += ["--project", project]
    args += [item["path"]]
    return args


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=Path("/Volumes/Genomics/D2C/_manifests/biosample_map.tsv"))
    ap.add_argument("--d2c-root", type=Path, default=Path("/Volumes/Genomics/D2C"))
    ap.add_argument("--db", type=Path, default=Path.home() / ".decodingus" / "navigator-rs.db")
    ap.add_argument("--navigator", type=Path, help="navigator binary (default: the repo's release build)")
    ap.add_argument("--project", default="D2C", help='project to file subjects under ("" to disable; default: D2C)')
    ap.add_argument("--apply", action="store_true", help="write to the workspace (default: dry run)")
    ap.add_argument("--no-backup", action="store_true", help="skip the pre-flight backup (NOT recommended)")
    ap.add_argument("--limit", type=int, help="ingest at most N samples (a bounded smoke test)")
    ap.add_argument("--lab", action="append", help="only this lab (repeatable) — for a staged rollout")
    ap.add_argument("--mode", choices=["dir", "file"], help="only samples ingested as a directory / single file")
    ap.add_argument("--match", choices=["ext", "donor", "new"], help="only matched-by-vendor-id / matched-by-donor / new subjects")
    ap.add_argument("--state", type=Path, help="resumable ledger of ingested subjects (default: alongside the report)")
    ap.add_argument("--report", type=Path, default=Path("d2c_load_report.tsv"))
    ap.add_argument("--index-cache", type=Path, default=Path("d2c_files.txt"))
    ap.add_argument("--refresh-index", action="store_true", help="re-walk the D2C tree instead of using the cache")
    args = ap.parse_args()

    navigator = args.navigator
    if not navigator:
        repo = Path(__file__).resolve().parents[2]
        for cand in (repo / "target/release/navigator", repo / "target/debug/navigator"):
            if cand.is_file() and os.access(cand, os.X_OK):
                navigator = cand
                break
    if not navigator or not Path(navigator).is_file():
        sys.exit("error: navigator binary not found — build it (cargo build --release -p navigator-ui) or pass --navigator")
    for p, what in ((args.manifest, "manifest"), (args.d2c_root, "D2C root"), (args.db, "workspace DB")):
        if not p.exists():
            sys.exit(f"error: {what} not found: {p} (is the volume mounted?)")

    state_path = args.state or args.report.with_suffix(".state")
    done = set()
    if state_path.is_file():
        done = {json.loads(ln)["subject"] for ln in state_path.open() if ln.strip()}

    print(f"manifest  : {args.manifest}")
    print(f"d2c root  : {args.d2c_root}")
    print(f"db        : {args.db}")
    print(f"navigator : {navigator}")
    print(f"project   : {args.project or '(none)'}")
    print(f"mode      : {'APPLY' if args.apply else 'DRY RUN (pass --apply to write)'}")
    if done:
        print(f"resuming  : {len(done)} subject(s) already ingested per {state_path}")
    print()

    rows = list(csv.DictReader(args.manifest.open(), delimiter="\t"))
    print(f"indexing {args.d2c_root} …", flush=True)
    index = build_index(args.d2c_root, args.index_cache, args.refresh_index)
    gvcfs: dict[str, list[str]] = defaultdict(list)
    for p in index:
        if p.endswith("/gatk4/chrY.g.vcf.gz"):
            gvcfs[Path(p).relative_to(args.d2c_root).parts[0]].append(p)
    print(f"  {len(index)} files, {len(rows)} manifest rows\n")

    ext, donors, linked, with_data = read_workspace(args.db)
    items = plan(rows, args.d2c_root, index, gvcfs, ext, donors, linked, with_data, done)

    actions = Counter(i["action"] for i in items)
    todo = [i for i in items if i["action"] == "ingest"]
    print("── plan ─────────────────────────────────")
    print(f"  ingest                : {actions['ingest']}")
    print(f"    ├─ matched existing : {sum(1 for i in todo if i['match'] != 'new')}")
    print(f"    └─ new [lab]-[kit]  : {sum(1 for i in todo if i['match'] == 'new')}")
    print(f"  skip: already present : {actions['skip-present']}")
    print(f"  skip: superseded      : {actions['skip-superseded']}")
    print(f"  skip: no data on disk : {actions['skip-no-data']}")
    print(f"  skip: unsupported fmt : {actions['skip-unsupported']}")
    print(f"  skip: done (resume)   : {actions['skip-done']}")
    print(f"  ── by lab ──")
    for lab, n in Counter(i["row"]["lab"] for i in todo).most_common():
        print(f"    {lab:<8} {n}")
    print(f"  ── by ingest mode ──")
    for mode, n in Counter(f'{i["mode"]} ({i["how"]})' for i in todo).most_common():
        print(f"    {mode:<24} {n}")
    print()

    with args.report.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["action", "why", "lab", "kit", "subject", "target", "match", "mode", "path"])
        for i in items:
            r = i["row"]
            w.writerow([i["action"], i["why"], r["lab"], r["kit"], r["subject"], i["target"], i["match"], i["mode"], i["path"] or ""])
    print(f"wrote per-row plan to {args.report}\n")

    if args.lab:
        labs = {l.lower() for l in args.lab}
        todo = [i for i in todo if i["row"]["lab"].lower() in labs]
    if args.mode:
        todo = [i for i in todo if i["mode"] == args.mode]
    if args.match:
        todo = [i for i in todo if i["match"] == args.match]
    if args.lab or args.mode or args.match:
        picked = ", ".join(f"--{k} {v}" for k, v in (("lab", args.lab), ("mode", args.mode), ("match", args.match)) if v)
        print(f"filtered ({picked}): {len(todo)} of {actions['ingest']} sample(s)\n")
    if args.limit:
        todo = todo[: args.limit]
        print(f"--limit {args.limit}: ingesting only the first {len(todo)}\n")

    if not args.apply:
        print(f"dry run — re-run with --apply to ingest the {actions['ingest']} sample(s) above")
        return 0
    if not todo:
        print("nothing to ingest")
        return 0

    if not args.no_backup:
        backup_db(args.db)

    ok = failed = 0
    started = time.time()
    with state_path.open("a") as ledger:
        for n, item in enumerate(todo, 1):
            cmd = [str(navigator)] + ingest_args(item, args.project or None)
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode == 0:
                ok += 1
                ledger.write(json.dumps({"subject": item["row"]["subject"], "target": item["target"]}) + "\n")
                ledger.flush()
            else:
                failed += 1
                # `ingest` created the subject before it hit the bad file; don't leave it behind.
                pruned = item["match"] == "new" and prune_empty_subject(args.db, item["row"]["name"])
                print(
                    f"FAIL {item['target']}{' (empty subject pruned)' if pruned else ''}\n"
                    f"     {' '.join(cmd)}\n     {proc.stderr.strip()[:400]}",
                    file=sys.stderr,
                    flush=True,
                )

            if n % 25 == 0 or n == len(todo):
                rate = n / max(time.time() - started, 1e-9)
                eta = (len(todo) - n) / rate if rate else 0
                print(f"  {n}/{len(todo)}  ok={ok} failed={failed}  {rate*60:.0f}/min  eta {eta/60:.0f}m", flush=True)

    print(f"\n── done ─────────────────────────────────")
    print(f"  ingested : {ok}")
    print(f"  failed   : {failed}")
    print(f"  elapsed  : {(time.time()-started)/60:.1f}m")
    print(f"  ledger   : {state_path} (re-run to resume; delete to start over)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
