# Migrations — do not edit an applied file, *including its comments*

`sqlx` records a **checksum of each migration file's full contents** in the `_sqlx_migrations`
table when it applies it. On the next open it re-checksums the files on disk and refuses to proceed
if one has changed:

```
error: could not open workspace ~/.decodingus/navigator-rs.db:
       migration failed: migration 33 was previously applied but has been modified
```

That is a **hard failure of the whole app** for anyone with an existing workspace — not a warning,
and not something the user can clear without deleting their database.

The checksum covers the entire file, so a comment-only edit breaks it exactly as thoroughly as a
schema change. This has already happened once: a repo-wide documentation-path rewrite
(`docs/` → `documents/`) updated a `-- See …` line inside `0033_retire_pca_ancient_ancestry.up.sql`,
which silently invalidated every existing database until it was restored byte-for-byte.

**Consequences for routine maintenance:**

- Never run a bulk find-and-replace across `crates/navigator-store/migrations/`. Exclude this
  directory from repo-wide rewrites (path renames, formatting sweeps, link fixers, licence headers).
- `0033_retire_pca_ancient_ancestry.up.sql` deliberately still says `docs/design/…`, a path that no
  longer exists. **Leave it.** The design docs moved to `documents/design/`, but correcting the
  comment would break every deployed database to fix a stale pointer in a file nobody reads at
  runtime. The same applies to any future stale reference in an applied migration.
- To change behaviour, add a **new** migration. Applied files are history, not source.
