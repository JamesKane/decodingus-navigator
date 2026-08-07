# Ancestral Origin Records

The genealogical context behind a tree placement: where a lineage's most distant known ancestor
came from, and roughly when. The AppView mirrors these into `fed.ancestral_origin` and renders them
as the genealogical-era origins icicle (`proposals/ancestral-origin-icicle.md` in the AppView repo).

**Status:** implemented — publisher `navigator publish-origins`, consumer `du-jobs` jetstream.

---

## Why this record is different

Every other record in this namespace is anonymized analysis. This one is **genealogy about a named
family**, and it is the only record type that carries a surname. It exists because the AppView had
no locality data at all for the genealogical era: of 9,642 placed Y samples only 1,380 carry a donor
coordinate, and the 7,882 direct-to-consumer tips — the entire genealogical era — have **three**
between them.

The posture that makes it publishable is the project's existing one
(`proposals/biosample-identifier-dedup.md`): an MDKA is the earliest *documented* ancestor on a
line, long dead, and their surname, parish and dates are genealogical context rather than
living-donor PII. What is new is that the rule is now **enforced in code at both ends** rather than
asserted in prose.

### What may never appear in this record

- **A given name.** The record has a `surname` field and no `ancestorName` field, and
  `navigator_domain::identity::surname_of` reduces the stored full name before serialization. The
  full name never leaves the workspace.
- **Anything about the living tester.** No donor identifier, no kit number on any public surface, no
  notes field.
- **An ancestor who might still be alive.** `birthYear` must be ≤ 1900.

---

## Ancestral Origin Record

**NSID:** `com.decodingus.atmosphere.ancestralOrigin`
**Key:** `origin-<lineage>-<biosample uuid>` — deterministic, so correcting an MDKA and
re-publishing overwrites that ancestor's record via `putRecord` rather than accumulating duplicates.

```json
{
  "$type": "com.decodingus.atmosphere.ancestralOrigin",
  "biosampleRef": "at://did:plc:.../com.decodingus.atmosphere.biosample/bio-...",
  "externalIds": [{ "namespace": "FTDNA", "value": "B5163" }],
  "lineage": "Y_DNA",
  "surname": "Kane",
  "originPlace": "Creegh South, Co. Clare, Ireland",
  "originCountry": "Ireland",
  "birthYear": 1830,
  "deathYear": 1908,
  "lat": "52.75",
  "lon": "-9.43",
  "createdAt": "2026-08-07T00:00:00Z"
}
```

| field | notes |
|---|---|
| `biosampleRef` | at-uri of the parent biosample record, when the sample is federated. |
| `externalIds` | **the join that actually resolves.** No placed sample on the AppView's tree carries an at-uri — the tips were bulk-loaded — while 7,548 carry an `FTDNA` row in `core.biosample_identifier`. Vendor namespaces stay off every public surface. |
| `lineage` | `Y_DNA` \| `MT_DNA`. `Auto` is not published: it has no tree to hang from. |
| `surname` | family name only. Particles are preserved (`van der Berg`, `Ó Súilleabháin`). |
| `originPlace` | **as recorded.** Normalization is the AppView's job (`du_db::place`) — one implementation, fixable without a client release, re-runnable over records already ingested. |
| `lat` / `lon` | strings, because DAG-CBOR has no float type. Coarsened to 2dp (~1 km). |

---

## The gates

Applied at the edge by `AncestralOriginRecord::build`, and applied **again independently** by the
AppView on ingest — a buggy or hostile client must not be able to widen them. A record failing any
of them is refused outright rather than stored and hidden.

1. **Surname only.** The AppView rejects a `surname` that is anything but one name token plus
   particles, so a given name cannot arrive through a field labelled `surname`.
2. **`birthYear` ≤ 1900**, and ≥ 1000 — a year outside that is corrupt, not merely old. This is the
   check that makes "not living-donor PII" verifiable rather than asserted.
3. **No birth year → country only.** Nothing then establishes that the ancestor is long dead, so
   `originPlace` and the coordinate are withheld. The record is still published; the detail is not.
4. **Coordinates coarsened to 2dp**, at publish *and* re-applied at ingest. A rooftop coordinate
   plus a surname narrows to one family; a county-scale view cannot use the precision anyway.
5. **A record with no join key is refused** — it could never be rendered, and could still be read.

---

## Who may publish

Publishing is not offered for every MDKA the workspace holds. The predicate
(`navigator_store::mdka::publishable`) requires **both**:

- **the workspace holds the subject's primary data** — a `variant_set`, or a `sequence_run` with an
  `alignment`. A roster row alone is someone else's kit that this workspace merely knows of;
  publishing its genealogy would be republishing what that tester gave a vendor, not us.
- **and the tester has not opted out.** `ftdna_member.publicly_shares` is the member's own FTDNA
  sharing setting, imported with the roster. A `0` is an explicit "do not show me" and overrides
  everything else. No roster row means there is no opt-out to honour.

Measured on the reference workspace: 583 Y MDKA rows sit on subjects with primary data; **558**
publish and **25 never do**.

---

## Publishing

```
navigator publish-origins --lineage y          # dry run: reports, sends nothing
navigator publish-origins --lineage y --apply  # queue for the next sync
```

Dry-run by default — genealogy is not something to publish by accident. Records are enqueued to the
sync outbox rather than posted directly, so the batch retries, survives being offline, and is
resumable: the deterministic rkey maps onto `putRecord`, making a second run a no-op on the AppView.

Reference workspace, Y lineage: 558 considered → 556 pass the field gates (2 refused), of which 207
carry a place and 306 are country-only.

---

## Still open

- **No GUI trigger.** Same gap as `private-y --project` and `rebuild-signatures --stale-tree`
  (`design/project-block-tree.md` §11.4) and it wants the same answer, not a bespoke button.
- **Retraction is unspecified.** Withdrawing consent needs the PDS record deleted *and* the AppView
  mirror tombstoned. The jetstream `delete` path handles the mechanism; nothing yet handles the
  workflow.
- **The importer leaves HTML entities** in `mdka.ancestor_name` (`mac L&#225;ire`). `surname_of`
  decodes them at the boundary so nothing mojibake is published, but the root cause is in
  `navigator-domain::ftdna` and 61 of 6,218 names are affected.
