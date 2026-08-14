-- Provenance for a realigned alignment: which alignment it was made from, and how.
--
-- Realignment (design/realignment-module.md) re-maps a vendor GRCh37/38 alignment's reads to
-- CHM13v2 and registers the result as a *new* alignment row under the same sequence_run — the same
-- physical library, only the mapping changed. Until now alignments were independent rows keyed
-- only to a run, so nothing recorded that one was derived from another.
--
-- This matters beyond bookkeeping. A subject can end up with two alignments of one library on two
-- builds, and every consumer that picks "the" alignment for an analysis needs to be able to tell
-- which is the vendor's original and which Navigator produced, and to say so in the UI. It is also
-- what makes the realigned row safe to delete and rebuild without touching the source.
--
-- Both columns are NULL for every existing row, which is correct: they were not derived from
-- anything. `derived_from_alignment_id` NULL means "this is an original".
ALTER TABLE alignment ADD COLUMN derived_from_alignment_id INTEGER REFERENCES alignment(id);

-- How the derivation was done, as `realign:<backend>-<preset>` (e.g. `realign:minimap2-sr`).
-- Free text rather than an enum because the backend and preset are recorded together and both may
-- gain values; `Alignment.aligner` continues to carry the mapper alone.
ALTER TABLE alignment ADD COLUMN derivation TEXT;

-- Finding a source alignment's derivatives is the common query — the UI asks "has this been
-- realigned already?" before offering the action, and the delete path asks before removing a row.
CREATE INDEX IF NOT EXISTS ix_alignment_derived_from ON alignment (derived_from_alignment_id);
