-- Per-call evidence from the source VCF, and a schema tag saying whether a set has it.
--
-- The VCF import previously kept only contig/position/ref/alt/rsID/genotype and threw QUAL, FILTER
-- and the whole FORMAT column away. That left every downstream consumer with nothing to gate on —
-- a private-Y engine over imported VCFs could not distinguish a 40x hom-alt call from a 2-read
-- artefact, which is the difference between a real new branch and mapping noise.
--
-- Existing rows keep NULL evidence, which is honest: the data was never captured, and NULL must not
-- be read as zero. `call_schema` on the set records which import wrote it, so a consumer can
-- *require* evidence (schema 2) rather than silently accept sets that can never satisfy a quality
-- gate. Re-importing a set upgrades it.

ALTER TABLE variant_call ADD COLUMN qual   REAL;
ALTER TABLE variant_call ADD COLUMN filter TEXT;     -- NULL when PASS/absent; only failures stored
ALTER TABLE variant_call ADD COLUMN dp     INTEGER;
ALTER TABLE variant_call ADD COLUMN gq     INTEGER;
ALTER TABLE variant_call ADD COLUMN ad_ref INTEGER;
ALTER TABLE variant_call ADD COLUMN ad_alt INTEGER;  -- depth of the *called* ALT, not simply ALT[0]

-- 1 = basic (pre-evidence). New VCF imports write 2; see variants::CALL_SCHEMA_*.
ALTER TABLE variant_set ADD COLUMN call_schema INTEGER NOT NULL DEFAULT 1;
