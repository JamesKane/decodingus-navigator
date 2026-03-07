-- ============================================
-- Enrich Haplogroup Reconciliation with Atmosphere Lexicon fields
-- Adds: heteroplasmy observations, identity verification,
--        manual override, and audit log
-- ============================================

ALTER TABLE haplogroup_reconciliation ADD COLUMN heteroplasmy_observations JSON DEFAULT '[]';
ALTER TABLE haplogroup_reconciliation ADD COLUMN identity_verification JSON;
ALTER TABLE haplogroup_reconciliation ADD COLUMN manual_override JSON;
ALTER TABLE haplogroup_reconciliation ADD COLUMN audit_log JSON DEFAULT '[]';
