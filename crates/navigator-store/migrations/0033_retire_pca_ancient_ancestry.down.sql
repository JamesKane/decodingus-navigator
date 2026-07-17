-- Irreversible by design: the deleted rows were fabricated output from estimators that no longer
-- exist in the codebase, so there is nothing to restore and no code left that could reproduce them.
-- Re-running the consensus ancestry estimate repopulates deep ancestry with `ANCIENT_ADMIXTURE`.
SELECT 1;
