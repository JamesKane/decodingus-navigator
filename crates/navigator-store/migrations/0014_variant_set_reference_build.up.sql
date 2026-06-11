-- Reference build the set's positions are on (e.g. 'hs1', 'GRCh38'), when known. Lets Y-SNP
-- panel placement read the build directly instead of re-deriving it from the subject's
-- alignments. NULL for sources of unknown build (generic VCF/CSV imports).
ALTER TABLE variant_set ADD COLUMN reference_build TEXT;
