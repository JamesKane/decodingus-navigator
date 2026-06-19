-- Fit residual for distance-minimizing ancestry models (G25_NMONTE): the Euclidean
-- distance between the sample and its fitted mixture in PC space. NULL for non-fit
-- methods (ADMIXTURE / AF_LIKELIHOOD / PCA_PROJECTION_GMM).
ALTER TABLE ancestry_result ADD COLUMN fit_distance REAL;
