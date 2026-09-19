# Calibration results

Generated on 2026-09-16 by inst/scripts/10-calibration.R with 3 workers,
300 replicates per null setting, in 83 minutes.

* type1_error.csv: family-wise error rate of each test under the global
  null (no effect anywhere), by network.
* null_pvalues_*.csv: smallest p-value per null replicate, for QQ plots.
* power.csv: probability of at least one significant result when the
  hub module of the network is shifted by `effect` in 4 vs 4 samples.
* layout_sensitivity*.csv: the same data on twenty Fruchterman-Reingold
  layouts of the 300-node graph; pairwise Jaccard of the genes in
  significant regions, versus the layout-free TFCE result.

levi 1.99.0, R 4.6.1.
