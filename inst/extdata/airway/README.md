# airway preprocessed files

Generated on 2026-09-13 by inst/scripts/11-airway-preprocess.R.

Source: Himes et al. (2014), airway smooth muscle cells treated with
dexamethasone, Bioconductor package `airway`. DESeq2 model ~ cell + dex;
contrast trt vs untrt; genes with >= 10 counts in >= 4 samples.

* airway_dex_genes.tsv: the 78 strongest responders (by adjusted
  p-value, then |Wald statistic|) that STRING recognised, with DESeq2
  statistics.
* airway_dex_logcpm.tsv: log2 CPM (edgeR, prior count 1) of those genes
  in the 8 samples.
* airway_dex_samples.tsv: sample, cell line and treatment.
* airway_string_nodes.tsv / airway_string_edges.tsv: STRING v11.5
  interactions (combined score >= 400) among those genes, with a
  Kamada-Kawai layout (set.seed(42)).

Versions: DESeq2 1.52.0, STRINGdb 2.24.0, org.Hs.eg.db 3.23.1, levi 3.11.1.
