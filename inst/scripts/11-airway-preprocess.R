# =============================================================================
# levi -- 11. Preprocess the airway dataset for the validation vignette
#
# Goal:     run the heavy, network-dependent steps ONCE, offline, and save
#           small text files that the vignette can load in seconds. This
#           script is not run when the package is built.
#
# Steps:    DESeq2 on airway (dexamethasone vs untreated, blocked by cell
#           line), Ensembl -> symbol mapping, the 80 strongest responders,
#           their STRING network with a Kamada-Kawai layout, and the
#           per-sample log2 CPM of those genes for the sample-label tests.
#
# Output:   inst/extdata/airway/airway_dex_genes.tsv      (80 rows)
#           inst/extdata/airway/airway_dex_samples.tsv    (8 rows)
#           inst/extdata/airway/airway_string_nodes.tsv   (<= 80 rows)
#           inst/extdata/airway/airway_string_edges.tsv
#           inst/extdata/airway/README.md                 (provenance)
#
# Requires: airway, DESeq2, org.Hs.eg.db, AnnotationDbi, STRINGdb, edgeR
# Runtime:  two to three minutes, plus the STRING download on first use.
# =============================================================================

suppressPackageStartupMessages({
    library(airway); library(DESeq2); library(levi)
})
options(timeout = 900)  # the STRING alias file is 21 MB
stopifnot(requireNamespace("org.Hs.eg.db"), requireNamespace("STRINGdb"),
          requireNamespace("edgeR"))

out_dir <- file.path("inst", "extdata", "airway")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- 1. Differential expression --------------------------------------------
data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("dex", "trt", "untrt"))

symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(res), column = "SYMBOL", keytype = "ENSEMBL",
    multiVals = "first")
tab <- data.frame(Ensembl = rownames(res), Symbol = unname(symbols),
    baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue,
    padj = res$padj, stringsAsFactors = FALSE)
tab <- tab[!is.na(tab$Symbol) & !is.na(tab$padj), ]
tab <- tab[!duplicated(tab$Symbol), ]

# The 80 genes with the strongest evidence of a dexamethasone response.
top <- head(tab[order(tab$padj, -abs(tab$stat)), ], 80)

# -- 2. STRING network with a fixed layout ---------------------------------
set.seed(42)
net <- leviFromSTRING(top$Symbol, species = 9606, score_threshold = 400,
                      layout = "kk")
nodes <- net$nodes
edges <- net$edges
top <- top[top$Symbol %in% nodes$name, ]

# -- 3. Per-sample log2 CPM for the sample-label tests ---------------------
cpm <- edgeR::cpm(counts(dds)[top$Ensembl, ], log = TRUE, prior.count = 1)
rownames(cpm) <- top$Symbol
colnames(cpm) <- colnames(dds)
samples <- data.frame(Sample = colnames(dds), cell = colData(dds)$cell,
                      dex = colData(dds)$dex, stringsAsFactors = FALSE)

# -- 4. Write ---------------------------------------------------------------
write.table(top, file.path(out_dir, "airway_dex_genes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(round(cpm, 4), file.path(out_dir, "airway_dex_logcpm.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(samples, file.path(out_dir, "airway_dex_samples.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(nodes, file.path(out_dir, "airway_string_nodes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(edges, file.path(out_dir, "airway_string_edges.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

writeLines(c(
    "# airway preprocessed files",
    "",
    sprintf("Generated on %s by inst/scripts/11-airway-preprocess.R.",
            format(Sys.Date())),
    "",
    "Source: Himes et al. (2014), airway smooth muscle cells treated with",
    "dexamethasone, Bioconductor package `airway`. DESeq2 model ~ cell + dex;",
    "contrast trt vs untrt; genes with >= 10 counts in >= 4 samples.",
    "",
    sprintf("* airway_dex_genes.tsv: the %d strongest responders (by adjusted",
            nrow(top)),
    "  p-value, then |Wald statistic|) that STRING recognised, with DESeq2",
    "  statistics.",
    "* airway_dex_logcpm.tsv: log2 CPM (edgeR, prior count 1) of those genes",
    "  in the 8 samples.",
    "* airway_dex_samples.tsv: sample, cell line and treatment.",
    sprintf("* airway_string_nodes.tsv / airway_string_edges.tsv: STRING v11.5"),
    "  interactions (combined score >= 400) among those genes, with a",
    "  Kamada-Kawai layout (set.seed(42)).",
    "",
    sprintf("Versions: DESeq2 %s, STRINGdb %s, org.Hs.eg.db %s, levi %s.",
            packageVersion("DESeq2"), packageVersion("STRINGdb"),
            packageVersion("org.Hs.eg.db"), packageVersion("levi"))),
    file.path(out_dir, "README.md"))

cat("Wrote", nrow(top), "genes,", nrow(nodes), "nodes,", nrow(edges),
    "edges to", out_dir, "\n")
