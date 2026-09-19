# =============================================================================
# levi -- 07. Starting from Bioconductor objects
#
# Goal:     feed levi from SummarizedExperiment, ExpressionSet and from the
#           differential expression tables of DESeq2, edgeR, limma and Seurat,
#           without going through an intermediate file.
#
# Note:     the blocks below are guarded by requireNamespace(), so the script
#           runs in full even without the optional packages installed.
#
# Runtime:  less than a minute.
# =============================================================================

library(levi)

ex      <- function(f) system.file("extdata", f, package = "levi")
network <- ex("hub_network.dat")


# -----------------------------------------------------------------------------
# 1. expressionInput accepts a data.frame
# -----------------------------------------------------------------------------
# This is the basis of everything that follows: there is no need to write a
# file to disk. All the leviFrom* adapters return a data.frame ready to go
# straight into expressionInput.

df <- data.frame(
    ID      = c("HUB", paste0("N", 1:8)),
    Test    = c(200, 200, 200, 200, 200, 5, 5, 5, 5),
    Control = c(10, 10, 10, 10, 10, 200, 200, 200, 200)
)

r <- levi(expressionInput         = df,
          networkCoordinatesInput = network,
          fileTypeInput           = "dat",
          geneSymbolInput         = "ID",
          readExpColumn           = readExpColumn("Test-Control"),
          resolutionValueInput    = 40,
          smoothValueInput        = 50)
cat("--- from a data.frame ---\n")
print(head(r$scores, 3), row.names = FALSE)


# -----------------------------------------------------------------------------
# 2. SummarizedExperiment
# -----------------------------------------------------------------------------
# leviFromSE() extracts an assay and computes the mean per group. There are
# two ways to indicate the groups: by label in colData, or by sample name.

if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {

    counts <- matrix(
        c(200, 200, 200, 200, 200,   5,   5,   5,   5,     # tumor1
          190, 210, 195, 205, 200,   6,   4,   5,   5,     # tumor2
           10,  10,  10,  10,  10, 200, 200, 200, 200),    # normal1
        nrow = 9,
        dimnames = list(c("HUB", paste0("N", 1:8)),
                        c("tumor1", "tumor2", "normal1")))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(counts = counts),
        colData = data.frame(condition = c("Tumor", "Tumor", "Normal"),
                             row.names = colnames(counts)))

    # form 1: by condition label
    expr_se <- leviFromSE(se,
                          assay_name    = "counts",
                          condition_col = "condition",
                          test_level    = "Tumor",
                          ctrl_level    = "Normal",
                          gene_col      = "ID")

    cat("\n--- leviFromSE ---\n")
    print(head(expr_se, 4), row.names = FALSE)

    r_se <- levi(expressionInput         = expr_se,
                 networkCoordinatesInput = network,
                 fileTypeInput           = "dat",
                 geneSymbolInput         = "ID",
                 readExpColumn           = readExpColumn("Test-Control"),
                 resolutionValueInput    = 40,
                 smoothValueInput        = 50)
    cat("range of the scores:",
        sprintf("%.3f to %.3f", min(r_se$scores$LandscapeScore),
                                max(r_se$scores$LandscapeScore)), "\n")

    # form 2: by sample name, when there is no condition column
    #   leviFromSE(se, test_col = c("tumor1","tumor2"), ctrl_col = "normal1")
    #
    # For raw counts, log_transform = TRUE applies log2(x + 1) before
    # aggregating -- and in that case use signal_mode = "logfc".

    # leviFromBioc() detects the class of the object and dispatches on its
    # own, which helps when the type comes from outside the function.
    expr_auto <- leviFromBioc(se, condition_col = "condition",
                              test_level = "Tumor", ctrl_level = "Normal",
                              gene_col = "ID")
    cat("leviFromBioc returns the same:", identical(expr_se, expr_auto), "\n")

} else {
    cat("\nSummarizedExperiment not installed; block skipped.\n")
}


# -----------------------------------------------------------------------------
# 3. ExpressionSet (microarray)
# -----------------------------------------------------------------------------
if (requireNamespace("Biobase", quietly = TRUE)) {

    # microarray intensities already come in log2
    mat <- matrix(
        c(9.5, 9.4, 9.6, 9.5, 9.5, 5.1, 5.0, 5.2, 5.1,
          9.6, 9.5, 9.5, 9.4, 9.6, 5.0, 5.1, 5.1, 5.0,
          5.2, 5.1, 5.0, 5.1, 5.2, 9.5, 9.4, 9.6, 9.5),
        nrow = 9,
        dimnames = list(c("HUB", paste0("N", 1:8)),
                        c("t1", "t2", "n1")))

    pd <- Biobase::AnnotatedDataFrame(
        data.frame(condition = c("Tumor", "Tumor", "Normal"),
                   row.names = colnames(mat)))
    eset <- Biobase::ExpressionSet(assayData = mat, phenoData = pd)

    expr_es <- leviFromExpressionSet(eset,
                                     condition_col = "condition",
                                     test_level    = "Tumor",
                                     ctrl_level    = "Normal",
                                     gene_col      = "ID")
    cat("\n--- leviFromExpressionSet ---\n")
    print(head(expr_es, 4), row.names = FALSE)

    # The values are in log2, so the correct mode is "logfc", not "ratio".
    r_es <- levi(expressionInput         = expr_es,
                 networkCoordinatesInput = network,
                 fileTypeInput           = "dat",
                 geneSymbolInput         = "ID",
                 readExpColumn           = readExpColumn("Test-Control"),
                 resolutionValueInput    = 40,
                 smoothValueInput        = 50,
                 signal_mode             = "logfc")
    cat("range of the scores (signal_mode = 'logfc'):",
        sprintf("%.3f to %.3f", min(r_es$scores$LandscapeScore),
                                max(r_es$scores$LandscapeScore)), "\n")

} else {
    cat("\nBiobase not installed; block skipped.\n")
}


# -----------------------------------------------------------------------------
# 4. Differential expression tables
# -----------------------------------------------------------------------------
# leviFromDESeq2(), leviFromEdgeR(), leviFromLimma() and leviFromSeurat()
# convert the output of those tools into a data.frame for levi.
#
# MIND the signal mode. These tables carry an already computed logFC beside a
# measure of mean abundance (baseMean, logCPM, AveExpr, pct.2). The two
# columns do NOT form a comparable test/control pair: baseMean is in the
# hundreds and log2FoldChange in the units, and a ratio between them has no
# biological meaning.
#
# The correct approach is to use only the logFC column, in single-column mode:
#
#     readExpColumn("log2FoldChange-log2FoldChange")  +  signal_mode = "logfc"

res_de <- data.frame(
    baseMean       = c(1200, 1100, 1300, 1150, 1250, 900, 950, 880, 910),
    log2FoldChange = c(4.3, 4.1, 4.4, 4.2, 4.3, -5.3, -5.1, -5.4, -5.2),
    row.names      = c("HUB", paste0("N", 1:8)))

expr_de <- leviFromDESeq2(res_de, gene_col = "ID")
cat("\n--- leviFromDESeq2 ---\n")
print(head(expr_de, 4), row.names = FALSE)

r_de <- levi(expressionInput         = expr_de,
             networkCoordinatesInput = network,
             fileTypeInput           = "dat",
             geneSymbolInput         = "ID",
             readExpColumn           = readExpColumn("log2FoldChange-log2FoldChange"),
             resolutionValueInput    = 40,
             smoothValueInput        = 50,
             signal_mode             = "logfc")

cat("scores from the logFC:\n")
print(r_de$scores[, c("Gene", "LandscapeScore", "Rank")], row.names = FALSE)

# The others follow the same pattern, changing only the column name:
#
#   leviFromEdgeR(fit)   -> readExpColumn("logFC-logFC")
#   leviFromLimma(fit)   -> readExpColumn("logFC-logFC")
#   leviFromSeurat(mk)   -> readExpColumn("Test-Test")
#
# All of them accept either the tool's object or a data.frame with the
# expected columns, which makes testing easier without installing the whole
# package.


# -----------------------------------------------------------------------------
# 5. Networks from STRING
# -----------------------------------------------------------------------------
# leviFromSTRING() removes the need to download a network file: it takes gene
# symbols, queries STRING, computes a layout with igraph and returns $nodes
# and $edges ready for fileTypeInput = "stg".
#
#   set.seed(42)   # force-directed layouts are stochastic
#   net <- leviFromSTRING(c("TP53","BRCA1","EGFR","MYC","PTEN"),
#                         species = 9606, score_threshold = 700)
#   levi(expressionInput          = expr_de,
#        networkCoordinatesInput  = net$nodes,
#        networkInteractionsInput = net$edges,
#        fileTypeInput            = "stg",
#        geneSymbolInput          = "ID",
#        readExpColumn            = readExpColumn("log2FoldChange-log2FoldChange"),
#        signal_mode              = "logfc")
#
# It requires STRINGdb and network access, which is why it stays commented
# out here.


# -----------------------------------------------------------------------------
# 6. Summary
# -----------------------------------------------------------------------------
#   - expressionInput accepts a data.frame, which removes the need for a
#     temporary file;
#   - leviFromSE / leviFromExpressionSet aggregate by group and return Test
#     and Control; leviFromBioc dispatches on the class of the object;
#   - DE tables carry a ready-made logFC: use single-column mode and
#     signal_mode "logfc", do not pair the logFC with the mean abundance
#     column;
#   - log2 data (microarray, VST/rlog, proteomics) call for "logfc"; counts
#     and TPM call for "ratio".
#
# Next: 08-visual-parameters.R tunes the appearance of the figure.
