# =============================================================================
# levi -- 04. Signal modes: from raw value to score
#
# Goal:     choose between signal_mode = "ratio", "logfc" and "zscore"
#           according to the nature of the data, and understand the effect of
#           logfc_k.
#
# Data:     logfc_network.dat / logfc_expression.dat -- six genes in a chain
#           with values ALREADY on the log2 scale (Control = 8 for all;
#           Test ranges from 3 to 10).
#
# Runtime:  less than a minute.
# =============================================================================

library(levi)

ex <- function(f) system.file("extdata", f, package = "levi")

run <- function(...) {
    levi(expressionInput         = ex("logfc_expression.dat"),
         networkCoordinatesInput = ex("logfc_network.dat"),
         fileTypeInput           = "dat",
         geneSymbolInput         = "ID",
         readExpColumn           = readExpColumn("Test-Control"),
         resolutionValueInput    = 40,
         smoothValueInput        = 50, ...)
}

sc <- function(r) {
    s <- setNames(round(r$scores$LandscapeScore, 3), r$scores$Gene)
    s[c("GA", "GB", "GC", "GD", "GE", "GF")]
}


# -----------------------------------------------------------------------------
# 1. The data
# -----------------------------------------------------------------------------
cat("--- expression (log2 scale) ---\n")
d <- read.delim(ex("logfc_expression.dat"))
d$logFC <- d$Test - d$Control        # in log2, the difference is the log-ratio
print(d, row.names = FALSE)

# GA to GE are below the control; only GF is above. In log2, the logFC ranges
# from -5 to +2.


# -----------------------------------------------------------------------------
# 2. The three modes
# -----------------------------------------------------------------------------
cat("\n--- ratio (default) ---\n"); print(sc(run(signal_mode = "ratio")))
cat("\n--- logfc ---\n");           print(sc(run(signal_mode = "logfc")))
cat("\n--- zscore ---\n");          print(sc(run(signal_mode = "zscore")))

# All three preserve the order, but distribute the values differently:
#
#   ratio   Test / (Test + Control), then rescaled. Assumes a LINEAR scale.
#           Applied to log2 values, as here, it treats 3 and 8 as if they were
#           counts, which compresses the real differences.
#
#   logfc   sigmoid 1 / (1 + exp(-k * logFC)). Assumes a LOG scale and works
#           with Test - Control. This is the correct mode for these data.
#
#   zscore  standardises the logFC and applies pnorm(). The score becomes the
#           gene's relative position in the distribution, not its absolute
#           magnitude.


# -----------------------------------------------------------------------------
# 3. How to choose
# -----------------------------------------------------------------------------
# The decisive question is: are the values on a linear or a logarithmic scale?
#
#   LINEAR SCALE  -> signal_mode = "ratio"
#     raw counts, TPM, FPKM, linear proteomics LFQ.
#
#   LOG SCALE     -> signal_mode = "logfc"
#     microarray RMA, VST/rlog from DESeq2, log2 proteomics, and any table
#     that already carries a ready-made logFC.
#
#   HETEROGENEOUS COMPARISONS -> signal_mode = "zscore"
#     several conditions on different scales, or when the gene's relative
#     position matters more than the size of the change.
#
# One pitfall: expressionLog = TRUE exists to convert log2 data back to the
# linear scale BEFORE the calculation in ratio mode. With signal_mode "logfc"
# or "zscore" it is ignored (levi emits a warning), because those modes need
# precisely the log scale.

cat("\n--- expressionLog with signal_mode = 'logfc' ---\n")
invisible(run(signal_mode = "logfc", expressionLog = TRUE))


# -----------------------------------------------------------------------------
# 4. logfc_k: the steepness of the sigmoid
# -----------------------------------------------------------------------------
# k controls how quickly the score leaves the neutral point as the logFC grows.

cat("\n--- effect of logfc_k ---\n")
for (k in c(0.5, 1, 2, 3)) {
    cat(sprintf("k = %-3g  ", k)); print(sc(run(signal_mode = "logfc", logfc_k = k)))
}

# A low k flattens the scale and keeps gradation between genes; a high k
# saturates, separating "changed" from "unchanged" almost like a threshold.
#
#   k = 0.5  large fold-changes (scRNA-seq, where |logFC| exceeds 5)
#   k = 1    default, suitable for most cases
#   k = 2-3  small fold-changes (microarray, |logFC| near 1)


# -----------------------------------------------------------------------------
# 5. Single-column input
# -----------------------------------------------------------------------------
# When the table already carries the computed logFC -- log2FoldChange from
# DESeq2, logFC from edgeR, avg_log2FC from Seurat -- there is no test/control
# pair to compare. In that case the column name is repeated on both sides, and
# levi switches to single-column mode automatically.

dfc <- data.frame(ID = d$ID, logFC = d$Test - d$Control)
cat("\n--- single logFC column ---\n")
print(dfc, row.names = FALSE)

r_single <- levi(
    expressionInput         = dfc,                # data.frame directly
    networkCoordinatesInput = ex("logfc_network.dat"),
    fileTypeInput           = "dat",
    geneSymbolInput         = "ID",
    readExpColumn           = readExpColumn("logFC-logFC"),   # same name twice
    resolutionValueInput    = 40,
    smoothValueInput        = 50,
    signal_mode             = "logfc"
)
cat("\nscores in single-column mode:\n")
print(sc(r_single))

# Note that expressionInput accepts a data.frame directly, with no file
# involved. That is how the leviFrom* adapters connect to levi -- see script
# 07.


# -----------------------------------------------------------------------------
# 6. Summary
# -----------------------------------------------------------------------------
#   - linear scale -> ratio; log scale -> logfc; mixed scales -> zscore;
#   - using ratio on log2 data compresses the differences and is the most
#     common mistake;
#   - expressionLog acts only in ratio mode;
#   - logfc_k tunes the contrast of the sigmoid, not the direction;
#   - ready-made logFC -> readExpColumn("column-column") + signal_mode
#     "logfc".
#
# Next: 05-batch-and-comparison.R generates several landscapes at once and
# compares two of them.
