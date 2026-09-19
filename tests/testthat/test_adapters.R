library(levi)

# Tests for the Bioconductor entry points. They build the objects on the fly
# instead of downloading anything, so the whole file runs offline and without
# the optional packages: what needs DESeq2, edgeR, limma or Seurat is exercised
# through the data.frame form those adapters also accept.

genes <- c("HUB", paste0("N", 1:8))

hub_net <- system.file("extdata", "hub_network.dat", package = "levi")

# A SummarizedExperiment with two tumour samples and one normal, matching the
# nodes of hub_network.dat.
make_se <- function() {
    counts <- matrix(
        c(200, 200, 200, 200, 200,   5,   5,   5,   5,
          190, 210, 195, 205, 200,   6,   4,   5,   5,
           10,  10,  10,  10,  10, 200, 200, 200, 200),
        nrow = 9,
        dimnames = list(genes, c("tumor1", "tumor2", "normal1")))

    SummarizedExperiment::SummarizedExperiment(
        assays  = list(counts = counts),
        colData = data.frame(condition = c("Tumor", "Tumor", "Normal"),
                             row.names = colnames(counts)))
}

# levi() on the hub network, at the smallest resolution that still works.
run_levi <- function(expr, comparison = "Test-Control", ...) {
    levi(expressionInput         = expr,
         networkCoordinatesInput = hub_net,
         fileTypeInput           = "dat",
         geneSymbolInput         = "ID",
         readExpColumn           = readExpColumn(comparison),
         resolutionValueInput    = 10,
         smoothValueInput        = 5, ...)
}


# --- leviFromSE --------------------------------------------------------------

test_that("leviFromSE: aggregates by condition label", {
    out <- leviFromSE(make_se(), assay_name = "counts",
                      condition_col = "condition",
                      test_level = "Tumor", ctrl_level = "Normal",
                      gene_col = "ID")

    expect_s3_class(out, "data.frame")
    expect_identical(colnames(out), c("ID", "Test", "Control"))
    expect_identical(out$ID, genes)

    # Test is the mean of the two tumour samples, Control the single normal.
    expect_equal(out$Test[out$ID == "HUB"], mean(c(200, 190)))
    expect_equal(out$Control[out$ID == "HUB"], 10)
})

test_that("leviFromSE: selects samples by name and by index", {
    se <- make_se()

    by_name <- leviFromSE(se, test_col = c("tumor1", "tumor2"),
                          ctrl_col = "normal1", gene_col = "ID")
    by_index <- leviFromSE(se, test_col = 1:2, ctrl_col = 3L, gene_col = "ID")

    expect_equal(by_name, by_index)
    expect_equal(by_name$Control[by_name$ID == "N5"], 200)
})

test_that("leviFromSE: log_transform applies log2(x + 1)", {
    out <- leviFromSE(make_se(), condition_col = "condition",
                      test_level = "Tumor", ctrl_level = "Normal",
                      gene_col = "ID", log_transform = TRUE)

    expect_equal(out$Test[out$ID == "HUB"], mean(log2(c(200, 190) + 1)))
})

test_that("leviFromSE: gene_col renames the identifier column", {
    out <- leviFromSE(make_se(), condition_col = "condition",
                      test_level = "Tumor", ctrl_level = "Normal",
                      gene_col = "Symbol")
    expect_identical(colnames(out)[1], "Symbol")
})

test_that("leviFromSE: falls back to the first assay with a message", {
    expect_message(
        out <- leviFromSE(make_se(), assay_name = "vst",
                          condition_col = "condition",
                          test_level = "Tumor", ctrl_level = "Normal",
                          gene_col = "ID"),
        "not found")
    expect_identical(colnames(out), c("ID", "Test", "Control"))
})

test_that("leviFromSE: rejects incomplete or unknown group specifications", {
    se <- make_se()

    expect_error(leviFromSE(se, gene_col = "ID"),
                 "Provide either")
    expect_error(leviFromSE(se, condition_col = "group",
                            test_level = "Tumor", ctrl_level = "Normal"),
                 "not found in colData")
    expect_error(leviFromSE(se, condition_col = "condition"),
                 "test_level and ctrl_level")
    expect_error(leviFromSE(se, condition_col = "condition",
                            test_level = "Treated", ctrl_level = "Normal"),
                 "Treated")
    expect_error(leviFromSE(se, test_col = "tumor9", ctrl_col = "normal1"),
                 "tumor9")
})


# --- leviFromBioc ------------------------------------------------------------

test_that("leviFromBioc: dispatches SummarizedExperiment to leviFromSE", {
    se <- make_se()
    args <- list(condition_col = "condition", test_level = "Tumor",
                 ctrl_level = "Normal", gene_col = "ID")

    expect_identical(do.call(leviFromBioc, c(list(se), args)),
                     do.call(leviFromSE,   c(list(se), args)))
})

test_that("leviFromBioc: names the unsupported class in the error", {
    expect_error(leviFromBioc(data.frame(a = 1)), "Unsupported class")
    expect_error(leviFromBioc(1:10), "integer")
})


# --- leviFromExpressionSet ---------------------------------------------------

test_that("leviFromExpressionSet: aggregates a microarray ExpressionSet", {
    skip_if_not_installed("Biobase")

    mat <- matrix(
        c(9.5, 9.4, 9.6, 9.5, 9.5, 5.1, 5.0, 5.2, 5.1,
          9.6, 9.5, 9.5, 9.4, 9.6, 5.0, 5.1, 5.1, 5.0,
          5.2, 5.1, 5.0, 5.1, 5.2, 9.5, 9.4, 9.6, 9.5),
        nrow = 9, dimnames = list(genes, c("t1", "t2", "n1")))

    pd <- Biobase::AnnotatedDataFrame(
        data.frame(condition = c("Tumor", "Tumor", "Normal"),
                   row.names = colnames(mat)))
    eset <- Biobase::ExpressionSet(assayData = mat, phenoData = pd)

    out <- leviFromExpressionSet(eset, condition_col = "condition",
                                 test_level = "Tumor", ctrl_level = "Normal",
                                 gene_col = "ID")

    expect_identical(colnames(out), c("ID", "Test", "Control"))
    expect_identical(out$ID, genes)
    expect_equal(out$Test[out$ID == "HUB"], mean(c(9.5, 9.6)))

    # leviFromBioc must reach the same adapter for this class.
    expect_identical(
        leviFromBioc(eset, condition_col = "condition", test_level = "Tumor",
                     ctrl_level = "Normal", gene_col = "ID"),
        out)
})


# --- differential expression tables ------------------------------------------

test_that("leviFromDESeq2: keeps baseMean and log2FoldChange, drops NA", {
    res <- data.frame(
        baseMean       = c(1200, 1100, NA),
        log2FoldChange = c(4.3, -5.1, 2.0),
        padj           = c(0.001, 0.002, 0.5),
        row.names      = c("HUB", "N5", "N6"))

    out <- leviFromDESeq2(res, gene_col = "ID")

    expect_identical(colnames(out), c("ID", "baseMean", "log2FoldChange"))
    expect_identical(out$ID, c("HUB", "N5"))   # the NA row is dropped
    expect_equal(out$log2FoldChange, c(4.3, -5.1))
})

test_that("leviFromDESeq2: requires the DESeq2 result columns", {
    expect_error(leviFromDESeq2(data.frame(logFC = 1, row.names = "HUB")),
                 "baseMean")
})

test_that("leviFromEdgeR: accepts a topTags-shaped data.frame", {
    tt <- data.frame(logFC = c(4.3, -5.1), logCPM = c(10.2, 9.8),
                     FDR = c(0.001, 0.002), row.names = c("HUB", "N5"))

    out <- leviFromEdgeR(tt, gene_col = "ID")

    expect_identical(colnames(out), c("ID", "logCPM", "logFC"))
    expect_equal(out$logFC, c(4.3, -5.1))
})

test_that("leviFromLimma: accepts a topTable-shaped data.frame", {
    tt <- data.frame(logFC = c(4.3, -5.1), AveExpr = c(8.1, 7.7),
                     adj.P.Val = c(0.001, 0.002), row.names = c("HUB", "N5"))

    out <- leviFromLimma(tt, gene_col = "ID")

    expect_identical(colnames(out), c("ID", "AveExpr", "logFC"))
    expect_equal(out$AveExpr, c(8.1, 7.7))
})

test_that("leviFromSeurat: maps pct.2 to Control and avg_log2FC to Test", {
    mk <- data.frame(avg_log2FC = c(2.5, -1.8), pct.1 = c(0.9, 0.2),
                     pct.2 = c(0.3, 0.8), row.names = c("HUB", "N5"))

    out <- leviFromSeurat(mk, gene_col = "ID")

    expect_identical(colnames(out), c("ID", "Control", "Test"))
    expect_equal(out$Test, c(2.5, -1.8))     # avg_log2FC
    expect_equal(out$Control, c(0.3, 0.8))   # pct.2
})

test_that("leviFromSeurat: reports every missing column at once", {
    expect_error(leviFromSeurat(data.frame(x = 1, row.names = "HUB")),
                 "avg_log2FC")
    expect_error(leviFromSeurat(data.frame(avg_log2FC = 1, row.names = "HUB")),
                 "pct.2")
})


# --- adapter output feeds levi() ---------------------------------------------

test_that("leviFromSE output is accepted by levi()", {
    expr <- leviFromSE(make_se(), condition_col = "condition",
                       test_level = "Tumor", ctrl_level = "Normal",
                       gene_col = "ID")

    pdf(NULL); on.exit(dev.off(), add = TRUE)
    res <- run_levi(expr)

    expect_named(res, c("comparison", "landscape", "scores", "peaks", "regions",
                        "pvalues", "plot", "plot3d", "raw_pvalues", "metadata"))
    expect_setequal(res$scores$Gene, genes)
    expect_true(all(res$scores$LandscapeScore >= 0 &
                    res$scores$LandscapeScore <= 1))

    # The core is over-expressed and the corners repressed in this dataset.
    sc <- setNames(res$scores$LandscapeScore, res$scores$Gene)
    expect_gt(sc[["HUB"]], 0.5)
    expect_lt(sc[["N5"]], 0.5)
})

test_that("leviFromDESeq2 output feeds levi() in single-column mode", {
    res_de <- data.frame(
        baseMean       = rep(1000, 9),
        log2FoldChange = c(4.3, 4.1, 4.4, 4.2, 4.3, -5.3, -5.1, -5.4, -5.2),
        row.names      = genes)

    expr <- leviFromDESeq2(res_de, gene_col = "ID")

    pdf(NULL); on.exit(dev.off(), add = TRUE)
    res <- run_levi(expr, comparison = "log2FoldChange-log2FoldChange",
                    signal_mode = "logfc")

    sc <- setNames(res$scores$LandscapeScore, res$scores$Gene)
    expect_gt(sc[["HUB"]], 0.5)
    expect_lt(sc[["N5"]], 0.5)
})


# --- leviDiff and leviGrid ---------------------------------------------------

test_that("leviDiff: a result against itself is flat zero", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    res <- run_levi(system.file("extdata", "hub_expression.dat",
                                package = "levi"))

    d <- leviDiff(res, res)

    expect_named(d, c("comparison", "diff", "plot"))
    expect_identical(colnames(d$diff), c("X", "Y", "Diff"))
    expect_true(all(d$diff$Diff == 0, na.rm = TRUE))
})

test_that("leviDiff: rejects landscapes of different sizes", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    small <- run_levi(expr)
    big <- levi(expressionInput         = expr,
                networkCoordinatesInput = hub_net,
                fileTypeInput           = "dat",
                geneSymbolInput         = "ID",
                readExpColumn           = readExpColumn("Test-Control"),
                resolutionValueInput    = 20,
                smoothValueInput        = 5)

    expect_error(leviDiff(small, big), "different sizes")
})

test_that("leviDiff: requires the landscape field", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    res <- run_levi(system.file("extdata", "hub_expression.dat",
                                package = "levi"))
    stripped <- res
    stripped$landscape <- NULL

    expect_error(leviDiff(stripped, res), "landscape")
})

test_that("leviGrid: accepts a single result and a list of results", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    res <- run_levi(system.file("extdata", "hub_expression.dat",
                                package = "levi"))

    expect_silent(g1 <- leviGrid(res))
    expect_false(is.null(g1))

    g2 <- leviGrid(list(res, res), ncol = 2,
                   titles = c("first", "second"))
    expect_false(is.null(g2))
})


# --- leviEnrich and leviFromSTRING -------------------------------------------

test_that("leviEnrich: refuses a result with no scores", {
    skip_if_not_installed("clusterProfiler")
    expect_error(leviEnrich(list(plot = NULL)), "result\\$scores not found")
})

test_that("leviFromSTRING: validates arguments before touching the network", {
    skip_if_not_installed("STRINGdb")

    # A data.frame of genes needs id_col to say which column holds them.
    expect_error(leviFromSTRING(data.frame(gene = "TP53")), "id_col")

    # match.arg rejects unknown layouts and network types offline.
    expect_error(leviFromSTRING("TP53", layout = "spiral"), "'arg'")
    expect_error(leviFromSTRING("TP53", network_type = "genetic"), "'arg'")
})
