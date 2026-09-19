# =============================================================================
# test_readExpColumn.R
# Tests for readExpColumn() — the comparison-specification helper.
# =============================================================================

library(levi)

# =============================================================================
# 1. Basic return type and structure
# =============================================================================

test_that("readExpColumn: returns a list", {
    expect_type(readExpColumn("Test-Control"), "list")
})

test_that("readExpColumn: single arg — list has 2 elements (call + string)", {
    res <- readExpColumn("Test-Control")
    # match.call() captures: [[1]] = readExpColumn (call), [[2]] = first arg
    expect_length(res, 2L)
})

test_that("readExpColumn: first element is the function call", {
    res <- readExpColumn("Test-Control")
    expect_equal(as.character(res[[1]]), "readExpColumn")
})

test_that("readExpColumn: second element matches the provided string", {
    res <- readExpColumn("TumorSmoker-NormalNeverSmoker")
    expect_equal(as.character(res[[2]]), "TumorSmoker-NormalNeverSmoker")
})

# =============================================================================
# 2. Multiple comparisons (batch mode)
# =============================================================================

test_that("readExpColumn: two args — list has 3 elements", {
    res <- readExpColumn("A-B", "A-C")
    expect_length(res, 3L)
})

test_that("readExpColumn: three args — list has 4 elements", {
    res <- readExpColumn("A-B", "A-C", "B-C")
    expect_length(res, 4L)
})

test_that("readExpColumn: args are captured in order", {
    res <- readExpColumn("Cond_A-Cond_B", "Cond_A-Cond_C")
    expect_equal(as.character(res[[2]]), "Cond_A-Cond_B")
    expect_equal(as.character(res[[3]]), "Cond_A-Cond_C")
})

# =============================================================================
# 3. Special patterns
# =============================================================================

test_that("readExpColumn: single-col logFC pattern accepted", {
    expect_no_error(readExpColumn("log2FoldChange-log2FoldChange"))
})

test_that("readExpColumn: self-comparison (null landscape) accepted", {
    expect_no_error(readExpColumn("Sample1-Sample1"))
})

test_that("readExpColumn: spaces in column names are accepted", {
    expect_no_error(readExpColumn("Gene expression-Control expression"))
})

# =============================================================================
# 4. Integration: readExpColumn drives correct column selection in levi()
# =============================================================================

test_that("readExpColumn: wrong column name causes levi() to stop", {
    hub_n <- system.file("extdata", "hub_network.dat",    package = "levi")
    hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")
    expect_error(
        levi(networkCoordinatesInput = hub_n,
             expressionInput         = hub_e,
             fileTypeInput           = "dat",
             geneSymbolInput         = "ID",
             readExpColumn           = readExpColumn("NONEXISTENT-Control"),
             contrastValueInput      = 50, resolutionValueInput = 10,
             zoomValueInput          = 50, smoothValueInput     = 5)
    )
})
