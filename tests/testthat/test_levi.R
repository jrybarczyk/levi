library(levi)

hub_n <- system.file("extdata", "hub_network.dat",    package = "levi")
hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")

# One shared levi() call — reused by all tests below (resolution=10, smooth=5)
.res <- levi(
    networkCoordinatesInput = hub_n,
    expressionInput         = hub_e,
    fileTypeInput           = "dat",
    geneSymbolInput         = "ID",
    readExpColumn           = readExpColumn("Test-Control"),
    contrastValueInput      = 50,
    resolutionValueInput    = 10,
    zoomValueInput          = 50,
    smoothValueInput        = 5
)

test_that("levi() returns a list with scores, landscape and comparison", {
    expect_type(.res, "list")
    expect_true(all(c("scores", "landscape", "comparison") %in% names(.res)))
})

test_that("scores is a data.frame with LandscapeScore in [0, 1]", {
    expect_s3_class(.res$scores, "data.frame")
    sc <- .res$scores$LandscapeScore
    expect_true(all(sc >= 0 & sc <= 1))
})

test_that("comparison field matches readExpColumn argument", {
    expect_equal(.res$comparison, "Test-Control")
})

test_that("landscape is a data.frame with numeric columns", {
    expect_s3_class(.res$landscape, "data.frame")
    expect_true(nrow(.res$landscape) > 0)
})

# --- nodes without expression value -----------------------------------------

# Until 2.0.0 this detection read fixed column indices of the edge table,
# which had drifted onto coordinate columns: the check never fired, so a total
# identifier mismatch produced a flat neutral landscape in silence.

test_that("levi(): reports how many nodes have no expression value", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expect_message(
        levi(networkCoordinatesInput = system.file("extdata",
                 "sparse_network.dat", package = "levi"),
             expressionInput = system.file("extdata",
                 "sparse_expression.dat", package = "levi"),
             fileTypeInput  = "dat",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "10 nodes without expression value")
})

test_that("levi(): writes the unmatched node names to the log it names", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    suppressMessages(
        levi(networkCoordinatesInput = system.file("extdata",
                 "sparse_network.dat", package = "levi"),
             expressionInput = system.file("extdata",
                 "sparse_expression.dat", package = "levi"),
             fileTypeInput  = "dat",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5))

    logPath <- file.path(tempdir(), "Test-Control", "levi.log")
    expect_true(file.exists(logPath))
    expect_setequal(readLines(logPath),
                    c("G2", "G3", "G4", "G6", "G7", "G9",
                      "G10", "G11", "G13", "G14"))
})

test_that("levi(): warns when no identifier matches at all", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    # Ensembl-style ids against a network named with symbols: the landscape
    # comes back uniformly neutral, which looks exactly like "nothing changed".
    expect_warning(
        res <- levi(networkCoordinatesInput = system.file("extdata",
                        "hub_network.dat", package = "levi"),
                    expressionInput = data.frame(
                        ID      = c("ENSG001", "ENSG002"),
                        Test    = c(200, 5),
                        Control = c(10, 200)),
                    fileTypeInput  = "dat",
                    geneSymbolInput = "ID",
                    readExpColumn  = readExpColumn("Test-Control"),
                    resolutionValueInput = 10,
                    smoothValueInput     = 5),
        "None of the 9 network nodes matched")

    expect_true(all(res$scores$LandscapeScore == 0.5))
})

test_that("levi(): stays quiet when every node has a value", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expect_no_message(
        expect_no_warning(
            levi(networkCoordinatesInput = system.file("extdata",
                     "hub_network.dat", package = "levi"),
                 expressionInput = system.file("extdata",
                     "hub_expression.dat", package = "levi"),
                 fileTypeInput  = "dat",
                 geneSymbolInput = "ID",
                 readExpColumn  = readExpColumn("Test-Control"),
                 resolutionValueInput = 10,
                 smoothValueInput     = 5)))
})

# --- early validation of the inputs ------------------------------------------

# Both checks used to surface as internal errors: indexing a missing column
# gave "undefined columns selected", and a format mismatch escaped as a dplyr
# error about seq(delimiter_edge + 1, ...). Neither named the argument to fix.

test_that("levi(): names the bad column and lists the available ones", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    err <- expect_error(
        levi(networkCoordinatesInput = system.file("extdata",
                 "hub_network.dat", package = "levi"),
             expressionInput = system.file("extdata",
                 "hub_expression.dat", package = "levi"),
             fileTypeInput  = "dat",
             geneSymbolInput = "GENE",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "geneSymbolInput")

    expect_match(conditionMessage(err), "GENE", fixed = TRUE)
    expect_match(conditionMessage(err), "ID, Test, Control", fixed = TRUE)
})

test_that("levi(): reports a network format mismatch as such", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    hub_expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    # A Medusa file declared as Pajek.
    expect_error(
        levi(networkCoordinatesInput = system.file("extdata",
                 "hub_network.dat", package = "levi"),
             expressionInput = hub_expr,
             fileTypeInput  = "net",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "no '\\*Edges' section")

    # A Pajek file declared as Medusa.
    pajek <- tempfile(fileext = ".net")
    writeLines(c("*Vertices 2", "1 \"A\" 0.1 0.2", "2 \"B\" 0.3 0.4",
                 "*Edges", "1 2 1"), pajek)

    expect_error(
        levi(networkCoordinatesInput = pajek,
             expressionInput = hub_expr,
             fileTypeInput  = "dat",
             geneSymbolInput = "ID",
             readExpColumn  = readExpColumn("Test-Control"),
             resolutionValueInput = 10,
             smoothValueInput     = 5),
        "no '\\*nodes' section")
})
