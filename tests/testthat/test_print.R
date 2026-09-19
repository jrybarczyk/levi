library(levi)

# levi() used to return a bare list: the user got no way to inspect the result
# without knowing its six fields by heart.

hub_result <- function() {
    levi(networkCoordinatesInput = system.file("extdata", "hub_network.dat",
             package = "levi"),
         expressionInput = system.file("extdata", "hub_expression.dat",
             package = "levi"),
         fileTypeInput  = "dat",
         geneSymbolInput = "ID",
         readExpColumn  = readExpColumn("Test-Control"),
         resolutionValueInput = 10,
         smoothValueInput     = 5)
}

test_that("levi(): the result carries the levi_result class", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- hub_result()
    expect_s3_class(res, "levi_result")
    expect_true(is.list(res))

    # The class must not get in the way of the fields.
    expect_named(res, c("comparison", "landscape", "scores", "peaks", "regions",
                        "pvalues", "plot", "plot3d", "raw_pvalues", "metadata"))
    expect_s3_class(res$scores, "data.frame")
})

test_that("print(): summarises the result without redrawing the plot", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- hub_result()
    out <- paste(capture.output(print(res)), collapse = "\n")

    expect_match(out, "Test-Control", fixed = TRUE)
    expect_match(out, "nodes scored: 9", fixed = TRUE)
    expect_match(out, "0.5 = no change", fixed = TRUE)
    expect_match(out, "HUB", fixed = TRUE)
    expect_match(out, "regions:", fixed = TRUE)
    expect_match(out, "not run (n_perm = 0)", fixed = TRUE)
})

test_that("print(): honours n and returns the object invisibly", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- hub_result()

    one <- paste(capture.output(print(res, n = 1)), collapse = "\n")
    expect_match(one, "highest: HUB")
    # With n = 1 only a single gene is listed on the "highest" line.
    expect_false(grepl("highest: HUB \\(1.000\\), ", one))

    # capture.output keeps the summary out of the test log; the point here is
    # the return value, not the text.
    invisible(capture.output(vis <- withVisible(print(res))))
    expect_false(vis$visible)
    expect_identical(vis$value, res)
})

test_that("print(): reports the permutation test when it ran", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    set.seed(1)
    res <- levi(networkCoordinatesInput = system.file("extdata",
                    "bimodal_network.dat", package = "levi"),
                expressionInput = system.file("extdata",
                    "bimodal_expression.dat", package = "levi"),
                fileTypeInput  = "dat",
                geneSymbolInput = "ID",
                readExpColumn  = readExpColumn("Test-Control"),
                resolutionValueInput = 10,
                smoothValueInput     = 5,
                n_perm               = 19,
                inference_unit       = "cell")

    out <- paste(capture.output(print(res)), collapse = "\n")
    expect_match(out, "masked cells significant over")
    expect_false(grepl("not run", out, fixed = TRUE))
})

test_that("leviGrid(): still recognises a single classed result", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- hub_result()
    expect_silent(g <- leviGrid(res))
    expect_false(is.null(g))
})

# --- the 3D surface ----------------------------------------------------------

# plot3d had no test at all: it depends on an optional package and the figure
# used to be shown and discarded, so nothing could reach it.

hub_args <- list(
    networkCoordinatesInput = system.file("extdata", "hub_network.dat",
        package = "levi"),
    expressionInput = system.file("extdata", "hub_expression.dat",
        package = "levi"),
    fileTypeInput   = "dat",
    geneSymbolInput = "ID",
    readExpColumn   = readExpColumn("Test-Control"),
    resolutionValueInput = 10,
    smoothValueInput     = 5)

test_that("plot3d = FALSE leaves the field empty", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- do.call(levi, hub_args)
    expect_null(res$plot3d)
})

test_that("plot3d = TRUE returns the plotly surface", {
    skip_if_not_installed("plotly")
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- do.call(levi, c(hub_args, list(plot3d = TRUE)))

    expect_false(is.null(res$plot3d))
    expect_s3_class(res$plot3d, "plotly")

    # The surface carries the landscape as a z matrix, on the same 0-1 scale
    # as the 2D figure.
    built <- plotly::plotly_build(res$plot3d)
    z <- built$x$data[[1]]$z
    expect_true(is.matrix(z) || is.list(z))
    expect_equal(built$x$data[[1]]$type, "surface")
})

test_that("plot3d carries the significance boundary when the test ran", {
    skip_if_not_installed("plotly")
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    # The 2D figure has always drawn the p = sig_level contours. The 3D view
    # used to show the bare surface, dropping the only information about which
    # part of the landscape holds up against the null.
    set.seed(1)
    res <- levi(networkCoordinatesInput = system.file("extdata",
                    "medusa.dat", package = "levi"),
                expressionInput = system.file("extdata",
                    "expression.dat", package = "levi"),
                fileTypeInput   = "dat",
                geneSymbolInput = "ID",
                readExpColumn   = readExpColumn(
                    "TumorCurrentSmoker-NormalNeverSmoker"),
                resolutionValueInput = 20,
                smoothValueInput     = 50,
                n_perm               = 49,
                inference_unit       = "cell",
                p_adjust_method      = "none",
                plot3d               = TRUE)

    built <- plotly::plotly_build(res$plot3d)
    types <- vapply(built$x$data, function(d) d$type, character(1))

    expect_equal(sum(types == "surface"), 1L)
    expect_gt(sum(types == "scatter3d"), 0)

    names <- unlist(lapply(built$x$data, function(d) d$name))
    expect_true(any(grepl("p <=", names, fixed = TRUE)))
})

test_that("plot3d stays a bare surface when no permutation ran", {
    skip_if_not_installed("plotly")
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- do.call(levi, c(hub_args, list(plot3d = TRUE)))
    built <- plotly::plotly_build(res$plot3d)

    expect_equal(length(built$x$data), 1L)
    expect_equal(built$x$data[[1]]$type, "surface")
})
