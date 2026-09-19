test_that("compiled null extraction agrees with full eight-connected regions", {
    set.seed(62)
    for (i in 1:20) {
        z <- matrix(sample(c(NA, .2, .4, .5, .6, .8), 80, TRUE), 8)
        for (min_cells in c(1, 3, 10)) {
            r <- levi:::.extractLandscapeRegions(z, -.1, .02, .1, min_cells)$summary
            expected <- vapply(c("over", "under"), function(d)
                max(c(0, r$Mass[r$Direction == d])), numeric(1))
            expect_equal(levi:::.maximumRegionMass(z, .02, .1, min_cells), expected)
        }
    }
})

test_that("regional mass rewards distributed support rather than only peaks", {
    z <- matrix(NA_real_, 8, 8); z[1:3, 1:3] <- .65; z[8, 8] <- .95
    r <- levi:::.extractLandscapeRegions(z, 0, 1, .1, 1)
    expect_equal(r$summary$Cells, c(9L, 1L))
    expect_equal(r$summary$Mass, c(.45, .35))
    expect_equal(r$summary$PeakScore, c(.65, .95))
})

test_that("regional correction uses both signs and counts ties", {
    z <- matrix(.5, 4, 4); z[1, 1] <- .8; z[4, 4] <- .2
    r <- levi:::.extractLandscapeRegions(z, 0, 1, .1, 1)
    null <- cbind(over = c(.1, .2, .05), under = c(.3, .1, .05))
    p <- levi:::.regionalPvalues(r, null)
    expect_equal(p$summary$PSpatial, c(.75, .75))
    expect_equal(p$null_max_mass$both, c(.3, .2, .05))
    r$summary$Mass[] <- 2
    expect_equal(levi:::.regionalPvalues(r, null)$summary$PSpatial, c(.25, .25))
})

test_that("permutations redetect areas and preserve missing positions", {
    draws <- 0L
    values <- cbind(c(-1, 1, NA), c(0, 0, NA)); edge <- matrix(c(1, 2), 1)
    z <- matrix(.5, 3, 3); z[1:2, 1:2] <- .7
    r <- levi:::.extractLandscapeRegions(z, 0, 1, .1, 1)
    testthat::local_mocked_bindings(landscape_gauss = function(..., SignalOut) {
        draws <<- draws + 1L
        expect_equal(SignalOut[3, 1], .5)
        expect_equal(SignalOut[4, 1], .5)
        x <- matrix(.5, 3, 3); x[3, 3] <- if (draws == 1L) .9 else .1
        list(m1 = x)
    }, .package = "levi")
    sig <- levi:::.networkSignals(values, edge, FALSE, "logfc", 1)
    p <- levi:::.permutationPvalues(matrix(0, 4, 2), sig$signal, sig$test,
        sig$control, z, 3, 0, 1, 1, .05, 2, node_values = values,
        edge_index = edge, signal_mode = "logfc", regions = r)
    expect_equal(p$null_max_mass$over, c(.3, 0))
    expect_equal(p$null_max_mass$under, c(0, .3))
    expect_equal(p$summary$PSpatial, 1/3)
})

test_that("regional API is reproducible and handles absence of effect", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    run <- function(n_perm, unit = "region") levi(
        expressionInput = data.frame(ID = c("HUB", paste0("N", 1:8)), Test = 1, Control = 1),
        networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
        fileTypeInput = "dat", geneSymbolInput = "ID", readExpColumn = readExpColumn("Test-Control"),
        resolutionValueInput = 1, n_perm = n_perm, inference_unit = unit)
    set.seed(3); a <- run(3)
    set.seed(3); b <- run(3)
    expect_equal(a$regions, b$regions)
    expect_equal(nrow(a$regions$summary), 0L)
    expect_equal(a$regions$null_max_mass$both, c(0, 0, 0))
    expect_null(a$pvalues); expect_null(a$raw_pvalues)
    expect_output(print(a), "regional permutation test")
    expect_null(run(0)$regions$inference)
    expect_true(is.list(run(3, "cell")$pvalues))
    expect_error(run(0, "gene"), "arg")
})

test_that("boundaries follow grid cells even on NA edges", {
    z <- matrix(NA_real_, 4, 4); z[1, 1] <- .8; z[4, 4] <- .2
    r <- levi:::.extractLandscapeRegions(z, 0, 1, .1, 1)
    b <- levi:::.regionBoundaries(r, 4, "over_01")
    expect_equal(nrow(b), 4L)
    expect_equal(range(c(b$x, b$xend)), c(.5, 1.5))
    expect_equal(range(c(b$y, b$yend)), c(3.5, 4.5))
    expect_equal(nrow(levi:::.regionBoundaries(r, 4, character())), 0L)
})

test_that("public regional results label areas and preserve supplied logFC", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    set.seed(12)
    x <- levi(expressionInput = data.frame(ID = c("HUB", paste0("N", 1:8)), logFC = .4),
        networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
        fileTypeInput = "dat", geneSymbolInput = "ID", readExpColumn = readExpColumn("logFC-logFC"),
        signal_mode = "logfc", inference_unit = "region", n_perm = 3,
        resolutionValueInput = 1, region_threshold = .02)
    expect_equal(range(x$landscape$z, na.rm = TRUE), rep(plogis(.4), 2))
    expect_true(nrow(x$regions$summary) > 0)
    expect_true(all(x$regions$summary$PSpatial == 1))
    expect_false(any(x$regions$summary$Significant))
    # No region passes sig_level, so none is labelled: a "p = 1.000" box on
    # every region is what the figure used to show.
    labels <- Filter(function(layer) is.data.frame(layer$data) &&
                      "Label" %in% names(layer$data), x$plot$layers)
    expect_length(labels, 0)
    expect_equal(x$metadata$inference_unit, "region")

    # Without a permutation test the regions are descriptive and keep names.
    y <- levi(expressionInput = data.frame(ID = c("HUB", paste0("N", 1:8)), logFC = .4),
        networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
        fileTypeInput = "dat", geneSymbolInput = "ID", readExpColumn = readExpColumn("logFC-logFC"),
        signal_mode = "logfc", inference_unit = "region", n_perm = 0,
        resolutionValueInput = 1, region_threshold = .02)
    labels <- Filter(function(layer) is.data.frame(layer$data) &&
                      "Label" %in% names(layer$data), y$plot$layers)
    expect_length(labels, 1)
    expect_true(all(grepl("^over_", labels[[1]]$data$Label)))
})
