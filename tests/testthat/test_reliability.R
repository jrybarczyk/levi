# Scientific contracts reproduced from the technical review.
review_run <- function(values, control = values, mode = "logfc", log = FALSE,
                       single = TRUE, ...) {
    levi(expressionInput = data.frame(ID = c("HUB", paste0("N", 1:8)),
                                     Test = values, Control = control),
         networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
         fileTypeInput = "dat", geneSymbolInput = "ID",
         readExpColumn = readExpColumn(if (single) "Test-Test" else "Test-Control"),
         resolutionValueInput = 10, smoothValueInput = 5,
         signal_mode = mode, expressionLog = log, ...)
}

test_that("ratio preserves the absolute neutral point and is context invariant", {
    signal <- function(x) as.numeric(levi:::.computeSignalOut(
        matrix(x, ncol = 1), matrix(1, length(x), 1)))
    expect_equal(signal(c(1, 2, 3)), c(0.5, 2/3, 3/4))
    expect_equal(signal(c(1, 2, 3, 100))[1:3], signal(1:3))
})

test_that("one-sided logFC agrees across surface, node scores, peaks and 3D", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    for (direction in c(-1, 1)) {
        res <- review_run(direction * seq_len(9))
        expect_true(all(direction * (res$landscape$z - 0.5) > 0, na.rm = TRUE))
        expect_true(all(direction * (res$scores$LandscapeScore - 0.5) > 0))
        if (nrow(res$peaks)) {
            n <- sqrt(nrow(res$landscape))
            surface <- matrix(res$landscape$z, n)
            expect_equal(res$peaks$Score, round(surface[cbind(
                res$peaks$MatrixRow, n + 1L - res$peaks$MatrixCol)], 4))
        }
    }
    skip_if_not_installed("plotly")
    res <- review_run(1:9, plot3d = TRUE)
    built <- plotly::plotly_build(res$plot3d)$x$data
    surf <- Filter(function(x) identical(x$type, "surface"), built)[[1]]
    expect_equal(as.numeric(surf$z), res$landscape$z)
})

test_that("expressionLog is ignored identically by scripts and GUI outside ratio", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    app <- system.file("shiny", "app.R", package = "levi")
    skip_if(!nzchar(app))
    env <- new.env()
    for (expr in parse(app)) {
        if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
            identical(expr[[2]], as.name("levi_shiny"))) eval(expr, env)
    }
    net <- levi:::.parseNetwork(system.file("extdata", "hub_network.dat",
                                           package = "levi"), NA, "dat")
    for (mode in c("logfc", "zscore")) for (single in c(FALSE, TRUE)) {
        args <- list(values = 8:16, control = rep(7, 9), mode = mode, single = single)
        a <- do.call(review_run, args)
        b <- suppressMessages(do.call(review_run, c(args, list(log = TRUE))))
        expect_equal(a$landscape, b$landscape)
        expect_equal(a$scores, b$scores)
        g <- suppressMessages(env$levi_shiny(
            expression = data.frame(ID = c("HUB", paste0("N", 1:8)),
                                    Test = 8:16, Control = 7),
            fileType = "dat", networkCoord = net$edges, networkInterac = net$nodes,
            geneSymbol = "ID", baseTest = "Test",
            baseControl = if (single) "Test" else "Control",
            alphaValue = 50, betaValue = 50, backValue = 10, smoothValue = 5,
            expressionLog = TRUE, signal_mode = mode))
        expect_equal(as.numeric(g$surface), a$landscape$z)
        expect_equal(g$scores, a$scores)
    }
    linear <- review_run(2^(1:9), control = rep(2^5, 9), mode = "ratio", single = FALSE)
    logged <- review_run(1:9, control = rep(5, 9), mode = "ratio", single = FALSE, log = TRUE)
    expect_equal(linear$landscape, logged$landscape)
})

test_that("uniform absence of effect stays neutral under permutation", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    for (mode in c("ratio", "logfc", "zscore")) {
        res <- review_run(rep(7, 9), mode = mode, single = FALSE, n_perm = 9,
                          inference_unit = "cell")
        expect_true(all(abs(res$landscape$z - 0.5) < 1e-12, na.rm = TRUE))
        expect_true(all(res$raw_pvalues$over == 1, na.rm = TRUE))
        expect_true(all(res$raw_pvalues$under == 1, na.rm = TRUE))
    }
})

test_that("node-label permutations recalculate edge signals and preserve missing positions", {
    captured <- list()
    testthat::local_mocked_bindings(landscape_gauss = function(..., SignalOut) {
        captured[[length(captured) + 1L]] <<- SignalOut
        list(m1 = matrix(0.5, 2, 2))
    }, .package = "levi")
    values <- cbind(c(1, 3, 7, NA), c(1, 1, 1, NA))
    edge <- matrix(c(1, 2), nrow = 1)
    sig <- matrix(rep(0.5, 5), ncol = 1)
    set.seed(42)
    levi:::.permutationPvalues(matrix(0, 5, 2), sig, sig, sig,
        matrix(0.5, 2, 2), 2, 0, 1, 1, 0.05, 12,
        node_values = values, edge_index = edge, signal_mode = "ratio")
    for (signal in captured) {
        # Recover the linear abundance from ratio with Control = 1.
        node_abundance <- signal[1:3] / (1 - signal[1:3])
        expect_equal(sort(node_abundance), c(1, 3, 7))
        expect_equal(signal[4], 0.5)
        midpoint <- mean(node_abundance[1:2])
        expect_equal(signal[5], midpoint / (midpoint + 1))
    }
})

test_that("both directional cell families share one adjustment and keep the mask", {
    p <- list(over = matrix(c(.001, .2, NA, .4), 2),
              under = matrix(c(.99, .8, NA, .6), 2))
    adjusted <- levi:::.adjustLandscapePvalues(p, "BY")
    raw <- unlist(p, use.names = FALSE)
    actual <- unlist(adjusted, use.names = FALSE)
    expect_identical(is.na(actual), is.na(raw))
    expect_equal(actual[!is.na(raw)], p.adjust(raw[!is.na(raw)], "BY"))
})

test_that("regions retain broad coherent components instead of only peaks", {
    z <- matrix(0.5, 6, 6)
    z[2:3, 2:3] <- c(0.7, 0.8, 0.75, 0.9)
    z[5:6, 5:6] <- c(0.3, 0.25, 0.2, 0.35)
    z[1, 6] <- 0.95 # singleton must be removed by min_cells
    regions <- levi:::.extractLandscapeRegions(z, zoomValue = 0,
        increase = 1, threshold = 0.1, min_cells = 3)
    expect_equal(nrow(regions$summary), 2)
    expect_setequal(regions$summary$Direction, c("over", "under"))
    expect_true(all(regions$summary$Cells == 4))
    expect_equal(nrow(regions$cells), 8)
    expect_true(all(regions$cells$Score != 0.95))
})

test_that("levi returns regions and records their fixed definition", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    res <- review_run(1:9, region_threshold = 0.05, region_min_cells = 1)
    expect_named(res$regions, c("summary", "cells", "threshold", "min_cells"))
    expect_equal(res$regions$threshold, 0.05)
    expect_equal(res$metadata$region_definition,
                 "8-connected grid cells beyond neutral score")
    expect_output(print(res), "regions:")
})

test_that("provenance prevents equal-sized incompatible differences", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    set.seed(53)
    res <- review_run(1:9, n_perm = 3, inference_unit = "cell")
    expect_identical(res$metadata$p_adjust_method, "BY")
    assign(".Random.seed", res$metadata$rng_state, envir = .GlobalEnv)
    again <- review_run(1:9, n_perm = 3, inference_unit = "cell")
    expect_equal(res$raw_pvalues, again$raw_pvalues)
    for (field in c("signal_mode", "nodes", "grid", "missing_genes")) {
        other <- res; other$metadata[[field]] <- "incompatible"
        expect_error(leviDiff(res, other), field)
    }
    other <- res; other$metadata <- NULL
    expect_error(leviDiff(res, other), "metadata")
    other <- res; other$landscape$Var1[1] <- -100
    expect_error(leviDiff(res, other), "grid coordinates")
})

test_that("edgeR adapters accept real exact and GLM test objects", {
    skip_if_not_installed("edgeR")
    set.seed(4)
    y <- edgeR::DGEList(matrix(rnbinom(40 * 6, mu = 50, size = 10), 40, 6),
                        group = rep(c("A", "B"), each = 3))
    rownames(y) <- paste0("gene", seq_len(40))
    exact <- edgeR::exactTest(y, dispersion = 0.1)
    design <- model.matrix(~ y$samples$group)
    fit <- edgeR::glmFit(y, design, dispersion = 0.1)
    lrt <- edgeR::glmLRT(fit, coef = 2)
    for (obj in list(exact, lrt)) {
        tt <- edgeR::topTags(obj, n = Inf)$table
        out <- leviFromEdgeR(obj)
        expect_equal(out$GeneID, rownames(tt))
        expect_equal(out$logFC, tt$logFC)
        expect_error(leviFromEdgeR(obj, coef = 2), "Select coef")
    }
    expect_error(leviFromEdgeR(fit), "test result")
})

test_that("each batch comparison records its own RNG starting state", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    expr <- data.frame(ID = c("HUB", paste0("N", 1:8)), Test = 1:9, Control = 1)
    args <- list(expressionInput = expr, fileTypeInput = "dat",
        networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
        geneSymbolInput = "ID", resolutionValueInput = 10, smoothValueInput = 5,
        n_perm = 3, signal_mode = "logfc", inference_unit = "cell")
    set.seed(321)
    batch <- do.call(levi, c(args, list(readExpColumn =
        readExpColumn("Test-Control", "Test-Control"))))
    expect_false(identical(batch[[1]]$metadata$rng_state, batch[[2]]$metadata$rng_state))
    assign(".Random.seed", batch[[2]]$metadata$rng_state, envir = .GlobalEnv)
    replay <- do.call(levi, c(args, list(readExpColumn = readExpColumn("Test-Control"))))
    expect_equal(replay$raw_pvalues, batch[[2]]$raw_pvalues)
})
