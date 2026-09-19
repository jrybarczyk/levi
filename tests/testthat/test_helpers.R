# =============================================================================
# test_helpers.R
# Unit tests for internal helper functions:
#   .computeSignalOut(), .rescale_for_cpp(), %||%
# These are accessed via levi::: which is acceptable in tests.
# =============================================================================

library(levi)

# ── Helper constructor ────────────────────────────────────────────────────────
mat <- function(x) matrix(x, ncol = 1)

# =============================================================================
# .computeSignalOut — ratio mode
# =============================================================================

test_that(".computeSignalOut ratio: output is in [0, 1]", {
    exp_ <- mat(c(100, 50, 200, 10))
    ctrl <- mat(c(10,  50, 100, 200))
    out  <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "ratio")
    expect_true(all(out >= 0 & out <= 1))
})

test_that(".computeSignalOut ratio: higher test → higher score", {
    exp_  <- mat(c(200, 100, 10))
    ctrl  <- mat(c(10,  10,  10))
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "ratio")
    expect_true(all(diff(as.numeric(out)) <= 0))  # scores should be decreasing
})

test_that(".computeSignalOut ratio: equal test and ctrl → score near 0.5", {
    exp_  <- mat(c(50, 50, 50))
    ctrl  <- mat(c(50, 50, 50))
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "ratio")
    expect_true(all(abs(out - 0.5) < 1e-9))
})

test_that(".computeSignalOut ratio: output stays within [0,1]", {
    exp_  <- mat(c(200, 100, 10))
    ctrl  <- mat(c(10,  10,  200))
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "ratio")
    expect_true(max(out) <= 1 + 1e-9)
    expect_true(min(out) >= 0 - 1e-9)
})

# =============================================================================
# .computeSignalOut — logfc mode (dual column)
# =============================================================================

test_that(".computeSignalOut logfc: neutral logFC=0 → score near 0.5", {
    exp_  <- mat(c(5, 5, 5))
    ctrl  <- mat(c(5, 5, 5))
    out   <- levi:::.computeSignalOut(exp_, ctrl,
                                       expressionLog = FALSE,
                                       signal_mode   = "logfc",
                                       logfc_k       = 1,
                                       single_col    = FALSE)
    expect_true(all(abs(out - 0.5) < 1e-9))
})

test_that(".computeSignalOut logfc: positive logFC → score > 0.5", {
    exp_  <- mat(c(8, 8, 8))  # log2 test
    ctrl  <- mat(c(5, 5, 5))  # log2 ctrl   → logFC = +3
    out   <- levi:::.computeSignalOut(exp_, ctrl,
                                       signal_mode = "logfc",
                                       logfc_k     = 1,
                                       single_col  = FALSE)
    expect_true(all(out > 0.5))
})

test_that(".computeSignalOut logfc: negative logFC → score < 0.5", {
    exp_  <- mat(c(5, 5, 5))   # log2 test
    ctrl  <- mat(c(8, 8, 8))   # log2 ctrl   → logFC = -3
    out   <- levi:::.computeSignalOut(exp_, ctrl,
                                       signal_mode = "logfc",
                                       logfc_k     = 1,
                                       single_col  = FALSE)
    expect_true(all(out < 0.5))
})

test_that(".computeSignalOut logfc: higher k → steeper response", {
    exp_  <- mat(c(7))
    ctrl  <- mat(c(5))
    out1  <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "logfc",
                                       logfc_k = 0.5, single_col = FALSE)
    out2  <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "logfc",
                                       logfc_k = 2.0, single_col = FALSE)
    # Higher k → score further from 0.5
    expect_gt(out2 - 0.5, out1 - 0.5)
})

test_that(".computeSignalOut logfc: single_col mode treats signalExp as logFC", {
    # signalExp IS the logFC → value 3 should map to sigmoid(3) ≈ 0.953
    exp_  <- mat(c(3))
    ctrl  <- mat(c(1))  # ignored in single_col mode
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "logfc",
                                       logfc_k = 1, single_col = TRUE)
    expected <- 1 / (1 + exp(-1 * 3))
    expect_equal(as.numeric(out), expected, tolerance = 1e-9)
})

test_that(".computeSignalOut logfc: expressionLog=TRUE leaves log-scale inputs unchanged", {
    exp_  <- mat(8)
    ctrl  <- mat(5)
    out   <- levi:::.computeSignalOut(exp_, ctrl, expressionLog = TRUE,
                                       signal_mode = "logfc",
                                       logfc_k = 1, single_col = FALSE)
    expected <- 1 / (1 + exp(-1 * 3))
    expect_equal(as.numeric(out), expected, tolerance = 1e-9)
})

# =============================================================================
# .computeSignalOut — zscore mode
# =============================================================================

test_that(".computeSignalOut zscore: output is in [0, 1]", {
    exp_  <- mat(c(8, 6, 7, 5, 6))
    ctrl  <- mat(c(5, 5, 5, 5, 5))
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "zscore",
                                       single_col = FALSE)
    expect_true(all(out >= 0 & out <= 1))
})

test_that(".computeSignalOut zscore: median logFC gene has score near 0.5", {
    # Symmetric distribution → median → z=0 → pnorm(0)=0.5
    lfc   <- c(-3, -1, 0, 1, 3)
    exp_  <- mat(5 + lfc)
    ctrl  <- mat(rep(5, 5))
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "zscore",
                                       single_col = FALSE)
    # The middle gene (lfc=0) should be near 0.5
    expect_true(abs(out[3] - 0.5) < 0.05)
})

test_that(".computeSignalOut zscore: single constant input → all scores ~0.5", {
    exp_  <- mat(c(5, 5, 5, 5))
    ctrl  <- mat(c(5, 5, 5, 5))
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "zscore",
                                       single_col = FALSE)
    # sd ≈ 0 → z = 0 for all → pnorm(0) = 0.5
    expect_true(all(abs(out - 0.5) < 0.01))
})

test_that(".computeSignalOut zscore: higher logFC → higher score than lower logFC", {
    exp_  <- mat(c(8, 5))
    ctrl  <- mat(c(5, 5))
    out   <- levi:::.computeSignalOut(exp_, ctrl, signal_mode = "zscore",
                                       single_col = FALSE)
    expect_gt(out[1], out[2])
})

# =============================================================================
# .rescale_for_cpp
# =============================================================================

test_that(".rescale_for_cpp: output is in [0.05, 0.95]", {
    x   <- mat(c(-5, -1, 0, 2, 8))
    out <- levi:::.rescale_for_cpp(x)
    expect_true(min(out) >= 0.05 - 1e-9)
    expect_true(max(out) <= 0.95 + 1e-9)
})

test_that(".rescale_for_cpp: min maps to 0.05, max maps to 0.95", {
    x   <- mat(c(0, 5, 10))
    out <- levi:::.rescale_for_cpp(x)
    expect_equal(min(out), 0.05, tolerance = 1e-9)
    expect_equal(max(out), 0.95, tolerance = 1e-9)
})

test_that(".rescale_for_cpp: constant input returns matrix of 0.5", {
    x   <- mat(c(3, 3, 3, 3))
    out <- levi:::.rescale_for_cpp(x)
    expect_true(all(out == 0.5))
})

test_that(".rescale_for_cpp: output preserves row count", {
    x   <- mat(seq_len(20))
    out <- levi:::.rescale_for_cpp(x)
    expect_equal(nrow(out), 20L)
    expect_equal(ncol(out), 1L)
})

test_that(".rescale_for_cpp: monotonic input produces monotonic output", {
    x   <- mat(1:10)
    out <- levi:::.rescale_for_cpp(x)
    expect_true(all(diff(as.numeric(out)) > 0))
})

# =============================================================================
# %||% null-coalescing operator
# =============================================================================

test_that("%||%: returns left when left is non-NULL", {
    expect_equal("left" %||% "right", "left")
})

test_that("%||%: returns right when left is NULL", {
    expect_equal(NULL %||% "right", "right")
})

test_that("%||%: chain: first non-NULL wins", {
    expect_equal(NULL %||% NULL %||% "found", "found")
})

test_that("%||%: numeric types are preserved", {
    expect_equal(42L %||% 99L, 42L)
    expect_equal(NULL %||% 99L, 99L)
})

# =============================================================================
# landscape_gauss / nearest_node_grid  (nucleo C++)
# =============================================================================

# support grid: 4 points, one per quadrant
.pts <- matrix(c(0.25, 0.75, 0.25, 0.75,
                 0.25, 0.25, 0.75, 0.75), ncol = 2)
.lg <- function(sig, sigma = 2, occ = 0.05) {
    m <- matrix(sig, ncol = 1)
    levi:::landscape_gauss(.pts, m, m, m,
                           resolutionValue = 20, zoomValue = 0,
                           increase = 1/19, sigma = sigma, occFrac = occ)
}

test_that("landscape_gauss: returns m1, m2, m3 and occ with the shape of the grid", {
    L <- .lg(c(0.1, 0.4, 0.6, 0.9))
    expect_true(all(c("m1", "m2", "m3", "occ") %in% names(L)))
    expect_equal(dim(L$m1), c(20L, 20L))
    expect_equal(dim(L$occ), c(20L, 20L))
})

test_that("landscape_gauss: the result stays within the range of the signals", {
    sig <- c(0.1, 0.4, 0.6, 0.9)
    v <- as.numeric(.lg(sig)$m1)
    v <- v[!is.na(v)]
    expect_true(all(v >= min(sig) - 1e-9))
    expect_true(all(v <= max(sig) + 1e-9))
})

test_that("landscape_gauss: a constant signal returns that value across the network", {
    # central property of the normalised convolution: dividing by sum(w) makes
    # the result a convex average, so data with no variation gives the neutral
    # value.
    for (sg in c(1, 2, 5)) {
        v <- as.numeric(.lg(rep(0.5, 4), sigma = sg)$m1)
        v <- v[!is.na(v)]
        expect_true(length(v) > 0)
        expect_true(all(abs(v - 0.5) < 1e-9))
    }
})

test_that("landscape_gauss: the scale does not depend on sigma", {
    sig <- c(0, 0, 1, 1)
    mx <- vapply(c(1, 2, 4, 8), function(sg) {
        max(.lg(sig, sigma = sg)$m1, na.rm = TRUE)
    }, numeric(1))
    # with a wide kernel the maximum drops a little, but does not collapse by 1/k
    expect_true(all(mx > 0.5))
})

test_that("landscape_gauss: outside the silhouette the value is NA, not zero", {
    L <- .lg(c(0.5, 0.5, 0.5, 0.5), sigma = 1, occ = 0.2)
    expect_true(anyNA(L$m1))
    expect_true(all(L$occ >= 0))
})

test_that("landscape_gauss: higher value near the point with the strongest signal", {
    L <- .lg(c(0, 0, 0, 1))          # only the fourth point, at (0.75, 0.75)
    hi <- L$m1[16, 16]               # cell near (0.75, 0.75)
    lo <- L$m1[5, 5]                 # cell near (0.25, 0.25)
    expect_true(!is.na(hi) && !is.na(lo))
    expect_gt(hi, lo)
})

test_that("landscape_gauss: rejects inconsistent inputs", {
    m <- matrix(c(0.5, 0.5), ncol = 1)
    expect_error(levi:::landscape_gauss(.pts, m, m, m, 20, 0, 1/19, 2, 0.05))
    m4 <- matrix(rep(0.5, 4), ncol = 1)
    expect_error(levi:::landscape_gauss(.pts, m4, m4, m4, 0, 0, 1/19, 2, 0.05))
})

test_that("nearest_node_grid: labels each cell with the nearest node", {
    g <- levi:::nearest_node_grid(.pts, resolutionValue = 20,
                                  zoomValue = 0, increase = 1/19)
    expect_equal(dim(g), c(20L, 20L))
    expect_true(all(g >= 1 & g <= 4))
    expect_setequal(unique(as.integer(g)), 1:4)
    # each quadrant belongs to its own point
    expect_equal(g[3, 3],  1L)
    expect_equal(g[18, 3], 2L)
    expect_equal(g[3, 18], 3L)
    expect_equal(g[18, 18], 4L)
})

test_that(".labelPermutations counts large unblocked designs without enumerating", {
    # 58 vs 49 unblocked: choose(107, 49) ~ 1e31 arrangements. Enumerating
    # them with combn() overflowed the integer range and aborted; the count
    # must come from choose() and the scheme must fall back to Monte Carlo.
    groups <- rep(c("Tumor", "Normal"), c(58, 49))
    set.seed(1)
    perms <- levi:::.labelPermutations(groups, n_perm = 9L)
    expect_false(perms$exact)
    expect_equal(perms$possible, choose(107, 49))
    expect_length(perms$labels, 9L)
    expect_true(all(vapply(perms$labels, function(x)
        sum(x == "Tumor") == 58, logical(1))))
    # The exact scheme is unchanged for a small design.
    ex <- levi:::.labelPermutations(rep(c("a", "b"), each = 3), n_perm = 9L)
    expect_true(ex$exact)
    expect_equal(ex$possible, 20)
    expect_length(ex$labels, 19L)
})

# -----------------------------------------------------------------------------
# edge_weighting: support weights and their effect on the landscape
# -----------------------------------------------------------------------------
test_that(".supportWeights: degree weights sum to one around every node", {
    # star: node 1 is the hub, 2..9 are leaves
    ei <- cbind(1L, 2:9)
    w <- levi:::.supportWeights(9L, ei, "degree")
    expect_length(w, 9 + 8)
    expect_equal(w[1:9], rep(1, 9))
    expect_equal(sum(w[-(1:9)]), 8 * (1 / 8 + 1) / 2)   # 4.5: each leaf-side 1/2, hub-side 1/16
    expect_equal(levi:::.supportWeights(9L, ei, "midpoint"), rep(1, 17))
    expect_equal(levi:::.supportWeights(9L, ei, "none"), c(rep(1, 9), rep(0, 8)))
})

test_that("landscape_gauss: a weight of zero removes the point, ones reproduce the default", {
    pts <- matrix(c(.2, .2, .8, .8, .5, .5), ncol = 2, byrow = TRUE)
    m <- function(v) matrix(v, ncol = 1)
    args <- list(resolutionValue = 30L, zoomValue = 0, increase = 1 / 30,
                 sigma = 2, occFrac = 0.01)
    full <- do.call(levi:::landscape_gauss, c(list(pts, m(c(0, 1, 1)), m(c(0, 1, 1)), m(c(0, 1, 1))), args))
    ones <- do.call(levi:::landscape_gauss, c(list(pts, m(c(0, 1, 1)), m(c(0, 1, 1)), m(c(0, 1, 1))), args, list(weights = rep(1, 3))))
    expect_equal(ones$m1, full$m1)
    drop3 <- do.call(levi:::landscape_gauss, c(list(pts, m(c(0, 1, 1)), m(c(0, 1, 1)), m(c(0, 1, 1))), args, list(weights = c(1, 1, 0))))
    two <- do.call(levi:::landscape_gauss, c(list(pts[1:2, ], m(c(0, 1)), m(c(0, 1)), m(c(0, 1))), args))
    expect_equal(drop3$m1, two$m1)
    expect_error(do.call(levi:::landscape_gauss, c(list(pts, m(c(0, 1, 1)), m(c(0, 1, 1)), m(c(0, 1, 1))), args, list(weights = 1:2))),
                 "one value per point")
})

test_that("edge_weighting: the hub emphasis follows the option", {
    # Hub in the centre, leaves on a circle: the midpoints sit halfway.
    ang <- seq(0, 2 * pi, length.out = 9)[-9]
    net <- list(nodes = data.frame(name = c("HUB", paste0("L", 1:8)),
                                   x = c(50, 50 + 40 * cos(ang)),
                                   y = c(50, 50 + 40 * sin(ang))),
                edges = data.frame(V1 = "HUB", V2 = paste0("L", 1:8)))
    expr <- data.frame(ID = c("HUB", paste0("L", 1:8)), logFC = c(3, rep(0, 8)))
    run <- function(mode) suppressMessages(levi(expressionInput = expr,
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        fileTypeInput = "stg", geneSymbolInput = "ID",
        readExpColumn = readExpColumn("logFC-logFC"), signal_mode = "logfc",
        resolutionValueInput = 20, smoothValueInput = 60,
        edge_weighting = mode, .draw = FALSE))
    sc <- lapply(c("midpoint", "degree", "none"), function(mode) {
        r <- run(mode)
        expect_identical(r$metadata$edge_weighting, mode)
        s <- r$scores; setNames(s$LandscapeScore, s$Gene)
    })
    names(sc) <- c("midpoint", "degree", "none")
    # Midpoints carry half the hub signal towards the leaves: the more weight
    # they get, the lower the hub reads and the higher the leaves read.
    expect_lt(sc$midpoint["HUB"], sc$degree["HUB"])
    expect_lt(sc$degree["HUB"], sc$none["HUB"])
    leaves <- paste0("L", 1:8)
    expect_gt(mean(sc$midpoint[leaves]), mean(sc$degree[leaves]))
    expect_gt(mean(sc$degree[leaves]), mean(sc$none[leaves]))
    # A constant signal is neutral under every weighting.
    flat <- data.frame(ID = expr$ID, logFC = 0)
    for (mode in c("midpoint", "degree", "none")) {
        r <- suppressMessages(levi(expressionInput = flat,
            networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
            fileTypeInput = "stg", geneSymbolInput = "ID",
            readExpColumn = readExpColumn("logFC-logFC"), signal_mode = "logfc",
            edge_weighting = mode, .draw = FALSE))
        z <- r$landscape$z
        expect_true(all(abs(z[!is.na(z)] - 0.5) < 1e-9))
    }
})

test_that("edge_weighting is carried into the sample-label landscape null", {
    ang <- seq(0, 2 * pi, length.out = 9)[-9]
    net <- list(nodes = data.frame(name = c("HUB", paste0("L", 1:8)),
                                   x = c(50, 50 + 40 * cos(ang)), y = c(50, 50 + 40 * sin(ang))),
                edges = data.frame(V1 = "HUB", V2 = paste0("L", 1:8)))
    genes <- c("HUB", paste0("L", 1:8))
    set.seed(11)
    x <- matrix(rnorm(9 * 8, 6, .2), 9, 8, dimnames = list(genes, NULL))
    g <- rep(c("C", "T"), each = 4); x["HUB", g == "T"] <- x["HUB", g == "T"] + 2
    r <- suppressMessages(leviReplicateInference(x, g, test = "T", control = "C",
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        fileTypeInput = "stg", resolutionValueInput = 15, smoothValueInput = 30,
        n_perm = 20, seed = 1, edge_weighting = "degree"))
    expect_identical(r$metadata$edge_weighting, "degree")
    expect_length(r$metadata$support_weights, 9 + 8)
})

test_that("the grid margin adapts to the kernel so the silhouette is never clipped", {
    net <- system.file("extdata", "medusa.dat", package = "levi")
    expr <- system.file("extdata", "expression.dat", package = "levi")
    border_occupied <- function(smooth, contrast, zoom, resolution = 30) {
        r <- suppressMessages(levi(expressionInput = expr, networkCoordinatesInput = net,
            fileTypeInput = "dat", geneSymbolInput = "ID",
            readExpColumn = readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
            smoothValueInput = smooth, contrastValueInput = contrast,
            zoomValueInput = zoom, resolutionValueInput = resolution, .draw = FALSE))
        z <- r$landscape; n <- max(z$Var1)
        border <- z$Var1 %in% c(1, n) | z$Var2 %in% c(1, n)
        c(occupied = sum(!is.na(z$z[border])), inside = sum(!is.na(z$z)))
    }
    for (smooth in c(10, 50, 100)) for (contrast in c(0, 50, 100)) for (zoom in c(0, 50, 100)) {
        b <- border_occupied(smooth, contrast, zoom)
        expect_equal(unname(b["occupied"]), 0,
            info = sprintf("smooth %d contrast %d zoom %d", smooth, contrast, zoom))
        expect_gt(unname(b["inside"]), 0)
    }
    # zoom 0 frames wider than zoom 100: fewer occupied cells on the same grid
    expect_lt(border_occupied(50, 50, 0)["inside"], border_occupied(50, 50, 100)["inside"])
})
