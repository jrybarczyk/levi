# Graph-native statistics: TFCE, autocorrelation, spectral smoothness and
# rewiring nulls. The permutation recipe shared by all of them lives in
# leviPermutationEngine.R.

# -- Threshold-free cluster enhancement on a graph ---------------------------

# Connected components of the subgraph induced by `active` nodes, on a graph
# built once by the caller.
.componentsAt <- function(active, g) {
    if (!length(active)) return(list())
    sub <- igraph::induced_subgraph(g, active)
    split(active, igraph::components(sub)$membership)
}

# Smith & Nichols (2009) TFCE with network components as clusters. Returns a
# nodes x 2 matrix with the "over" (positive) and "under" (negative) scores.
#
# The graph is built once per call and the components are recomputed only
# when the active set changes between consecutive thresholds: with a few
# nodes most of the n_steps thresholds share the same active set, and this
# routine runs once per permutation draw.
.tfce <- function(statistic, nodes, edges, E = .5, H = 2, n_steps = 100L) {
    x <- as.numeric(statistic[nodes])
    ans <- matrix(0, length(nodes), 2,
                  dimnames = list(nodes, c("over", "under")))
    g <- .graphFromEdges(length(nodes), edges)
    for (side in c("over", "under")) {
        v <- if (side == "over") pmax(x, 0) else pmax(-x, 0)
        top <- max(v, na.rm = TRUE)
        if (!is.finite(top) || top == 0) next
        dh <- top / n_steps
        active_prev <- NULL
        members <- list()
        for (h in seq(dh, top, length.out = n_steps)) {
            active <- which(v >= h)
            if (!identical(active, active_prev)) {
                members <- .componentsAt(active, g)
                active_prev <- active
            }
            for (cc in members)
                ans[cc, side] <- ans[cc, side] + length(cc)^E * h^H * dh
        }
    }
    ans
}

# Column maxima of a TFCE matrix: the per-draw statistic of every TFCE test.
.tfceMaxima <- function(statistic, nodes, edges, E, H, n_steps) {
    apply(.tfce(statistic, nodes, edges, E, H, n_steps), 2, max)
}

#' Threshold-free graph-cluster enhancement with sample-label inference
#'
#' Fits a moderated `limma` t-statistic for every gene, integrates the signed
#' statistics over all thresholds with threshold-free cluster enhancement
#' (TFCE) using biological network components as clusters, and assigns
#' FWER-controlled P-values by permuting sample labels. `PGlobal` uses the
#' joint maximum across genes and both directions.
#' @param expression Log-scale genes-by-samples expression matrix with gene
#'   identifiers as row names.
#' @param groups Two-level sample-group vector, one value per column.
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput LEVI
#'   network inputs: node coordinates, optional interactions and the file type
#'   accepted by [levi()].
#' @param test,control Group labels defining the `test - control` contrast.
#' @param blocks Optional exchangeability blocks (e.g. donor identifiers).
#'   Labels are permuted only within blocks and the block enters the linear
#'   model as a fixed effect.
#' @param n_perm Number of Monte Carlo sample-label permutations.
#' @param permutation_method `"auto"`, `"monte_carlo"` or `"exact"`.
#'   `"exact"` enumerates every distinct label assignment; `"monte_carlo"`
#'   draws `n_perm` random permutations; `"auto"` uses the exact scheme when
#'   the number of assignments is at most `max_exact` and Monte Carlo otherwise.
#' @param max_exact Maximum number of label assignments enumerated by the exact
#'   permutation scheme.
#' @param E,H TFCE extent and height exponents applied to the component size
#'   and to the threshold height at every integration step.
#' @param n_steps Number of threshold steps used to integrate the TFCE score
#'   between zero and the maximum statistic.
#' @param seed Optional RNG seed.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling how
#'   the permutations are evaluated. The default runs them serially.
#' @return A list with signed gene statistics, TFCE scores, FWER-adjusted
#'   P-values and the permutation null maxima.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C"), X = 1:3, Y = 1:3)
#' edges <- data.frame(V1 = c("A", "B"), V2 = c("B", "C"))
#' expression <- matrix(rnorm(18), 3, dimnames = list(c("A", "B", "C"), NULL))
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   leviGraphTFCEInference(expression, rep(c("control", "case"), each = 3),
#'     network, edges, fileTypeInput = "stg", n_perm = 3)
#' }
#' @export
leviGraphTFCEInference <- function(expression, groups, networkCoordinatesInput,
    networkInteractionsInput = NA, fileTypeInput = "dat",
    test = unique(groups)[2], control = unique(groups)[1], blocks = NULL,
    n_perm = 999L, permutation_method = c("auto", "monte_carlo", "exact"),
    max_exact = 50000L, E = .5, H = 2, n_steps = 100L, seed = NULL,
    BPPARAM = SerialParam()) {
    if (!requireNamespace("limma", quietly = TRUE))
        stop("'limma' is required.", call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    levels <- c(control, test)
    statistic <- function(g) .moderatedT(expression, g, levels, blocks)

    observed_t <- statistic(groups)
    observed <- .tfce(observed_t, ni$nodes, ni$edges, E, H, n_steps)

    perms <- .labelPermutations(groups, blocks, n_perm, permutation_method,
                                max_exact)
    null <- .permutationNull(perms$labels, function(g)
        .tfceMaxima(statistic(g), ni$nodes, ni$edges, E, H, n_steps),
        BPPARAM = BPPARAM, names = c("over", "under"))
    null_global <- apply(null, 1L, max)
    p_global <- .maxPvalue(null_global, pmax(observed[, 1], observed[, 2]))

    list(
        statistic = data.frame(
            Gene = ni$nodes,
            T = as.numeric(observed_t[ni$nodes]),
            TFCE = observed[, 1] - observed[, 2],
            POver = .maxPvalue(null[, 1], observed[, 1]),
            PUnder = .maxPvalue(null[, 2], observed[, 2]),
            PGlobal = p_global,
            GlobalSignificant = p_global <= .05),
        null_max = null, null_global = null_global,
        exact = perms$exact, possible_permutations = perms$possible,
        method = paste("limma t; graph TFCE; joint maximum across network",
                       "and direction sample-label permutation"))
}

# -- Spatial autocorrelation --------------------------------------------------

#' Moran global and local autocorrelation on a biological graph
#'
#' Computes global Moran's I and local Moran's I for gene scores using the
#' binary adjacency of the biological network, with hotspot classes and
#' node-label permutation P-values.
#' @param scores Named numeric vector of gene scores (for example moderated
#'   t-statistics or log fold-changes) indexed by gene identifier. Network
#'   nodes without a finite score are dropped before testing.
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput LEVI
#'   network inputs: node coordinates, optional interactions and the file type
#'   accepted by [levi()].
#' @param n_perm Number of node-label permutations used to build the null.
#' @param seed Optional RNG seed.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling how
#'   the permutations are evaluated. The default runs them serially.
#' @return A list with one global Moran's I test and a table of local Moran's I
#'   scores, hotspot classes, two-sided permutation P-values and
#'   Benjamini-Hochberg adjusted P-values across nodes.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C"), X = 1:3, Y = 1:3)
#' edges <- data.frame(V1 = c("A", "B"), V2 = c("B", "C"))
#' leviGraphMoran(c(A = 1, B = 2, C = 3), network, edges,
#'   fileTypeInput = "stg", n_perm = 9)
#' @export
leviGraphMoran <- function(scores, networkCoordinatesInput,
    networkInteractionsInput = NA, fileTypeInput = "dat", n_perm = 999L,
    seed = NULL, BPPARAM = SerialParam()) {
    if (!is.null(seed)) withr::local_seed(seed)
    a <- .alignedGraph(scores, networkCoordinatesInput,
                       networkInteractionsInput, fileTypeInput)
    n <- length(a$x)
    W <- .adjacencyMatrix(n, a$edges)

    # Global I and the local I of every node, in one vector.
    total_weight <- sum(W)
    statistic <- function(x) {
        z <- x - mean(x)
        global <- n / total_weight * .quadraticForm(z, W) / sum(z^2)
        zs <- as.numeric(scale(x))
        c(global, zs * as.numeric(W %*% zs))
    }
    observed <- statistic(a$x)
    null <- .permutationNull(.nodeShuffles(a$x, n_perm), statistic,
                             BPPARAM = BPPARAM)

    # Two-sided p-values: |I| under the null at least as large as observed.
    p <- (1 + colSums(abs(null) >= rep(abs(observed), each = nrow(null)))) /
        (nrow(null) + 1)

    degree <- Matrix::rowSums(W)
    neighbour_mean <- ifelse(degree > 0, as.numeric(W %*% a$x) / degree, NA)
    centre <- mean(a$x)
    class <- ifelse(a$x >= centre,
        ifelse(neighbour_mean >= centre, "High-High", "High-Low"),
        ifelse(neighbour_mean <= centre, "Low-Low", "Low-High"))

    # Local p-values form one family per network; adjust them together.
    list(
        global = data.frame(MoranI = observed[1], P = p[1]),
        local = data.frame(Gene = a$nodes, I = observed[-1], Class = class,
                           P = p[-1],
                           PAdjusted = stats::p.adjust(p[-1], "BH")),
        method = "node-label permutation")
}

#' Getis-Ord hotspots and coldspots on a biological graph
#'
#' Computes the local Getis--Ord Gi* statistic for every gene, including the
#' node itself in its neighbourhood, and assigns one-sided node-label
#' permutation P-values for hotspots and coldspots. A node adjacent to every
#' other node has no variance to compare against and receives `NA`.
#' @inheritParams leviGraphMoran
#' @return A data frame of local Getis--Ord Gi* scores with one-sided
#'   permutation P-values for each direction, the two-sided P-value (twice the
#'   smaller one-sided value), its Benjamini-Hochberg adjustment across nodes
#'   and a `Significant` flag at adjusted P <= 0.05.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C"), X = 1:3, Y = 1:3)
#' edges <- data.frame(V1 = c("A", "B"), V2 = c("B", "C"))
#' leviGraphGetisOrd(c(A = 1, B = 2, C = 3), network, edges,
#'   fileTypeInput = "stg", n_perm = 9)
#' @export
leviGraphGetisOrd <- function(scores, networkCoordinatesInput,
    networkInteractionsInput = NA, fileTypeInput = "dat", n_perm = 999L,
    seed = NULL, BPPARAM = SerialParam()) {
    if (!is.null(seed)) withr::local_seed(seed)
    a <- .alignedGraph(scores, networkCoordinatesInput,
                       networkInteractionsInput, fileTypeInput)
    n <- length(a$x)
    W <- .adjacencyMatrix(n, a$edges, self = TRUE)
    wi <- Matrix::rowSums(W)
    wi2 <- Matrix::rowSums(W^2)

    # A node whose neighbourhood covers the whole graph (the centre of a
    # star, say) has n * wi2 == wi^2 and an undefined Gi*, reported as NA.
    statistic <- function(x) {
        denominator <- stats::sd(x) * sqrt((n * wi2 - wi^2) / pmax(n - 1, 1))
        gi <- as.numeric((W %*% x - mean(x) * wi) / denominator)
        gi[!is.finite(gi)] <- NA_real_
        gi
    }
    observed <- statistic(a$x)
    null <- .permutationNull(.nodeShuffles(a$x, n_perm), statistic,
                             BPPARAM = BPPARAM)

    obs_rows <- rep(observed, each = nrow(null))
    p_over  <- (1 + colSums(null >= obs_rows, na.rm = TRUE)) / (nrow(null) + 1)
    p_under <- (1 + colSums(null <= obs_rows, na.rm = TRUE)) / (nrow(null) + 1)
    p_over[!is.finite(observed)] <- NA_real_
    p_under[!is.finite(observed)] <- NA_real_

    # Taking the smaller one-sided p-value doubles the size of the test, so
    # the two-sided p-value is twice that minimum (capped at one). Nodes are
    # then adjusted as one family.
    p_two_sided <- pmin(1, 2 * pmin(p_over, p_under))
    p_adjusted <- stats::p.adjust(p_two_sided, "BH")

    data.frame(Gene = a$nodes, GiStar = observed, POver = p_over,
               PUnder = p_under, PTwoSided = p_two_sided,
               PAdjusted = p_adjusted,
               Class = ifelse(observed > 0, "hotspot", "coldspot"),
               Significant = p_adjusted <= .05,
               method = "node-label permutation")
}

#' Weighted Moran and Laplacian statistics for scored biological interactions
#'
#' `edge_weights` must have one value per interaction row (or be the name/index
#' of a score column in an interaction data frame). It lets STRING confidence
#' scores contribute to neighbourhood and smoothness calculations.
#'
#' @inheritParams leviGraphMoran
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput LEVI
#'   network inputs. `networkInteractionsInput` must be supplied because the
#'   edge weights are read from its rows.
#' @param edge_weights Positive edge weights, one per interaction row, or the
#'   name or index of the column in `networkInteractionsInput` holding them
#'   (default: the third column).
#' @return A one-row data frame with weighted Moran's I, Laplacian energy and
#'   node-label permutation P-values.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C"), X = 1:3, Y = 1:3)
#' edges <- data.frame(V1 = c("A", "B"), V2 = c("B", "C"), score = c(1, 2))
#' leviGraphWeightedTopology(c(A = 1, B = 2, C = 3), network, edges,
#'   edge_weights = "score", fileTypeInput = "stg", n_perm = 9)
#' @export
leviGraphWeightedTopology <- function(scores, networkCoordinatesInput,
    networkInteractionsInput, edge_weights = 3L, fileTypeInput = "stg",
    n_perm = 999L, seed = NULL, BPPARAM = SerialParam()) {
    if (!is.null(seed)) withr::local_seed(seed)
    p <- .parseNetwork(networkCoordinatesInput, networkInteractionsInput,
                       fileTypeInput)
    nodes <- as.character(p$nodes[, 1])
    x <- as.numeric(scores[nodes])
    keep <- is.finite(x)
    nodes <- nodes[keep]
    x <- x[keep]

    # The weights come from the raw interaction table, not from the parsed
    # edges, so that any score column can be used.
    raw <- if (is.data.frame(networkInteractionsInput) ||
               is.matrix(networkInteractionsInput)) {
        networkInteractionsInput
    } else {
        utils::read.table(networkInteractionsInput, header = TRUE,
                          sep = "\t", check.names = FALSE)
    }
    is_column <- length(edge_weights) == 1L &&
        (is.character(edge_weights) || is.numeric(edge_weights))
    w <- if (is_column) as.numeric(raw[[edge_weights]]) else
        as.numeric(edge_weights)
    if (length(w) != nrow(raw))
        stop("'edge_weights' must align to interaction rows.", call. = FALSE)
    i <- match(as.character(raw[, 1]), nodes)
    j <- match(as.character(raw[, 2]), nodes)
    ok <- is.finite(i) & is.finite(j) & i != j & is.finite(w) & w > 0
    W <- .adjacencyMatrix(length(nodes), cbind(i[ok], j[ok]), weights = w[ok])
    L <- .laplacian(W)
    total_weight <- sum(W)

    statistic <- function(z) {
        zc <- z - mean(z)
        c(MoranI = length(z) / total_weight * .quadraticForm(zc, W) /
              sum(zc^2),
          LaplacianEnergy = .quadraticForm(z, L))
    }
    observed <- statistic(x)
    null <- .permutationNull(.nodeShuffles(x, n_perm), statistic,
                             BPPARAM = BPPARAM,
                             names = c("MoranI", "LaplacianEnergy"))

    data.frame(
        MoranI = observed[1],
        MoranP = (1 + sum(abs(null[, 1]) >= abs(observed[1]))) /
            (nrow(null) + 1),
        LaplacianEnergy = observed[2],
        EnergyP = (1 + sum(null[, 2] <= observed[2])) / (nrow(null) + 1),
        Nodes = length(nodes), Edges = Matrix::nnzero(W) / 2,
        method = "STRING-weighted node-label permutation")
}

# -- Spectral smoothness ------------------------------------------------------

#' Laplacian smoothness and graph-Fourier low-frequency energy
#'
#' Measures how smoothly gene scores vary over the biological network through
#' the Laplacian quadratic form and the share of signal energy in the lowest
#' quarter of the graph-Fourier spectrum, with node-label permutation P-values.
#' The spectral basis requires a full eigendecomposition of the Laplacian,
#' computed once; the cost grows with the cube of the node count, so the test
#' is intended for pathway-scale networks (up to a few thousand nodes).
#' @inheritParams leviGraphMoran
#' @return A one-row data frame of Laplacian energy, low-frequency energy ratio
#'   and node-label permutation P-values.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C"), X = 1:3, Y = 1:3)
#' edges <- data.frame(V1 = c("A", "B"), V2 = c("B", "C"))
#' leviGraphSpectrum(c(A = 1, B = 2, C = 3), network, edges,
#'   fileTypeInput = "stg", n_perm = 9)
#' @export
leviGraphSpectrum <- function(scores, networkCoordinatesInput,
    networkInteractionsInput = NA, fileTypeInput = "dat", n_perm = 999L,
    seed = NULL, BPPARAM = SerialParam()) {
    if (!is.null(seed)) withr::local_seed(seed)
    a <- .alignedGraph(scores, networkCoordinatesInput,
                       networkInteractionsInput, fileTypeInput)
    n <- length(a$x)
    W <- .adjacencyMatrix(n, a$edges)
    L <- .laplacian(W)

    # The Laplacian does not change between permutations, so the basis of
    # the lowest quarter of the spectrum is computed once. The graph-Fourier
    # basis is orthonormal, so the total energy sum(c^2) equals sum(x^2) and
    # the remaining eigenvectors are never needed.
    n_low <- max(1, floor(n / 4) + 1)
    basis <- .lowFrequencyBasis(L, n_low)

    statistic <- function(x) {
        coefficients <- crossprod(basis, x)
        c(energy = .quadraticForm(x, L),
          low_frequency = sum(coefficients^2) / sum(x^2))
    }
    observed <- statistic(a$x)
    null <- .permutationNull(.nodeShuffles(a$x, n_perm), statistic,
                             BPPARAM = BPPARAM,
                             names = c("energy", "low_frequency"))

    data.frame(
        LaplacianEnergy = observed[1],
        EnergyP = (1 + sum(null[, 1] <= observed[1])) / (nrow(null) + 1),
        LowFrequency = observed[2],
        LowFrequencyP = (1 + sum(null[, 2] >= observed[2])) / (nrow(null) + 1))
}

# -- Degree-preserving rewiring nulls -----------------------------------------

#' Degree-preserving rewiring test for graph clusters
#'
#' Thresholds gene scores, forms connected components on the biological
#' network and compares the observed maximum cluster mass with the mass
#' obtained on degree-preserving rewired networks. This tests whether the
#' clustering depends on the specific wiring rather than on the degree
#' sequence.
#' @inheritParams leviGraphMoran
#' @param threshold Absolute score threshold for cluster formation.
#' @param n_perm Number of rewired null networks.
#' @param rewire_steps Number of degree-preserving edge swaps passed to
#'   [igraph::keeping_degseq()] for every null network. Defaults to ten times
#'   the number of edges (at least 10).
#' @return A list of observed graph regions, degree-preserving rewiring null
#'   maxima and topology P-values.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C", "D"), X = 1:4, Y = 1:4)
#' edges <- data.frame(V1 = c("A", "B", "C", "D"), V2 = c("B", "C", "D", "A"))
#' leviGraphRewiringInference(c(A = 3, B = 2.5, C = -1, D = -2), network, edges,
#'   fileTypeInput = "stg", n_perm = 9)
#' @export
leviGraphRewiringInference <- function(scores, networkCoordinatesInput,
    networkInteractionsInput = NA, fileTypeInput = "dat", threshold = 2,
    n_perm = 999L, rewire_steps = NULL, seed = NULL,
    BPPARAM = SerialParam()) {
    if (!is.null(seed)) withr::local_seed(seed)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    if (is.null(rewire_steps)) rewire_steps <- max(10, 10 * nrow(ni$edges))

    observed <- list(summary = .graphRegions(scores, ni$nodes, ni$edges,
                                             threshold))
    draws <- .rewiredEdgeLists(ni$edges, length(ni$nodes), n_perm,
                               rewire_steps)
    null <- .permutationNull(draws, function(e)
        .directionMaxima(.graphRegions(scores, ni$nodes, e, threshold)),
        BPPARAM = BPPARAM, names = c("over", "under"))

    regions <- .regionalPvalues(observed, null)
    regions$summary$Significant <- regions$summary$PSpatial <= .05
    list(regions = regions, null_max = null,
         method = "degree-preserving graph rewiring; maximum cluster mass")
}

#' Degree-preserving rewiring sensitivity test for graph TFCE
#'
#' Computes graph TFCE scores for fixed gene scores and compares them with the
#' TFCE maxima obtained on degree-preserving rewired networks.
#' @inheritParams leviGraphRewiringInference
#' @inheritParams leviGraphTFCEInference
#' @return A data frame with TFCE scores and degree-preserving rewiring
#'   P-values.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C", "D"), X = 1:4, Y = 1:4)
#' edges <- data.frame(V1 = c("A", "B", "C", "D"), V2 = c("B", "C", "D", "A"))
#' leviGraphTFCERewiring(c(A = 3, B = 2.5, C = -1, D = -2), network, edges,
#'   fileTypeInput = "stg", n_perm = 9)
#' @export
leviGraphTFCERewiring <- function(scores, networkCoordinatesInput,
    networkInteractionsInput = NA, fileTypeInput = "dat", n_perm = 999L,
    rewire_steps = NULL, E = .5, H = 2, n_steps = 100L, seed = NULL,
    BPPARAM = SerialParam()) {
    if (!is.null(seed)) withr::local_seed(seed)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    if (is.null(rewire_steps)) rewire_steps <- max(10, 10 * nrow(ni$edges))

    observed <- .tfce(scores, ni$nodes, ni$edges, E, H, n_steps)
    draws <- .rewiredEdgeLists(ni$edges, length(ni$nodes), n_perm,
                               rewire_steps)
    null <- .permutationNull(draws, function(e)
        .tfceMaxima(scores, ni$nodes, e, E, H, n_steps),
        BPPARAM = BPPARAM, names = c("over", "under"))
    p_global <- .maxPvalue(apply(null, 1L, max),
                           pmax(observed[, 1], observed[, 2]))

    data.frame(Gene = ni$nodes, TFCEOver = observed[, 1],
               TFCEUnder = observed[, 2], PRewireGlobal = p_global,
               Significant = p_global <= .05)
}

# -- Freedman-Lane ------------------------------------------------------------

#' Graph TFCE with Freedman--Lane residual permutation
#'
#' Tests a two-level condition while retaining nuisance covariates. Residuals
#' from the nuisance-only model are permuted, then added back to fitted values
#' before refitting the full model. `blocks` restrict residual exchanges.
#'
#' @inheritParams leviGraphTFCEInference
#' @param covariates Data frame of nuisance covariates with one row per
#'   sample. They define the reduced model whose residuals are permuted.
#' @param blocks Optional exchangeability blocks; residuals are exchanged only
#'   within blocks.
#' @param n_perm Number of Freedman--Lane residual permutations.
#' @return A list with gene-level statistics, TFCE scores and a global
#'   Freedman--Lane permutation null distribution.
#' @examples
#' network <- data.frame(ID = c("A", "B", "C"), X = 1:3, Y = 1:3)
#' edges <- data.frame(V1 = c("A", "B"), V2 = c("B", "C"))
#' expression <- matrix(rnorm(18), 3, dimnames = list(c("A", "B", "C"), NULL))
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   leviGraphTFCEFreedmanLane(expression, rep(c("control", "case"), each = 3),
#'     data.frame(batch = rep(1:2, 3)), network, edges, fileTypeInput = "stg",
#'     n_perm = 19)
#' }
#' @export
leviGraphTFCEFreedmanLane <- function(expression, groups, covariates,
    networkCoordinatesInput, networkInteractionsInput = NA,
    fileTypeInput = "dat", test = unique(groups)[2],
    control = unique(groups)[1], blocks = NULL, n_perm = 999L, E = .5, H = 2,
    n_steps = 100L, seed = NULL, BPPARAM = SerialParam()) {
    if (!requireNamespace("limma", quietly = TRUE))
        stop("'limma' is required.", call. = FALSE)
    expression <- as.matrix(expression)
    groups <- as.character(groups)
    if (ncol(expression) != length(groups) ||
        nrow(covariates) != length(groups))
        stop("Expression, groups and covariates have incompatible sample ",
             "counts.", call. = FALSE)
    if (!setequal(unique(groups), c(control, test)))
        stop("Exactly test and control are required.", call. = FALSE)
    if (is.null(blocks)) blocks <- rep("all", length(groups))
    if (length(blocks) != length(groups))
        stop("Invalid 'blocks'.", call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)

    # Reduced (nuisance-only) and full designs.
    x0 <- stats::model.matrix(~ ., data = as.data.frame(covariates))
    condition <- as.numeric(factor(groups, levels = c(control, test))) - 1
    x1 <- cbind(x0, condition = condition)
    fitted0 <- limma::lmFit(expression, x0)$coefficients %*% t(x0)
    residuals0 <- expression - fitted0

    statistic <- function(y) {
        fit <- limma::eBayes(limma::lmFit(y, x1))
        stats::setNames(fit$t[, ncol(x1)], rownames(y))
    }
    observed_t <- statistic(expression)
    observed <- .tfce(observed_t, ni$nodes, ni$edges, E, H, n_steps)

    # Draws are residual orderings, exchanged within blocks.
    by_block <- split(seq_along(groups), blocks)
    draws <- replicate(as.integer(n_perm), {
        ord <- seq_along(groups)
        for (i in by_block) ord[i] <- sample(i)
        ord
    }, simplify = FALSE)
    null <- .permutationNull(draws, function(ord) {
        y <- fitted0 + residuals0[, ord, drop = FALSE]
        .tfceMaxima(statistic(y), ni$nodes, ni$edges, E, H, n_steps)
    }, BPPARAM = BPPARAM, names = c("over", "under"))
    null_global <- apply(null, 1L, max)
    p_global <- .maxPvalue(null_global, pmax(observed[, 1], observed[, 2]))

    list(
        statistic = data.frame(
            Gene = ni$nodes, T = as.numeric(observed_t[ni$nodes]),
            TFCEOver = observed[, 1], TFCEUnder = observed[, 2],
            PGlobal = p_global, GlobalSignificant = p_global <= .05),
        null_global = null_global,
        method = "limma t; graph TFCE; Freedman-Lane residual permutation")
}
