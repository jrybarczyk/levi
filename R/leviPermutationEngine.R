# Shared machinery for the permutation tests in leviNetworkStatistics.R,
# leviSingleCellStatistics.R, leviInference.R and leviRNAseqInference.R.
#
# Every one of those tests follows the same recipe: compute a statistic on the
# observed data, recompute it on a list of null draws (permuted sample labels,
# shuffled node scores or rewired networks), keep the maximum per draw and
# compare. This file holds the recipe once. The public functions only define
# which statistic and which kind of draw they use.
#
# All randomness lives in the functions that build the draws, which run in the
# calling process. The statistic applied to each draw is deterministic, so the
# result does not depend on the BiocParallel back-end and set.seed() before a
# call reproduces it with any BPPARAM.

# -- Null draws ---------------------------------------------------------------

# Sample-label permutations, optionally within exchangeability blocks.
#
# Returns a list with $labels (list of label vectors), $exact (whether every
# distinct assignment was enumerated) and $possible (their number). In the
# exact scheme the observed assignment is left out of $labels; the (k + 1) /
# (n + 1) estimator adds it back.
.labelPermutations <- function(groups, blocks = NULL, n_perm = 999L,
                               method = c("auto", "monte_carlo", "exact"),
                               max_exact = 50000L, alpha = 0.05) {
    method <- match.arg(method)
    groups <- as.character(groups)
    if (length(unique(groups)) != 2L)
        stop("Exactly two groups are required.", call. = FALSE)
    if (is.null(blocks)) blocks <- rep("all", length(groups))
    by_block <- split(seq_along(groups), blocks)

    # Count the distinct arrangements before enumerating anything: a block
    # that holds a single condition (a donor with only controls, say) has
    # nothing to exchange and contributes one fixed arrangement. Counting
    # with choose() keeps an unblocked 58-vs-49 design (about 1e31
    # arrangements) from being expanded with combn(), which used to overflow
    # the integer range and abort before the Monte Carlo branch was reached.
    n_possible <- prod(vapply(by_block, function(i) {
        lev <- unique(groups[i])
        if (length(lev) == 1L) return(1)
        choose(length(i), sum(groups[i] == lev[2]))
    }, numeric(1)))
    exact <- method == "exact" ||
        (method == "auto" && n_possible <= max_exact)
    if (exact && n_possible > max_exact)
        stop("Exact permutation space exceeds 'max_exact'.", call. = FALSE)

    # The smallest p-value a test can report is 1 / (draws + 1). When that
    # floor sits above alpha no result can be declared significant, which is
    # a property of the design and not of the data, so say so up front.
    n_draws <- if (exact) n_possible - 1 else as.integer(n_perm)
    floor_p <- 1 / (n_draws + 1)
    if (floor_p > alpha) {
        cause <- if (exact) {
            sprintf("Only %s distinct label arrangements are possible",
                    format(n_possible, big.mark = ","))
        } else {
            sprintf("With n_perm = %d Monte Carlo draws", n_draws)
        }
        text <- sprintf("%s, so the smallest attainable p-value is %.3f, ",
                        cause, floor_p)
        text <- sprintf("%sabove alpha = %g. No result of this test can be ",
                        text, alpha)
        warning(text, "significant at that level.", call. = FALSE)
    }

    if (!exact) {
        labels <- replicate(as.integer(n_perm), {
            x <- groups
            for (i in by_block) x[i] <- sample(x[i])
            x
        }, simplify = FALSE)
        return(list(labels = labels, exact = FALSE, possible = n_possible))
    }

    # Every way of placing the test labels inside one block.
    choices <- lapply(by_block, function(i) {
        lev <- unique(groups[i])
        if (length(lev) == 1L) return(list(groups[i]))
        test <- lev[2]
        k <- sum(groups[i] == test)
        lapply(utils::combn(seq_along(i), k, simplify = FALSE), function(pos) {
            x <- rep(lev[1], length(i))
            x[pos] <- test
            x
        })
    })
    grid <- expand.grid(lapply(choices, seq_along), KEEP.OUT.ATTRS = FALSE)
    labels <- lapply(seq_len(nrow(grid)), function(r) {
        x <- groups
        for (j in seq_along(by_block))
            x[by_block[[j]]] <- choices[[j]][[grid[r, j]]]
        x
    })
    labels <- labels[!vapply(labels, identical, logical(1), groups)]
    list(labels = labels, exact = TRUE, possible = n_possible)
}

# Node-label draws: the same scores in a random order.
.nodeShuffles <- function(x, n_perm) {
    replicate(as.integer(n_perm), sample(x), simplify = FALSE)
}

# Degree-preserving rewiring draws: one edge index matrix per null network.
.rewiredEdgeLists <- function(edges, n_nodes, n_perm, steps) {
    g <- igraph::simplify(.graphFromEdges(n_nodes, edges))
    replicate(as.integer(n_perm), {
        rewired <- igraph::rewire(g,
            with = igraph::keeping_degseq(niter = steps))
        igraph::as_edgelist(rewired, names = FALSE)
    }, simplify = FALSE)
}

# -- Null distribution and p-values ------------------------------------------

# Applies `statistic` to every draw and stacks the results as rows. Runs
# through BiocParallel; the default SerialParam keeps examples and tests
# single-threaded, and MulticoreParam() or SnowParam() parallelise the loop
# without any other change.
.permutationNull <- function(draws, statistic, BPPARAM = SerialParam(),
                             names = NULL) {
    rows <- bplapply(draws, statistic, BPPARAM = BPPARAM)
    null <- do.call(rbind, lapply(rows, as.numeric))
    if (is.null(null)) null <- matrix(numeric(), 0L, length(names))
    if (!is.null(names)) colnames(null) <- names
    null
}

# (k + 1) / (n + 1) upper-tail p-values of each observed value against one
# null vector of maxima. The +1 counts the observed data as a draw, so no
# finite number of permutations reports zero.
.maxPvalue <- function(null_max, observed) {
    (1 + colSums(outer(null_max, observed, `>=`))) / (length(null_max) + 1)
}

# Largest cluster mass in each direction from a .graphRegions() summary.
.directionMaxima <- function(regions) {
    c(over  = max(c(0, regions$Mass[regions$Direction == "over"])),
      under = max(c(0, regions$Mass[regions$Direction == "under"])))
}

# -- Network helpers ----------------------------------------------------------

# Node names and a de-duplicated edge index matrix, self-loops removed.
.networkIndex <- function(networkCoordinatesInput,
                          networkInteractionsInput = NA,
                          fileTypeInput = "dat") {
    p <- .parseNetwork(networkCoordinatesInput, networkInteractionsInput,
                       fileTypeInput)
    nodes <- unique(as.character(p$nodes[, 1]))
    e <- cbind(match(as.character(p$edges[, 1]), nodes),
               match(as.character(p$edges[, 2]), nodes))
    keep <- stats::complete.cases(e) & e[, 1] != e[, 2]
    list(nodes = nodes, edges = unique(e[keep, , drop = FALSE]))
}

# Restricts the network to the nodes that carry a finite score and renumbers
# the edges accordingly.
.alignedGraph <- function(scores, networkCoordinatesInput,
                          networkInteractionsInput, fileTypeInput) {
    ni <- .networkIndex(networkCoordinatesInput, networkInteractionsInput,
                        fileTypeInput)
    x <- as.numeric(scores[ni$nodes])
    keep <- is.finite(x)
    map <- cumsum(keep)
    e <- ni$edges[keep[ni$edges[, 1]] & keep[ni$edges[, 2]], , drop = FALSE]
    if (nrow(e)) e <- cbind(map[e[, 1]], map[e[, 2]])
    list(nodes = ni$nodes[keep], x = x[keep], edges = e)
}

.graphFromEdges <- function(n_nodes, edges) {
    g <- igraph::make_empty_graph(n = n_nodes, directed = FALSE)
    if (nrow(edges)) g <- igraph::add_edges(g, as.vector(t(edges)))
    g
}

# Symmetric sparse adjacency matrix, so that the autocorrelation statistics
# cost O(edges) rather than O(nodes^2). `weights` gives one value per edge
# row; `self` puts ones on the diagonal (Getis-Ord Gi* includes the node
# itself). A repeated edge keeps its last weight, as an assignment would.
.adjacencyMatrix <- function(n_nodes, edges, weights = NULL, self = FALSE) {
    i <- c(edges[, 1], edges[, 2])
    j <- c(edges[, 2], edges[, 1])
    x <- if (is.null(weights)) rep(1, length(i)) else rep(weights, 2)
    if (self) {
        i <- c(i, seq_len(n_nodes))
        j <- c(j, seq_len(n_nodes))
        x <- c(x, rep(1, n_nodes))
    }
    Matrix::sparseMatrix(i = as.integer(i), j = as.integer(j), x = x,
                         dims = c(n_nodes, n_nodes), use.last.ij = TRUE)
}

# Combinatorial Laplacian D - W of a (weighted) adjacency matrix.
.laplacian <- function(W) {
    Matrix::Diagonal(x = Matrix::rowSums(W)) - W
}

# Quadratic form t(x) M x for a sparse symmetric M.
.quadraticForm <- function(x, M) {
    as.numeric(Matrix::crossprod(x, M %*% x))
}

# Eigenvectors of the `k` smallest Laplacian eigenvalues. The statistic
# needs the lowest quarter of the spectrum, and partial solvers (RSpectra,
# with or without shift-invert) were measured slower than a full symmetric
# eigen() for k = n/4, so the dense decomposition is used throughout. It is
# O(n^3): about 25 s at 2500 nodes and 3.5 min at 5000 on one core.
.lowFrequencyBasis <- function(L, k, cost_limit = 4000L) {
    n <- nrow(L)
    k <- min(k, n)
    if (n > cost_limit)
        message("Eigendecomposition of a ", n, "-node Laplacian; this takes ",
                "minutes. The spectrum test is meant for pathway-scale ",
                "networks; Moran, Getis-Ord and rewiring tests scale to ",
                "large graphs.")
    vectors <- eigen(as.matrix(L), symmetric = TRUE)$vectors
    vectors[, seq.int(n - k + 1L, n), drop = FALSE]
}

# -- Expression helpers -------------------------------------------------------

# Moderated t for the last design column: `test` versus `control`, with an
# optional block fixed effect. `trend = TRUE` is the limma-trend variant for
# logCPM-type data, where the variance depends on the mean.
.moderatedT <- function(expression, groups, levels, blocks = NULL,
                        trend = FALSE) {
    g <- factor(groups, levels = levels)
    design <- if (is.null(blocks)) stats::model.matrix(~ g) else
        stats::model.matrix(~ factor(blocks) + g)
    fit <- limma::eBayes(limma::lmFit(expression, design), trend = trend)
    stats::setNames(fit$t[, ncol(design)], rownames(expression))
}

# log2 counts per million with a pseudocount of one. With `normalize =
# "TMM"` the library sizes are scaled by edgeR's trimmed mean of M-values so
# that a few highly expressed genes do not drive the composition.
.logCPM <- function(counts, normalize = c("TMM", "none")) {
    normalize <- match.arg(normalize)
    lib <- colSums(counts)
    if (normalize == "TMM") {
        if (!requireNamespace("edgeR", quietly = TRUE))
            stop("'edgeR' is required for TMM normalisation; install it or ",
                 "use normalize = \"none\".", call. = FALSE)
        # normLibSizes() replaced calcNormFactors() in edgeR 4.0.
        tmm <- if (exists("normLibSizes", asNamespace("edgeR")))
            edgeR::normLibSizes else edgeR::calcNormFactors
        lib <- lib * tmm(counts, method = "TMM")
    }
    expr <- log2(t(t(counts) / lib * 1e6) + 1)
    rownames(expr) <- rownames(counts)
    expr
}

# Keeps the cell types whose pseudobulks cover every required column key
# (donor, or donor::condition in the paired design).
.completePseudobulkTypes <- function(pb, keys, types, paired) {
    have <- if (paired) paste(pb$donor, pb$condition, sep = "::") else pb$donor
    types[vapply(types, function(tp) all(keys %in% have[pb$cell_type == tp]),
                 logical(1))]
}

# One logCPM matrix per cell type, columns in the order of `keys`.
#
# Genes that are not expressed are dropped before normalisation and the
# limma fit: a droplet matrix carries tens of thousands of all-zero rows,
# and fitting them gives more than half of the residual variances exactly
# zero, which breaks the eBayes prior and the variance trend ("eBayes
# unreliable") for the genes that matter. A gene is kept when it has at
# least `min_count` counts in `min_samples` pseudobulks, or when it belongs
# to the network (`keep_genes`), so every node keeps a statistic.
.pseudobulkExpression <- function(pb, keys, types, paired,
                                  normalize = "TMM", keep_genes = NULL,
                                  min_count = 10L, min_samples = 2L) {
    have <- if (paired) paste(pb$donor, pb$condition, sep = "::") else pb$donor
    out <- lapply(types, function(tp) {
        at <- which(pb$cell_type == tp)
        cols <- at[match(keys, have[at])]
        counts <- pb$counts[, cols, drop = FALSE]
        expressed <- rowSums(counts >= min_count) >= min_samples
        if (!is.null(keep_genes))
            expressed <- expressed | rownames(counts) %in% keep_genes
        .logCPM(counts[expressed, , drop = FALSE], normalize)
    })
    names(out) <- types
    out
}
