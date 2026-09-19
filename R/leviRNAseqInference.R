#' Bulk RNA-seq graph-cluster inference
#'
#' Fits a DESeq2 negative-binomial model to raw counts, projects the signed Wald
#' statistic onto a biological graph, and assigns maximum-cluster-mass P-values
#' by permuting sample labels. It is intended for a two-group design without
#' additional covariates. When `blocks` is supplied, the block is included as
#' a fixed effect in the DESeq2 model and labels are permuted only within
#' blocks. This supports paired two-condition studies.
#' @param counts Integer genes-by-samples count matrix with gene row names.
#' @param groups Two-level group vector.
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput Network inputs.
#' @param test,control Group labels defining test minus control.
#' @param threshold Absolute Wald-statistic threshold for graph components.
#' @param n_perm Number of sample-label permutations.
#' @param blocks Optional exchangeability blocks.
#' @param seed Optional RNG seed.
#' @param refit_dispersions Logical. `FALSE` (default) keeps the size factors
#'   and dispersions estimated on the observed labels and refits only the
#'   negative-binomial GLM in each permutation, a test conditional on those
#'   estimates and considerably faster. `TRUE` re-estimates dispersions in
#'   every permutation.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling how
#'   the permutations are evaluated. The default runs them serially.
#' @param permutation_method `"auto"`, `"monte_carlo"` or `"exact"`.
#'   `"exact"` enumerates every distinct label assignment; `"monte_carlo"`
#'   draws `n_perm` random permutations; `"auto"` uses the exact scheme when
#'   the number of assignments is at most `max_exact` and Monte Carlo otherwise.
#' @param max_exact Maximum number of label assignments enumerated by the exact
#'   permutation scheme.
#' @details For a fast alternative on the same design, transform the counts
#'   with `limma::voom()` or `edgeR::cpm(log = TRUE)` and use
#'   [leviGraphClusterInference()], whose moderated t refits in milliseconds.
#' @return List containing the DESeq2 statistic, graph regions, null maxima
#'   and the dispersion estimator used (`"trend"` for the standard `DESeq()`
#'   pipeline, `"genewise"` when the trend could not be fitted).
#' @examples
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#' \donttest{
#'   set.seed(3)
#'   genes <- c("HUB", paste0("N", 1:8))
#'   counts <- matrix(rpois(length(genes) * 6, 30), length(genes), 6,
#'     dimnames = list(genes, paste0("s", 1:6)))
#'   counts[, 4:6] <- counts[, 4:6] + 30L
#'   out <- leviBulkGraphInference(counts, c("C", "C", "C", "T", "T", "T"),
#'     system.file("extdata", "hub_network.dat", package = "levi"),
#'     test = "T", control = "C", threshold = 1, n_perm = 3, seed = 1)
#'   out$regions$summary
#' }
#' }
#' @export
leviBulkGraphInference <- function(counts, groups, networkCoordinatesInput,
    networkInteractionsInput = NA, fileTypeInput = "dat",
    test = unique(groups)[2], control = unique(groups)[1], threshold = 2,
    n_perm = 999L, blocks = NULL, seed = NULL,
    permutation_method = c("auto", "monte_carlo", "exact"), max_exact = 50000L,
    refit_dispersions = FALSE, BPPARAM = SerialParam()) {
    if (!requireNamespace("DESeq2", quietly = TRUE))
        stop("'DESeq2' is required for bulk RNA-seq inference.", call. = FALSE)
    counts <- as.matrix(counts)
    if (is.null(rownames(counts))) stop("'counts' needs gene row names.", call. = FALSE)
    if (ncol(counts) != length(groups)) stop("Invalid 'groups'.", call. = FALSE)
    if (any(!is.finite(counts)) || any(counts < 0))
        stop("'counts' must be finite and non-negative.", call. = FALSE)
    counts <- round(counts)
    groups <- as.character(groups)
    if (!setequal(unique(groups), c(as.character(test), as.character(control))))
        stop("'groups' must contain exactly test and control.", call. = FALSE)
    has_blocks <- !is.null(blocks)
    if (is.null(blocks)) blocks <- rep("all", length(groups))
    if (length(blocks) != length(groups)) stop("Invalid 'blocks'.", call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)
    parsed <- .parseNetwork(networkCoordinatesInput, networkInteractionsInput, fileTypeInput)
    nodes <- unique(as.character(parsed$nodes[, 1]))
    edge_index <- cbind(match(as.character(parsed$edges[, 1]), nodes),
                        match(as.character(parsed$edges[, 2]), nodes))
    edge_index <- edge_index[stats::complete.cases(edge_index), , drop = FALSE]
    sample_names <- colnames(counts) %||% paste0("sample", seq_len(ncol(counts)))
    colnames(counts) <- sample_names
    make_dds <- function(g) {
        coldata <- data.frame(
            group = factor(g, levels = c(control, test)),
            block = factor(blocks), row.names = sample_names)
        DESeq2::DESeqDataSetFromMatrix(
            countData = counts, colData = coldata,
            design = if (has_blocks) ~ block + group else ~ group)
    }
    # Two dispersion estimators. "trend" is the standard DESeq() pipeline;
    # "genewise" skips the parametric trend fit, which small pilot studies
    # with few genes cannot always support.
    fit_dds <- function(dds, estimator) {
        if (estimator == "trend") return(DESeq2::DESeq(dds, quiet = TRUE))
        dds <- DESeq2::estimateSizeFactors(dds)
        dds <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
        DESeq2::dispersions(dds) <- S4Vectors::mcols(dds)$dispGeneEst
        DESeq2::nbinomWaldTest(dds, quiet = TRUE)
    }
    wald_of <- function(dds) {
        r <- DESeq2::results(dds, contrast = c("group", test, control))
        setNames(as.numeric(r$log2FoldChange / r$lfcSE), rownames(r))
    }
    # The estimator is chosen on the observed data. The same estimator must
    # produce every null draw; mixing the two would compare values on
    # different scales.
    fit_observed <- function(estimator) list(
        estimator = estimator, dds = fit_dds(make_dds(groups), estimator))
    observed_fit <- tryCatch(fit_observed("trend"),
                             error = function(e) fit_observed("genewise"))
    perms <- .labelPermutations(groups, blocks, n_perm, permutation_method, max_exact)
    # By default the size factors and dispersions estimated on the observed
    # labels are kept fixed and only the GLM is refitted per permutation: a
    # permutation test conditional on those estimates, about an order of
    # magnitude cheaper. refit_dispersions = TRUE re-estimates everything.
    permuted_fit <- function(g, fit) {
        if (refit_dispersions) return(fit_dds(make_dds(g), fit$estimator))
        dds <- make_dds(g)
        DESeq2::sizeFactors(dds) <- DESeq2::sizeFactors(fit$dds)
        DESeq2::dispersions(dds) <- DESeq2::dispersions(fit$dds)
        DESeq2::nbinomWaldTest(dds, quiet = TRUE)
    }
    if (refit_dispersions)
        message("Refitting DESeq2 dispersions for each of ",
                length(perms$labels), " label permutations. This is slow on ",
                "real data; pass BPPARAM = BiocParallel::MulticoreParam() ",
                "to parallelise or use refit_dispersions = FALSE.")
    attempt <- function(fit) list(fit = fit, null = .permutationNull(
        perms$labels,
        function(perm) .directionMaxima(.graphRegions(
            wald_of(permuted_fit(perm, fit)), nodes, edge_index, threshold)),
        BPPARAM = BPPARAM, names = c("over", "under")))
    run <- tryCatch(attempt(observed_fit), error = function(e) {
        if (observed_fit$estimator == "genewise" || !refit_dispersions)
            stop(e)
        message("DESeq2's dispersion trend could not be fitted in every ",
                "permutation; recomputing observed and null statistics ",
                "with gene-wise dispersions.")
        attempt(fit_observed("genewise"))
    })
    estimator <- run$fit$estimator
    observed_stat <- wald_of(run$fit$dds)
    null <- run$null
    observed <- list(summary = .graphRegions(observed_stat, nodes, edge_index, threshold))
    regions <- .regionalPvalues(observed, null)
    regions$summary$Significant <- regions$summary$PSpatial <= .05
    list(statistic = data.frame(Gene = names(observed_stat), Wald = as.numeric(observed_stat)),
         regions = regions, threshold = threshold, n_perm = nrow(null), exact = perms$exact,
         dispersion_estimator = estimator,
         refit_dispersions = refit_dispersions,
         method = "DESeq2 Wald; graph-connected maximum cluster mass; sample labels")
}

#' Aggregate single-cell counts to pseudobulks
#'
#' @param counts Genes-by-cells count matrix.
#' @param donor Donor identifier per cell.
#' @param cell_type Cell-type label per cell.
#' @param condition Optional condition label per cell. When supplied, cells are
#'   aggregated separately for each donor-by-condition-by-cell-type combination.
#' @param min_cells Minimum cells required for a pseudobulk.
#' @return List with pseudobulk counts and donor, cell-type and (when supplied)
#'   condition metadata.
#' @examples
#' counts <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), NULL))
#' leviPseudobulk(counts, donor = c("d1", "d1", "d2", "d2"),
#'   cell_type = rep("T", 4), min_cells = 1)
#' @export
leviPseudobulk <- function(counts, donor, cell_type, condition = NULL,
    min_cells = 10L) {
    if (ncol(counts) != length(donor) || length(donor) != length(cell_type))
        stop("donor and cell_type must have one value per cell.", call. = FALSE)
    if (!is.null(condition) && length(condition) != ncol(counts))
        stop("'condition' must have one value per cell.", call. = FALSE)
    if (is.null(rownames(counts))) stop("'counts' needs gene row names.", call. = FALSE)
    if (is.null(condition)) condition <- rep(NA_character_, ncol(counts))
    has_condition <- !all(is.na(condition))
    key <- if (has_condition) paste(donor, condition, cell_type, sep = "::") else
        paste(donor, cell_type, sep = "::")
    groups <- split(seq_along(key), key)
    groups <- groups[lengths(groups) >= min_cells]
    if (!length(groups)) stop("No donor-cell-type pseudobulk passes min_cells.", call. = FALSE)
    pb <- vapply(groups, function(i) Matrix::rowSums(counts[, i, drop = FALSE]),
                 numeric(nrow(counts)))
    if (is.null(dim(pb))) pb <- matrix(pb, ncol = 1L)
    rownames(pb) <- rownames(counts)
    info <- do.call(rbind, strsplit(names(groups), "::", fixed = TRUE))
    if (has_condition) {
        colnames(pb) <- paste0(info[, 1], "_", info[, 2], "_", info[, 3])
        return(list(counts = pb, donor = info[, 1], condition = info[, 2],
                    cell_type = info[, 3], n_cells = lengths(groups)))
    }
    colnames(pb) <- paste0(info[, 1], "_", info[, 2])
    list(counts = pb, donor = info[, 1], cell_type = info[, 2],
         n_cells = lengths(groups))
}

#' Single-cell pseudobulk graph-cluster inference
#'
#' Aggregates cells by donor and cell type, then applies [leviBulkGraphInference]
#' separately to each cell type. Donors, not cells, are permuted.
#' @param counts Genes-by-cells raw counts.
#' @param donor,cell_type Per-cell donor and cell-type labels.
#' @param condition Named condition vector indexed by donor.
#' @param min_cells Minimum cells required for each donor-cell-type pseudobulk;
#'   passed to [leviPseudobulk()].
#' @param ... Arguments forwarded to [leviBulkGraphInference].
#' @return Named list of cell-type-specific graph inference results.
#' @examples
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#' \donttest{
#'   set.seed(1)
#'   genes <- c("HUB", paste0("N", 1:8))
#'   donor <- rep(paste0("d", 1:6), each = 4)
#'   cell_type <- rep(rep(c("A", "B"), each = 2), 6)
#'   counts <- matrix(rpois(length(genes) * length(donor), 30), length(genes),
#'     dimnames = list(genes, NULL))
#'   treated <- donor %in% c("d4", "d5", "d6")
#'   counts[, treated] <- counts[, treated] + 30L
#'   condition <- c(d1 = "C", d2 = "C", d3 = "C", d4 = "T", d5 = "T", d6 = "T")
#'   out <- leviSingleCellGraphInference(counts, donor, cell_type,
#'     condition = condition, min_cells = 2,
#'     networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
#'     fileTypeInput = "dat", test = "T", control = "C", threshold = 1,
#'     n_perm = 3, seed = 1)
#'   names(out)
#' }
#' }
#' @export
leviSingleCellGraphInference <- function(counts, donor, cell_type, condition,
    min_cells = 10L, ...) {
    pb <- leviPseudobulk(counts, donor, cell_type, min_cells = min_cells)
    if (is.null(names(condition))) stop("'condition' must be named by donor.", call. = FALSE)
    types <- unique(pb$cell_type)
    out <- lapply(types, function(type) {
        keep <- pb$cell_type == type
        groups <- unname(condition[pb$donor[keep]])
        if (anyNA(groups) || length(unique(groups)) != 2L) return(NULL)
        leviBulkGraphInference(pb$counts[, keep, drop = FALSE], groups, ...)
    })
    names(out) <- types
    out[!vapply(out, is.null, logical(1))]
}
