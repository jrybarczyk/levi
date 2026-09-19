#' Adjust regional LEVI results across pathways or networks
#'
#' Applies a second multiplicity correction to regional spatial P-values from
#' several independently fitted LEVI networks. Each input result must use
#' regional inference. The within-network P-values already use the maximum
#' regional mass; this function controls the additional family of networks.
#' @param results Named list of `levi_result` objects.
#' @param method Method accepted by [stats::p.adjust()]. Default `"holm"`.
#' @param level Significance threshold.
#' @return The input list with `PFamily` and `FamilySignificant` added to every
#' regional summary.
#' @examples
#' set.seed(1)
#' net <- system.file("extdata", "hub_network.dat", package = "levi")
#' genes <- c("HUB", paste0("N", 1:8))
#' run <- function(lfc) levi(expressionInput = data.frame(ID = genes, logFC = lfc),
#'   networkCoordinatesInput = net, fileTypeInput = "dat", geneSymbolInput = "ID",
#'   readExpColumn = readExpColumn("logFC-logFC"), signal_mode = "logfc",
#'   inference_unit = "region", n_perm = 3, resolutionValueInput = 1,
#'   region_threshold = .02)
#' results <- list(pathway1 = run(rnorm(9, 1)), pathway2 = run(rnorm(9, -1)))
#' adjusted <- leviAdjustPathways(results)
#' adjusted$pathway1$regions$summary
#' @export
leviAdjustPathways <- function(results, method = "holm", level = 0.05) {
    if (!is.list(results) || !length(results))
        stop("'results' must be a non-empty list of LEVI results.", call. = FALSE)
    method <- match.arg(method, stats::p.adjust.methods)
    where <- lapply(results, function(x) {
        if (is.null(x$regions$summary$PSpatial))
            stop("Every result must contain regional PSpatial values.", call. = FALSE)
        seq_len(nrow(x$regions$summary))
    })
    p <- unlist(lapply(results, function(x) x$regions$summary$PSpatial),
                use.names = FALSE)
    adjusted <- stats::p.adjust(p, method = method)
    at <- 0L
    for (i in seq_along(results)) {
        n <- length(where[[i]])
        values <- if (n) adjusted[at + seq_len(n)] else numeric()
        results[[i]]$regions$summary$PFamily <- values
        results[[i]]$regions$summary$FamilySignificant <- values <= level
        results[[i]]$regions$family_inference <- list(method = method,
            level = level, networks = length(results), regions = length(p))
        at <- at + n
    }
    results
}

#' Attribute a LEVI region to nearby network genes
#'
#' Ranks nodes by their Gaussian support inside each detected region, weighted
#' by the node's departure from the neutral score. This is an interpretation
#' score and not a gene-level P-value.
#' @param result A `levi_result`.
#' @param top_n Maximum genes returned per region; `Inf` returns all.
#' @return Data frame with regional gene contribution scores and ranks.
#' @examples
#' set.seed(1)
#' genes <- c("HUB", paste0("N", 1:8))
#' x <- levi(expressionInput = data.frame(ID = genes, logFC = rnorm(9, 1)),
#'   networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
#'   fileTypeInput = "dat", geneSymbolInput = "ID",
#'   readExpColumn = readExpColumn("logFC-logFC"), signal_mode = "logfc",
#'   inference_unit = "region", n_perm = 3, resolutionValueInput = 1,
#'   region_threshold = .02)
#' leviRegionGenes(x, top_n = 3)
#' @export
leviRegionGenes <- function(result, top_n = 20L) {
    cells <- result$regions$cells
    coords <- result$metadata$node_coordinates
    signal <- result$metadata$node_signal
    nodes <- as.character(result$metadata$nodes[, 1])
    grid <- result$metadata$grid
    if (is.null(signal)) stop("The result does not store node signals; rerun levi().",
                              call. = FALSE)
    if (!nrow(cells)) return(data.frame())
    sigma_coord <- grid$sigma * grid$increase
    pieces <- lapply(split(cells, cells$Region), function(region) {
        xy <- as.matrix(region[, c("X", "Y")])
        support <- vapply(seq_len(nrow(coords)), function(i) {
            d2 <- (xy[, 1] - coords[i, 1])^2 + (xy[, 2] - coords[i, 2])^2
            sum(exp(-d2 / (2 * sigma_coord^2)))
        }, numeric(1))
        contribution <- support * abs(signal - 0.5)
        out <- data.frame(Region = region$Region[1], Gene = nodes,
            NodeSignal = signal, Contribution = contribution,
            stringsAsFactors = FALSE)
        out <- out[order(out$Contribution, decreasing = TRUE), ]
        out$Rank <- seq_len(nrow(out))
        if (is.finite(top_n)) out <- head(out, as.integer(top_n))
        out
    })
    ans <- do.call(rbind, pieces)
    rownames(ans) <- NULL
    ans
}

#' Regional inference from biological replicates
#'
#' Permutes sample labels, recalculates gene log fold-changes, reconstructs the
#' complete LEVI landscape and redetects regions at every iteration. This tests
#' treatment association across biological replicates rather than node-label
#' localisation on one fixed effect vector.
#' @param expression Numeric genes-by-samples matrix on a log scale.
#' @param groups Two-level condition vector, one value per sample.
#' @param gene_ids Gene identifiers; defaults to row names.
#' @param test,control Levels defining the `test - control` contrast.
#' @param blocks Optional permutation blocks (e.g. donor identifiers).
#' @param n_perm Number of label permutations.
#' @param sig_level Significance level.
#' @param seed Optional RNG seed.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling how
#'   the permutations are evaluated. The default runs them serially.
#' @param permutation_method `"auto"`, `"monte_carlo"` or `"exact"`.
#'   `"exact"` enumerates every distinct label assignment; `"monte_carlo"`
#'   draws `n_perm` random permutations; `"auto"` uses the exact scheme when
#'   the number of assignments is at most `max_exact` and Monte Carlo otherwise.
#' @param max_exact Maximum number of label assignments enumerated by the exact
#'   permutation scheme.
#' @param ... Network and landscape arguments passed to [levi()]. Do not supply
#' expression, comparison, signal mode, permutation, or inference arguments.
#' @return A `levi_result` whose regional P-values are based on sample labels.
#' @examples
#' set.seed(1)
#' genes <- c("HUB", paste0("N", 1:8))
#' # Four replicates per group: 70 label arrangements, enumerated exactly.
#' expression <- matrix(rnorm(length(genes) * 8, sd = .1), length(genes), 8,
#'   dimnames = list(genes, NULL))
#' expression[, 5:8] <- expression[, 5:8] + 1
#' out <- leviReplicateInference(expression, rep(c("C", "T"), each = 4),
#'   test = "T", control = "C", seed = 2,
#'   networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
#'   fileTypeInput = "dat", resolutionValueInput = 1, region_threshold = .02)
#' out$regions$summary
#' @export
leviReplicateInference <- function(expression, groups, gene_ids = rownames(expression),
    test = unique(groups)[2], control = unique(groups)[1], blocks = NULL,
    n_perm = 999L, sig_level = 0.05, seed = NULL,
    permutation_method = c("auto", "monte_carlo", "exact"),
    max_exact = 50000L, BPPARAM = SerialParam(), ...) {
    expression <- as.matrix(expression)
    storage.mode(expression) <- "double"
    if (ncol(expression) != length(groups))
        stop("'groups' must have one value per expression column.", call. = FALSE)
    if (is.null(gene_ids) || length(gene_ids) != nrow(expression))
        stop("'gene_ids' must have one value per expression row.", call. = FALSE)
    groups <- as.character(groups)
    if (!setequal(unique(groups), c(as.character(test), as.character(control))))
        stop("'groups' must contain exactly the requested test and control levels.", call. = FALSE)
    if (is.null(blocks)) blocks <- rep("all", length(groups))
    if (length(blocks) != length(groups))
        stop("'blocks' must have one value per sample.", call. = FALSE)
    n_perm <- as.integer(n_perm)
    if (!is.finite(n_perm) || n_perm < 1L)
        stop("'n_perm' must be a positive integer.", call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)
    effect <- function(g) rowMeans(expression[, g == test, drop = FALSE], na.rm = TRUE) -
        rowMeans(expression[, g == control, drop = FALSE], na.rm = TRUE)
    dots <- list(...)
    required <- c("networkCoordinatesInput", "fileTypeInput")
    if (!all(required %in% names(dots)))
        stop("Supply 'networkCoordinatesInput' and 'fileTypeInput' in ...",
             call. = FALSE)
    network_edges <- dots$networkInteractionsInput
    if (is.null(network_edges)) network_edges <- NA
    parsed <- .parseNetwork(dots$networkCoordinatesInput, network_edges,
                            dots$fileTypeInput)
    run <- function(lfc) do.call(levi, c(list(
        expressionInput = data.frame(.Gene = gene_ids, logFC = lfc),
        geneSymbolInput = ".Gene", readExpColumn = readExpColumn("logFC-logFC"),
        signal_mode = "logfc", n_perm = 0L, inference_unit = "region",
        sig_level = sig_level, plot3d = FALSE, .parsed_network = parsed,
        .draw = FALSE), dots))
    observed <- run(effect(groups))
    perms <- .labelPermutations(groups, blocks, n_perm, permutation_method,
                                max_exact, alpha = sig_level)
    node_names <- as.character(observed$metadata$nodes[, 1])
    edge_index <- observed$metadata$edges
    node_coord <- observed$metadata$node_coordinates
    coord <- rbind(node_coord,
        (node_coord[edge_index[, 1], , drop = FALSE] +
         node_coord[edge_index[, 2], , drop = FALSE]) / 2)
    grid <- observed$metadata$grid
    k <- if (is.null(dots$logfc_k)) 1 else dots$logfc_k
    project <- function(lfc) {
        collapsed <- tapply(lfc, gene_ids, mean, na.rm = TRUE)
        node_lfc <- unname(collapsed[node_names])
        .networkSignals(cbind(node_lfc, node_lfc), edge_index, TRUE,
                        "logfc", k)$signal
    }
    null <- .permutationNull(perms$labels, function(perm) {
        signal <- project(effect(perm))
        z <- landscape_gauss(coord = coord, SignalOut = signal,
            signalExp = signal, signalCtrl = signal,
            resolutionValue = grid$resolution, zoomValue = grid$zoom,
            increase = grid$increase, sigma = grid$sigma,
            occFrac = grid$occupancy,
            weights = observed$metadata$support_weights %||% numeric(0))$m1
        .maximumRegionMass(z, grid$increase,
            observed$metadata$region_threshold,
            observed$metadata$region_min_cells)
    }, BPPARAM = BPPARAM, names = c("over", "under"))
    observed$regions <- .regionalPvalues(observed$regions, null)
    observed$regions$summary$Significant <-
        observed$regions$summary$PSpatial <= sig_level
    observed$regions$inference$method <- "sample_label"
    observed$regions$inference$scope <- "biological replicate label permutation"
    observed$metadata$inference <- "sample-label randomisation; maximum regional mass"
    observed$metadata$n_perm <- nrow(null)
    observed$metadata$permutation_exact <- perms$exact
    observed$metadata$possible_permutations <- perms$possible
    observed <- .annotateReplicateRegions(observed, sig_level)
    observed
}

.annotateReplicateRegions <- function(result, sig_level = 0.05) {
    summary <- result$regions$summary
    if (!nrow(summary)) return(result)
    # Remove descriptive labels created before the replicate P-values existed.
    is_label <- vapply(result$plot$layers, function(layer)
        is.data.frame(layer$data) && "Label" %in% names(layer$data), logical(1))
    result$plot$layers <- result$plot$layers[!is_label]
    n <- result$metadata$grid$resolution
    selected <- summary$Region[summary$PSpatial <= sig_level]
    labels <- summary[summary$Region %in% selected, , drop = FALSE]
    labels$Label <- sprintf("%s\np = %.3f", labels$Region, labels$PSpatial)
    boundaries <- .regionBoundaries(result$regions, n, selected)
    if (nrow(boundaries)) result$plot <- result$plot +
        ggplot2::geom_segment(data = boundaries,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            colour = "white", linewidth = .9, inherit.aes = FALSE)
    if (nrow(labels)) result$plot <- result$plot + ggplot2::geom_label(
        data = labels,
        ggplot2::aes(x = PeakRow, y = n + 1L - PeakCol, label = Label),
        size = 3, inherit.aes = FALSE)
    result
}

.graphRegions <- function(statistic, node_names, edge_index, threshold) {
    statistic <- as.numeric(statistic[node_names])
    graph <- igraph::graph_from_edgelist(edge_index, directed = FALSE)
    graph <- igraph::add_vertices(graph, max(0, length(node_names) - igraph::vcount(graph)))
    one_side <- function(direction) {
        keep <- if (direction == "over") statistic >= threshold else
            statistic <= -threshold
        ids <- which(is.finite(statistic) & keep)
        if (!length(ids)) return(NULL)
        sub <- igraph::induced_subgraph(graph, ids)
        cmp <- igraph::components(sub)
        lapply(split(seq_along(ids), cmp$membership), function(pos) {
            members <- ids[pos]
            value <- statistic[members]
            peak <- if (direction == "over") which.max(value) else which.min(value)
            data.frame(Direction = direction, Nodes = length(members),
                Mass = sum(abs(value)), PeakStatistic = value[peak],
                PeakGene = node_names[members[peak]],
                Genes = paste(node_names[members], collapse = ";"),
                stringsAsFactors = FALSE)
        })
    }
    pieces <- c(one_side("over"), one_side("under"))
    if (!length(pieces)) return(data.frame(Region = character(), Direction = character(),
        Nodes = integer(), Mass = numeric(), PeakStatistic = numeric(),
        PeakGene = character(), Genes = character()))
    out <- do.call(rbind, pieces)
    out$Region <- sprintf("%s_%02d", out$Direction,
                          ave(seq_len(nrow(out)), out$Direction, FUN = seq_along))
    out[, c("Region", setdiff(names(out), "Region"))]
}

#' Graph-connected cluster-mass inference from biological replicates
#'
#' Fits a moderated `limma` t-statistic for every gene, thresholds the signed
#' statistics, forms components using biological network edges, and assigns
#' FWER-controlled P-values by permuting sample labels. The graph components,
#' rather than adjacency on the 2D layout grid, are the inferential regions.
#' @param expression Log-scale genes-by-samples expression matrix.
#' @param groups Two-level sample-group vector.
#' @param networkCoordinatesInput,networkInteractionsInput,fileTypeInput LEVI
#' network inputs.
#' @param test,control Group labels defining test minus control.
#' @param threshold Absolute moderated-t threshold for cluster formation.
#' @param n_perm Number of sample-label permutations.
#' @param blocks Optional exchangeability blocks.
#' @param seed Optional RNG seed.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling how
#'   the permutations are evaluated. The default runs them serially.
#' @param permutation_method `"auto"`, `"monte_carlo"` or `"exact"`.
#'   `"exact"` enumerates every distinct label assignment; `"monte_carlo"`
#'   draws `n_perm` random permutations; `"auto"` uses the exact scheme when
#'   the number of assignments is at most `max_exact` and Monte Carlo otherwise.
#' @param max_exact Maximum number of label assignments enumerated by the exact
#'   permutation scheme.
#' @return List with observed moderated statistics, graph regions and null maxima.
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   genes <- c("HUB", paste0("N", 1:8))
#'   # Five paired donors: 32 within-donor label swaps, enumerated exactly.
#'   expression <- matrix(rnorm(length(genes) * 10, sd = .1), length(genes), 10,
#'     dimnames = list(genes, NULL))
#'   expression[, 6:10] <- expression[, 6:10] + 2
#'   out <- leviGraphClusterInference(expression, rep(c("C", "T"), each = 5),
#'     system.file("extdata", "hub_network.dat", package = "levi"),
#'     test = "T", control = "C", threshold = 1,
#'     blocks = rep(paste0("d", 1:5), 2), seed = 1)
#'   out$regions$summary
#' }
#' @export
leviGraphClusterInference <- function(expression, groups,
    networkCoordinatesInput, networkInteractionsInput = NA,
    fileTypeInput = "dat", test = unique(groups)[2], control = unique(groups)[1],
    threshold = 2, n_perm = 999L, blocks = NULL, seed = NULL,
    permutation_method = c("auto", "monte_carlo", "exact"), max_exact = 50000L,
    BPPARAM = SerialParam()) {
    if (!requireNamespace("limma", quietly = TRUE))
        stop("'limma' is required for moderated-t graph inference.", call. = FALSE)
    expression <- as.matrix(expression)
    if (is.null(rownames(expression)))
        stop("'expression' must have gene identifiers as row names.", call. = FALSE)
    groups <- as.character(groups)
    if (ncol(expression) != length(groups))
        stop("'groups' must have one value per expression column.", call. = FALSE)
    if (!setequal(unique(groups), c(as.character(test), as.character(control))))
        stop("'groups' must contain exactly test and control.", call. = FALSE)
    has_blocks <- !is.null(blocks)
    if (is.null(blocks)) blocks <- rep("all", length(groups))
    if (length(blocks) != length(groups)) stop("Invalid 'blocks'.", call. = FALSE)
    if (!is.null(seed)) withr::local_seed(seed)
    parsed <- .parseNetwork(networkCoordinatesInput, networkInteractionsInput,
                            fileTypeInput)
    nodes <- unique(as.character(parsed$nodes[, 1]))
    edge_index <- cbind(match(as.character(parsed$edges[, 1]), nodes),
                        match(as.character(parsed$edges[, 2]), nodes))
    edge_index <- edge_index[stats::complete.cases(edge_index), , drop = FALSE]
    statistic <- function(g) {
        if (has_blocks) {
            design <- stats::model.matrix(~ factor(blocks) +
                factor(g, levels = c(control, test)))
            fit <- limma::eBayes(limma::lmFit(expression, design))
            return(setNames(fit$t[, ncol(design)], rownames(expression)))
        }
        design <- stats::model.matrix(~ 0 + factor(g, levels = c(control, test)))
        colnames(design) <- c("control", "test")
        fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expression, design),
            contrasts = cbind(test_minus_control = c(-1, 1))))
        setNames(fit$t[, 1], rownames(expression))
    }
    observed_stat <- statistic(groups)
    summary <- .graphRegions(observed_stat, nodes, edge_index, threshold)
    observed <- list(summary = summary)
    perms <- .labelPermutations(groups, blocks, n_perm, permutation_method, max_exact)
    null <- .permutationNull(perms$labels, function(perm)
        .directionMaxima(.graphRegions(statistic(perm), nodes, edge_index,
                                       threshold)),
        BPPARAM = BPPARAM, names = c("over", "under"))
    observed <- .regionalPvalues(observed, null)
    observed$summary$Significant <- observed$summary$PSpatial <= .05
    list(statistic = data.frame(Gene = names(observed_stat), T = as.numeric(observed_stat)),
         regions = observed, threshold = threshold, n_perm = nrow(null), exact = perms$exact,
         method = "limma moderated t; graph-connected maximum cluster mass; sample labels")
}
