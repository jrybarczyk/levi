# Null-coalescing helper (backport for R < 4.4)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# -- Signal transformation engine ---------------------------------------------
#
# .computeSignalOut(): converts raw expression values to a landscape score
# in [0, 1] according to the chosen signal_mode.
#
# Arguments
#   signalExp    matrix (n x 1) - test expression (possibly log-scale)
#   signalCtrl   matrix (n x 1) - control expression, OR matrix(1, n, 1) for
#                single-column mode
#   expressionLog retained for internal compatibility; back-transformation
#                is performed upstream, and only for ratio mode.
#   signal_mode  one of "ratio", "logfc", "zscore"
#   logfc_k      steepness of the sigmoid (logfc mode); ignored otherwise
#   single_col   logical - TRUE when Test == Control (one-column mode)
#
# Supported data-type guide:
#   Raw counts / TPM / FPKM / linear LFQ -> signal_mode = "ratio"
#   Log2-normalised microarray (RMA), VST/rlog DESeq2, log2 Proteomics
#       with separate test & control columns -> signal_mode = "logfc",
#       expressionLog = FALSE  (keep log-scale; function computes Test-Control)
#   DESeq2 log2FoldChange only, edgeR logFC only, Seurat avg_log2FC,
#       scRNA-seq FindMarkers single-column -> signal_mode = "logfc",
#       single-column mode (readExpColumn("logFC-logFC"))
#   Relative position within a logFC distribution ->
#       signal_mode = "zscore"
#
.computeSignalOut <- function(signalExp, signalCtrl,
                               expressionLog = FALSE,
                               signal_mode   = "ratio",
                               logfc_k       = 1,
                               single_col    = FALSE) {
    # -- ratio (current default) -------------------------------------------
    if (signal_mode == "ratio") {
        # Test/(Test + Control) is only meaningful on a non-negative linear
        # scale. With signed values it maps a strongly down-regulated gene onto
        # the same score as a strongly up-regulated one, so a table of logFCs
        # run through the default mode comes back inverted, and nothing about
        # the figure says so.
        if (any(signalExp < 0, na.rm = TRUE) ||
            any(signalCtrl < 0, na.rm = TRUE)) {
            warning("signal_mode = \"ratio\" expects non-negative values ",
                "(raw counts, TPM, FPKM), but the expression data contain ",
                "negative ones. Test/(Test + Control) is not meaningful on a ",
                "signed scale and can invert the landscape. For log-scale ",
                "data or a ready-made logFC, use signal_mode = \"logfc\".",
                call. = FALSE)
        }

        out <- signalExp / (signalExp + signalCtrl)

        # Test + Control == 0 gives +/-Inf, which would leave the score outside
        # the [0, 1] the whole package is built on. Treat those as undefined.
        nonFinite <- !is.finite(out)
        if (any(nonFinite)) {
            warning(sum(nonFinite), " value(s) have Test + Control == 0, ",
                "leaving Test/(Test + Control) undefined. They were set to ",
                "the neutral score 0.5.", call. = FALSE)
            out[nonFinite] <- NA_real_
        }

        # Keep the absolute ratio: min-max scaling moves unchanged genes
        # away from 0.5 and makes scores depend on unrelated genes.
        out[is.na(out)] <- 0.5
        return(.clampScore(out))
    }

    # -- derive log fold-change --------------------------------------------
    # A single column already contains logFC. Two columns contain log-scale
    # expression, regardless of expressionLog (which only affects ratio).
    if (single_col) {
        lfc <- as.numeric(signalExp)
    } else {
        lfc <- as.numeric(signalExp) - as.numeric(signalCtrl)
    }

    if (signal_mode == "logfc") {
        # sigmoid: maps 0 -> 0.5 (neutral), +Inf -> 1, -Inf -> 0
        out <- matrix(1 / (1 + exp(-logfc_k * lfc)), ncol = 1)
        return(.clampScore(out))
    }

    # -- zscore (pnorm) ----------------------------------------------------
    # Maps the distribution of logFC values to [0,1] via the normal CDF.
    # Average gene -> 0.5; over-expressed tail -> near 1; under-expressed -> 0.
    # This is mean/SD standardisation, not an outlier-robust estimator.
    z   <- (lfc - mean(lfc, na.rm = TRUE)) /
           (stats::sd(lfc,  na.rm = TRUE) + 1e-10)
    .clampScore(matrix(stats::pnorm(z), ncol = 1))
}

# The landscape score is documented as living in [0, 1] with 0.5 meaning "no
# change", and the normalised convolution only preserves that if the signals
# handed to it already respect it. This is the single place where every mode
# passes through, so the contract is enforced here: anything undefined becomes
# the neutral point and the rest is clamped.
.clampScore <- function(x) {
    x[!is.finite(x)] <- 0.5
    pmin(pmax(x, 0), 1)
}

# .rescale_for_cpp(): normalises an expression matrix to [0.05, 0.95] so it
# is safe to pass to the C++ routines regardless of the original scale.
.rescale_for_cpp <- function(x) {
    rng <- range(x, finite = TRUE)
    if (diff(rng) > 0) {
        (as.matrix(x) - rng[1]) / diff(rng) * 0.90 + 0.05
    } else {
        matrix(0.5, nrow = nrow(x), ncol = ncol(x))
    }
}

#' @title leviGrid
#' @description Arrange multiple \code{levi()} results side by side for visual
#' comparison. Each panel shows one comparison's landscape at the same colour
#' scale.
#' @param results A \code{levi} result or a \strong{list} of \code{levi}
#' results (e.g., the output of a batch call with multiple
#' \code{readExpColumn()} comparisons).
#' @param ncol Integer. Number of columns in the grid. Default is \code{2}.
#' @param titles Character vector of panel titles. If \code{NULL} (default),
#' uses the \code{comparison} field from each result.
#' @param ... Additional arguments passed to the layout engine.
#' @return Invisibly returns the arranged plot object (class depends on the
#' available package: \code{patchwork}, \code{cowplot}, or \code{gridExtra}).
#' @details
#' Requires one of the following packages (checked in order):
#' \code{patchwork}, \code{cowplot}, or \code{gridExtra}. Install any one with
#' \code{install.packages()}.
#' @examples
#' hub_n <- system.file("extdata", "hub_network.dat", package = "levi")
#' hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")
#' res <- levi(expressionInput         = hub_e,
#'             networkCoordinatesInput = hub_n,
#'             fileTypeInput           = "dat",
#'             geneSymbolInput         = "ID",
#'             readExpColumn           = readExpColumn("Test-Control"),
#'             resolutionValueInput    = 10,
#'             smoothValueInput        = 5)
#' leviGrid(res)
#' \donttest{
#' # Batch mode: one landscape per comparison, arranged in a single column
#' multi_e <- system.file("extdata", "hub_multicomp_expression.dat",
#'                        package = "levi")
#' res2 <- levi(expressionInput         = multi_e,
#'              networkCoordinatesInput = hub_n,
#'              fileTypeInput           = "dat",
#'              geneSymbolInput         = "ID",
#'              readExpColumn           = readExpColumn("Cond_A-Cond_C",
#'                                                      "Cond_B-Cond_C"),
#'              resolutionValueInput    = 10,
#'              smoothValueInput        = 5)
#' leviGrid(res2, ncol = 1)
#' }
#' @export
leviGrid <- function(results, ncol = 2L, titles = NULL, ...) {
    # Normalise: single result -> list of 1. The class covers results from
    # levi >= 2.0.0; the $plot test keeps older objects working.
    if (inherits(results, "levi_result") || !is.null(results$plot))
        results <- list(results)

    plots <- lapply(seq_along(results), function(i) {
        r <- results[[i]]
        p <- r$plot
        ttl <- if (!is.null(titles) && length(titles) >= i) {
            titles[[i]]
        } else if (!is.null(r$comparison)) {
            r$comparison
        } else {
            paste0("Result ", i)
        }
        p + ggplot2::ggtitle(ttl)
    })

    if (requireNamespace("patchwork", quietly = TRUE)) {
        out <- patchwork::wrap_plots(plots, ncol = ncol, ...)
        methods::show(out)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        out <- cowplot::plot_grid(plotlist = plots, ncol = ncol, ...)
        methods::show(out)
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        out <- gridExtra::grid.arrange(grobs = plots, ncol = ncol, ...)
    } else {
        message("Install patchwork, cowplot, or gridExtra for side-by-side layout.",
                " Printing plots sequentially.")
        lapply(plots, methods::show)
        out <- plots
    }
    invisible(out)
}

#' @title leviDiff
#' @description Compute a differential landscape by subtracting two
#' \code{levi()} results cell by cell. The resulting map highlights regions
#' where expression changed \strong{between} conditions A and B.
#' @param result_a A single \code{levi()} result (list with
#' \code{$landscape}).
#' @param result_b A single \code{levi()} result to subtract from
#' \code{result_a}.
#' @param label_a Character. Label for condition A. Default: uses
#' \code{result_a$comparison} or \code{"A"}.
#' @param label_b Character. Label for condition B. Default: uses
#' \code{result_b$comparison} or \code{"B"}.
#' @param setcolor Character palette name accepted by \code{levi()}.
#' Default \code{"default"} uses blue-to-red.
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{\code{diff}}{data.frame with columns X, Y, Diff (B minus A, range -1 to 1)}
#'   \item{\code{plot}}{ggplot2 differential landscape object}
#'   \item{\code{comparison}}{character string summarising the subtraction}
#' }
#' @details
#' Both results must have been generated with the \strong{same} network and
#' \strong{same} resolution (\code{resolutionValueInput}). The \code{Diff}
#' value is positive where condition B is more expressed and negative where
#' condition A is more expressed.
#' @examples
#' hub_n <- system.file("extdata", "hub_network.dat", package = "levi")
#' hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")
#' res <- levi(expressionInput         = hub_e,
#'             networkCoordinatesInput = hub_n,
#'             fileTypeInput           = "dat",
#'             geneSymbolInput         = "ID",
#'             readExpColumn           = readExpColumn("Test-Control"),
#'             resolutionValueInput    = 10,
#'             smoothValueInput        = 5)
#' diff_res <- leviDiff(res, res)
#' head(diff_res$diff)
#' \donttest{
#' # Two comparisons from the same network, then their differential landscape
#' multi_e <- system.file("extdata", "hub_multicomp_expression.dat",
#'                        package = "levi")
#' res2 <- levi(expressionInput         = multi_e,
#'              networkCoordinatesInput = hub_n,
#'              fileTypeInput           = "dat",
#'              geneSymbolInput         = "ID",
#'              readExpColumn           = readExpColumn("Cond_A-Cond_C",
#'                                                      "Cond_B-Cond_C"),
#'              resolutionValueInput    = 10,
#'              smoothValueInput        = 5)
#' diff_res <- leviDiff(res2[[1]], res2[[2]])
#' diff_res$plot
#' }
#' @details Both results must contain compatible metadata: signal mode,
#' input scale, network/layout, grid/smoothing, missing genes and versions.
#' Legacy results without metadata must be recalculated. zscore differences
#' compare separately standardised relative positions and issue a warning.
#' @export
leviDiff <- function(result_a, result_b,
                     label_a  = NULL,
                     label_b  = NULL,
                     setcolor = "default") {
    if (is.null(result_a$landscape) || is.null(result_b$landscape))
        stop("Both results must contain a $landscape data.frame. ",
             "Re-run levi() with levi >= 2.0.0.")

    la <- result_a$landscape
    lb <- result_b$landscape

    if (nrow(la) != nrow(lb))
        stop("Landscapes have different sizes (",
             nrow(la), " vs ", nrow(lb), "). ",
             "Use the same resolutionValueInput for both calls.")

    if (is.null(result_a$metadata) || is.null(result_b$metadata))
        stop("Both results need metadata; re-run levi() before leviDiff().")
    compatible <- c("signal_mode", "single_col", "logfc_k", "expressionLog",
                    "nodes", "edges", "grid", "missing_genes", "versions")
    for (field in compatible) {
        if (!isTRUE(all.equal(result_a$metadata[[field]],
                              result_b$metadata[[field]])))
            stop("Incompatible landscapes: metadata field '", field, "' differs.")
    }
    if (result_a$metadata$signal_mode == "zscore")
        warning("zscore differences compare relative positions within separately ",
                "standardised distributions, not absolute effects.", call. = FALSE)
    if (anyDuplicated(la[, c("Var1", "Var2")]) ||
        anyDuplicated(lb[, c("Var1", "Var2")]) ||
        !setequal(paste(la$Var1, la$Var2), paste(lb$Var1, lb$Var2)))
        stop("Landscapes have incompatible grid coordinates.")

    lbl_a <- label_a %||% result_a$comparison %||% "A"
    lbl_b <- label_b %||% result_b$comparison %||% "B"

    # align on X,Y coordinates
    merged <- merge(la[, c("Var1","Var2","z")],
                    lb[, c("Var1","Var2","z")],
                    by = c("Var1","Var2"), suffixes = c("_a","_b"))
    merged$Diff <- merged$z_b - merged$z_a

    title_str <- paste0(lbl_b, " - ", lbl_a)

    diff_colors <- c("#053061","#2166ac","#4393c3","#92c5de",
                     "#d1e5f0","#f7f7f7",
                     "#fddbc7","#f4a582","#d6604d","#b2182b","#67001f")

    p <- ggplot2::ggplot(merged, ggplot2::aes(x = Var1, y = Var2)) +
        ggplot2::geom_raster(ggplot2::aes(fill = Diff),
                             interpolate = TRUE, hjust = 0.5, vjust = 0.5) +
        ggplot2::scale_fill_gradientn(
            colours = diff_colors,
            limits  = c(-1, 1),
            breaks  = seq(-1, 1, 0.5),
            guide   = ggplot2::guide_colorbar(
                title          = "Diff",
                title.position = "right",
                title.hjust    = 0.5,
                title.theme    = ggplot2::element_text(angle = 270, size = 9),
                barwidth = 1, barheight = 10)) +
        ggplot2::theme_void() +
        ggplot2::ggtitle(title_str) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5, face = "bold",
            margin = ggplot2::margin(t = 10, b = -10)),
            legend.margin = ggplot2::margin(0, 0, 0, -20)) +
        ggplot2::coord_fixed(ratio = 1)

    methods::show(p)

    out <- list(
        comparison = title_str,
        diff       = data.frame(X = merged$Var1, Y = merged$Var2,
                                Diff = round(merged$Diff, 4)),
        plot       = p
    )
    invisible(out)
}

# Internal helper: detect 2D peaks and valleys in the landscape matrix
.detectPeaks2D <- function(matrixOut, coord, nnodes, nodesCoord,
                            zoomValue, increase, resolutionValue,
                            radius = 5, min_score = 0.65) {
    n <- resolutionValue
    is_peak   <- matrix(FALSE, n, n)
    is_valley <- matrix(FALSE, n, n)

    # Cells outside the network silhouette are NA and do not compete for extremes.
    for (i in seq(radius + 1, n - radius)) {
        for (j in seq(radius + 1, n - radius)) {
            val <- matrixOut[i, j]
            if (is.na(val)) next
            nb  <- matrixOut[(i - radius):(i + radius),
                             (j - radius):(j + radius)]
            nb  <- nb[!is.na(nb)]
            if (!length(nb)) next
            if (val >= min_score       && val == max(nb)) is_peak[i, j]   <- TRUE
            if (val <= (1 - min_score) && val == min(nb)) is_valley[i, j] <- TRUE
        }
    }

    nodeX     <- coord[seq_len(nnodes), 1]
    nodeY     <- coord[seq_len(nnodes), 2]
    geneNames <- as.character(nodesCoord[, 1])

    .build_table <- function(mask, type) {
        idx <- which(mask, arr.ind = TRUE)
        if (nrow(idx) == 0L) return(NULL)

        px <- zoomValue + (idx[, 1] - 1) * increase
        py <- zoomValue + (idx[, 2] - 1) * increase

        nearest <- vapply(seq_len(nrow(idx)), function(p) {
            geneNames[which.min((nodeX - px[p])^2 + (nodeY - py[p])^2)]
        }, character(1))

        data.frame(
            Type        = type,
            NearestGene = nearest,
            MatrixRow   = idx[, 1],
            MatrixCol   = idx[, 2],
            Score       = round(matrixOut[idx], 4),
            stringsAsFactors = FALSE
        )
    }

    result <- rbind(.build_table(is_peak, "peak"),
                    .build_table(is_valley, "valley"))
    if (is.null(result)) return(data.frame())
    result <- result[order(result$Score, decreasing = TRUE), ]
    # One label per gene (highest-scoring entry wins — no duplicate labels)
    result <- result[!duplicated(result$NearestGene), ]
    rownames(result) <- NULL
    result
}

# Extract coherent spatial territories from a landscape.  This is deliberately
# separate from peak detection: a broad component can be a region even when it
# has no unique local maximum.
.extractLandscapeRegions <- function(matrixOut, zoomValue, increase,
                                     threshold = 0.1, min_cells = 3L) {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold <= 0 || threshold >= 0.5)
        stop("'region_threshold' must be a number strictly between 0 and 0.5.",
             call. = FALSE)
    if (!is.numeric(min_cells) || length(min_cells) != 1L ||
        !is.finite(min_cells) || min_cells < 1 || min_cells != floor(min_cells))
        stop("'region_min_cells' must be a positive integer.", call. = FALSE)

    empty_summary <- data.frame(
        Region = character(), Direction = character(), Cells = integer(),
        Area = numeric(), Mass = numeric(), PeakScore = numeric(),
        PeakRow = integer(), PeakCol = integer(), CentroidX = numeric(),
        CentroidY = numeric(), stringsAsFactors = FALSE)
    empty_cells <- data.frame(
        Region = character(), Direction = character(), MatrixRow = integer(),
        MatrixCol = integer(), X = numeric(), Y = numeric(), Score = numeric(),
        stringsAsFactors = FALSE)

    walk_components <- function(mask, direction) {
        nr <- nrow(mask)
        nc <- ncol(mask)
        visited <- matrix(FALSE, nr, nc)
        summaries <- list()
        cells <- list()
        component <- 0L
        neighbours <- expand.grid(row = seq.int(-1L, 1L), col = seq.int(-1L, 1L))
        neighbours <- neighbours[!(neighbours$row == 0L & neighbours$col == 0L), ]

        starts <- which(mask, arr.ind = TRUE)
        if (!nrow(starts)) return(list(summary = summaries, cells = cells))
        max_cells <- sum(mask)
        queue_r <- integer(max_cells)
        queue_c <- integer(max_cells)
        member_r <- integer(max_cells)
        member_c <- integer(max_cells)

        for (s in seq_len(nrow(starts))) {
            sr <- starts[s, 1]
            sc <- starts[s, 2]
            if (visited[sr, sc]) next

            head <- 1L
            tail <- 1L
            queue_r[tail] <- sr
            queue_c[tail] <- sc
            visited[sr, sc] <- TRUE
            members <- 0L

            while (head <= tail) {
                r <- queue_r[head]
                c <- queue_c[head]
                head <- head + 1L
                members <- members + 1L
                member_r[members] <- r
                member_c[members] <- c
                for (d in seq_len(nrow(neighbours))) {
                    rr <- r + neighbours$row[d]
                    cc <- c + neighbours$col[d]
                    if (rr >= 1L && rr <= nr && cc >= 1L && cc <= nc &&
                        mask[rr, cc] && !visited[rr, cc]) {
                        tail <- tail + 1L
                        queue_r[tail] <- rr
                        queue_c[tail] <- cc
                        visited[rr, cc] <- TRUE
                    }
                }
            }

            if (members < min_cells) next
            component <- component + 1L
            region <- sprintf("%s_%02d", direction, component)
            rows <- member_r[seq_len(members)]
            cols <- member_c[seq_len(members)]
            scores <- matrixOut[cbind(rows, cols)]
            excess <- if (direction == "over") scores - 0.5 - threshold else
                0.5 - scores - threshold
            excess <- pmax(excess, 0)
            weight <- if (sum(excess) > 0) excess else rep(1, length(excess))
            x <- zoomValue + (rows - 1L) * increase
            y <- zoomValue + (cols - 1L) * increase
            peak <- if (direction == "over") which.max(scores) else which.min(scores)

            summaries[[length(summaries) + 1L]] <- data.frame(
                Region = region, Direction = direction, Cells = length(rows),
                Area = length(rows) * increase^2, Mass = sum(excess) * increase^2,
                PeakScore = scores[peak], PeakRow = rows[peak], PeakCol = cols[peak],
                CentroidX = stats::weighted.mean(x, weight),
                CentroidY = stats::weighted.mean(y, weight),
                stringsAsFactors = FALSE)
            cells[[length(cells) + 1L]] <- data.frame(
                Region = region, Direction = direction, MatrixRow = rows,
                MatrixCol = cols, X = x, Y = y, Score = scores,
                stringsAsFactors = FALSE)
        }
        list(summary = summaries, cells = cells)
    }

    over <- walk_components(is.finite(matrixOut) & matrixOut >= 0.5 + threshold,
                            "over")
    under <- walk_components(is.finite(matrixOut) & matrixOut <= 0.5 - threshold,
                             "under")
    summary <- do.call(rbind, c(over$summary, under$summary))
    cells <- do.call(rbind, c(over$cells, under$cells))
    if (is.null(summary)) summary <- empty_summary
    if (is.null(cells)) cells <- empty_cells
    if (nrow(summary)) {
        summary <- summary[order(summary$Mass, decreasing = TRUE), , drop = FALSE]
        rownames(summary) <- NULL
    }
    if (nrow(cells)) rownames(cells) <- NULL
    list(summary = summary, cells = cells, threshold = threshold,
         min_cells = as.integer(min_cells))
}
