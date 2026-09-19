#' @title Print a levi result
#'
#' @description Compact summary of the object returned by \code{\link{levi}}:
#' the comparison, how many nodes were scored, the strongest and weakest genes,
#' the peaks and valleys detected, and whether a permutation test was run.
#'
#' @param x A \code{levi_result} object, as returned by \code{\link{levi}}.
#' @param n Integer. How many genes to list at each end of the ranking.
#'   Default is 3.
#' @param ... Ignored, present for compatibility with the generic.
#'
#' @return \code{x}, invisibly. Called for the summary it prints.
#'
#' @details The landscape plot is drawn by \code{\link{levi}} itself, so this
#' method prints text only and never redraws the figure. Use \code{x$plot} to
#' get the ggplot object back.
#'
#' @seealso \code{\link{levi}}
#'
#' @examples
#' hub_n <- system.file("extdata", "hub_network.dat", package = "levi")
#' hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")
#'
#' res <- levi(networkCoordinatesInput = hub_n,
#'             expressionInput         = hub_e,
#'             fileTypeInput           = "dat",
#'             geneSymbolInput         = "ID",
#'             readExpColumn           = readExpColumn("Test-Control"),
#'             resolutionValueInput    = 10,
#'             smoothValueInput        = 5)
#' res
#'
#' @export
print.levi_result <- function(x, n = 3L, ...) {
    cat("levi landscape:", x$comparison, "\n")

    scores <- x$scores
    n <- min(as.integer(n), nrow(scores))

    cat("  nodes scored: ", nrow(scores), sep = "")
    if (!is.null(x$landscape)) {
        cells <- nrow(x$landscape)
        bg <- sum(is.na(x$landscape$z))
        cat(sprintf("  |  grid: %d cells, %.0f%% background",
            cells, 100 * bg / cells))
    }
    cat("\n")

    if (nrow(scores) > 0) {
        rng <- range(scores$LandscapeScore)
        cat(sprintf("  score range: %.3f to %.3f  (%s)\n",
            rng[1], rng[2], x$metadata$meaning %||% "signal interpretation unavailable"))

        fmt <- function(idx) paste(sprintf("%s (%.3f)",
            scores$Gene[idx], scores$LandscapeScore[idx]), collapse = ", ")

        cat("  highest: ", fmt(seq_len(n)), "\n", sep = "")
        cat("  lowest:  ", fmt(rev(seq(nrow(scores) - n + 1, nrow(scores)))),
            "\n", sep = "")
    }

    peaks <- x$peaks
    if (is.null(peaks) || nrow(peaks) == 0) {
        cat("  peaks: none detected\n")
    } else {
        cat(sprintf("  peaks: %d peak(s), %d valley(s)\n",
            sum(peaks$Type == "peak"), sum(peaks$Type == "valley")))
    }

    regions <- x$regions
    if (is.null(regions) || is.null(regions$summary)) {
        cat("  regions: unavailable (re-run with levi >= 2.0.0)\n")
    } else {
        summary <- regions$summary
        cat(sprintf("  regions: %d over, %d under (threshold %.3f; min %d cells)\n",
            sum(summary$Direction == "over"), sum(summary$Direction == "under"),
            regions$threshold, regions$min_cells))
    }

    if (!is.null(regions$inference)) {
        cat("  regional permutation test: ", sum(regions$summary$Significant),
            " region(s) pass the joint maximum-mass threshold; ",
            regions$inference$n_perm, " permutations\n", sep = "")
    } else if (is.null(x$pvalues)) {
        cat("  permutation test: not run (n_perm = 0)\n")
    } else {
        inMask <- sum(!is.na(x$pvalues$over))
        level <- x$metadata$sig_level %||% 0.05
        adjustment <- x$metadata$p_adjust_method %||% "unknown"
        cat(sprintf(
            "  permutation test: %d of %d masked cells significant over, %d under\n",
            sum(x$pvalues$over <= level, na.rm = TRUE), inMask,
            sum(x$pvalues$under <= level, na.rm = TRUE)))
    }

    if (!is.null(x$pvalues))
        cat("  threshold: ", level, "; adjustment: ", adjustment, "\n", sep = "")
    cat("  fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
    invisible(x)
}
