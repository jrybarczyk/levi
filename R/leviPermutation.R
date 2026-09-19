# Permutation significance test, kept apart from levi_function() so that the
# calculation can be exercised without building a plot.
#
# The null it samples from: expression values are shuffled among the network
# positions with the layout fixed. This tests spatial association conditional
# on this network and layout; it is not a test between biological replicates.

# Returns a list with the $over and $under p-value matrices, each
# resolutionValue x resolutionValue.
# progress: optional function(i, n) called once per iteration. It exists so the
# Shiny interface can drive its progress bar without keeping its own copy of
# this loop; NULL in script mode.
.permutationPvalues <- function(coord, SignalOut, signalExp, signalCtrl,
                                matrixOut, resolutionValue, zoomValue,
                                increase, sigma, occFrac, n_perm,
                                progress = NULL, node_values = NULL, edge_index = NULL,
                                single_col = FALSE, signal_mode = "ratio",
                                logfc_k = 1, regions = NULL, perm_strata = NULL,
                                weights = numeric(0)) {
    numberCoord <- nrow(coord)

    # Only two counts per cell are ever needed, so they are accumulated as the
    # permutations are drawn. Keeping the whole res x res x n_perm stack (plus
    # a second one replicating matrixOut just to compare against it) cost
    # roughly 880 MB at the settings the documentation recommends for
    # publication -- resolution 100 with n_perm = 1000 -- for a result that
    # fits in two matrices. NA propagates through the sums exactly as it did
    # through rowSums(), so cells outside the silhouette stay NA.
    countOver  <- matrix(0L, resolutionValue, resolutionValue)
    countUnder <- matrix(0L, resolutionValue, resolutionValue)
    regional <- !is.null(regions)
    if (regional) null_max <- matrix(0, n_perm, 2L,
                                     dimnames = list(NULL, c("over", "under")))

    # Only the signals move; the coordinates are the same on every iteration,
    # so the silhouette and the convolution denominator do not need recomputing.
    for (p in seq_len(n_perm)) {
        if (!is.null(progress)) progress(p, n_perm)
        if (is.null(node_values)) {
            permIdx <- sample(seq_len(numberCoord))
            signals <- list(signal = SignalOut[permIdx, , drop = FALSE],
                            test = signalExp[permIdx, , drop = FALSE],
                            control = signalCtrl[permIdx, , drop = FALSE])
        } else {
            # Preserve the missingness pattern; only measured genes exchange labels.
            measured <- which(rowSums(!is.finite(node_values)) == 0L)
            permuted <- node_values
            if (is.null(perm_strata)) {
                permuted[measured, ] <- node_values[
                    measured[sample.int(length(measured))], , drop = FALSE]
            } else {
                for (idx in split(measured, perm_strata[measured]))
                    permuted[idx, ] <- node_values[sample(idx), , drop = FALSE]
            }
            signals <- .networkSignals(permuted, edge_index, single_col,
                                        signal_mode, logfc_k)
        }
        mf_p <- landscape_gauss(
            coord           = coord,
            SignalOut       = signals$signal,
            signalExp       = signals$test,
            signalCtrl      = signals$control,
            resolutionValue = resolutionValue,
            zoomValue       = zoomValue,
            increase        = increase,
            sigma           = sigma,
            occFrac         = occFrac,
            weights         = weights)

        if (regional) {
            # Redetect components in every randomisation: the observed areas
            # must not be treated as regions specified before seeing the data.
            null_max[p, ] <- .maximumRegionMass(mf_p$m1, increase,
                regions$threshold, regions$min_cells)
        } else {
            countOver  <- countOver  + (mf_p$m1 >= matrixOut)
            countUnder <- countUnder + (mf_p$m1 <= matrixOut)
        }
    }

    if (regional) return(.regionalPvalues(regions, null_max))

    # (k + 1)/(n + 1): avoids a p-value of exactly 0, which no finite number of
    # permutations can support.
    list(
        over  = (countOver  + 1) / (n_perm + 1),
        under = (countUnder + 1) / (n_perm + 1)
    )
}

# Same eight-connectivity and excess mass as .extractLandscapeRegions, using
# igraph's compiled component search; null draws need no cell/centroid tables.
.maximumRegionMass <- function(z, increase, threshold, min_cells) {
    vapply(c("over", "under"), function(side) {
        excess <- if (side == "over") z - .5 - threshold else .5 - z - threshold
        mask <- is.finite(z) & if (side == "over") z >= .5 + threshold else
            z <= .5 - threshold
        at <- which(mask, arr.ind = TRUE)
        if (!nrow(at)) return(0)
        ids <- matrix(0L, nrow(z), ncol(z)); ids[mask] <- seq_len(nrow(at))
        edges <- list()
        offsets <- rbind(c(1, 0), c(0, 1), c(1, 1), c(1, -1))
        for (k in seq_len(4L)) {
            rr <- at[, 1] + offsets[k, 1]; cc <- at[, 2] + offsets[k, 2]
            keep <- rr >= 1 & rr <= nrow(z) & cc >= 1 & cc <= ncol(z)
            target <- ids[cbind(rr[keep], cc[keep])]
            edges[[k]] <- cbind(which(keep)[target > 0], target[target > 0])
        }
        graph <- igraph::make_empty_graph(nrow(at), directed = FALSE)
        el <- do.call(rbind, edges)
        if (nrow(el)) graph <- igraph::add_edges(graph, as.vector(t(el)))
        cmp <- igraph::components(graph)
        mass <- rowsum(pmax(excess[mask], 0), cmp$membership, reorder = FALSE)
        # Align explicitly; do not depend on component traversal order.
        sizes <- cmp$csize[as.integer(rownames(mass))]
        max(c(0, mass[sizes >= min_cells, 1])) * increase^2
    }, numeric(1))
}

# Maximum mass over BOTH directions controls the search across the grid under
# the global exchangeable node-label null. It does not certify each pixel,
# strong control under arbitrary partial nulls, or biological replication.
.regionalPvalues <- function(regions, null_max) {
    joint <- pmax(null_max[, "over"], null_max[, "under"])
    regions$summary$PSpatial <- vapply(regions$summary$Mass, function(mass) {
        tolerance <- 1e-12 * max(1, abs(mass))
        (1 + sum(joint >= mass - tolerance)) / (length(joint) + 1)
    }, numeric(1))
    regions$null_max_mass <- data.frame(null_max, both = joint)
    regions$inference <- list(method = "node_label", statistic = "maximum regional mass",
        correction = "maximum over all regions and both directions",
        scope = "global spatial null conditional on network/layout; not biological replication",
        n_perm = length(joint))
    regions
}

# Exact grid-cell boundary segments, including regions touching the silhouette.
# Plot-space y is reversed in the native landscape.
.regionBoundaries <- function(regions, n, selected = regions$summary$Region) {
    cells <- regions$cells[regions$cells$Region %in% selected, , drop = FALSE]
    empty <- data.frame(x = numeric(), y = numeric(), xend = numeric(),
                        yend = numeric(), Region = character())
    if (!nrow(cells)) return(empty)
    pieces <- lapply(split(cells, cells$Region), function(d) {
        r <- d$MatrixRow; c <- d$MatrixCol
        keys <- paste(r, c)
        x <- r; y <- n + 1 - c
        segments <- list()
        offsets <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
        for (k in seq_len(4L)) {
            keep <- !paste(r + offsets[k, 1], c + offsets[k, 2]) %in% keys
            if (k <= 2) {
                xx <- x[keep] + offsets[k, 1] / 2
                segments[[k]] <- data.frame(x = xx, xend = xx,
                    y = y[keep] - .5, yend = y[keep] + .5)
            } else {
                yy <- y[keep] - offsets[k, 2] / 2
                segments[[k]] <- data.frame(x = x[keep] - .5,
                    xend = x[keep] + .5, y = yy, yend = yy)
            }
        }
        out <- do.call(rbind, segments)
        out$Region <- rep(d$Region[1], nrow(out))
        out
    })
    do.call(rbind, pieces)
}

# Adds one significance contour to the landscape. Called once per side, since
# the two differ only in the matrix, the line style and the wording.
#
# Returns the chart, unchanged when there is no contour to draw.
.significanceContour <- function(chart, pvals, i, sig_level, linetype, side) {
    df <- as.data.frame(melt(pvals[i, rev(i)], value.name = "pval"))

    # Outside the silhouette the p-value is NA; without dropping those rows
    # stat_contour warns on every call that it removed them.
    df <- df[!is.na(df$pval), , drop = FALSE]

    # With no values on both sides of sig_level there is no contour to draw;
    # warn instead of letting ggplot complain.
    #
    # The comparison has to be strict. p-values live on the lattice of
    # multiples of 1/(n_perm + 1), so sig_level often coincides with a value
    # exactly: with n_perm = 39 and sig_level = 0.05 the whole under side sits
    # at or above 0.05 with nothing below it. min() > sig_level was FALSE, the
    # guard let it through, and stat_contour generated nothing -- which reached
    # the user as "Zero contours were generated" plus two range warnings from
    # inside ggplot, none of them naming the cause.
    if (nrow(df) == 0 ||
        !any(df$pval < sig_level) || !any(df$pval > sig_level)) {
        message("No ", side, " contour at sig_level = ", sig_level,
                ": all p-values fall on one side of the threshold. ",
                "No significant boundary is drawn.")
        return(chart)
    }

    chart + geom_contour(
        data        = df,
        aes(x = Var1, y = Var2, z = pval),
        breaks      = sig_level,
        colour      = "white",
        linetype    = linetype,
        linewidth   = 0.6,
        inherit.aes = FALSE
    )
}

# Traces the significance boundary onto the 3D surface.
#
# The 2D figure gets the boundary as contour lines; without this the 3D view
# drops the only information saying which part of the landscape holds up
# against the null, and nothing on the figure admits the omission.
#
# contourLines() gives the p = sig_level isolines in matrix coordinates. Each
# vertex is then lifted to the height of the surface underneath, so the
# boundary lies on the landscape instead of floating over it.
#
# Returns the figure, unchanged when there is no boundary to draw.
.addSignificance3D <- function(fig, zmat, pvals, i, sig_level, perm_side) {
    sides <- list()
    if (perm_side %in% c("both", "over"))
        sides[["over-expressed"]] <- list(m = pvals$over[i, rev(i)],
                                          dash = "solid")
    if (perm_side %in% c("both", "under"))
        sides[["under-expressed"]] <- list(m = pvals$under[i, rev(i)],
                                           dash = "dot")

    for (nm in names(sides)) {
        m <- sides[[nm]]$m

        # Same condition the 2D contour uses: with every p-value on one side of
        # the threshold there is no isoline to draw.
        if (all(is.na(m)) || min(m, na.rm = TRUE) > sig_level ||
            max(m, na.rm = TRUE) < sig_level) next

        lines <- grDevices::contourLines(x = seq_len(nrow(m)),
                                         y = seq_len(ncol(m)),
                                         z = m, levels = sig_level)
        first <- TRUE
        for (ln in lines) {
            r <- pmin(pmax(round(ln$x), 1L), nrow(zmat))
            cc <- pmin(pmax(round(ln$y), 1L), ncol(zmat))

            # A small lift keeps the line from being swallowed by the surface
            # it sits on; the axis is capped at 1, so clamp there.
            zv <- pmin(zmat[cbind(r, cc)] + 0.015, 1)

            # inherit = FALSE: without it the trace picks up the surface's
            # colorscale and colorbar, which scatter3d does not accept.
            fig <- plotly::add_trace(fig, inherit = FALSE,
                x = cc - 1L, y = r - 1L, z = zv,
                type = "scatter3d", mode = "lines",
                line = list(color = "white", width = 4,
                            dash = sides[[nm]]$dash),
                name = paste0("p <= ", sig_level, ", ", nm),
                legendgroup = nm, showlegend = first,
                hoverinfo = "name")
            first <- FALSE
        }
    }
    fig
}
