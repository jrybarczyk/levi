# Builds the landscape figure. Split out of levi_function() so that the plot
# recipe -- raster, colour scale, the increase/decrease arrows and the optional
# contour lines -- sits in one place instead of in the middle of the pipeline.
#
# The Shiny interface calls it too, which is why the column names, the palette
# and the title are parameters: it keeps the surface under different column
# names, offers a 20-tone "Multicolor" option beside the named palettes, and
# puts the comparison in the page rather than over the figure. Before this was
# shared the two drawings had already drifted -- the scale arrows sat at
# different heights on each side.
#
# Significance contours are added afterwards by .significanceContour(), because
# they only exist when a permutation test ran.
#
# Arguments:
#   landgraphFinal  data.frame with the surface in long form
#   setcolor        palette name for .colorSet(); ignored when colours is given
#   titleChart      title over the figure, or NULL for none
#   contourLevi     TRUE to add the plain contour lines
#   colours         explicit colour vector, overriding setcolor
#   cols            names of the x, y and value columns, in that order
#
# Returns the ggplot object.
.buildLandscapeChart <- function(landgraphFinal, setcolor, titleChart,
                                 contourLevi, colours = NULL,
                                 cols = c("Var1", "Var2", "z")) {
    matrixSize <- sqrt(NROW(landgraphFinal))

    xs <- as.name(cols[1])
    ys <- as.name(cols[2])
    zs <- as.name(cols[3])

    if (is.null(colours)) colours <- .colorSet(setcolor)

    landgraphChart <- ggplot(data = landgraphFinal,
        aes(x = !!xs, y = !!ys)) +
        geom_raster(aes(fill = !!zs), interpolate = TRUE,
                    hjust = 0.5, vjust = 0.5) +
        scale_fill_gradientn(colours = colours,
        values=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
        breaks=seq(0,1,0.2), limits=c(0,1), na.value = "grey92",
        guide = guide_colorbar(title="Landscape score",
        title.position = "right", title.hjust = 0.5,
        title.theme = element_text(angle = 270, size = 9),
        barwidth= 1, barheight = 10)) +
        theme_void() +
        ggtitle(titleChart) +
        theme(plot.title = element_text(
        margin = margin(t = 10, b = -10), hjust = 0.5, lineheight=.8,
        face="bold"), legend.margin=margin(0,0,0,-20)) +
        annotate("text", x = c(matrixSize+2, matrixSize+2),
        y = c(matrixSize*0.42, matrixSize*0.6),
        label = c("lower score", "higher score"), size=3, angle=90) +
        annotate("segment", x = matrixSize*1.07, xend = matrixSize*1.07,
        y = matrixSize*0.49, yend = matrixSize*0.35, colour = "black",
        linewidth=0.2, alpha=0.6, arrow=arrow(type = "closed",
        length = unit(x = c(0.2), units = "cm"))) +
        annotate("segment", x = matrixSize*1.07, xend = matrixSize*1.07,
        y = matrixSize*0.51, yend = matrixSize*0.67, colour = "black",
        linewidth=0.2, alpha=0.6, arrow=arrow(type = "closed",
        length = unit(x = c(0.2), units = "cm"))) +
        coord_fixed(ratio = 1)

    if (isTRUE(contourLevi)) {
        # Drop the background: without this stat_contour warns that it
        # removed the NA cells every time the plot is drawn.
        keep <- !is.na(landgraphFinal[[cols[3]]])
        contour_df <- landgraphFinal[keep, , drop = FALSE]
        if (nrow(contour_df) > 0)
            landgraphChart <- landgraphChart +
                geom_contour(data = contour_df,
                             aes(x = !!xs, y = !!ys, z = !!zs),
                             inherit.aes = FALSE)
    }

    landgraphChart
}

# Builds the interactive 3D surface. Shared with the Shiny interface, which is
# why the palette arrives as a colour vector and the title can be NULL: the
# interface offers a "Multicolor" option that is not one of the named palettes,
# and shows the comparison in the page rather than over the figure.
#
# Significance boundaries are added by .addSignificance3D() when a permutation
# test ran, so both views carry the same information.
#
# Arguments:
#   zmat        the landscape as a matrix, already in [0, 1]
#   colours     colour vector for the surface scale
#   titleChart  title over the figure, or NULL for none
#   pvals       $over and $under matrices, or NULL when no test ran
#   i           reorientation index, as used by the 2D contours
#   sig_level   significance threshold for the boundary
#   perm_side   which sides to draw
#
# Returns the plotly figure, or NULL when plotly is not installed.
.buildSurface3D <- function(zmat, colours, titleChart = NULL, pvals = NULL,
                            i = NULL, sig_level = 0.05, perm_side = "both",
                            camera = NULL) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        warning("plotly is required for 3D visualization. ",
                "Install with: install.packages('plotly')", call. = FALSE)
        return(NULL)
    }

    n_colors   <- length(colours)
    colorscale <- lapply(seq_along(colours), function(idx) {
        list((idx - 1) / (n_colors - 1), colours[idx])
    })

    fig <- plotly::layout(
        plotly::plot_ly(
            z          = zmat,
            type       = "surface",
            colorscale = colorscale,
            cmin = 0, cmax = 1,
            colorbar   = list(title = "Landscape\nscore", len = 0.5)
        ),
        title = titleChart,
        scene = list(
            xaxis = list(title = "", showticklabels = FALSE,
                         showgrid = FALSE, zeroline = FALSE),
            yaxis = list(title = "", showticklabels = FALSE,
                         showgrid = FALSE, zeroline = FALSE),
            zaxis = list(title = "Landscape score", range = c(0, 1)),
            camera = camera
        )
    )
    # The toolbar camera icon exports the surface exactly as it is oriented
    # on screen; these options only raise the resolution of that PNG.
    fig <- plotly::config(fig,
        toImageButtonOptions = list(format = "png", filename = "levi_surface3D",
                                    width = 1600, height = 1200, scale = 2))

    if (!is.null(pvals) && !is.null(i))
        fig <- .addSignificance3D(fig, zmat, pvals, i, sig_level, perm_side)

    fig
}
