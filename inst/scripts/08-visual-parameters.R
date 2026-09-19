# =============================================================================
# levi -- 08. Visual parameters
#
# Goal:     understand what each image parameter does, which ones affect only
#           the appearance and which ones change the values, and how to export
#           the figure.
#
# Data:     medusa.dat with expression.dat -- a real network of 30 nodes and
#           325 interactions, where the differences are visible.
#
# Runtime:  about two minutes.
# =============================================================================

library(levi)

ex <- function(f) system.file("extdata", f, package = "levi")

run <- function(resolution = 40, smoothing = 50, contrast = 50, zoom = 50, ...) {
    levi(expressionInput         = ex("expression.dat"),
         networkCoordinatesInput = ex("medusa.dat"),
         fileTypeInput           = "dat",
         geneSymbolInput         = "ID",
         readExpColumn           = readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
         resolutionValueInput    = resolution,
         smoothValueInput        = smoothing,
         contrastValueInput      = contrast,
         zoomValueInput          = zoom, ...)
}

profile <- function(r) sprintf("%.3f to %.3f | background %2.0f%% | peaks %d valleys %d",
    min(r$landscape$z, na.rm = TRUE), max(r$landscape$z, na.rm = TRUE),
    100 * mean(is.na(r$landscape$z)),
    sum(r$peaks$Type == "peak"), sum(r$peaks$Type == "valley"))


# -----------------------------------------------------------------------------
# 1. resolutionValueInput -- the size of the grid
# -----------------------------------------------------------------------------
# Ranges from 1 to 100 and becomes a grid of 30 to 240 cells per side. The cost
# grows with the square of that number.

cat("--- resolution ---\n")
for (rv in c(10, 40, 80)) {
    t <- system.time(r <- run(resolution = rv))
    side <- as.integer((rv / 100) * 210 + 30)
    cat(sprintf("  %3d -> grid %3dx%-3d %s | %.1f s\n",
                rv, side, side, profile(r), t[["elapsed"]]))
}

# The resolution does not change the interpretation, only the detail. It is
# worth exploring at 30-40 and going up to 70-100 only for the final figure.


# -----------------------------------------------------------------------------
# 2. smoothValueInput -- the kernel width
# -----------------------------------------------------------------------------
# Sets sigma, the width of the Gaussian, in grid cells. Since sigma is
# proportional to the resolution, changing the resolution does NOT alter the
# apparent smoothing.

cat("\n--- smoothing ---\n")
for (sv in c(10, 30, 50, 80)) {
    cat(sprintf("  %3d -> %s\n", sv, profile(run(smoothing = sv))))
}

# Little smoothing preserves individual nodes and edges and produces a grainy
# surface; heavy smoothing reveals the regional trend and erases the detail.
#
# The important point: smoothing does NOT change the scale. The extremes move
# towards the neutral point because each cell now averages more neighbours,
# but the value 0.5 still means "no change". Choose sigma by the spatial scale
# of interest, not to "improve the contrast".


# -----------------------------------------------------------------------------
# 3. contrastValueInput -- the clipping of the silhouette
# -----------------------------------------------------------------------------
# Sets the occupancy threshold below which a cell becomes background (NA).

cat("\n--- contrast (silhouette) ---\n")
for (ct in c(5, 50, 95)) {
    cat(sprintf("  %3d -> %s\n", ct, profile(run(contrast = ct))))
}

# High contrast glues the silhouette to the network and can fragment sparse
# regions into islands; low contrast produces a continuous patch extending
# beyond the nodes. Note that the range of values barely changes: the
# parameter clips, it does not rescale.


# -----------------------------------------------------------------------------
# 4. zoomValueInput -- the framing
# -----------------------------------------------------------------------------
cat("\n--- zoom ---\n")
for (zv in c(0, 50, 100)) {
    cat(sprintf("  %3d -> %s\n", zv, profile(run(zoom = zv))))
}

# Shifts the window over the plane. Useful when the network sits off-centre or
# when a sub-region needs to be cropped.


# -----------------------------------------------------------------------------
# 5. Palettes
# -----------------------------------------------------------------------------
# setcolor = "default" uses the 20-tone scale, from indigo to red. There are
# also two-tone palettes, useful when the audience perceives colour
# differently or when the figure will be printed in black and white.

palettes <- c("default", "terrain", "rainbow", "heat", "topo", "cm",
              "purple_pink", "green_blue", "blue_yellow",
              "pink_green", "orange_purple", "green_marine")
cat("\navailable palettes:\n  ", paste(palettes, collapse = ", "), "\n")

invisible(run(setcolor = "blue_yellow"))

# contourLevi = TRUE adds contour lines over the surface, which helps delimit
# regions in monochrome printing.
invisible(run(setcolor = "green_blue", contourLevi = TRUE))


# -----------------------------------------------------------------------------
# 6. 3D surface
# -----------------------------------------------------------------------------
# plot3d = TRUE generates, besides the 2D map, an interactive surface with
# plotly. It serves exploration; for publication the 2D version usually
# communicates better.

if (requireNamespace("plotly", quietly = TRUE)) {
    cat("\nplotly available: use plot3d = TRUE for the interactive surface.\n")
} else {
    cat("\nplotly not installed; plot3d = TRUE would emit a warning.\n")
}


# -----------------------------------------------------------------------------
# 7. Exporting
# -----------------------------------------------------------------------------
# $plot is an ordinary ggplot object, so anything you do with ggplot2 applies.

r <- run(resolution = 60)

target <- file.path(tempdir(), "landscape.png")
ggplot2::ggsave(target, r$plot, width = 6, height = 5, dpi = 300)
cat("\nfigure saved at:", target, "\n")

# For publication, lossless TIFF:
#   ggplot2::ggsave("fig.tiff", r$plot, width = 6, height = 5, dpi = 300,
#                   compression = "lzw")

# The tables export like any data.frame:
#   write.csv(r$scores, "scores.csv", row.names = FALSE)
#   write.csv(r$peaks,  "peaks.csv",  row.names = FALSE)

# And the plot can be adjusted after it is generated:
#   r$plot + ggplot2::labs(title = "Tumour against Normal",
#                          caption = "STRING network, score >= 700")


# -----------------------------------------------------------------------------
# 8. Summary
# -----------------------------------------------------------------------------
#   - resolution: detail and cost, cost grows quadratically;
#   - smooth: spatial scale of the pattern; does not alter the value scale;
#   - contrast: clipping of the silhouette; does not alter the values;
#   - zoom: framing;
#   - setcolor and contourLevi: legibility of the figure.
#
# None of them changes the MEANING of the result. The ones that do are
# signal_mode, logfc_k and expressionLog, covered in script 04 -- if a result
# looks wrong, start there, not with the visual parameters.
#
# End of the series. For the graphical interface, see LEVIui(browser = TRUE).
