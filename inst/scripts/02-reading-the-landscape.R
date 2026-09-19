# =============================================================================
# levi -- 02. Reading the landscape
#
# Goal:     understand the scale, the neutral point, the silhouette and peak
#           detection, using datasets whose correct result is known.
#
# Data:     flat (null control), hub (asymmetric) and gradient (monotonic).
#
# Runtime:  less than a minute.
# =============================================================================

library(levi)

ex <- function(f) system.file("extdata", f, package = "levi")

run <- function(name, ...) {
    levi(expressionInput         = ex(paste0(name, "_expression.dat")),
         networkCoordinatesInput = ex(paste0(name, "_network.dat")),
         fileTypeInput           = "dat",
         geneSymbolInput         = "ID",
         readExpColumn           = readExpColumn("Test-Control"),
         resolutionValueInput    = 40,
         smoothValueInput        = 50, ...)
}


# -----------------------------------------------------------------------------
# 1. The neutral point: 0.5 means "no change"
# -----------------------------------------------------------------------------
# The flat dataset is the null control: nine genes with Test = Control = 100.
# Nothing changed between the conditions, so the correct answer is 0.5
# everywhere.

flat <- run("flat")

cat("--- flat: no variation between the conditions ---\n")
print(flat$scores[, c("Gene", "LandscapeScore")], row.names = FALSE)
cat("standard deviation of the scores:", sd(flat$scores$LandscapeScore), "\n")
cat("landscape range:",
    sprintf("%.3f to %.3f", min(flat$landscape$z, na.rm = TRUE),
                            max(flat$landscape$z, na.rm = TRUE)), "\n")

# Every value comes out at exactly 0.5. This is not an accident of
# parameterisation: the landscape is a weighted average of the signals, and
# the average of values all equal to 0.5 is 0.5, whatever the smoothing or the
# network density.
#
# It is worth using flat as a sanity check whenever a result looks odd: if a
# dataset with no variation does not give 0.5, something is wrong in the input
# data.


# -----------------------------------------------------------------------------
# 2. The scale is absolute, not relative
# -----------------------------------------------------------------------------
# A value of 0.8 means "over-expressed", not merely "higher than the rest of
# this figure". That is what makes two different landscapes comparable.

hub <- run("hub")

cat("\n--- hub: over-expressed core, repressed periphery ---\n")
print(hub$scores[, c("Gene", "LandscapeScore", "Rank")], row.names = FALSE)

above <- hub$scores$Gene[hub$scores$LandscapeScore > 0.5]
below <- hub$scores$Gene[hub$scores$LandscapeScore < 0.5]
cat("above 0.5:", paste(above, collapse = ", "), "\n")
cat("below 0.5:", paste(below, collapse = ", "), "\n")


# -----------------------------------------------------------------------------
# 3. The silhouette: NA is not zero
# -----------------------------------------------------------------------------
# Cells far from any node or edge receive no value: they are NA and are drawn
# in the background colour. The distinction matters, because 0 is a legitimate
# value meaning "strongly repressed gene".

cat("\n--- silhouette ---\n")
n_na <- sum(is.na(hub$landscape$z))
cat(sprintf("background: %d of %d cells (%.0f%%)\n",
            n_na, nrow(hub$landscape), 100 * n_na / nrow(hub$landscape)))

# The contrastValueInput parameter controls how tightly the silhouette hugs
# the network.
for (ct in c(10, 50, 90)) {
    r <- run("hub", contrastValueInput = ct)
    cat(sprintf("  contrastValueInput = %2d  ->  background %2.0f%%\n",
                ct, 100 * mean(is.na(r$landscape$z))))
}

# High contrast tightens the silhouette around the network; low contrast lets
# the patch spread out. It does not change the values, only the clipping.


# -----------------------------------------------------------------------------
# 4. Peaks and valleys
# -----------------------------------------------------------------------------
# $peaks lists local maxima and minima: cells that dominate their neighbourhood
# and cross a threshold (0.65 for a peak, 0.35 for a valley). Each one is
# labelled with the nearest node.

cat("\n--- peaks and valleys of hub ---\n")
print(hub$peaks, row.names = FALSE)

cat("\ncount by type:\n")
print(table(hub$peaks$Type))

# A landscape with no peaks or valleys is not an error: it means no region
# reached sufficient intensity. flat is the extreme example.
cat("\npeaks/valleys in flat:", nrow(flat$peaks), "\n")


# -----------------------------------------------------------------------------
# 5. Ordering and monotonicity
# -----------------------------------------------------------------------------
# The gradient dataset has six genes in a chain, with expression rising
# monotonically from GA to GF. The landscape should preserve that order.

grad <- run("gradient")
s <- setNames(grad$scores$LandscapeScore, grad$scores$Gene)
ordering <- s[c("GA", "GB", "GC", "GD", "GE", "GF")]

cat("\n--- gradient: is the order preserved? ---\n")
print(round(ordering, 3))
cat("strictly increasing:", !is.unsorted(ordering, strictly = TRUE), "\n")

# Smoothing pulls the extremes towards the centre, but does not invert the
# order. The stronger the smoothing, the more the values converge to the
# middle -- see script 08.


# -----------------------------------------------------------------------------
# 6. Practical summary
# -----------------------------------------------------------------------------
# When interpreting a landscape:
#
#   - 0.5 is the neutral point; above it is over-expression, below it is
#     repression;
#   - the scale is absolute, so two landscapes with the same parameters are
#     comparable side by side (see script 05);
#   - the grey background is not a low value, it is absence of network;
#   - $peaks points to where to look, but ties produce several peaks;
#   - broad regions weigh more than isolated cells, because they reflect
#     several neighbouring nodes agreeing.
#
# Next: 03-topology-and-landscape.R shows how the shape of the network moulds
# the result.
