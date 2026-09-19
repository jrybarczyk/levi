# =============================================================================
# levi -- 06. Significance by permutation
#
# Goal:     inspect spatial association conditional on a fixed network
#           and layout, and read adjusted directional cell p-values.
#
# Note:     this script uses inference_unit = "cell", the legacy per-cell
#           test, because it is about reading p-value matrices. The default
#           since 2.0.0 is inference_unit = "region" (maximum regional
#           mass); see vignette("levi_inference") for the difference.
#
# Data:     bimodal (two opposing modules), medusa.dat with real expression
#           data, and flat (null control).
#
# Runtime:  about two minutes.
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
         smoothValueInput        = 50, inference_unit = "cell", ...)
}


# -----------------------------------------------------------------------------
# 1. The question the test answers
# -----------------------------------------------------------------------------
# Shuffle measured node expression pairs, preserving missing positions, and
# recalculate derived edge signals. This tests association between expression
# and network positions conditional on this layout. It does not remove layout
# sensitivity or test differential expression between biological replicates.
# Both tails over occupied cells form one family, adjusted with BY by default.

set.seed(42)   # the test is random; fixing the seed makes the result reproducible

bim <- run("bimodal", n_perm = 199, sig_level = 0.05, perm_side = "both")


# -----------------------------------------------------------------------------
# 2. The p-value matrices
# -----------------------------------------------------------------------------
cat("\n--- structure of $pvalues ---\n")
cat("components:", paste(names(bim$pvalues), collapse = ", "), "\n")
cat("dimensions:", paste(dim(bim$pvalues$over), collapse = " x "), "\n")

po <- bim$pvalues$over
pu <- bim$pvalues$under
cat(sprintf("p(over)  minimum %.4f | cells with p <= 0.05: %d\n",
            min(po, na.rm = TRUE), sum(po <= 0.05, na.rm = TRUE)))
cat(sprintf("p(under) minimum %.4f | cells with p <= 0.05: %d\n",
            min(pu, na.rm = TRUE), sum(pu <= 0.05, na.rm = TRUE)))

# $over  answers "is this region higher than expected by chance?"
# $under answers "is this region lower than expected by chance?"
#
# Both are needed because a repressed region is as informative as an
# over-expressed one, and a one-sided test would miss half the picture.
#
# Cells outside the network are NA in both matrices.


# -----------------------------------------------------------------------------
# 3. The p-value floor depends on n_perm
# -----------------------------------------------------------------------------
# The unadjusted p-value (stored in $raw_pvalues) is computed as (k + 1) / (n_perm + 1), where k counts how many
# permutations reached or exceeded the observed value. The "+1" prevents
# claiming p = 0, which no finite number of permutations can support.
#
# The raw floor is 1 / (n_perm + 1); adjustment can raise it substantially.

for (n in c(99, 199, 999)) {
    cat(sprintf("  n_perm = %3d  ->  smallest possible p-value = %.4f\n",
                n, 1 / (n + 1)))
}

# With n_perm = 19 the floor is 0.05, and no region can fall BELOW
# sig_level = 0.05: the test would have no resolution to detect anything. That
# is why n_perm must be comfortably larger than 1 / sig_level.
#
# More permutations improve numerical resolution; no fixed number establishes
# inferential validity. Inspect sensitivity and validate the chosen null model.


# -----------------------------------------------------------------------------
# 4. The contours on the plot
# -----------------------------------------------------------------------------
# With perm_side = "both", levi draws two white contour lines over the
# landscape, both at the value sig_level:
#
#   dashed  boundary of the significantly over-expressed region
#   dotted  boundary of the significantly repressed region
#
# perm_side selects which tails are displayed, without changing the adjustment
# family. No contour is drawn unless values occur on both sides of the adjusted
# threshold. With few permutations, BY may leave no significant cells.

set.seed(42)
med <- levi(
    expressionInput         = ex("expression.dat"),
    networkCoordinatesInput = ex("medusa.dat"),
    fileTypeInput           = "dat",
    geneSymbolInput         = "ID",
    readExpColumn           = readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
    resolutionValueInput    = 40,
    smoothValueInput        = 50,
    n_perm                  = 199)

in_mask <- sum(!is.na(med$pvalues$over))
cat("\n--- medusa.dat, real data ---\n")
cat(sprintf("cells in the mask: %d\n", in_mask))
cat(sprintf("  significant upwards:   %4d (%.0f%%)\n",
            sum(med$pvalues$over  <= 0.05, na.rm = TRUE),
            100 * sum(med$pvalues$over  <= 0.05, na.rm = TRUE) / in_mask))
cat(sprintf("  significant downwards: %4d (%.0f%%)\n",
            sum(med$pvalues$under <= 0.05, na.rm = TRUE),
            100 * sum(med$pvalues$under <= 0.05, na.rm = TRUE) / in_mask))

# Report the computed adjusted values; a visible colour pattern need not
# cross the significance threshold.


# -----------------------------------------------------------------------------
# 5. The null control
# -----------------------------------------------------------------------------
# In flat there is no variation at all, so shuffling the values changes
# nothing and no region can stand out. This checks a limiting case; it does
# not establish calibration under general random null signals.

set.seed(42)
flat <- run("flat", n_perm = 199)

cat("\n--- flat under permutation ---\n")
cat("cells with p(over) <= 0.05:",
    sum(flat$pvalues$over <= 0.05, na.rm = TRUE), "\n")
cat("cells with p(under) <= 0.05:",
    sum(flat$pvalues$under <= 0.05, na.rm = TRUE), "\n")

# When no cell crosses the threshold there is no contour to draw and levi
# warns instead of leaving the plot silent. This is a result, not an error.


# -----------------------------------------------------------------------------
# 6. Cost
# -----------------------------------------------------------------------------
# Each permutation recomputes the landscape, including its derived edge
# signals and Gaussian convolution, on the same coordinates.

t <- system.time(invisible(run("bimodal", n_perm = 199)))
cat(sprintf("\n199 permutations on the bimodal network: %.1f s\n", t[["elapsed"]]))

# The cost grows linearly with n_perm and with the square of the resolution,
# but hardly depends on the size of the network. It is worth calibrating the
# resolution before raising n_perm.


# -----------------------------------------------------------------------------
# 7. Summary
# -----------------------------------------------------------------------------
#   - permutation shuffles the expression while keeping the topology, which
#     conditions the test on the fixed layout;
#   - $pvalues$over and $under are jointly adjusted directional cell tests;
#   - the smallest possible raw p-value is 1 / (n_perm + 1): n_perm must be well
#     above 1 / sig_level;
#   - the contours mark where the signal holds up, not where there is colour.
#     If every cell falls on the same side of the threshold there is no curve
#     to draw, and levi warns;
#   - set a seed with set.seed() so the figure can be reproduced.
#
# Next: 07-bioconductor-data.R starts from real Bioconductor objects.
