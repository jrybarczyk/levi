# =============================================================================
# levi -- 05. Batch mode and comparison between landscapes
#
# Goal:     generate several landscapes in a single call, lay them out side by
#           side and subtract one from another to see what changed between
#           conditions.
#
# Data:     hub_network.dat with hub_multicomp_expression.dat, which carries
#           three conditions (Cond_A, Cond_B, Cond_C) for the same nine nodes.
#
# Runtime:  about a minute.
# =============================================================================

library(levi)

ex <- function(f) system.file("extdata", f, package = "levi")

cat("--- expression with three conditions ---\n")
print(read.delim(ex("hub_multicomp_expression.dat")), row.names = FALSE)

# Cond_A: high core, low periphery.
# Cond_B: the reverse.
# Cond_C: intermediate values, serving as a common reference.


# -----------------------------------------------------------------------------
# 1. Several comparisons in one call
# -----------------------------------------------------------------------------
# readExpColumn() accepts as many comparisons as needed. The network is read
# and processed only once, not once per comparison.

results <- levi(
    expressionInput         = ex("hub_multicomp_expression.dat"),
    networkCoordinatesInput = ex("hub_network.dat"),
    fileTypeInput           = "dat",
    geneSymbolInput         = "ID",
    readExpColumn           = readExpColumn("Cond_A-Cond_C",
                                            "Cond_B-Cond_C"),
    resolutionValueInput    = 40,
    smoothValueInput        = 50
)

# With more than one comparison, levi() returns a LIST of results. With a
# single one, it returns the result directly -- worth checking before
# indexing.
cat("\nnumber of results:", length(results), "\n")
for (r in results)
    cat(sprintf("  %-14s scores from %.3f to %.3f\n", r$comparison,
                min(r$scores$LandscapeScore), max(r$scores$LandscapeScore)))


# -----------------------------------------------------------------------------
# 2. Side by side with leviGrid()
# -----------------------------------------------------------------------------
# Since the scale is absolute and the parameters are the same, the panels can
# be read together: the same colour means the same value in both.

leviGrid(results, ncol = 2)

# leviGrid uses patchwork, cowplot or gridExtra, whichever is installed. With
# none of them, it prints the plots in sequence.


# -----------------------------------------------------------------------------
# 3. The differential landscape
# -----------------------------------------------------------------------------
# leviDiff() subtracts two surfaces cell by cell. The result answers a
# different question: not "where is expression altered", but "where do the two
# conditions DIVERGE".

dif <- leviDiff(results[[1]], results[[2]])

cat("\n--- differential landscape ---\n")
cat("comparison:", dif$comparison, "\n")
cat(sprintf("range of Diff: %.3f to %.3f\n",
            min(dif$diff$Diff, na.rm = TRUE), max(dif$diff$Diff, na.rm = TRUE)))
print(head(dif$diff), row.names = FALSE)

# Diff = B - A, therefore:
#
#   Diff > 0  (red)    the region is more expressed in B than in A
#   Diff ~ 0  (white)  the two conditions agree there
#   Diff < 0  (blue)   the region is more expressed in A
#
# Here Cond_A and Cond_B were built as opposites, so the differential map
# saturates at both extremes: the core comes out strongly negative and the
# periphery strongly positive.

cat("\ncells with |Diff| > 0.5:",
    sum(abs(dif$diff$Diff) > 0.5, na.rm = TRUE), "of", nrow(dif$diff), "\n")


# -----------------------------------------------------------------------------
# 4. Requirements for comparing
# -----------------------------------------------------------------------------
# For two landscapes to be comparable, whether in the grid or in the
# differential:
#
#   - SAME network. Cells only correspond to one another if the coordinates
#     are the same; leviDiff() stops with an error if the sizes differ;
#   - SAME resolutionValueInput. Different resolutions produce grids of
#     different sizes;
#   - SAME signal_mode. Mixing ratio with logfc compares distinct scales;
#   - preferably the same smoothValueInput, so that the level of detail is the
#     same in both.
#
# The safest route is to generate the comparisons in a single call, as in
# step 1: that way all parameters are necessarily identical.


# -----------------------------------------------------------------------------
# 5. Accumulating results from separate calls
# -----------------------------------------------------------------------------
# When the comparisons come from different files, simply gather the outputs
# into a list before passing them to leviGrid.

one <- levi(expressionInput         = ex("hub_expression.dat"),
            networkCoordinatesInput = ex("hub_network.dat"),
            fileTypeInput           = "dat",
            geneSymbolInput         = "ID",
            readExpColumn           = readExpColumn("Test-Control"),
            resolutionValueInput    = 40,
            smoothValueInput        = 50)

leviGrid(list(one, results[[1]]), ncol = 2,
         titles = c("hub_expression", "multicomp Cond_A-Cond_C"))

# titles replaces the automatic labels, useful when the name of the comparison
# does not say much on its own.


# -----------------------------------------------------------------------------
# 6. Summary
# -----------------------------------------------------------------------------
#   - several comparisons in one call guarantee identical parameters and avoid
#     reprocessing the network;
#   - with one comparison the return is a result; with several, a list;
#   - leviGrid compares absolute values; leviDiff highlights divergence;
#   - Diff is always B minus A, in the order the arguments were passed.
#
# Next: 06-significance.R separates signal from chance.
