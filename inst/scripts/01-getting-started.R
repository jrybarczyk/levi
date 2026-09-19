# =============================================================================
# levi -- 01. Getting started
#
# Goal:     import a network and an expression file, generate the first
#           landscape and understand what the function returns.
#
# Data:     hub_network.dat / hub_expression.dat, shipped with the package.
#           Star-shaped network with 9 nodes: a central HUB, four inner
#           neighbours (N1-N4) and four outer ones (N5-N8).
#
# Runtime:  a few seconds.
# =============================================================================

library(levi)


# -----------------------------------------------------------------------------
# 1. The two input files
# -----------------------------------------------------------------------------
# levi always needs two files:
#
#   a) the NETWORK, with the coordinates of each node and the list of
#      interactions;
#   b) the EXPRESSION, with an identifier column and one or more value
#      columns.
#
# The identifier in the expression file must match the node name in the
# network.

network    <- system.file("extdata", "hub_network.dat",    package = "levi")
expression <- system.file("extdata", "hub_expression.dat", package = "levi")

# The network is in Medusa (DAT) format: the *edges section lists the
# interactions and the *nodes section lists "name  x  y".
cat("--- network ---\n")
writeLines(readLines(network))

cat("\n--- expression ---\n")
print(read.delim(expression), row.names = FALSE)

# In this dataset, HUB and N1-N4 are far higher in the test (200 against 10)
# and N5-N8 far lower (5 against 200). The expected result is therefore an
# over-expressed core surrounded by four repressed corners.


# -----------------------------------------------------------------------------
# 2. The minimal call
# -----------------------------------------------------------------------------
# readExpColumn() states which columns to compare, in "Test-Control" form.

result <- levi(
    expressionInput         = expression,
    networkCoordinatesInput = network,
    fileTypeInput           = "dat",          # dat, dyn, net or stg
    geneSymbolInput         = "ID",           # identifier column
    readExpColumn           = readExpColumn("Test-Control"),
    resolutionValueInput    = 40,
    smoothValueInput        = 50
)

# The plot is drawn automatically. The ggplot object lives in $plot, so it can
# be saved or modified afterwards:
#
#   ggplot2::ggsave("landscape.png", result$plot, width = 6, height = 5)


# -----------------------------------------------------------------------------
# 3. What levi() returns
# -----------------------------------------------------------------------------
cat("\n--- fields of the result ---\n")
str(result, max.level = 1, give.attr = FALSE)

# comparison -- the label of the comparison, useful in batch mode
cat("\ncomparison:", result$comparison, "\n")

# scores -- one value per network node, ordered from highest to lowest
cat("\n--- scores ---\n")
print(result$scores, row.names = FALSE)

# LandscapeScore lives in [0, 1] and has a well-defined neutral point:
#
#     1.0  maximum over-expression in the test
#     0.5  no change between test and control
#     0.0  maximum repression in the test
#
# Here HUB and N1-N4 reach 1.0 and N5-N8 fall to 0.0, exactly the asymmetry
# we put into the data.

# peaks -- automatically detected local maxima and minima
cat("\n--- peaks and valleys ---\n")
print(result$peaks, row.names = FALSE)

# HUB and N1-N4 share the same maximum value, so several peaks are reported;
# there is no single winner when there is a tie.

# landscape -- the plot matrix in long format (Var1, Var2, z)
cat("\n--- landscape ---\n")
cat("dimensions:", nrow(result$landscape), "rows\n")
print(head(result$landscape), row.names = FALSE)

# Cells outside the network are NA, not 0. That distinction matters: 0 means
# "there is network here, strongly repressed", while NA means "there is no
# network here".
cat("background cells (NA):",
    sum(is.na(result$landscape$z)), "of", nrow(result$landscape), "\n")

# pvalues -- NULL while n_perm = 0; see script 06.
cat("pvalues:", if (is.null(result$pvalues)) "NULL (n_perm = 0)" else "present", "\n")


# -----------------------------------------------------------------------------
# 4. How to read the result
# -----------------------------------------------------------------------------
# The landscape is a continuous surface built from the nodes and the midpoints
# of the edges. Each cell receives the average of the neighbouring signals,
# weighted by a Gaussian, which produces regions rather than isolated points:
# that regional reading is what levi adds to an ordinary heatmap.
#
# In practice, three questions guide the reading:
#
#   1. Are there red regions (above 0.5) and blue ones (below)? They point to
#      territories of the network coherently altered in the same direction.
#   2. Do the genes at the top and bottom of $scores make biological sense?
#   3. Are the regions broad or pointwise? A broad region suggests a
#      functional module; an isolated peak, a lone gene.
#
# Next: 02-reading-the-landscape.R goes deeper into the scale and into peak
# detection.
