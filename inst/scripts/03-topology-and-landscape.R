# =============================================================================
# levi -- 03. How network topology moulds the landscape
#
# Goal:     compare five topologies under the same parameters, to separate
#           what comes from the EXPRESSION from what comes from the SHAPE of
#           the network.
#
# Data:     the five toy datasets of the package, each with an expected
#           result documented in ?hub_dataset and friends.
#
# Runtime:  about a minute.
# =============================================================================

library(levi)

ex <- function(f) system.file("extdata", f, package = "levi")

run <- function(name) {
    levi(expressionInput         = ex(paste0(name, "_expression.dat")),
         networkCoordinatesInput = ex(paste0(name, "_network.dat")),
         fileTypeInput           = "dat",
         geneSymbolInput         = "ID",
         readExpColumn           = readExpColumn("Test-Control"),
         resolutionValueInput    = 40,
         smoothValueInput        = 50)
}

summarise <- function(r, name, comment) {
    s <- r$scores$LandscapeScore
    cat(sprintf("\n== %-9s %s\n", name, comment))
    cat(sprintf("   nodes: %2d | scores from %.3f to %.3f | median %.3f\n",
                nrow(r$scores), min(s), max(s), median(s)))
    cat(sprintf("   peaks: %d | valleys: %d | background: %.0f%%\n",
                sum(r$peaks$Type == "peak"), sum(r$peaks$Type == "valley"),
                100 * mean(is.na(r$landscape$z))))
    invisible(r)
}


# -----------------------------------------------------------------------------
# 1. Star: a centre and a periphery
# -----------------------------------------------------------------------------
# Nine nodes: HUB at the centre, N1-N4 nearby, N5-N8 at the corners. The
# centre and the inner neighbours are over-expressed; the corners, repressed.

hub <- summarise(run("hub"), "hub",
    "star, high core against low periphery")

# The landscape shows a compact central patch and four separate islands.
# Since N5-N8 have no neighbours among themselves, each corner becomes an
# isolated valley.


# -----------------------------------------------------------------------------
# 2. Chain: a spatial gradient
# -----------------------------------------------------------------------------
# Six nodes in a line, with expression rising from one end to the other.

grad <- summarise(run("gradient"), "gradient",
    "linear chain, expression rising from GA to GF")

s <- setNames(round(grad$scores$LandscapeScore, 3), grad$scores$Gene)
cat("   order: ", paste(names(s[order(s)]), collapse = " < "), "\n")

# Here each node's value depends on its neighbours in the chain, so the
# landscape becomes a continuous ramp. No isolated node dominates, which is
# why peak detection finds little: there is no local maximum in the middle of
# a ramp.


# -----------------------------------------------------------------------------
# 3. Bimodal: two modules in opposite directions
# -----------------------------------------------------------------------------
# Two clusters joined by a bridge. Cluster A is entirely over-expressed and
# cluster B entirely repressed.

bim <- summarise(run("bimodal"), "bimodal",
    "two opposing clusters, joined by a bridge")

sb <- setNames(bim$scores$LandscapeScore, bim$scores$Gene)
cat(sprintf("   cluster A (mean %.3f) against cluster B (mean %.3f)\n",
            mean(sb[grep("^A", names(sb))]), mean(sb[grep("^B", names(sb))])))

# This is the case levi handles better than a heatmap: the two territories
# appear as contiguous patches, and the boundary between them is visible. A
# ranked gene list would show the same values, but would not show that they
# cluster in space.


# -----------------------------------------------------------------------------
# 4. Flat: the null control
# -----------------------------------------------------------------------------
flat <- summarise(run("flat"), "flat",
    "no variation between the conditions")

cat("   all scores equal to", unique(round(flat$scores$LandscapeScore, 3)), "\n")

# No variation, no structure: a uniformly neutral surface and no peaks. This
# is the result that validates the reading of all the others.


# -----------------------------------------------------------------------------
# 5. Sparse: network larger than the data
# -----------------------------------------------------------------------------
# Fifteen nodes in the network, but the expression file carries only five. A
# common situation in practice: the network comes from an interaction database
# and covers more genes than the experiment measured.

sparse <- summarise(run("sparse"), "sparse",
    "15 nodes in the network, only 5 with a measured value")

cat("\n   scores:\n")
print(sparse$scores[, c("Gene", "LandscapeScore", "Rank")], row.names = FALSE)

# The ten unmeasured genes receive the neutral value and therefore pull the
# surrounding region towards 0.5 instead of distorting it to an extreme. levi
# also writes a log with the names of those nodes; the path appears in the
# message emitted during execution.


# -----------------------------------------------------------------------------
# 6. What the comparison teaches
# -----------------------------------------------------------------------------
# Same parameters, five very different results. Worth keeping in mind:
#
#   - EXPRESSION sets the height, TOPOLOGY sets the drawing. The same values
#     in a different layout produce a different figure, so the network layout
#     is an analytical choice, not merely an aesthetic one;
#   - neighbouring nodes reinforce each other: an altered gene surrounded by
#     genes altered in the same direction produces a broad region, while an
#     isolated gene barely shows. That is the information levi adds;
#   - the absence of peaks is informative (gradient and flat), not a failure;
#   - unmeasured genes do not invent signal: they stay at the neutral point.
#
# Next: 04-signal-modes.R covers how raw values become a score.
