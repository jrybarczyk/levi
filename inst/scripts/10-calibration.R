# =============================================================================
# levi -- 10. Calibration study: type I error, power and layout sensitivity
#
# Goal:     measure, by simulation, whether the permutation tests behave as
#           their documentation claims. Unit tests prove the code does what
#           it says; this script checks that the p-values are uniform under
#           the null, that the family-wise error rate stays at alpha, how
#           power grows with the effect size, and how much the landscape
#           test depends on the layout.
#
# This script is NOT run when the package is built. It takes roughly an hour
# on four cores. It writes small CSV summaries to inst/extdata/validation/,
# which the vignette "levi_validation" loads and plots in seconds.
#
# Run with: Rscript inst/scripts/10-calibration.R [workers] [replicates]
# =============================================================================

suppressPackageStartupMessages({ library(levi); library(BiocParallel) })
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 3L
n_rep   <- if (length(args) >= 2) as.integer(args[2]) else 300L
BPPARAM <- if (workers > 1) MulticoreParam(workers) else SerialParam()
out_dir <- file.path("inst", "extdata", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
alpha <- 0.05
t_start <- Sys.time()
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# -- Networks ----------------------------------------------------------------
# medusa: the 30-node real network shipped with levi.
# ba300: a 300-node scale-free graph (Barabasi-Albert, m = 2) laid out with
# Kamada-Kawai, closer in size and degree distribution to a STRING module.
medusa <- system.file("extdata", "medusa.dat", package = "levi")
set.seed(2026)
ba_graph <- igraph::sample_pa(300, m = 2, directed = FALSE)
igraph::V(ba_graph)$name <- paste0("g", seq_len(300))
ba <- leviFromEdges(ba_graph, layout = "kk")
networks <- list(
    medusa = list(coord = medusa, edges = NA, type = "dat",
                  genes = unique(levi:::.parseNetwork(medusa, NA, "dat")$nodes[, 1])),
    ba300  = list(coord = ba$nodes, edges = ba$edges, type = "stg",
                  genes = ba$nodes$name))

hub_module <- function(net, size = 8) {
    # The highest-degree node and its nearest neighbours: the "responding
    # module" in the power simulations.
    ni <- levi:::.networkIndex(net$coord, net$edges, net$type)
    g <- levi:::.graphFromEdges(length(ni$nodes), ni$edges)
    hub <- which.max(igraph::degree(g))
    order <- igraph::distances(g, v = hub)[1, ]
    ni$nodes[order(order)][seq_len(size)]
}

# -- Simulation helpers ------------------------------------------------------
sim_expression <- function(genes, n_per_group = 4, effect = 0, module = NULL) {
    x <- matrix(rnorm(length(genes) * 2 * n_per_group, sd = 1), length(genes),
                dimnames = list(genes, NULL))
    groups <- rep(c("C", "T"), each = n_per_group)
    if (effect != 0 && length(module))
        x[module, groups == "T"] <- x[module, groups == "T"] + effect
    list(x = x, groups = groups)
}

regional_min_p <- function(net, lfc, n_perm = 199, smooth = 50) {
    r <- levi(expressionInput = data.frame(ID = names(lfc), logFC = lfc),
              networkCoordinatesInput = net$coord,
              networkInteractionsInput = net$edges, fileTypeInput = net$type,
              geneSymbolInput = "ID", readExpColumn = readExpColumn("logFC-logFC"),
              signal_mode = "logfc", resolutionValueInput = 20,
              smoothValueInput = smooth, n_perm = n_perm, .draw = FALSE)
    p <- r$regions$summary$PSpatial
    significant <- r$regions$summary$Region[r$regions$summary$Significant]
    genes <- if (length(significant)) {
        attribution <- leviRegionGenes(r, top_n = 5)
        unique(attribution$Gene[attribution$Region %in% significant])
    } else character()
    list(min_p = if (length(p)) min(p) else NA_real_, n_regions = length(p),
         sig_genes = genes)
}

replicate_min_p <- function(net, sim) {
    r <- suppressWarnings(leviReplicateInference(sim$x, sim$groups, test = "T",
        control = "C", networkCoordinatesInput = net$coord,
        networkInteractionsInput = net$edges, fileTypeInput = net$type,
        resolutionValueInput = 20, smoothValueInput = 50))
    p <- r$regions$summary$PSpatial
    if (length(p)) min(p) else NA_real_
}

tfce_min_p <- function(net, sim) {
    r <- suppressWarnings(leviGraphTFCEInference(sim$x, sim$groups, net$coord,
        net$edges, fileTypeInput = net$type, test = "T", control = "C",
        n_steps = 50))
    min(r$statistic$PGlobal)
}

rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(...)

# -- 1. Type I error under the global null ------------------------------------
# One replicate = one dataset with no effect anywhere. The family-wise error
# rate is the fraction of replicates whose smallest p-value is <= alpha; with
# regions redetected in every permutation and the maximum taken over both
# directions, it should be at most alpha.
for (nm in names(networks)) {
    net <- networks[[nm]]
    say("Type I error, node-label regional test,", nm)
    node <- bplapply(seq_len(n_rep), function(i) {
        set.seed(1000 + i)
        regional_min_p(net, setNames(rnorm(length(net$genes)), net$genes))$min_p
    }, BPPARAM = BPPARAM)
    node <- unlist(node)
    add(test = "node-label regional", network = nm, design = "iid logFC",
        replicates = n_rep, alpha = alpha,
        fwer = mean(node <= alpha, na.rm = TRUE),
        fwer_se = sqrt(mean(node <= alpha, na.rm = TRUE) *
                       (1 - mean(node <= alpha, na.rm = TRUE)) / sum(!is.na(node))),
        no_region = mean(is.na(node)))
    write.csv(data.frame(test = "node-label regional", network = nm, min_p = node),
              file.path(out_dir, sprintf("null_pvalues_regional_%s.csv", nm)),
              row.names = FALSE)

    say("Type I error, sample-label replicate test,", nm)
    rep_p <- unlist(bplapply(seq_len(n_rep), function(i) {
        set.seed(2000 + i)
        replicate_min_p(net, sim_expression(net$genes))
    }, BPPARAM = BPPARAM))
    add(test = "sample-label replicate", network = nm, design = "4 vs 4, exact",
        replicates = n_rep, alpha = alpha,
        fwer = mean(rep_p <= alpha, na.rm = TRUE),
        fwer_se = sqrt(mean(rep_p <= alpha, na.rm = TRUE) *
                       (1 - mean(rep_p <= alpha, na.rm = TRUE)) / sum(!is.na(rep_p))),
        no_region = mean(is.na(rep_p)))
    write.csv(data.frame(test = "sample-label replicate", network = nm, min_p = rep_p),
              file.path(out_dir, sprintf("null_pvalues_replicate_%s.csv", nm)),
              row.names = FALSE)

    say("Type I error, graph TFCE,", nm)
    n_tfce <- min(n_rep, 200L)
    tf_p <- unlist(bplapply(seq_len(n_tfce), function(i) {
        set.seed(3000 + i)
        tfce_min_p(net, sim_expression(net$genes))
    }, BPPARAM = BPPARAM))
    add(test = "graph TFCE", network = nm, design = "4 vs 4, exact",
        replicates = n_tfce, alpha = alpha,
        fwer = mean(tf_p <= alpha),
        fwer_se = sqrt(mean(tf_p <= alpha) * (1 - mean(tf_p <= alpha)) / n_tfce),
        no_region = 0)
    write.csv(data.frame(test = "graph TFCE", network = nm, min_p = tf_p),
              file.path(out_dir, sprintf("null_pvalues_tfce_%s.csv", nm)),
              row.names = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "type1_error.csv"),
          row.names = FALSE)

# -- 2. Power against effect size on the hub module ---------------------------
power_rows <- list()
n_power <- min(n_rep, 100L)
for (nm in names(networks)) {
    net <- networks[[nm]]
    module <- hub_module(net)
    for (effect in c(0.5, 1, 1.5, 2)) {
        say("Power,", nm, "effect", effect)
        res <- bplapply(seq_len(n_power), function(i) {
            set.seed(4000 + 100 * effect + i)
            sim <- sim_expression(net$genes, effect = effect, module = module)
            c(replicate = replicate_min_p(net, sim), tfce = tfce_min_p(net, sim))
        }, BPPARAM = BPPARAM)
        res <- do.call(rbind, res)
        for (test in colnames(res))
            power_rows[[length(power_rows) + 1L]] <- data.frame(
                test = c(replicate = "sample-label replicate",
                         tfce = "graph TFCE")[[test]],
                network = nm, module_size = length(module), effect = effect,
                replicates = n_power,
                power = mean(res[, test] <= alpha, na.rm = TRUE))
    }
}
write.csv(do.call(rbind, power_rows), file.path(out_dir, "power.csv"),
          row.names = FALSE)

# -- 3. Layout sensitivity ----------------------------------------------------
# Same graph, same expression, twenty force-directed layouts. How stable is
# the set of genes attributed to significant regions? Compared with TFCE,
# which never looks at coordinates and therefore gives one answer. The effect
# (a 20-gene module shifted by 3 SD, smoothing 20) is deliberately large:
# sensitivity to the layout can only be measured where the test rejects.
say("Layout sensitivity")
set.seed(5000)
module <- hub_module(networks$ba300, size = 20)
sim <- sim_expression(networks$ba300$genes, effect = 3, module = module)
lfc <- rowMeans(sim$x[, sim$groups == "T"]) - rowMeans(sim$x[, sim$groups == "C"])
layouts <- bplapply(seq_len(20), function(i) {
    set.seed(6000 + i)
    lay <- leviFromEdges(ba_graph, layout = "fr")
    net <- list(coord = lay$nodes, edges = lay$edges, type = "stg")
    out <- regional_min_p(net, lfc, smooth = 20)
    list(min_p = out$min_p, n_regions = out$n_regions, genes = out$sig_genes)
}, BPPARAM = BPPARAM)
sets <- lapply(layouts, `[[`, "genes")
jaccard <- function(a, b) if (!length(a) && !length(b)) 1 else
    length(intersect(a, b)) / length(union(a, b))
pairs <- combn(seq_along(sets), 2)
layout_rows <- data.frame(
    layout = seq_along(layouts),
    min_p = vapply(layouts, `[[`, numeric(1), "min_p"),
    n_regions = vapply(layouts, `[[`, integer(1), "n_regions"),
    n_sig_genes = lengths(sets),
    module_recovered = vapply(sets, function(s) mean(module %in% s), numeric(1)))
write.csv(layout_rows, file.path(out_dir, "layout_sensitivity.csv"),
          row.names = FALSE)
write.csv(data.frame(
    mean_pairwise_jaccard = mean(apply(pairs, 2, function(p)
        jaccard(sets[[p[1]]], sets[[p[2]]]))),
    layouts = length(sets),
    tfce_min_p = tfce_min_p(networks$ba300, sim)),
    file.path(out_dir, "layout_sensitivity_summary.csv"), row.names = FALSE)

writeLines(c(
    "# Calibration results",
    "",
    sprintf("Generated on %s by inst/scripts/10-calibration.R with %d workers,",
            format(Sys.Date()), workers),
    sprintf("%d replicates per null setting, in %.0f minutes.", n_rep,
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
    "",
    "* type1_error.csv: family-wise error rate of each test under the global",
    "  null (no effect anywhere), by network.",
    "* null_pvalues_*.csv: smallest p-value per null replicate, for QQ plots.",
    "* power.csv: probability of at least one significant result when the",
    "  hub module of the network is shifted by `effect` in 4 vs 4 samples.",
    "* layout_sensitivity*.csv: the same data on twenty Fruchterman-Reingold",
    "  layouts of the 300-node graph; pairwise Jaccard of the genes in",
    "  significant regions, versus the layout-free TFCE result.",
    "",
    sprintf("levi %s, R %s.", packageVersion("levi"), getRversion())),
    file.path(out_dir, "README.md"))
say("Done in", round(as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
    "minutes")
