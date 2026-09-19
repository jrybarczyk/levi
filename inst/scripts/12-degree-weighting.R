# =============================================================================
# levi -- 12. Degree weighting study: does edge_weighting change calibration,
#           power or layout sensitivity of the landscape tests?
#
# Compares edge_weighting = "midpoint" (historical), "degree" and "none" on
# the 300-node Barabasi-Albert network used by 10-calibration.R:
#   1. family-wise error under the global null (node-label regional and
#      sample-label replicate tests);
#   2. power of the replicate test when the hub module responds;
#   3. layout sensitivity: pairwise Jaccard of significant genes over twenty
#      Fruchterman-Reingold layouts.
# Not run at build time. Run with:
#   Rscript inst/scripts/12-degree-weighting.R [workers] [replicates] [out_dir]
# =============================================================================
suppressPackageStartupMessages({ library(levi); library(BiocParallel) })
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 3L
n_rep   <- if (length(args) >= 2) as.integer(args[2]) else 100L
out_dir <- if (length(args) >= 3) args[3] else file.path("inst", "extdata", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
BPPARAM <- if (workers > 1) MulticoreParam(workers) else SerialParam()
alpha <- 0.05; modes <- c("midpoint", "degree", "none")
t_start <- Sys.time(); say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

set.seed(2026)
ba_graph <- igraph::sample_pa(300, m = 2, directed = FALSE)
igraph::V(ba_graph)$name <- paste0("g", seq_len(300))
ba <- leviFromEdges(ba_graph, layout = "kk")
net <- list(coord = ba$nodes, edges = ba$edges, type = "stg", genes = ba$nodes$name)
deg <- igraph::degree(ba_graph); cat("graus: max", max(deg), " mediana", median(deg), "\n")

hub_module <- function(size) {
    hub <- which.max(deg); d <- igraph::distances(ba_graph, v = hub)[1, ]
    net$genes[order(d)][seq_len(size)]
}
sim_expression <- function(n_per_group = 4, effect = 0, module = NULL) {
    x <- matrix(rnorm(300 * 2 * n_per_group), 300, dimnames = list(net$genes, NULL))
    groups <- rep(c("C", "T"), each = n_per_group)
    if (effect != 0) x[module, groups == "T"] <- x[module, groups == "T"] + effect
    list(x = x, groups = groups)
}
regional <- function(lfc, mode, n_perm = 199, smooth = 50, coord = net$coord, edges = net$edges) {
    r <- levi(expressionInput = data.frame(ID = names(lfc), logFC = lfc),
              networkCoordinatesInput = coord, networkInteractionsInput = edges,
              fileTypeInput = "stg", geneSymbolInput = "ID",
              readExpColumn = readExpColumn("logFC-logFC"), signal_mode = "logfc",
              resolutionValueInput = 20, smoothValueInput = smooth, n_perm = n_perm,
              edge_weighting = mode, .draw = FALSE)
    p <- r$regions$summary$PSpatial
    sig <- r$regions$summary$Region[r$regions$summary$Significant]
    genes <- if (length(sig)) { a <- leviRegionGenes(r, top_n = 5); unique(a$Gene[a$Region %in% sig]) } else character()
    list(min_p = if (length(p)) min(p) else NA_real_, genes = genes)
}
replicate_p <- function(sim, mode) {
    r <- suppressWarnings(leviReplicateInference(sim$x, sim$groups, test = "T", control = "C",
        networkCoordinatesInput = net$coord, networkInteractionsInput = net$edges,
        fileTypeInput = "stg", resolutionValueInput = 20, smoothValueInput = 50,
        edge_weighting = mode))
    p <- r$regions$summary$PSpatial; if (length(p)) min(p) else NA_real_
}

rows <- list()
for (mode in modes) {
    say("Type I, node-label regional,", mode)
    p_node <- unlist(bplapply(seq_len(n_rep), function(i) { set.seed(1000 + i)
        regional(setNames(rnorm(300), net$genes), mode)$min_p }, BPPARAM = BPPARAM))
    say("Type I, sample-label replicate,", mode)
    p_rep <- unlist(bplapply(seq_len(n_rep), function(i) { set.seed(2000 + i)
        replicate_p(sim_expression(), mode) }, BPPARAM = BPPARAM))
    for (nm in c("node", "rep")) { p <- get(paste0("p_", nm)); f <- mean(p <= alpha, na.rm = TRUE)
        rows[[length(rows) + 1]] <- data.frame(edge_weighting = mode,
            test = c(node = "node-label regional", rep = "sample-label replicate")[[nm]],
            replicates = n_rep, fwer = f, fwer_se = sqrt(f * (1 - f) / sum(!is.na(p))), no_region = mean(is.na(p))) }
}
type1 <- do.call(rbind, rows); print(type1)
write.csv(type1, file.path(out_dir, "degree_weighting_type1.csv"), row.names = FALSE)

module <- hub_module(8); power_rows <- list(); n_power <- min(n_rep, 100L)
for (effect in c(1, 1.5)) for (mode in modes) {
    say("Power, effect", effect, mode)
    p <- unlist(bplapply(seq_len(n_power), function(i) { set.seed(4000 + 100 * effect + i)
        replicate_p(sim_expression(effect = effect, module = module), mode) }, BPPARAM = BPPARAM))
    power_rows[[length(power_rows) + 1]] <- data.frame(edge_weighting = mode, test = "sample-label replicate",
        module_size = 8, effect = effect, replicates = n_power, power = mean(p <= alpha, na.rm = TRUE))
}
power <- do.call(rbind, power_rows); print(power)
write.csv(power, file.path(out_dir, "degree_weighting_power.csv"), row.names = FALSE)

say("Layout sensitivity")
set.seed(5000); module20 <- hub_module(20)
sim <- sim_expression(effect = 3, module = module20)
lfc <- rowMeans(sim$x[, sim$groups == "T"]) - rowMeans(sim$x[, sim$groups == "C"])
layouts <- lapply(seq_len(20), function(i) { set.seed(6000 + i); leviFromEdges(ba_graph, layout = "fr") })
jaccard <- function(a, b) if (!length(a) && !length(b)) 1 else length(intersect(a, b)) / length(union(a, b))
lay_rows <- list()
for (mode in modes) {
    out <- bplapply(layouts, function(lay) regional(lfc, mode, smooth = 20, coord = lay$nodes, edges = lay$edges), BPPARAM = BPPARAM)
    sets <- lapply(out, `[[`, "genes"); pairs <- combn(20, 2)
    lay_rows[[length(lay_rows) + 1]] <- data.frame(edge_weighting = mode,
        mean_pairwise_jaccard = mean(apply(pairs, 2, function(p) jaccard(sets[[p[1]]], sets[[p[2]]]))),
        layouts_significant = sum(vapply(out, function(o) isTRUE(o$min_p <= alpha), logical(1))),
        mean_sig_genes = mean(lengths(sets)),
        module_recovered = mean(vapply(sets, function(s) mean(module20 %in% s), numeric(1))))
}
lay <- do.call(rbind, lay_rows); print(lay)
write.csv(lay, file.path(out_dir, "degree_weighting_layout.csv"), row.names = FALSE)
say("DONE em", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min")
