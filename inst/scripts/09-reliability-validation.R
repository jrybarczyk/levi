# Reproducible numerical diagnostics after the 2026-09 technical review.
# Run with: Rscript inst/scripts/09-reliability-validation.R [output_directory]
# This small simulation is a regression diagnostic, not inferential certification.
library(levi)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else tempdir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260905)
xy <- as.matrix(expand.grid(x = c(.15, .5, .85), y = c(.15, .5, .85)))
edges <- cbind(1:8, 2:9)
coords <- rbind(xy, (xy[edges[, 1], ] + xy[edges[, 2], ]) / 2)
landscape <- function(values, mode, sigma = 1) {
    signals <- levi:::.networkSignals(values, edges, FALSE, mode, 1)
    matrix <- levi:::landscape_gauss(coords, signals$signal,
        signals$test, signals$control, 15, 0, 1/14, sigma, .05)$m1
    list(signals = signals, matrix = matrix)
}
rows <- list()
for (mode in c("ratio", "logfc", "zscore")) {
    for (missing in c(0L, 2L)) {
        rates <- numeric(40)
        discoveries <- logical(40)
        for (replicate in seq_len(40)) {
            values <- cbind(rnorm(9), rnorm(9))
            if (mode == "ratio") values <- exp(values)
            if (missing > 0L) values[seq_len(missing), ] <- NA_real_
            obs <- landscape(values, mode)
            s <- obs$signals
            p <- levi:::.permutationPvalues(coords, s$signal, s$test, s$control,
                obs$matrix, 15, 0, 1/14, 1, .05, 199,
                node_values = values, edge_index = edges, signal_mode = mode)
            q <- levi:::.adjustLandscapePvalues(p, "BY")
            rates[replicate] <- mean(unlist(p) <= .05, na.rm = TRUE)
            discoveries[replicate] <- any(unlist(q) <= .05, na.rm = TRUE)
        }
        rows[[length(rows) + 1L]] <- data.frame(mode, missing,
            replicates = 40, permutations = 199,
            raw_directional_cell_rejection_rate = mean(rates),
            raw_rate_mc_se = sd(rates) / sqrt(length(rates)),
            BY_any_discovery_fraction = mean(discoveries))
    }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "null_diagnostics.csv"), row.names = FALSE)

# Sensitivity of node scores to smoothing, layout, and missing coverage.
net <- levi:::.parseNetwork(system.file("extdata", "hub_network.dat", package = "levi"), NA, "dat")
net$nodes[, 2:3] <- lapply(net$nodes[, 2:3], as.numeric)
expr <- read.delim(system.file("extdata", "hub_expression.dat", package = "levi"))
run <- function(nodes, smooth = 5, expression = expr) {
    res <- suppressMessages(levi:::levi_function(expression, "dat", nodes, net$edges,
        "ID", readExpColumn("Test-Control"), 50, 50, 10, smooth, FALSE,
        FALSE, "default", .parsed_network = list(nodes = nodes, edges = net$edges),
        .draw = FALSE))
    setNames(res$scores$LandscapeScore, res$scores$Gene)
}
baseline <- run(net$nodes)
shifted <- net$nodes
shifted[, 2:3] <- shifted[, 2:3] + 100
changed <- net$nodes
changed[, 2:3] <- changed[rev(seq_len(nrow(changed))), 2:3]
sparse <- expr; sparse$Test[1:2] <- NA_real_
variants <- list(translation = run(shifted), smoothing_30 = run(net$nodes, 30),
                 reassigned_layout = run(changed), missing_two = run(net$nodes, expression = sparse))
sensitivity <- do.call(rbind, lapply(names(variants), function(name) {
    score <- variants[[name]][names(baseline)]
    data.frame(scenario = name, score_correlation = cor(baseline, score),
               max_absolute_score_change = max(abs(score - baseline)))
}))
write.csv(sensitivity, file.path(out_dir, "sensitivity_diagnostics.csv"), row.names = FALSE)
print(do.call(rbind, rows))
print(sensitivity)
cat("Small null sample and coarse Monte Carlo resolution limit interpretation.\n",
    "Absence of BY discoveries here does not establish calibration or power.\n")
