test_that("pathway adjustment adds a family-wise correction", {
    make_result <- function(p) list(regions = list(summary = data.frame(PSpatial = p)))
    out <- leviAdjustPathways(list(a = make_result(c(.01, .2)),
                                   b = make_result(.03)), method = "holm")
    expect_equal(out$a$regions$summary$PFamily, c(.03, .2))
    expect_equal(out$b$regions$summary$PFamily, .06)
    expect_false(out$b$regions$summary$FamilySignificant)
})

test_that("stratified node permutations do not cross strata", {
    values <- cbind(c(-2, -1, 1, 2), c(0, 0, 0, 0))
    edge <- matrix(c(1, 2), 1)
    observed <- matrix(.5, 2, 2)
    regions <- levi:::.extractLandscapeRegions(observed, 0, 1, .1, 1)
    testthat::local_mocked_bindings(landscape_gauss = function(..., SignalOut) {
        expect_true(all(SignalOut[1:2] < .5))
        expect_true(all(SignalOut[3:4] > .5))
        list(m1 = observed)
    }, .package = "levi")
    sig <- levi:::.networkSignals(values, edge, FALSE, "logfc", 1)
    levi:::.permutationPvalues(matrix(0, 5, 2), sig$signal, sig$test,
        sig$control, observed, 2, 0, 1, 1, .05, 4,
        node_values = values, edge_index = edge, signal_mode = "logfc",
        regions = regions, perm_strata = c("low", "low", "high", "high"))
})

test_that("region gene attribution ranks local non-neutral nodes", {
    x <- list(
        regions = list(cells = data.frame(Region = "over_01", X = 0, Y = 0)),
        metadata = list(node_coordinates = rbind(c(0, 0), c(1, 1)),
            node_signal = c(.8, .9), nodes = data.frame(V1 = c("A", "B")),
            grid = list(sigma = 1, increase = .1)))
    out <- leviRegionGenes(x)
    expect_equal(out$Gene[1], "A")
    expect_true(out$Contribution[1] > out$Contribution[2])
})

test_that("replicate inference permutes sample labels and returns regional p-values", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    genes <- c("HUB", paste0("N", 1:8))
    mat <- matrix(0, length(genes), 4, dimnames = list(genes, NULL))
    mat[, 3:4] <- 1
    set.seed(4)
    expect_warning(out <- leviReplicateInference(mat, c("C", "C", "T", "T"),
        test = "T", control = "C", n_perm = 3, seed = 2,
        networkCoordinatesInput = system.file("extdata", "hub_network.dat", package = "levi"),
        fileTypeInput = "dat", resolutionValueInput = 1,
        region_threshold = .02), "smallest attainable p-value")
    expect_equal(out$regions$inference$method, "sample_label")
    expect_true(all(out$regions$summary$PSpatial >= .25))
    expect_equal(out$metadata$n_perm, 5L)
    expect_true(out$metadata$permutation_exact)
})

test_that("graph cluster inference uses network-connected limma clusters", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8))
    x <- matrix(rnorm(length(genes) * 6, sd = .1), length(genes), 6,
                dimnames = list(genes, NULL))
    x[, 4:6] <- x[, 4:6] + 2
    expect_warning(out <- leviGraphClusterInference(x,
        c("C", "C", "C", "T", "T", "T"),
        system.file("extdata", "hub_network.dat", package = "levi"),
        test = "T", control = "C", threshold = 1, n_perm = 3,
        blocks = c("d1", "d2", "d3", "d1", "d2", "d3"), seed = 1),
        "smallest attainable p-value")
    expect_true(nrow(out$regions$summary) > 0)
    expect_true(all(out$regions$summary$Nodes >= 1))
    expect_true(all(out$regions$summary$PSpatial >= .05))
    expect_true(out$exact)
})

test_that("network statistics provide TFCE, Moran, spectrum and rewiring nulls", {
    skip_if_not_installed("limma")
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    genes <- c("HUB", paste0("N", 1:8))
    x <- matrix(rnorm(genes |> length() * 4), length(genes), 4,
                dimnames = list(genes, NULL))
    expect_warning(tfce <- leviGraphTFCEInference(x, c("C", "C", "T", "T"), net,
        test = "T", control = "C", n_perm = 2, seed = 1),
        "smallest attainable p-value")
    expect_equal(nrow(tfce$statistic), length(genes))
    expect_true(all(c("PGlobal", "GlobalSignificant") %in% names(tfce$statistic)))
    score <- setNames(rnorm(length(genes)), genes)
    expect_true(is.finite(leviGraphMoran(score, net, n_perm = 3)$global$MoranI))
    expect_true(is.finite(leviGraphSpectrum(score, net, n_perm = 3)$LaplacianEnergy))
    expect_true(is.list(leviGraphRewiringInference(score, net, n_perm = 2)))
})

test_that("single-cell TFCE uses joint donor permutations across cell types", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8)); donor <- rep(paste0("d", 1:4), each = 4)
    type <- rep(rep(c("A", "B"), each = 2), 4)
    counts <- matrix(rpois(length(genes) * length(donor), 10), length(genes),
        dimnames = list(genes, NULL))
    expect_warning(out <- leviSingleCellTFCEInference(counts, donor, type,
        c(d1 = "C", d2 = "C", d3 = "T", d4 = "T"),
        system.file("extdata", "hub_network.dat", package = "levi"), min_cells = 2,
        permutation_method = "exact", seed = 1), "smallest attainable p-value")
    expect_true(out$exact)
    expect_equal(length(out$null_global), 5L)
    expect_true(all(c("PGlobal", "GlobalSignificant") %in% names(out$results$A)))
})

test_that("single-cell TFCE respects paired donor-by-condition pseudobulks", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8))
    donor <- rep(paste0("d", 1:4), each = 4)
    condition <- rep(rep(c("ctrl", "stim"), each = 2), 4)
    type <- rep(c("A", "B"), times = 8)
    counts <- matrix(rpois(length(genes) * length(donor), 20), length(genes),
        dimnames = list(genes, NULL))
    expect_warning(out <- leviSingleCellTFCEInference(counts, donor, type,
        condition,
        system.file("extdata", "hub_network.dat", package = "levi"), min_cells = 1,
        permutation_method = "exact", seed = 1), "smallest attainable p-value")
    expect_true(out$exact)
    expect_equal(length(out$null_global), 15L)
    expect_match(out$method, "within-donor")
})

test_that("pseudobulk sums cells by donor and cell type", {
    x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
                dimnames = list(c("A", "B"), NULL))
    out <- leviPseudobulk(x, donor = c("d1", "d1", "d2"),
        cell_type = c("T", "T", "B"), min_cells = 1)
    expect_equal(unname(out$counts[, "d1_T"]), c(4, 6))
    expect_equal(unname(out$counts[, "d2_B"]), c(5, 6))
})

test_that("bulk RNA-seq graph inference returns regional P-values", {
    skip_if_not_installed("DESeq2")
    set.seed(3)
    genes <- c("HUB", paste0("N", 1:8))
    counts <- matrix(rpois(length(genes) * 6, 30), length(genes), 6,
        dimnames = list(genes, paste0("s", 1:6)))
    counts[, 4:6] <- counts[, 4:6] + 30L
    out <- leviBulkGraphInference(counts, c("C", "C", "C", "T", "T", "T"),
        system.file("extdata", "hub_network.dat", package = "levi"),
        test = "T", control = "C", threshold = 1, n_perm = 3, seed = 1)
    expect_true(all(out$regions$summary$PSpatial >= .25))
})

test_that("pseudobulk expression drops unexpressed genes but keeps network nodes", {
    skip_if_not_installed("limma")
    genes <- c("HUB", paste0("N", 1:8)); donor <- rep(paste0("d", 1:4), each = 4)
    condition <- rep(rep(c("ctrl", "stim"), each = 2), 4)
    type <- rep("A", 16)
    set.seed(3)
    counts <- matrix(rpois(length(genes) * length(donor), 20), length(genes),
        dimnames = list(genes, NULL))
    counts["N8", ] <- 0                           # a network node never detected
    zeros <- matrix(0L, 2000, ncol(counts),        # a droplet-style zero block
        dimnames = list(paste0("Z", 1:2000), NULL))
    counts <- rbind(counts, zeros)
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    # Fitting the 2000 all-zero rows made limma warn "eBayes unreliable" on
    # every permutation; they are filtered out now, so no such warning.
    warns <- character()
    out <- withCallingHandlers(
        leviSingleCellTFCEInference(counts, donor, type, condition, net,
            min_cells = 1, permutation_method = "exact", seed = 1),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
        })
    expect_false(any(grepl("eBayes unreliable", warns)))
    expect_true(any(grepl("smallest attainable p-value", warns)))
    # Every network node, including the undetected one, keeps a statistic.
    expect_setequal(out$results$A$Gene, genes)
    pb <- leviPseudobulk(counts, donor, type, condition, min_cells = 1)
    keys <- paste(pb$donor, pb$condition, sep = "::")
    expr <- levi:::.pseudobulkExpression(pb, keys, "A", paired = TRUE,
        normalize = "none", keep_genes = genes)$A
    expect_true(all(genes %in% rownames(expr)))
    expect_false(any(grepl("^Z", rownames(expr))))
})
