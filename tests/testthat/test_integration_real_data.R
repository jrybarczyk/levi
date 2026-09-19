# Integration tests: replay the real-data analyses of the supplement and
# check that the qualitative conclusions still hold.
#
# They are skipped unless LEVI_INTEGRATION is set, because the GEO and Kang
# 2018 cases need raw files that are not shipped (LEVI_REAL_DATA) and the
# whole file takes a few minutes. Run them from the package root with
#
#   LEVI_INTEGRATION=1 LEVI_REAL_DATA=../real_data_tests \
#       Rscript -e 'testthat::test_local(filter = "integration")'
#
# Reference values (n_perm = 999, seed = 1) are those recorded in the
# supplement; the bounds below leave room for Monte Carlo error so the
# tests document conclusions, not digits.

# -- airway (RNA-seq, shipped with the package) ------------------------------

test_that("airway: node-label, replicate and TFCE agree on the SPARCL1 module", {
    skip_unless_integration()
    aw <- function(f) system.file("extdata", "airway", f, package = "levi")
    genes <- read.delim(aw("airway_dex_genes.tsv"))
    logcpm <- as.matrix(read.delim(aw("airway_dex_logcpm.tsv"), row.names = 1))
    samples <- read.delim(aw("airway_dex_samples.tsv"))
    nodes <- read.delim(aw("airway_string_nodes.tsv"))
    edges <- read.delim(aw("airway_string_edges.tsv"))
    pdf(NULL); on.exit(dev.off())

    set.seed(1)
    land <- do.call(levi, c(list(expressionInput = genes,
        networkCoordinatesInput = nodes, networkInteractionsInput = edges,
        geneSymbolInput = "Symbol", signal_mode = "logfc",
        readExpColumn = readExpColumn("log2FoldChange-log2FoldChange")),
        landscape_args))
    s <- land$regions$summary
    # supplement: over_02, 3137 cells, p = 0.014
    expect_equal(s$Region[1], "over_02")
    expect_true(s$Significant[1])
    expect_lt(s$PSpatial[1], 0.04)
    expect_false(any(s$Significant[-1]))

    rep <- do.call(leviReplicateInference, c(list(logcpm, samples$dex,
        test = "trt", control = "untrt", networkCoordinatesInput = nodes,
        networkInteractionsInput = edges, seed = 1), landscape_args))
    expect_true(rep$metadata$permutation_exact)
    expect_equal(rep$metadata$possible_permutations, 70)
    r <- rep$regions$summary
    expect_equal(r$Region[1], "over_02")
    expect_equal(r$PSpatial[1], 2 / 70)    # exact test: deterministic
    expect_equal(sum(r$Significant), 1L)

    tf <- leviGraphTFCEInference(logcpm, samples$dex, nodes, edges,
        fileTypeInput = "stg", test = "trt", control = "untrt",
        n_perm = 999, seed = 1)
    expect_true(tf$exact)
    st <- tf$statistic
    expect_equal(st$Gene[which.min(st$PGlobal)], "SPARCL1")
    expect_equal(min(st$PGlobal), 1 / 70)
    expect_true(all(st$GlobalSignificant[st$Gene %in% c("SPARCL1", "DUSP1",
        "COL1A1", "SOX4")]))
    expect_gte(sum(st$GlobalSignificant), 8)
})

# -- GSE10072 (lung adenocarcinoma microarray) -------------------------------

test_that("GSE10072: node-label misses the tumour module that sample-label tests detect", {
    skip_unless_integration()
    g <- load_gse10072()
    net <- read_saved_network("gse10072")
    expect_equal(nrow(net$nodes), 72L)
    expect_true(all(net$nodes$name %in% g$tt$Symbol))
    pdf(NULL); on.exit(dev.off())

    sig <- g$tt[g$tt$Symbol %in% net$nodes$name, c("Symbol", "logFC")]
    set.seed(1)
    land <- do.call(levi, c(list(expressionInput = sig,
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        geneSymbolInput = "Symbol", signal_mode = "logfc",
        readExpColumn = readExpColumn("logFC-logFC")), landscape_args))
    s <- land$regions$summary
    # supplement: under_01 (2675 cells) p = 0.92 -- the whole network responds
    expect_equal(s$Region[1], "under_01")
    expect_gt(s$PSpatial[1], 0.5)
    expect_false(any(s$Significant))

    gm <- gene_matrix(g$expr, g$tt, net$nodes$name)
    expect_equal(dim(gm), c(72L, 107L))
    expect_equal(as.vector(table(g$group)[c("Tumor", "Normal")]), c(58L, 49L))
    rep <- do.call(leviReplicateInference, c(list(gm, g$group,
        test = "Tumor", control = "Normal", networkCoordinatesInput = net$nodes,
        networkInteractionsInput = net$edges, seed = 1), landscape_args))
    r <- rep$regions$summary
    # supplement: under_01 p = 0.001, under_02 p = 0.013
    expect_equal(r$Region[1], "under_01")
    expect_true(r$Significant[1])
    expect_lte(region_p(r, "under_01"), 0.005)
    expect_lt(region_p(r, "under_02"), 0.05)
    # the region masses come from the same observed landscape
    expect_equal(r$Mass[match(s$Region, r$Region)], s$Mass, tolerance = 1e-8)

    paired <- names(which(table(g$patient) == 2))
    keep <- g$patient %in% paired
    expect_equal(length(paired), 33L)
    blk <- do.call(leviReplicateInference, c(list(gm[, keep], g$group[keep],
        test = "Tumor", control = "Normal", blocks = g$patient[keep],
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        seed = 1), landscape_args))
    b <- blk$regions$summary
    # supplement: under_01 p = 0.001, under_02 p = 0.17 (paired by patient)
    expect_lte(region_p(b, "under_01"), 0.005)
    expect_gt(region_p(b, "under_02"), 0.05)

    cl <- leviGraphClusterInference(gm, g$group, net$nodes, net$edges,
        fileTypeInput = "stg", test = "Tumor", control = "Normal",
        n_perm = 999, seed = 1)
    c2 <- cl$regions$summary
    big <- c2[which.max(c2$Nodes), ]
    # supplement: one under cluster of 62 nodes (CD36 peak), p = 0.001
    expect_equal(big$Direction, "under")
    expect_equal(big$Nodes, 62L)
    expect_true(big$Significant)
})

test_that("GSE10072 smoking: the normal-lung signature network separates current from never smokers", {
    skip_unless_integration()
    g <- load_gse10072()
    net <- read_saved_network("gse10072_smoking")
    expect_equal(nrow(net$nodes), 36L)
    pdf(NULL); on.exit(dev.off())

    design <- model.matrix(~ 0 + factor(g$strata))
    colnames(design) <- levels(factor(g$strata))
    fit <- limma::lmFit(g$expr, design)
    con <- limma::makeContrasts(Normal_Current - Normal_Never, levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, con))
    tt <- g$collapse(limma::topTable(fit2, coef = 1, number = Inf, sort.by = "P"))
    expect_true(all(net$nodes$name %in% tt$Symbol))

    keep <- g$strata %in% c("Normal_Current", "Normal_Never")
    expect_equal(as.vector(table(g$strata[keep])), c(16L, 15L))
    gm <- gene_matrix(g$expr, tt, net$nodes$name)[, keep]
    rep <- do.call(leviReplicateInference, c(list(gm, g$strata[keep],
        test = "Normal_Current", control = "Normal_Never",
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        seed = 1), landscape_args))
    r <- rep$regions$summary
    # supplement: under_01 p = 0.001, over_02 p = 0.002, over_03 p = 0.017
    expect_true(all(r$Significant[r$Region %in% c("under_01", "over_02")]))
    expect_gte(sum(r$Significant), 2L)

    cl <- leviGraphClusterInference(gm, g$strata[keep], net$nodes, net$edges,
        fileTypeInput = "stg", test = "Normal_Current", control = "Normal_Never",
        n_perm = 999, seed = 1)
    c2 <- cl$regions$summary
    # supplement: CYP1B1 over-cluster (9 nodes) and ETV5 under-cluster (11
    # nodes), both p = 0.001
    cyp <- c2[c2$PeakGene == "CYP1B1", ]; etv <- c2[c2$PeakGene == "ETV5", ]
    expect_equal(c(cyp$Direction, etv$Direction), c("over", "under"))
    expect_equal(c(cyp$Nodes, etv$Nodes), c(9L, 11L))
    expect_true(cyp$Significant && etv$Significant)
})

test_that("GSE10072 focal adhesion (KEGG a priori): sample-label tests detect the pathway node-label misses", {
    skip_unless_integration()
    g <- load_gse10072()
    net <- read_saved_network("gse10072_hsa04510")
    expect_equal(c(nrow(net$nodes), nrow(net$edges)), c(186L, 1669L))
    pdf(NULL); on.exit(dev.off())

    de <- g$tt[g$tt$Symbol %in% net$nodes$name, c("Symbol", "logFC", "adj.P.Val")]
    expect_gte(sum(de$adj.P.Val < 0.05), 120)     # supplement: 134 of 186
    set.seed(1)
    land <- do.call(levi, c(list(expressionInput = de[, 1:2],
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        geneSymbolInput = "Symbol", signal_mode = "logfc",
        readExpColumn = readExpColumn("logFC-logFC")), landscape_args))
    s <- land$regions$summary
    expect_equal(s$Region[1], "under_01")
    expect_gt(s$PSpatial[1], 0.1)                  # supplement: 0.24
    expect_false(any(s$Significant))

    gm <- gene_matrix(g$expr, g$tt, net$nodes$name)
    rep <- do.call(leviReplicateInference, c(list(gm, g$group,
        test = "Tumor", control = "Normal", networkCoordinatesInput = net$nodes,
        networkInteractionsInput = net$edges, seed = 1), landscape_args))
    r <- rep$regions$summary
    expect_lte(region_p(r, "under_01"), 0.005)     # supplement: 0.001
    expect_lte(region_p(r, "under_02"), 0.01)      # supplement: 0.001

    tf <- leviGraphTFCEInference(gm, g$group, net$nodes, net$edges,
        fileTypeInput = "stg", test = "Tumor", control = "Normal",
        n_perm = 999, seed = 1)
    st <- tf$statistic
    # supplement: 127 of 186 genes with PGlobal <= 0.05 (48 over, 79 under);
    # VWF, SPP1 and CAV1 lead
    expect_gte(sum(st$GlobalSignificant), 110)
    expect_true(all(st$GlobalSignificant[st$Gene %in% c("VWF", "SPP1", "CAV1")]))
    expect_lt(st$TFCE[st$Gene == "VWF"], 0)
    expect_gt(st$TFCE[st$Gene == "SPP1"], 0)
    expect_gt(sum(st$GlobalSignificant & st$TFCE < 0),
              sum(st$GlobalSignificant & st$TFCE > 0))
})

# -- Kang 2018 (single cell, CD14+ monocytes, IFN-beta) ----------------------

test_that("Kang monocytes: node-label p ~ 1 on the STRING network, donor-level tests detect the response", {
    skip_unless_integration()
    k <- load_kang18()
    ct <- "CD14+ Monocytes"
    tt <- kang_pseudobulk_table(k, ct)
    expect_gt(sum(tt$adj.P.Val < 0.05), 1000)
    net <- read_saved_network("kang")
    expect_equal(nrow(net$nodes), 83L)
    expect_true(all(net$nodes$name %in% tt$Symbol))
    pdf(NULL); on.exit(dev.off())

    de <- tt[tt$Symbol %in% net$nodes$name, c("Symbol", "logFC")]
    set.seed(1)
    land <- do.call(levi, c(list(expressionInput = de,
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        geneSymbolInput = "Symbol", signal_mode = "logfc",
        readExpColumn = readExpColumn("logFC-logFC")), landscape_args))
    s <- land$regions$summary
    # supplement: a single over_01 region covering the network, p = 0.88
    expect_equal(s$Region[1], "over_01")
    expect_gt(s$PSpatial[1], 0.5)

    reg <- do.call(leviSingleCellRegionalInference, c(list(k$counts, k$donor,
        k$cell_type, k$condition, net$nodes, net$edges, cell_types = ct,
        seed = 1), landscape_args))
    expect_true(reg$exact)
    expect_equal(reg$possible_permutations, 256)
    cm <- reg$cluster_mass[[ct]]$summary
    # supplement: over_01 (71 nodes, NT5C3A peak) p = 2/256, nothing else
    expect_equal(cm$Region[1], "over_01")
    expect_equal(cm$Nodes[1], 71L)
    expect_equal(cm$PSpatial[1], 2 / 256)
    expect_equal(sum(cm$Significant), 1L)
    # the landscape region covering the whole network is not significant
    # under the same exact null (18/256): on this DE-derived network only
    # the graph-cluster statistic separates the response from the null
    ls <- reg$landscapes[[ct]]$regions$summary
    expect_equal(ls$Region[1], "over_01")
    expect_equal(ls$PSpatial[1], 18 / 256)

    tf <- leviSingleCellTFCEInference(k$counts, k$donor, k$cell_type,
        k$condition, net$nodes, net$edges, fileTypeInput = "stg",
        cell_types = ct, n_perm = 999, seed = 1)
    expect_true(tf$exact)
    st <- tf$results[[ct]]
    expect_gte(sum(st$GlobalSignificant), 60)      # supplement: 74 of 83
    expect_true(all(st$GlobalSignificant[st$Gene %in% c("ISG15", "IFIT1",
        "CXCL10")]))
})

test_that("Kang JAK-STAT (KEGG a priori): node-label n.s., regional and TFCE detect the pathway", {
    skip_unless_integration()
    k <- load_kang18()
    ct <- "CD14+ Monocytes"
    tt <- kang_pseudobulk_table(k, ct)
    net <- read_saved_network("kang_hsa04630")
    expect_equal(c(nrow(net$nodes), nrow(net$edges)), c(43L, 171L))
    pdf(NULL); on.exit(dev.off())

    de <- tt[tt$Symbol %in% net$nodes$name, c("Symbol", "logFC", "adj.P.Val")]
    expect_gte(mean(de$adj.P.Val < 0.05), 0.6)    # supplement: 70 %
    set.seed(1)
    land <- do.call(levi, c(list(expressionInput = de[, 1:2],
        networkCoordinatesInput = net$nodes, networkInteractionsInput = net$edges,
        geneSymbolInput = "Symbol", signal_mode = "logfc",
        readExpColumn = readExpColumn("logFC-logFC")), landscape_args))
    s <- land$regions$summary
    expect_equal(s$Region[1], "over_01")
    expect_gt(s$PSpatial[1], 0.05)                 # supplement: 0.18
    expect_false(any(s$Significant))

    reg <- do.call(leviSingleCellRegionalInference, c(list(k$counts, k$donor,
        k$cell_type, k$condition, net$nodes, net$edges, cell_types = ct,
        seed = 1), landscape_args))
    cm <- reg$cluster_mass[[ct]]$summary
    # supplement: over_01 with 23 genes (SOCS1, STAT1, STAT2, JAK2), p = 2/256
    expect_equal(cm$Region[1], "over_01")
    expect_equal(cm$PSpatial[1], 2 / 256)
    genes <- strsplit(cm$Genes[1], ";")[[1]]
    expect_equal(length(genes), 23L)
    expect_true(all(c("SOCS1", "STAT1", "STAT2", "JAK2") %in% genes))

    tf <- leviSingleCellTFCEInference(k$counts, k$donor, k$cell_type,
        k$condition, net$nodes, net$edges, fileTypeInput = "stg",
        cell_types = ct, n_perm = 999, seed = 1)
    st <- tf$results[[ct]]
    expect_gte(sum(st$PGlobal <= 0.05), 20)        # supplement: 23 of 43
    expect_true(all(st$PGlobal[st$Gene %in% c("STAT1", "STAT2", "SOCS1")] <= 0.05))
})
