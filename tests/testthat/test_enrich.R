library(levi)

# leviEnrich() was one of the two functions with only its offline validation
# under test, because enrichment needs an annotation database. With
# org.Hs.eg.db available the GO branch runs locally, with no network.

skip_without_annotation <- function() {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
}

# A result shaped like levi()'s, with a cell-cycle set at the top and a set of
# liver-specific genes at the bottom.
fake_result <- function() {
    up <- c("CDK1", "CCNB1", "CCNA2", "CDC20", "BUB1", "PLK1", "AURKA",
            "AURKB", "MAD2L1", "TTK", "CENPE", "KIF11", "TOP2A", "MKI67",
            "BIRC5")
    dn <- c("ALB", "APOA1", "APOB", "CYP3A4", "F2", "FGA", "FGB",
            "SERPINA1", "TF", "TTR")

    structure(
        list(comparison = "Tumor-Normal",
             scores = data.frame(
                 Gene = c(up, dn),
                 LandscapeScore = c(seq(1, 0.7, length.out = length(up)),
                                    seq(0.3, 0, length.out = length(dn))),
                 stringsAsFactors = FALSE)),
        class = "levi_result")
}


test_that("leviEnrich: refuses a result with no scores", {
    skip_if_not_installed("clusterProfiler")
    expect_error(leviEnrich(list(plot = NULL)), "result\\$scores not found")
})

# GO enrichment costs about half a minute, most of it loading the annotation
# database, so it runs once and both tests below read the same result.
enrich_cache <- local({
    cached <- NULL
    function() {
        if (is.null(cached))
            cached <<- suppressMessages(leviEnrich(fake_result(), top_n = 15,
                bottom_n = 10, orgdb = "org.Hs.eg.db", types = "GO_BP",
                pval_cutoff = 0.05))
        cached
    }
})

test_that("leviEnrich: returns one entry per tail", {
    skip_without_annotation()
    expect_named(enrich_cache(), c("over", "under"))
})

test_that("leviEnrich: finds the biology planted in the top of the ranking", {
    skip_without_annotation()

    res <- enrich_cache()
    expect_true("GO_BP" %in% names(res$over))

    terms <- as.data.frame(res$over$GO_BP)
    expect_gt(nrow(terms), 0)

    # The over-expressed set is made of cell-cycle genes, so the enrichment
    # has to say so somewhere in its terms.
    expect_true(any(grepl("cell cycle|mitotic|division", terms$Description,
                          ignore.case = TRUE)))
})

test_that("leviEnrich: warns instead of failing when the OrgDb is absent", {
    skip_if_not_installed("clusterProfiler")

    expect_warning(
        suppressMessages(leviEnrich(fake_result(), top_n = 5, bottom_n = 5,
            orgdb = "org.NoSuchOrganism.eg.db", types = "GO_BP")),
        "not installed")
})

test_that("KEGG receives the supplied background in the same identifier space", {
    skip_without_annotation()
    calls <- list()
    testthat::local_mocked_bindings(enrichKEGG = function(...) {
        calls[[length(calls) + 1L]] <<- list(...)
        NULL
    }, .package = "clusterProfiler")
    res <- fake_result()
    bg <- c(res$scores$Gene, "TP53", "EGFR")
    suppressMessages(leviEnrich(res, top_n = 3, bottom_n = 3, types = "KEGG",
        orgdb = "org.Hs.eg.db", universe = bg))
    mapping <- suppressMessages(clusterProfiler::bitr(bg, fromType = "SYMBOL",
        toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
    expect_length(calls, 2)
    for (call in calls) {
        expect_setequal(call$universe, mapping$ENTREZID)
        expect_identical(call$keyType, "ncbi-geneid")
        expect_true(all(call$gene %in% call$universe))
    }
    calls <- list()
    res$scores$Gene <- as.character(seq_len(nrow(res$scores)))
    suppressMessages(leviEnrich(res, types = "KEGG", keytype = "ENTREZID",
                               top_n = 3, bottom_n = 3, universe = as.character(1:30)))
    expect_length(calls, 2)
    expect_setequal(calls[[1]]$universe, as.character(1:30))
    expect_error(leviEnrich(fake_result(), types = "KEGG"), "orgdb")
})
