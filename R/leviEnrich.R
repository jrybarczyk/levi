#' @title leviEnrich
#' @description Run gene-set enrichment on the high-scoring (peaks) and
#' low-scoring (valleys) regions returned by \code{levi()}. Uses
#' \code{clusterProfiler} for GO and KEGG enrichment.
#' @param result A single \code{levi()} result (list with \code{$scores}).
#' @param top_n Integer. Number of top-ranked genes (by landscape score) to
#' test as the "over-expressed" gene set. Default is \code{50}.
#' @param bottom_n Integer. Number of bottom-ranked genes to test as the
#' "under-expressed" gene set. Default is \code{50}.
#' @param organism Character. KEGG organism code (e.g., \code{"hsa"} for
#' human, \code{"mmu"} for mouse). Also used to select the OrgDb for GO
#' enrichment: pass the \strong{OrgDb package name} as \code{orgdb}. See
#' \code{clusterProfiler::enrichGO()} for details.
#' @param orgdb Character. OrgDb package name for GO enrichment (e.g.,
#' \code{"org.Hs.eg.db"} for human). If \code{NULL}, GO enrichment is
#' skipped.
#' @param keytype Character. Key type of the gene identifiers in
#' \code{result$scores$Gene} (e.g., \code{"SYMBOL"}, \code{"ENTREZID"}).
#' Default is \code{"SYMBOL"}.
#' @param pval_cutoff Numeric. Adjusted p-value threshold. Default \code{0.05}.
#' @param types Character vector. Which analyses to run. Any subset of
#' \code{c("GO_BP","GO_MF","GO_CC","KEGG")}. Default runs GO_BP and KEGG.
#' @param universe Character vector. Background gene universe. If \code{NULL}
#' (default), all genes in \code{result$scores$Gene} are used.
#' For KEGG, both genes and universe are converted to ENTREZID using orgdb;
#' with keytype = "ENTREZID", no conversion is needed.
#' @return A named list with elements \code{over} and \code{under}, each
#' a named list of enrichment result objects (one per requested type).
#' Plots a dot-plot for each non-empty result automatically.
#' @details
#' Requires \code{clusterProfiler} (Bioconductor). Install with
#' \code{BiocManager::install("clusterProfiler")}. For GO enrichment also
#' install the relevant OrgDb, e.g.
#' \code{BiocManager::install("org.Hs.eg.db")}.
#' @examples
#' hub_n <- system.file("extdata", "hub_network.dat", package = "levi")
#' hub_e <- system.file("extdata", "hub_expression.dat", package = "levi")
#' res <- levi(expressionInput         = hub_e,
#'             networkCoordinatesInput = hub_n,
#'             fileTypeInput           = "dat",
#'             geneSymbolInput         = "ID",
#'             readExpColumn           = readExpColumn("Test-Control"),
#'             resolutionValueInput    = 10,
#'             smoothValueInput        = 5)
#' \dontrun{
#' # Requires clusterProfiler and the org.Hs.eg.db annotation package,
#' # plus real gene symbols in the network.
#' enrich <- leviEnrich(res,
#'                      top_n    = 50,
#'                      orgdb    = "org.Hs.eg.db",
#'                      organism = "hsa",
#'                      keytype  = "SYMBOL")
#' enrich$over$GO_BP
#' enrich$under$KEGG
#' }
#' @export
leviEnrich <- function(result,
                        top_n      = 50L,
                        bottom_n   = 50L,
                        organism   = "hsa",
                        orgdb      = NULL,
                        keytype    = "SYMBOL",
                        pval_cutoff = 0.05,
                        types      = c("GO_BP", "KEGG"),
                        universe   = NULL) {

    if (!requireNamespace("clusterProfiler", quietly = TRUE))
        stop("clusterProfiler is required. ",
             "Install with: BiocManager::install('clusterProfiler')")

    if (is.null(result$scores))
        stop("result$scores not found. Re-run levi() with levi >= 2.0.0.")

    scores <- result$scores
    all_genes <- scores$Gene
    bg <- if (is.null(universe)) all_genes else universe

    top_genes    <- head(scores$Gene, top_n)
    bottom_genes <- tail(scores$Gene, bottom_n)

    kegg_map <- NULL
    if ("KEGG" %in% types) {
        ids <- unique(as.character(c(all_genes, bg)))
        if (keytype == "ENTREZID") {
            kegg_map <- data.frame(input = ids, ENTREZID = ids)
        } else {
            if (is.null(orgdb) || !requireNamespace(orgdb, quietly = TRUE))
                stop("An installed orgdb is required to convert genes and ",
                     "universe to ENTREZID for KEGG.")
            kegg_map <- clusterProfiler::bitr(ids, fromType = keytype,
                toType = "ENTREZID", OrgDb = get(orgdb, asNamespace(orgdb)))
            names(kegg_map)[names(kegg_map) == keytype] <- "input"
        }
    }
    kegg_ids <- function(genes) unique(as.character(
        kegg_map$ENTREZID[kegg_map$input %in% genes & !is.na(kegg_map$ENTREZID)]))
    kegg_bg <- if (!is.null(kegg_map)) kegg_ids(bg) else NULL
    if ("KEGG" %in% types && !length(kegg_bg))
        stop("No universe identifiers could be mapped to ENTREZID for KEGG.")

    .run_one <- function(genes, label) {
        out <- list()
        message("\n[leviEnrich] Enriching ", label, " set (",
                length(genes), " genes)...")

        # GO enrichment ------------------------------------------------------------------
        if (!is.null(orgdb)) {
            if (!requireNamespace(orgdb, quietly = TRUE))
                warning("OrgDb package '", orgdb, "' not installed. ",
                        "Skipping GO enrichment.")
            else {
                db <- get(orgdb, envir = asNamespace(orgdb))
                go_subont <- c(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC")
                for (ont_key in intersect(types, c("GO_BP","GO_MF","GO_CC"))) {
                    res <- tryCatch(
                        clusterProfiler::enrichGO(
                            gene          = genes,
                            OrgDb         = db,
                            keyType       = keytype,
                            ont           = go_subont[[ont_key]],
                            pAdjustMethod = "BH",
                            pvalueCutoff  = pval_cutoff,
                            universe      = bg,
                            readable      = (keytype != "SYMBOL")),
                        error = function(e) { message(e$message); NULL })
                    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
                        methods::show(clusterProfiler::dotplot(res,
                            title = paste(label, ont_key, result$comparison)))
                        out[[ont_key]] <- res
                    } else {
                        message("[leviEnrich] No ", ont_key,
                                " enrichment for ", label, " set.")
                    }
                }
            }
        }

        # KEGG enrichment ---------------------------------------------------------------
        if ("KEGG" %in% types) {
            kegg_genes <- intersect(kegg_ids(genes), kegg_bg)
            if (!length(kegg_genes)) {
                message("[leviEnrich] No KEGG genes in the supplied universe.")
                return(out)
            }

            res <- tryCatch(
                clusterProfiler::enrichKEGG(
                    gene          = kegg_genes,
                    universe      = kegg_bg,
                    keyType       = "ncbi-geneid",
                    organism      = organism,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = pval_cutoff),
                error = function(e) { message(e$message); NULL })
            if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
                methods::show(clusterProfiler::dotplot(res,
                    title = paste(label, "KEGG", result$comparison)))
                out[["KEGG"]] <- res
            } else {
                message("[leviEnrich] No KEGG enrichment for ", label, " set.")
            }
        }
        out
    }

    list(
        over  = .run_one(top_genes,    "over-expressed"),
        under = .run_one(bottom_genes, "under-expressed")
    )
}
