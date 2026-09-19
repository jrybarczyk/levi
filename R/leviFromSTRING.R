#' @title leviFromSTRING
#' @description Build a levi-compatible network directly from the STRING
#' protein interaction database using the \pkg{STRINGdb} Bioconductor package.
#' Node coordinates are computed from the interaction graph via an igraph
#' layout algorithm, so no manual network file download is required.
#'
#' @param genes Character vector of gene symbols (or Entrez IDs / Ensembl
#'   IDs, matched by STRINGdb's mapping step). Genes not found in STRING are
#'   silently dropped and reported via \code{message()}.
#' @param species Integer. NCBI taxonomy ID. Common values:
#'   \code{9606} (human, default), \code{10090} (mouse),
#'   \code{10116} (rat), \code{7227} (Drosophila), \code{4932} (yeast),
#'   \code{6239} (C. elegans).
#' @param score_threshold Integer (0-1000). Minimum STRING combined score
#'   to retain an interaction. \code{400} = medium confidence (default),
#'   \code{700} = high confidence, \code{900} = highest confidence.
#' @param network_type Character. STRING network type:
#'   \code{"full"} (default, all channels) or \code{"physical"} (only
#'   physical protein-protein interactions).
#' @param layout Character. igraph layout algorithm used to assign (x, y)
#'   coordinates to nodes. Options: \code{"fr"} Fruchterman-Reingold
#'   (default, good for most sizes), \code{"kk"} Kamada-Kawai (small
#'   networks, <=100 nodes), \code{"lgl"} Large Graph Layout (>500 nodes),
#'   \code{"dh"} Davidson-Harel, \code{"circle"}.
#' @param version Character. STRING database version. Default \code{"11.5"}.
#' @param input_directory Character. Local cache directory for STRINGdb files.
#'   Defaults to a temporary directory (re-downloaded each session). Set a
#'   persistent path (e.g. \code{"~/.stringdb_cache"}) to avoid
#'   re-downloading.
#' @param id_col Character. Column in \code{genes} to map when \code{genes}
#'   is a data.frame. Ignored when \code{genes} is a character vector.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{$nodes}}{data.frame with columns \code{name} (gene symbol),
#'     \code{x}, \code{y} - ready for \code{networkCoordinatesInput}.}
#'   \item{\code{$edges}}{data.frame with columns \code{V1}, \code{V2}
#'     (interacting gene symbols) - ready for \code{networkInteractionsInput}.}
#'   \item{\code{$graph}}{The igraph object (invisible, for debugging or
#'     alternative layout experiments).}
#' }
#' Pass \code{$nodes} and \code{$edges} directly to \code{levi()} with
#' \code{fileTypeInput = "stg"}.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item STRINGdb maps gene symbols to STRING protein IDs.
#'   \item Interactions above \code{score_threshold} are retrieved.
#'   \item An igraph subnetwork is built and isolated nodes are retained.
#'   \item A 2-D layout is computed (coordinates scaled to \[0, 100\]).
#'   \item Gene symbols are recovered from the STRING node metadata.
#' }
#'
#' \strong{Layout choice guide:}
#' \tabular{ll}{
#'   \code{"fr"}     \tab Default. Force-directed, balanced for 50-500 nodes.\cr
#'   \code{"kk"}     \tab Kamada-Kawai. Better aesthetics for <100 nodes.\cr
#'   \code{"lgl"}    \tab Large Graph Layout. Use for >500 nodes.\cr
#'   \code{"dh"}     \tab Davidson-Harel. Slow but high quality.\cr
#'   \code{"circle"} \tab Ring layout. Use for pathway-like linear chains.\cr
#' }
#'
#' \strong{Tip - seed for reproducibility:}
#' Force-directed layouts are stochastic. Set \code{set.seed()} before
#' calling \code{leviFromSTRING()} when exact coordinates must be reproduced.
#'
#' @examples
#' # All examples below query the STRING database over the network and
#' # download a species-wide alias file (~20 MB), so they are not run
#' # during automated checks.
#' \dontrun{
#' genes <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN")
#' set.seed(42)
#' net <- leviFromSTRING(genes, species = 9606, score_threshold = 700)
#' head(net$nodes)
#'
#' # Human MAPK signalling genes
#' mapk_genes <- c("EGFR", "KRAS", "BRAF", "MAP2K1", "MAPK1", "MAPK3",
#'                 "RPS6KA1", "MYC", "JUN", "FOS")
#' set.seed(42)
#' net <- leviFromSTRING(mapk_genes, species = 9606, score_threshold = 400)
#' levi(
#'     expressionInput          = my_expression_df,   # your own data.frame
#'     networkCoordinatesInput  = net$nodes,
#'     networkInteractionsInput = net$edges,
#'     fileTypeInput            = "stg",
#'     geneSymbolInput          = "GeneID",
#'     readExpColumn            = readExpColumn("Tumor-Normal"),
#'     signal_mode              = "logfc"
#' )
#' }
#'
#' @author Jose Rafael Pilan (rafael.pilan@unesp.br)
#' @export
leviFromSTRING <- function(genes,
                            species         = 9606,
                            score_threshold = 400,
                            network_type    = c("full", "physical"),
                            layout          = c("fr", "kk", "lgl", "dh", "circle"),
                            version         = "11.5",
                            input_directory = tempdir(),
                            id_col          = NULL) {

    network_type <- match.arg(network_type)
    layout       <- match.arg(layout)

    if (!requireNamespace("STRINGdb", quietly = TRUE))
        stop("Package 'STRINGdb' is required. ",
             "Install with: BiocManager::install('STRINGdb')", call. = FALSE)
    if (!requireNamespace("igraph", quietly = TRUE))
        stop("Package 'igraph' is required. ",
             "Install with: install.packages('igraph')", call. = FALSE)
    if (!is.character(input_directory) || length(input_directory) != 1L ||
        is.na(input_directory) || !nzchar(input_directory))
        stop("'input_directory' must be one non-empty directory path.", call. = FALSE)
    if (!dir.exists(input_directory))
        dir.create(input_directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(input_directory))
        stop("Could not create STRING cache directory: ", input_directory,
             call. = FALSE)

    # --- 1. Accept data.frame or character vector ----------------------------
    if (is.data.frame(genes) || is.matrix(genes)) {
        if (is.null(id_col))
            stop("When 'genes' is a data.frame, supply 'id_col' with the ",
                 "column name containing gene identifiers.", call. = FALSE)
        gene_vec <- as.character(genes[[id_col]])
    } else {
        gene_vec <- as.character(genes)
    }
    gene_vec <- unique(gene_vec[!is.na(gene_vec) & nchar(gene_vec) > 0])

    # --- 2. Connect to STRINGdb ----------------------------------------------
    message("Connecting to STRING v", version, " (species = ", species,
            ", score_threshold = ", score_threshold, ") ...")
    sdb <- STRINGdb::STRINGdb$new(
        version         = version,
        species         = species,
        score_threshold = score_threshold,
        network_type    = network_type,
        input_directory = input_directory
    )

    # --- 3. Map gene symbols to STRING IDs -----------------------------------
    gene_df  <- data.frame(gene = gene_vec, stringsAsFactors = FALSE)
    mapped   <- sdb$map(gene_df, "gene", removeUnmappedRows = TRUE)
    n_mapped <- nrow(mapped)
    n_input  <- length(gene_vec)

    if (n_mapped == 0L)
        stop("None of the supplied genes could be mapped to STRING IDs. ",
             "Check gene symbols and species.", call. = FALSE)

    if (n_mapped < n_input) {
        not_found <- setdiff(gene_vec, mapped$gene)
        n_extra <- length(not_found) - 10L
        extra   <- if (n_extra > 0L) paste0("... (+", n_extra, " more)") else ""
        message(n_input - n_mapped,
                " gene(s) not found in STRING and removed: ",
                paste(head(not_found, 10), collapse = ", "), extra)
    }
    message(n_mapped, " genes successfully mapped.")

    # --- 4. Build igraph subnetwork ------------------------------------------
    g <- sdb$get_subnetwork(mapped$STRING_id)

    if (igraph::vcount(g) == 0L)
        stop("The STRING subnetwork is empty (no edges above score_threshold ",
             score_threshold, "). Try lowering score_threshold.", call. = FALSE)

    # Add gene symbols as vertex attribute
    # STRINGdb vertex names are STRING IDs; recover gene symbols via mapped df
    id_to_gene <- stats::setNames(mapped$gene, mapped$STRING_id)
    igraph::V(g)$gene_symbol <- id_to_gene[igraph::V(g)$name]

    # Nodes with no symbol mapping (edge-only neighbours): use STRING ID
    missing_sym <- is.na(igraph::V(g)$gene_symbol)
    igraph::V(g)$gene_symbol[missing_sym] <- igraph::V(g)$name[missing_sym]

    # --- 5. Compute layout ---------------------------------------------------
    # Shared with leviFromEdges(): the vertices here are named by STRING id,
    # so the gene symbols are passed in to label them.
    nodes <- .layoutNodes(g, layout, names = igraph::V(g)$gene_symbol)

    # --- 6. Build edges using gene symbols -----------------------------------
    el   <- igraph::as_edgelist(g)
    edges <- data.frame(
        V1 = id_to_gene[el[, 1]],
        V2 = id_to_gene[el[, 2]],
        stringsAsFactors = FALSE
    )
    # Replace any remaining unmapped IDs with the ID itself
    edges$V1[is.na(edges$V1)] <- el[is.na(edges$V1), 1]
    edges$V2[is.na(edges$V2)] <- el[is.na(edges$V2), 2]

    message("Network built: ", igraph::vcount(g), " nodes, ",
            igraph::ecount(g), " edges.")

    invisible(list(nodes = nodes, edges = edges, graph = g))
}
