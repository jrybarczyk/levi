#' @title leviFromDESeq2
#' @description Convert a DESeq2 result object into a data.frame compatible
#' with \code{levi}'s \code{expressionInput} parameter.
#' @param dds_result A \code{DESeqResults} object returned by
#' \code{DESeq2::results()}, or any \code{data.frame} with columns
#' \code{baseMean} and \code{log2FoldChange}.
#' @param gene_col Character. Name for the gene identifier column.
#' Default is \code{"GeneID"}.
#' @return A \code{data.frame} with columns: gene identifier, baseMean
#' (abundance annotation) and log2FoldChange, ready to pass
#' to \code{levi(expressionInput = ...)}.
#' @details
#' Use \code{signal_mode = "logfc"} and
#' \code{readExpColumn("log2FoldChange-log2FoldChange")}.
#' The \code{baseMean} column is annotation, not a control measurement.
#' Rows with \code{NA} in \code{log2FoldChange} are removed automatically.
#' @examples
#' mock_res <- data.frame(
#'     baseMean       = c(120, 45, 80),
#'     log2FoldChange = c(2.1, -1.4, 0.2),
#'     row.names      = c("GENE1", "GENE2", "GENE3")
#' )
#' expr_df <- leviFromDESeq2(mock_res)
#' head(expr_df)
#' \dontrun{
#' library(DESeq2)
#' res <- results(dds)
#' expr_df <- leviFromDESeq2(res, gene_col = "GeneSymbol")
#' }
#' @export
leviFromDESeq2 <- function(dds_result, gene_col = "GeneID") {
    result_df <- as.data.frame(dds_result)
    if (!all(c("baseMean", "log2FoldChange") %in% colnames(result_df))) {
        stop("dds_result must have 'baseMean' and 'log2FoldChange' columns.")
    }
    result_df[[gene_col]] <- rownames(result_df)
    out <- result_df[, c(gene_col, "baseMean", "log2FoldChange")]
    colnames(out) <- c(gene_col, "baseMean", "log2FoldChange")
    out <- na.omit(out)
    rownames(out) <- NULL
    return(out)
}

#' @title leviFromEdgeR
#' @description Convert an edgeR result object into a data.frame compatible
#' with \code{levi}'s \code{expressionInput} parameter.
#' @param fit A \code{DGELRT} or \code{DGEExact} object from edgeR, or any
#' \code{data.frame} with columns \code{logFC} and \code{logCPM}.
#' @param coef Retained for call compatibility. Omit for edgeR test results;
#' select the coefficient or contrast in \code{glmLRT()} or
#' \code{glmQLFTest()} before calling this adapter.
#' @param gene_col Character. Name for the gene identifier column.
#' Default is \code{"GeneID"}.
#' @return A \code{data.frame} with columns: gene identifier, logCPM
#' (abundance annotation) and logFC. Use \code{signal_mode = "logfc"}
#' with \code{readExpColumn("logFC-logFC")}; abundance is not a control.
#' @examples
#' mock_tt <- data.frame(
#'     logFC  = c(2.1, -1.4, 0.5),
#'     logCPM = c(5.2, 3.8, 4.1),
#'     row.names = c("GENE1", "GENE2", "GENE3")
#' )
#' expr_df <- leviFromEdgeR(mock_tt)
#' head(expr_df)
#' \dontrun{
#' library(edgeR)
#' fit <- glmQLFit(dge, design)
#' qlf <- glmQLFTest(fit, coef = 2)
#' expr_df <- leviFromEdgeR(qlf, gene_col = "GeneSymbol")
#' }
#' @export
leviFromEdgeR <- function(fit, coef = 1, gene_col = "GeneID") {
    if (is.data.frame(fit) && all(c("logFC", "logCPM") %in% colnames(fit))) {
        tt <- fit
    } else {
        if (!requireNamespace("edgeR", quietly = TRUE)) {
            stop("edgeR is required. ",
             "Install with: BiocManager::install('edgeR')")
        }
        if (!inherits(fit, "DGEExact") && !inherits(fit, "DGELRT"))
            stop("fit must be a DGEExact or DGELRT test result; run ",
                 "exactTest(), glmLRT() or glmQLFTest() first.")
        if (!missing(coef))
            stop("Select coef or contrast in glmLRT()/glmQLFTest() before ",
                 "calling leviFromEdgeR(); coef cannot change a test result.")
        tt <- edgeR::topTags(fit, n = Inf)$table
        if (!all(c("logFC", "logCPM") %in% names(tt)))
            stop("The edgeR result must contain one logFC column; ",
                 "select a single contrast before calling leviFromEdgeR().")
    }
    tt[[gene_col]] <- rownames(tt)
    out <- tt[, c(gene_col, "logCPM", "logFC")]
    colnames(out) <- c(gene_col, "logCPM", "logFC")
    out <- na.omit(out)
    rownames(out) <- NULL
    return(out)
}

#' @title leviFromLimma
#' @description Convert a limma result object into a data.frame compatible
#' with \code{levi}'s \code{expressionInput} parameter.
#' @param fit A \code{MArrayLM} object after \code{limma::eBayes()}, or any
#' \code{data.frame} with columns \code{logFC} and \code{AveExpr}.
#' @param coef Integer or character. Coefficient to extract (used only when
#' \code{fit} is a limma object). Default is \code{1}.
#' @param gene_col Character. Name for the gene identifier column.
#' Default is \code{"GeneID"}.
#' @return A \code{data.frame} with columns: gene identifier, AveExpr
#' (abundance annotation) and logFC. Use \code{signal_mode = "logfc"}
#' with \code{readExpColumn("logFC-logFC")}; abundance is not a control.
#' @examples
#' mock_top <- data.frame(
#'     logFC   = c(2.1, -1.4, 0.5),
#'     AveExpr = c(5.2, 3.8, 4.1),
#'     row.names = c("GENE1", "GENE2", "GENE3")
#' )
#' expr_df <- leviFromLimma(mock_top)
#' head(expr_df)
#' \dontrun{
#' library(limma)
#' fit2 <- eBayes(fit)
#' expr_df <- leviFromLimma(fit2, coef = 1, gene_col = "GeneSymbol")
#' }
#' @export
leviFromLimma <- function(fit, coef = 1, gene_col = "GeneID") {
    if (is.data.frame(fit) && all(c("logFC", "AveExpr") %in% colnames(fit))) {
        tt <- fit
    } else {
        if (!requireNamespace("limma", quietly = TRUE)) {
            stop("limma is required. ",
             "Install with: BiocManager::install('limma')")
        }
        tt <- limma::topTable(fit, coef = coef, n = Inf)
    }
    tt[[gene_col]] <- rownames(tt)
    out <- tt[, c(gene_col, "AveExpr", "logFC")]
    colnames(out) <- c(gene_col, "AveExpr", "logFC")
    out <- na.omit(out)
    rownames(out) <- NULL
    return(out)
}

#' @title leviFromSeurat
#' @description Convert a Seurat \code{FindMarkers} result into a data.frame
#' compatible with \code{levi}'s \code{expressionInput} parameter. Enables
#' single-cell RNA-seq data visualization on biological networks.
#' @param markers A \code{data.frame} returned by
#' \code{Seurat::FindMarkers()}, with columns \code{avg_log2FC} and
#' \code{pct.2}.
#' @param gene_col Character. Name for the gene identifier column.
#' Default is \code{"GeneID"}.
#' @return A \code{data.frame} with columns: gene identifier, pct.2
#' renamed to Control (annotation) and avg_log2FC renamed to Test.
#' @details
#' Use \code{signal_mode = "logfc"} and
#' \code{readExpColumn("Test-Test")}. The output Test column contains
#' avg_log2FC; Control contains pct.2 as annotation, not a comparable control.
#' @examples
#' mock_markers <- data.frame(
#'     avg_log2FC = c(1.5, -0.8, 0.3),
#'     pct.2      = c(0.1, 0.4, 0.6),
#'     row.names  = c("GENE1", "GENE2", "GENE3")
#' )
#' expr_df <- leviFromSeurat(mock_markers)
#' head(expr_df)
#' \dontrun{
#' library(Seurat)
#' markers <- FindMarkers(seurat_obj, ident.1 = "Tumor", ident.2 = "Normal")
#' expr_df <- leviFromSeurat(markers, gene_col = "GeneSymbol")
#' levi(expressionInput          = expr_df,
#'      networkCoordinatesInput  = my_network,
#'      fileTypeInput            = "dat",
#'      geneSymbolInput          = "GeneSymbol",
#'      readExpColumn            = readExpColumn("Test-Test"),
#'      signal_mode              = "logfc")
#' }
#' @export
leviFromSeurat <- function(markers, gene_col = "GeneID") {
    required_cols <- c("avg_log2FC", "pct.2")
    missing_cols <- setdiff(required_cols, colnames(markers))
    if (length(missing_cols) > 0) {
        stop("markers data.frame is missing columns: ",
             paste(missing_cols, collapse = ", "),
             ". Ensure it was produced by Seurat::FindMarkers().")
    }
    out <- markers[, c("pct.2", "avg_log2FC")]
    out[[gene_col]] <- rownames(markers)
    out <- out[, c(gene_col, "pct.2", "avg_log2FC")]
    colnames(out) <- c(gene_col, "Control", "Test")
    out <- na.omit(out)
    rownames(out) <- NULL
    return(out)
}
