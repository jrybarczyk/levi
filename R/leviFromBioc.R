#' @title leviFromExpressionSet
#' @description Convert an \code{ExpressionSet} (Biobase) object into a
#' data.frame compatible with \code{levi()}'s \code{expressionInput} parameter.
#' Typical use case: Bioconductor microarray datasets (Affymetrix, Agilent, etc.).
#' @param eset An \code{ExpressionSet} object.
#' @param condition_col Character. Column in \code{pData(eset)} that contains
#' the condition labels. Used together with \code{test_level} and
#' \code{ctrl_level}.
#' @param test_level Character. Value in \code{condition_col} identifying the
#' test samples.
#' @param ctrl_level Character. Value in \code{condition_col} identifying the
#' control samples.
#' @param test_col Character vector. Sample names to use as test group
#' (alternative to \code{condition_col} + \code{test_level}).
#' @param ctrl_col Character vector. Sample names to use as control group.
#' @param gene_col Character. Name for the gene identifier column in the
#' output. Default \code{"GeneID"}.
#' @return A \code{data.frame} with columns \code{gene_col}, \code{"Test"},
#' and \code{"Control"} containing per-gene mean intensity for each group.
#' Pass directly to \code{levi(expressionInput = ...)}.
#' @details
#' Requires the \code{Biobase} package (Bioconductor). Install with
#' \code{BiocManager::install("Biobase")}.
#'
#' Two usage patterns:
#' \enumerate{
#'   \item \strong{Condition labels}: specify \code{condition_col},
#'     \code{test_level}, and \code{ctrl_level}.
#'   \item \strong{Sample names}: provide \code{test_col} and \code{ctrl_col}
#'     directly.
#' }
#' @examples
#' if (requireNamespace("Biobase", quietly = TRUE)) {
#'     mat <- matrix(c(5.1, 6.2, 4.8,   5.4, 6.0, 5.1,   7.3, 5.9, 4.2),
#'                   nrow = 3,
#'                   dimnames = list(
#'                       c("GENE1", "GENE2", "GENE3"),
#'                       c("tumor1", "tumor2", "normal1")))
#'     pd <- Biobase::AnnotatedDataFrame(data.frame(
#'         condition = c("Tumor", "Tumor", "Normal"),
#'         row.names = colnames(mat)))
#'     eset <- Biobase::ExpressionSet(assayData = mat, phenoData = pd)
#'     expr_df <- leviFromExpressionSet(eset,
#'         condition_col = "condition",
#'         test_level    = "Tumor",
#'         ctrl_level    = "Normal")
#'     head(expr_df)
#' }
#' @export
leviFromExpressionSet <- function(eset,
                                   condition_col = NULL,
                                   test_level    = NULL,
                                   ctrl_level    = NULL,
                                   test_col      = NULL,
                                   ctrl_col      = NULL,
                                   gene_col      = "GeneID") {
    if (!requireNamespace("Biobase", quietly = TRUE)) {
        stop("Biobase is required. ",
             "Install with: BiocManager::install('Biobase')")
    }

    mat <- Biobase::exprs(eset)
    pd  <- Biobase::pData(eset)

    if (!is.null(condition_col)) {
        if (!condition_col %in% colnames(pd))
            stop("condition_col '", condition_col, "' not found in pData.")
        grp      <- as.character(pd[[condition_col]])
        test_idx <- which(grp == test_level)
        ctrl_idx <- which(grp == ctrl_level)
        if (length(test_idx) == 0)
            stop("No samples with test_level '", test_level, "' found.")
        if (length(ctrl_idx) == 0)
            stop("No samples with ctrl_level '", ctrl_level, "' found.")
    } else if (!is.null(test_col) && !is.null(ctrl_col)) {
        test_idx <- match(test_col, colnames(mat))
        ctrl_idx <- match(ctrl_col, colnames(mat))
        if (anyNA(test_idx))
            stop("Some test sample names not found in ExpressionSet.")
        if (anyNA(ctrl_idx))
            stop("Some control sample names not found in ExpressionSet.")
    } else {
        stop("Provide (condition_col + test_level + ctrl_level) ",
             "or (test_col + ctrl_col).")
    }

    test_mean <- rowMeans(mat[, test_idx, drop = FALSE], na.rm = TRUE)
    ctrl_mean <- rowMeans(mat[, ctrl_idx, drop = FALSE], na.rm = TRUE)

    out <- data.frame(
        Gene    = rownames(mat),
        Test    = test_mean,
        Control = ctrl_mean,
        stringsAsFactors = FALSE
    )
    colnames(out)[1] <- gene_col
    out <- stats::na.omit(out)
    rownames(out) <- NULL
    out
}


#' @title leviFromBioc
#' @description Universal adapter for Bioconductor expression containers.
#' Automatically detects the class of \code{obj} and calls the appropriate
#' \code{leviFrom*()} function.
#' @param obj A Bioconductor expression object. Supported classes:
#' \itemize{
#'   \item \code{SummarizedExperiment} (and subclasses such as
#'     \code{SingleCellExperiment} or \code{RangedSummarizedExperiment})
#'   \item \code{ExpressionSet} (Biobase)
#' }
#' @param ... Additional arguments forwarded to \code{\link{leviFromSE}} or
#' \code{\link{leviFromExpressionSet}}.
#' @return A \code{data.frame} with columns for gene identifiers, test, and
#' control expression values, ready to pass to \code{levi()}.
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     counts <- matrix(
#'         c(100, 50, 200,   120, 60, 210,   80, 150, 60),
#'         nrow = 3,
#'         dimnames = list(
#'             c("GENE1", "GENE2", "GENE3"),
#'             c("tumor1", "tumor2", "normal1"))
#'     )
#'     col_data <- data.frame(
#'         condition = c("Tumor", "Tumor", "Normal"),
#'         row.names = colnames(counts)
#'     )
#'     se <- SummarizedExperiment::SummarizedExperiment(
#'         assays  = list(counts = counts),
#'         colData = col_data
#'     )
#'     expr_df <- leviFromBioc(se,
#'         condition_col = "condition",
#'         test_level    = "Tumor",
#'         ctrl_level    = "Normal")
#'     head(expr_df)
#' }
#' @seealso \code{\link{leviFromSE}}, \code{\link{leviFromExpressionSet}}
#' @export
leviFromBioc <- function(obj, ...) {
    if (methods::is(obj, "SummarizedExperiment")) {
        leviFromSE(obj, ...)
    } else if (methods::is(obj, "ExpressionSet")) {
        leviFromExpressionSet(obj, ...)
    } else {
        stop("Unsupported class: ", paste(class(obj), collapse = ", "),
             ". Supported: SummarizedExperiment (and subclasses), ExpressionSet.")
    }
}
