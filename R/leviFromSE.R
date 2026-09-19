#' @title leviFromSE
#' @description Convert a \code{SummarizedExperiment} object into a
#' data.frame compatible with \code{levi}'s \code{expressionInput} parameter.
#' Handles both single-assay and multi-assay SE objects.
#' @param se A \code{SummarizedExperiment} (or subclass) object.
#' @param assay_name Character. Name or index of the assay to extract.
#' Default is \code{"counts"} (falls back to assay index 1 if not found).
#' @param test_col Character or integer. Column index/name in
#' \code{colData(se)} identifying which samples are the \strong{test}
#' condition. Alternatively, a character vector of sample names.
#' @param ctrl_col Character or integer. Column index/name in
#' \code{colData(se)} identifying which samples are the \strong{control}
#' condition. Same format as \code{test_col}.
#' @param condition_col Character. Column in \code{colData(se)} that contains
#' condition labels. Used together with \code{test_level} / \code{ctrl_level}
#' to select samples. Ignored when \code{test_col} / \code{ctrl_col} are
#' provided directly.
#' @param test_level Character. Value in \code{condition_col} identifying test
#' samples.
#' @param ctrl_level Character. Value in \code{condition_col} identifying
#' control samples.
#' @param gene_col Character. Name for the gene identifier column in the
#' output. Default is \code{"GeneID"}.
#' @param log_transform Logical. If \code{TRUE}, apply \code{log2(x + 1)}
#' to the extracted counts before aggregating. Useful for count data.
#' Default is \code{FALSE}.
#' @return A \code{data.frame} with columns \code{gene_col}, \code{"Test"},
#' and \code{"Control"} containing per-gene mean expression for each group.
#' Pass directly to \code{levi(expressionInput = ...)}.
#' @details
#' Two usage patterns:
#' \enumerate{
#'   \item \strong{Direct sample selection}: provide sample names or indices
#'     via \code{test_col} / \code{ctrl_col}.
#'   \item \strong{Condition-label selection}: provide a \code{colData}
#'     column name via \code{condition_col} and the levels via
#'     \code{test_level} / \code{ctrl_level}.
#' }
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     counts <- matrix(
#'         c(100, 50, 200,   120, 60, 210,   80, 150, 60),
#'         nrow = 3,
#'         dimnames = list(
#'             c("GENE1", "GENE2", "GENE3"),
#'             c("tumor1", "tumor2", "normal1")
#'         )
#'     )
#'     col_data <- data.frame(
#'         condition = c("Tumor", "Tumor", "Normal"),
#'         row.names = colnames(counts)
#'     )
#'     se <- SummarizedExperiment::SummarizedExperiment(
#'         assays  = list(counts = counts),
#'         colData = col_data
#'     )
#'     expr_df <- leviFromSE(se,
#'         condition_col = "condition",
#'         test_level    = "Tumor",
#'         ctrl_level    = "Normal")
#'     head(expr_df)
#' }
#' \dontrun{
#' # Pattern 2: direct sample names, on your own SummarizedExperiment
#' expr_df <- leviFromSE(se,
#'     test_col = c("tumor_1", "tumor_2"),
#'     ctrl_col = c("normal_1", "normal_2"))
#' }
#' @export
leviFromSE <- function(se,
                        assay_name    = "counts",
                        test_col      = NULL,
                        ctrl_col      = NULL,
                        condition_col = NULL,
                        test_level    = NULL,
                        ctrl_level    = NULL,
                        gene_col      = "GeneID",
                        log_transform = FALSE) {

    # Extract assay --
    avail <- SummarizedExperiment::assayNames(se)
    if (!is.null(assay_name) && assay_name %in% avail) {
        mat <- SummarizedExperiment::assay(se, assay_name)
    } else {
        if (!is.null(assay_name) && !(assay_name %in% avail)) {
            message("Assay '", assay_name, "' not found. Using assay index 1.")
        }
        mat <- SummarizedExperiment::assay(se, 1)
    }

    if (log_transform) mat <- log2(mat + 1)

    # Identify test / control samples --
    if (!is.null(condition_col)) {
        cd <- as.data.frame(SummarizedExperiment::colData(se))
        if (!condition_col %in% colnames(cd))
            stop("condition_col '", condition_col, "' not found in colData.")
        grp <- cd[[condition_col]]

        if (is.null(test_level) || is.null(ctrl_level))
            stop("Provide test_level and ctrl_level when using condition_col.")

        test_idx <- which(grp == test_level)
        ctrl_idx <- which(grp == ctrl_level)
        if (length(test_idx) == 0)
            stop("No samples with test_level '", test_level, "' found.")
        if (length(ctrl_idx) == 0)
            stop("No samples with ctrl_level '", ctrl_level, "' found.")

    } else if (!is.null(test_col) && !is.null(ctrl_col)) {
        # accept names or integer indices
        if (is.character(test_col)) {
            test_idx <- match(test_col, colnames(mat))
            if (anyNA(test_idx))
                stop("Some test sample names not found: ",
                     paste(test_col[is.na(test_idx)], collapse = ", "))
        } else {
            test_idx <- test_col
        }
        if (is.character(ctrl_col)) {
            ctrl_idx <- match(ctrl_col, colnames(mat))
            if (anyNA(ctrl_idx))
                stop("Some control sample names not found: ",
                     paste(ctrl_col[is.na(ctrl_idx)], collapse = ", "))
        } else {
            ctrl_idx <- ctrl_col
        }
    } else {
        stop("Provide either (condition_col + test_level + ctrl_level) ",
             "or (test_col + ctrl_col).")
    }

    # Compute per-gene means --
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
