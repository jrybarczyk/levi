if(getRversion() >= "3.4.0") utils::globalVariables(c("Var1", "Var2", "z",
    "SigCoordPiso", "matrix_entrada", "matrix_saida",
    "x", "y", "label", "Diff", "pval"))
#' @title readExpColumn - Define expression comparisons for levi
#'
#' @description Specify which columns of the expression file represent the
#' \emph{test} condition and the \emph{control} condition for the landscape
#' computation. Multiple comparisons can be provided in a single call to
#' enable batch-mode processing (one landscape per comparison).
#'
#' @usage readExpColumn(x, ...)
#'
#' @param x Character string of the form \code{"TestColumn-ControlColumn"},
#'   where \code{TestColumn} and \code{ControlColumn} are column names in the
#'   expression file. The hyphen (\code{-}) separates test from control.
#'   \strong{Special cases:}
#'   \itemize{
#'     \item Single-column logFC input (e.g., DESeq2 \code{log2FoldChange}
#'           only): repeat the same column name on both sides -
#'           \code{"log2FoldChange-log2FoldChange"}.
#'     \item Self-comparison (null landscape): use the same column twice -
#'           \code{"Sample-Sample"}.
#'   }
#' @param ... Additional comparisons, each following the same
#'   \code{"Test-Control"} format. When multiple comparisons are supplied,
#'   \code{levi()} generates one landscape per comparison and returns a list.
#'
#' @return A named list used internally by \code{\link{levi}} to identify
#'   the expression columns for each comparison. The first element is always
#'   the function call itself (used as a sentinel); subsequent elements are
#'   the comparison strings.
#'
#' @details
#' Column names must match exactly the column headers in the expression file
#' or data.frame (case-sensitive). The hyphen separator is required; column
#' names that themselves contain hyphens should be renamed before use.
#'
#' When using \code{signal_mode = "logfc"} with a single-column fold-change
#' input (e.g., from \code{\link{leviFromDESeq2}}), repeat the logFC column
#' name: \code{readExpColumn("log2FoldChange-log2FoldChange")}. The function
#' detects \code{baseTest == baseControl} and routes to single-column mode
#' automatically.
#'
#' @seealso \code{\link{levi}}, \code{\link{leviFromDESeq2}},
#'   \code{\link{leviFromEdgeR}}, \code{\link{leviFromLimma}},
#'   \code{\link{leviFromSeurat}}
#'
#' @author Jose Rafael Pilan \email{rafael.pilan@@unesp.br}
#'
#' @examples
#' # Single comparison
#' readExpColumn("TumorCurrentSmoker-NormalNeverSmoker")
#'
#' # Multiple comparisons (batch mode) - levi() returns a list
#' readExpColumn(
#'     "TumorCurrentSmoker-NormalNeverSmoker",
#'     "TumorFormerSmoker-NormalFormerSmoker"
#' )
#'
#' # Single-column logFC (repeat the column name)
#' readExpColumn("log2FoldChange-log2FoldChange")
#'
#' # Self-comparison (null landscape, all scores ~ 0.5)
#' readExpColumn("Sample1-Sample1")
#'
#' @export
readExpColumn <- function(x, ...) {
    c(list(quote(readExpColumn)), list(x, ...))
}
