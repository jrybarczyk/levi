#' @title LEVIui - Interactive Shiny GUI for levi
#'
#' @description Launch the \pkg{levi} graphical user interface (GUI) locally
#' via Shiny. The GUI exposes all parameters available in script mode and
#' additionally provides:
#' \itemize{
#'   \item Interactive brushable region selection on the landscape heatmap
#'         (returns gene scores for the selected area).
#'   \item Gene search - type a gene symbol to highlight its position with a
#'         yellow circle on the landscape.
#'   \item Peak label overlay (\emph{requires} \pkg{ggrepel}).
#'   \item Signal transformation mode selector (\code{ratio}, \code{logfc},
#'         \code{zscore}) and \code{logfc_k} steepness control.
#'   \item Interactive 3D surface (\emph{requires} \pkg{plotly}), carrying
#'         the same significance boundary as the 2D map.
#'   \item Permutation significance test with per-iteration progress bar.
#'   \item Download buttons for the landscape plot, node score table, and
#'         peak/valley table (CSV).
#' }
#'
#' @param browser Logical. \code{TRUE} opens the app in the system web browser;
#'   \code{FALSE} (default) opens it in the RStudio Viewer pane.
#'
#' @return Runs the Shiny application; does not return an R value.
#'
#' @details
#' The GUI is a full-featured interface to \code{\link{levi}}. All file
#' uploads, parameter sliders, and result downloads are handled interactively.
#' For reproducible, automated, or batch analyses use \code{\link{levi}}
#' directly in script mode.
#'
#' The app is located in \code{inst/shiny/} and can also be launched with
#' \code{shiny::runApp(system.file("shiny", package = "levi"))}.
#'
#' @seealso \code{\link{levi}}, \code{\link{readExpColumn}}
#'
#' @author Jose Rafael Pilan \email{rafael.pilan@@unesp.br},
#'   Isabelle Mira da Silva
#'
#' @examples
#' if (interactive()) {
#'     LEVIui(browser = FALSE)   # opens in RStudio Viewer
#'     LEVIui(browser = TRUE)    # opens in system browser
#' }
#'
#' @export

LEVIui <- function(browser = FALSE) {
    if (!is.logical(browser) || length(browser) != 1L || is.na(browser))
        stop("'browser' must be a single TRUE or FALSE value.", call. = FALSE)

    appDir <- system.file("shiny", package = "levi")
    if (!nzchar(appDir))
        stop("The Shiny application was not found in the installed package. ",
            "Reinstall levi with BiocManager::install('levi').", call. = FALSE)

    # launch.browser is only set when TRUE: leaving it unset lets RStudio open
    # the app in its Viewer pane, which is what browser = FALSE promises.
    if (browser) {
        shiny::runApp(appDir, launch.browser = TRUE)
    } else {
        shiny::runApp(appDir)
    }
}
