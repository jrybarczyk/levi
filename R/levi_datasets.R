# =============================================================================
# Documentation for bundled toy datasets (inst/extdata/)
#
# Each dataset pair ships as two tab-delimited .dat files:
#   *_network.dat   - network topology + node coordinates (Medusa DAT format)
#   *_expression.dat - expression values (one or more comparisons)
#
# Use file.path(system.file(package="levi"), "extdata", "<filename>") to
# access any of these files from R.
# =============================================================================


#' Hub-topology toy dataset
#'
#' A minimal 9-node star network with a known, asymmetric expression pattern.
#' The central hub and its four closest spokes are strongly over-expressed
#' (Test = 200, Control = 10); the four outer corner nodes are strongly
#' under-expressed (Test = 5, Control = 200). The expected landscape shows a
#' clear peak at the centre and valleys at the periphery.
#'
#' @format Two tab-delimited files in \code{inst/extdata/}:
#' \describe{
#'   \item{\code{hub_network.dat}}{Medusa DAT network file with 9 nodes
#'     (\code{HUB}, \code{N1}-\code{N8}) and 8 edges (all connecting to HUB).
#'     Node coordinates place the hub at (0.50, 0.50) and spokes at increasing
#'     radial distances.}
#'   \item{\code{hub_expression.dat}}{Tab-delimited expression file with
#'     columns \code{ID}, \code{Test}, \code{Control}. Inner nodes
#'     (HUB, N1-N4): Test = 200, Control = 10. Outer nodes (N5-N8):
#'     Test = 5, Control = 200.}
#' }
#'
#' @section Expected results:
#' With default parameters (\code{signal_mode = "ratio"}):
#' \itemize{
#'   \item \code{score(HUB)} > all spoke scores > 0.5
#'   \item \code{score(N5)} ~ \code{score(N6)} ~ \code{score(N7)} ~
#'         \code{score(N8)} < 0.5
#'   \item Peaks detected over the over-expressed core (HUB and N1-N4 all
#' reach the same maximum, so several are reported); valleys near
#' N5-N8.
#' }
#'
#' @source Simulated data generated for unit testing. No biological meaning.
#'
#' @seealso \code{\link{levi}},
#' \code{\link[=gradient_dataset]{gradient dataset}}
#'
#' @return NULL (this object documents data files, not a function)
#'
#' @examples
#' \donttest{
#' hub_net  <- file.path(system.file(package = "levi"),
#'                        "extdata", "hub_network.dat")
#' hub_expr <- file.path(system.file(package = "levi"),
#'                        "extdata", "hub_expression.dat")
#'
#' res <- levi(
#'     networkCoordinatesInput = hub_net,
#'     expressionInput         = hub_expr,
#'     fileTypeInput           = "dat",
#'     geneSymbolInput         = "ID",
#'     readExpColumn           = readExpColumn("Test-Control"),
#'     contrastValueInput      = 50,
#'     resolutionValueInput    = 50,
#'     zoomValueInput          = 50,
#'     smoothValueInput        = 50
#' )
#' head(res$scores)
#' }
#'
#' @name hub_dataset
#' @aliases hub_network hub_expression
NULL


#' Gradient-topology toy dataset
#'
#' A 6-node linear chain where expression increases monotonically from left
#' (\code{GA}) to right (\code{GF}). Designed to verify that landscape scores
#' preserve strict rank order along a gradient and that no spurious peaks
#' appear in intermediate positions.
#'
#' @format Two tab-delimited files in \code{inst/extdata/}:
#' \describe{
#'   \item{\code{gradient_network.dat}}{Linear chain GA-GB-GC-GD-GE-GF with
#'     nodes equally spaced from x = 0.10 to x = 0.90 at y = 0.50.}
#'   \item{\code{gradient_expression.dat}}{Columns \code{ID}, \code{Test},
#'     \code{Control}. Values increase from GA (Test = 5, Control = 100) to
#'     GF (Test = 100, Control = 5), giving a smooth gradient of signal
#'     ratios from ~ 0.09 to ~ 0.91.}
#' }
#'
#' @section Expected results:
#' \itemize{
#'   \item Strictly monotone scores: GA < GB < GC < GD < GE < GF.
#'   \item GF ranked 1st; GA ranked last.
#'   \item \code{score(GF)} > 0.5; \code{score(GA)} < 0.5.
#'   \item No peak detected at intermediate nodes.
#' }
#'
#' @source Simulated data generated for unit testing.
#'
#' @seealso \code{\link{levi}}, \code{\link[=hub_dataset]{hub dataset}}
#'
#' @return NULL (this object documents data files, not a function)
#'
#' @examples
#' \donttest{
#' grad_net  <- file.path(system.file(package = "levi"),
#'                         "extdata", "gradient_network.dat")
#' grad_expr <- file.path(system.file(package = "levi"),
#'                         "extdata", "gradient_expression.dat")
#'
#' res <- levi(
#'     networkCoordinatesInput = grad_net,
#'     expressionInput         = grad_expr,
#'     fileTypeInput           = "dat",
#'     geneSymbolInput         = "ID",
#'     readExpColumn           = readExpColumn("Test-Control"),
#'     contrastValueInput      = 50,
#'     resolutionValueInput    = 50,
#'     zoomValueInput          = 50,
#'     smoothValueInput        = 50
#' )
#' res$scores[order(res$scores$Rank), c("Gene", "LandscapeScore", "Rank")]
#' }
#'
#' @name gradient_dataset
#' @aliases gradient_network gradient_expression
NULL


#' Bimodal toy dataset
#'
#' A 13-node network with two star clusters connected by a neutral bridge gene.
#' Cluster A is strongly over-expressed; cluster B is strongly under-expressed.
#' The resulting landscape has a peak on the left and a valley on the right -
#' ideal for validating simultaneous peak-and-valley detection and
#' \code{\link{leviDiff}}.
#'
#' @format Two tab-delimited files in \code{inst/extdata/}:
#' \describe{
#'   \item{\code{bimodal_network.dat}}{Two star sub-networks:
#'     \code{A_HUB} + A1-A5 (left, x ~ 0.1-0.35) and \code{B_HUB} + B1-B5
#'     (right, x ~ 0.65-0.90), linked via BRIDGE (x = 0.50).}
#'   \item{\code{bimodal_expression.dat}}{Columns \code{ID}, \code{Test},
#'     \code{Control}. Cluster A: Test = 200, Control = 10. BRIDGE:
#'     Test = Control = 100 (neutral). Cluster B: Test = 10, Control = 200.}
#' }
#'
#' @section Expected results:
#' \itemize{
#'   \item \code{score(A_HUB)} > \code{score(BRIDGE)} > \code{score(B_HUB)}.
#'   \item All cluster-A nodes score > 0.5; all cluster-B nodes < 0.5.
#'   \item Peak(s) detected near A cluster; valley(s) near B cluster.
#'   \item With \code{perm_side = "both"}: dashed contour on left,
#'         dotted contour on right.
#' }
#'
#' @source Simulated data generated for unit testing.
#'
#' @seealso \code{\link{levi}}, \code{\link{leviDiff}}
#'
#' @return NULL (this object documents data files, not a function)
#'
#' @examples
#' \donttest{
#' bim_net  <- file.path(system.file(package = "levi"),
#'                        "extdata", "bimodal_network.dat")
#' bim_expr <- file.path(system.file(package = "levi"),
#'                        "extdata", "bimodal_expression.dat")
#'
#' res <- levi(
#'     networkCoordinatesInput = bim_net,
#'     expressionInput         = bim_expr,
#'     fileTypeInput           = "dat",
#'     geneSymbolInput         = "ID",
#'     readExpColumn           = readExpColumn("Test-Control"),
#'     contrastValueInput      = 50,
#'     resolutionValueInput    = 50,
#'     zoomValueInput          = 50,
#'     smoothValueInput        = 50,
#'     contourLevi             = TRUE
#' )
#' }
#'
#' @name bimodal_dataset
#' @aliases bimodal_network bimodal_expression
NULL


#' Flat (null-expression) toy dataset
#'
#' A 3x3 grid network where all genes have identical expression
#' (Test = Control = 100). The landscape is expected to be nearly flat
#' (all scores ~ 0.5) with no meaningful peaks or valleys.
#' Used to verify that \code{levi} does not produce false positives when
#' there is no differential expression.
#'
#' @format Two tab-delimited files in \code{inst/extdata/}:
#' \describe{
#'   \item{\code{flat_network.dat}}{Nine nodes F11-F33 arranged as a 3x3 grid
#'     with horizontal and vertical edges (12 edges total). Node positions
#'     span (0.20, 0.20) to (0.80, 0.80).}
#'   \item{\code{flat_expression.dat}}{All 9 genes with
#'     Test = Control = 100, yielding a uniform signal ratio of 0.5.}
#' }
#'
#' @section Expected results:
#' \itemize{
#'   \item All landscape scores within [0.2, 0.8] (no extreme regions).
#'   \item Standard deviation of scores < 0.15.
#'   \item No peaks with score > 0.8; no valleys with score < 0.2.
#' }
#'
#' @source Simulated data generated for unit testing.
#'
#' @seealso \code{\link{levi}}
#'
#' @return NULL (this object documents data files, not a function)
#'
#' @examples
#' \donttest{
#' flat_net  <- file.path(system.file(package = "levi"),
#'                         "extdata", "flat_network.dat")
#' flat_expr <- file.path(system.file(package = "levi"),
#'                         "extdata", "flat_expression.dat")
#'
#' res <- levi(
#'     networkCoordinatesInput = flat_net,
#'     expressionInput         = flat_expr,
#'     fileTypeInput           = "dat",
#'     geneSymbolInput         = "ID",
#'     readExpColumn           = readExpColumn("Test-Control"),
#'     contrastValueInput      = 50,
#'     resolutionValueInput    = 50,
#'     zoomValueInput          = 50,
#'     smoothValueInput        = 50
#' )
#' round(res$scores$LandscapeScore, 3)
#' }
#'
#' @name flat_dataset
#' @aliases flat_network flat_expression
NULL


#' Log2-intensity (logfc-mode) toy dataset
#'
#' A 6-node linear chain identical in topology to the
#' \code{\link[=gradient_dataset]{gradient dataset}}, but with expression
#' values representing log2 microarray intensities (range 5.0-8.0).
#' The logFC between test and control spans -3 to +3 log2 units.
#' Designed to show that \code{signal_mode = "logfc"} extracts wider contrast
#' than \code{signal_mode = "ratio"} on log-scale data.
#'
#' @format Two tab-delimited files in \code{inst/extdata/}:
#' \describe{
#'   \item{\code{logfc_network.dat}}{Linear chain GA-GF, identical layout to
#'     \code{gradient_network.dat}.}
#'   \item{\code{logfc_expression.dat}}{Columns \code{ID}, \code{Test},
#'     \code{Control} in log2 intensity units:
#'     \tabular{lll}{
#'       Gene \tab Test \tab Control \cr
#'       GA   \tab 5.0  \tab 8.0    \cr
#'       GB   \tab 6.0  \tab 7.5    \cr
#'       GC   \tab 7.0  \tab 7.5    \cr
#'       GD   \tab 7.5  \tab 7.0    \cr
#'       GE   \tab 8.0  \tab 6.5    \cr
#'       GF   \tab 8.0  \tab 5.0    \cr
#'     }
#'     LogFC = Test - Control ranges from -3 (GA) to +3 (GF).
#'   }
#' }
#'
#' @section Signal mode comparison:
#' \describe{
#'   \item{\code{signal_mode = "logfc"}}{
#'     sigmoid(-3) ~ 0.047 (GA) -> sigmoid(+3) ~ 0.953 (GF).
#'     Score range ~ 0.91 - strong biological contrast.}
#'   \item{\code{signal_mode = "ratio"}}{
#'     5/(5+8) ~ 0.38 (GA) -> 8/(8+5) ~ 0.62 (GF).
#'     Score range ~ 0.24 - weaker contrast from treating
#'     log2 values as linear.}
#' }
#' Use \code{signal_mode = "logfc"} whenever Test and Control represent
#' log2-transformed intensities (RMA microarray, VST/rlog, log2-proteomics).
#'
#' @source Simulated log2 microarray intensities for unit testing.
#'
#' @seealso \code{\link{levi}},
#' \code{\link[=gradient_dataset]{gradient dataset}}
#'
#' @return NULL (this object documents data files, not a function)
#'
#' @examples
#' \donttest{
#' logfc_net  <- file.path(system.file(package = "levi"),
#'                          "extdata", "logfc_network.dat")
#' logfc_expr <- file.path(system.file(package = "levi"),
#'                          "extdata", "logfc_expression.dat")
#'
#' # logfc mode - recommended for log2-scale data
#' res_lfc <- levi(
#'     networkCoordinatesInput = logfc_net,
#'     expressionInput         = logfc_expr,
#'     fileTypeInput           = "dat",
#'     geneSymbolInput         = "ID",
#'     readExpColumn           = readExpColumn("Test-Control"),
#'     signal_mode             = "logfc",
#'     contrastValueInput      = 50,
#'     resolutionValueInput    = 50,
#'     zoomValueInput          = 50,
#'     smoothValueInput        = 50
#' )
#' diff(range(res_lfc$scores$LandscapeScore))  # score range with logfc mode
#' }
#'
#' @name logfc_dataset
#' @aliases logfc_network logfc_expression
NULL


#' Multi-comparison toy dataset
#'
#' A three-condition expression file for the
#' \code{\link[=hub_dataset]{hub network}}, enabling batch-mode processing and
#' testing of \code{\link{leviGrid}} and \code{\link{leviDiff}}.
#' Condition A is the maximum-contrast state (hub up, corners down);
#' Condition B is neutral; Condition C is the reverse of A.
#'
#' @format One tab-delimited file in \code{inst/extdata/} with columns
#' \code{ID}, \code{Cond_A}, \code{Cond_B}, \code{Cond_C}:
#' \tabular{llll}{
#'   Gene  \tab Cond_A \tab Cond_B \tab Cond_C \cr
#'   HUB   \tab 200    \tab 100    \tab 10     \cr
#'   N1-N4 \tab 200    \tab 100    \tab 10     \cr
#'   N5-N8 \tab 10     \tab 100    \tab 200    \cr
#' }
#' Designed to be used with \code{hub_network.dat} (\code{fileTypeInput = "dat"}).
#'
#' @section Use cases:
#' \describe{
#'   \item{Batch mode}{
#'     \code{readExpColumn("Cond_A-Cond_B", "Cond_A-Cond_C")} returns two
#'     landscapes in a single \code{levi()} call.}
#'   \item{leviGrid}{Compare the two landscapes side by side.}
#'   \item{leviDiff}{Subtract Cond_A-Cond_B from Cond_A-Cond_C to reveal
#'     the additional contrast gained by Condition C.}
#' }
#'
#' @source Simulated data generated for unit testing.
#'
#' @seealso \code{\link{levi}}, \code{\link{leviGrid}}, \code{\link{leviDiff}}
#'
#' @return NULL (this object documents data files, not a function)
#'
#' @examples
#' \donttest{
#' hub_net <- file.path(system.file(package = "levi"),
#'                       "extdata", "hub_network.dat")
#' mc_expr <- file.path(system.file(package = "levi"),
#'                       "extdata", "hub_multicomp_expression.dat")
#'
#' res_list <- levi(
#'     networkCoordinatesInput = hub_net,
#'     expressionInput         = mc_expr,
#'     fileTypeInput           = "dat",
#'     geneSymbolInput         = "ID",
#'     readExpColumn           = readExpColumn("Cond_A-Cond_B",
#'                                             "Cond_A-Cond_C")
#' )
#' leviGrid(res_list, titles = c("A vs B", "A vs C"))
#' leviDiff(res_list[[1]], res_list[[2]])
#' }
#'
#' @name multicomp_dataset
#' @aliases hub_multicomp_expression
NULL


#' Sparse-expression toy dataset
#'
#' A 15-node linear chain where only 5 of the 15 network genes have measured
#' expression values. The remaining 10 genes are absent from the expression
#' file and receive the mode-appropriate neutral value
#' (0.5 for \code{signal_mode = "ratio"}; 0.0 for \code{"logfc"} and
#' \code{"zscore"}). Used to verify that missing genes are handled gracefully
#' and that scores are returned for all 15 network nodes.
#'
#' @format Two tab-delimited files in \code{inst/extdata/}:
#' \describe{
#'   \item{\code{sparse_network.dat}}{Linear chain S01-S15 with 14 edges,
#'     nodes equally spaced from x = 0.05 to x = 0.95 at y = 0.50.}
#'   \item{\code{sparse_expression.dat}}{Only genes S01-S05 are present.
#'     Expression spans from strongly up-regulated (S01: Test = 200,
#'     Control = 10) to strongly down-regulated (S05: Test = 10,
#'     Control = 200). Genes S06-S15 are not listed and receive neutral
#'     imputation.}
#' }
#'
#' @section Expected results:
#' \itemize{
#'   \item \code{nrow(result$scores)} = 15 (all network nodes are scored).
#'   \item \code{score(S01)} > \code{score(S05)}.
#'   \item No error or warning is raised for the 10 missing genes.
#'   \item Both \code{signal_mode = "ratio"} and \code{"zscore"} complete
#'         successfully; missing-gene neutral values differ between modes
#'         but the outcome is valid for both.
#' }
#'
#' @source Simulated data generated for unit testing.
#'
#' @seealso \code{\link{levi}}
#'
#' @return NULL (this object documents data files, not a function)
#'
#' @examples
#' \donttest{
#' sparse_net  <- file.path(system.file(package = "levi"),
#'                           "extdata", "sparse_network.dat")
#' sparse_expr <- file.path(system.file(package = "levi"),
#'                           "extdata", "sparse_expression.dat")
#'
#' res <- levi(
#'     networkCoordinatesInput = sparse_net,
#'     expressionInput         = sparse_expr,
#'     fileTypeInput           = "dat",
#'     geneSymbolInput         = "ID",
#'     readExpColumn           = readExpColumn("Test-Control"),
#'     contrastValueInput      = 50,
#'     resolutionValueInput    = 50,
#'     zoomValueInput          = 50,
#'     smoothValueInput        = 50
#' )
#' nrow(res$scores)   # 15 - all network nodes scored
#' }
#'
#' @name sparse_dataset
#' @aliases sparse_network sparse_expression
NULL
