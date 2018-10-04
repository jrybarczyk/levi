if(getRversion() >= "3.4.0") utils::globalVariables(c("Var1", "Var2", "z",
    "SigCoordPiso", "matrix_entrada", "matrix_saida"))
#' @import ggplot2 reshape2 colorspace grid dplyr Rcpp
#' @importFrom stats complete.cases
#' @importFrom igraph graph.edgelist
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph as_long_data_frame
#' @importFrom grDevices bmp cm.colors dev.off heat.colors jpeg png rainbow
#' @importFrom grDevices terrain.colors tiff topo.colors
#' @importFrom stats aggregate setNames
#' @importFrom utils count.fields read.delim read.table unzip write.table
#' @importFrom testthat expect_message test_check
#' @title readExpColumn
#' @usage readExpColumn(x,...)
#' @description This function helps to prepare the data in the script mode. It
#' also allows the obtention of dataset plot as a batch processing.
#' @param x Names of two expression datasets to be compared. They should be
#' separated by hyphen (-)
#' @param ... To add more comparisons, each combination must be separated by
#' comma (,).
#' @note To generate a plot from a single dataset, the name of the sample
#' must be informed twice (Ex. "CaseA-CaseA")
#' @details List the names of the expression datasets that will be used
#' for comparison
#' @return Returns the names of comparisons to be used by Levi
#' @author José Rafael Pilan (rafael.pilan@unesp.br)
#' @examples
#' base <- readExpColumn(a="NormalNeverSmoker-NormalNeverSmoker")
#' @export
readExpColumn <- function(x,...) {
    c <- as.list(match.call())
    return(c)
}
