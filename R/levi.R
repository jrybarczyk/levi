#' @title levi
#' @import RColorBrewer shinydashboard httr rmarkdown markdown
#' @import xml2 Rcpp knitr methods
#' @importFrom Rcpp evalCpp
#' @useDynLib levi
#' @importFrom shinyjs hide enable disable reset useShinyjs extendShinyjs
#'         js inlineCSS toggleState
#' @importFrom shiny actionButton actionLink addResourcePath column
#'         conditionalPanel downloadButton downloadHandler
#'         eventReactive fileInput fluidPage helpText isolate
#'         mainPanel need numericInput observe observeEvent
#'         outputOptions plotOutput radioButtons
#'         reactive reactiveValues renderPlot renderUI runApp
#'         selectInput shinyApp shinyServer shinyUI sidebarLayout
#'         sidebarPanel sliderInput stopApp tabPanel tabsetPanel
#'         textInput textOutput titlePanel uiOutput tags HTML
#'         h4 img icon updateTabsetPanel updateTextInput validate
#'         wellPanel checkboxInput br checkboxGroupInput a strong
#'         renderPrint fluidRow showNotification brushedPoints
#' @importFrom dplyr union as_data_frame groups select slice data_frame
#' @importFrom dplyr arrange %>% filter
#' @importFrom methods is
#' @importFrom DT datatable dataTableOutput renderDataTable
#' @importFrom igraph graph.edgelist as_long_data_frame
#' @usage levi(expressionInput, fileTypeInput, networkCoordinatesInput,
#' networkInteractionsInput, geneSymbolnput, readExpColumn,
#' contrastValueInput, zoomValueInput, resolutionValueInput,
#' smoothValueInput, expressionLog, contourLevi, setcolor)
#' @description This is the Levi script mode. It allows you to create the
#' integration of networks and gene expression levels as batch
#' processing
#' @param expressionInput Filename of gene expression data, which is a
#' numeric data.frame or matrix. The rows represent genes/proteins and the
#' columns represent the experiment (RNA-seq, microarray, etc).
#' @param fileTypeInput Filename of biological network. Levi can read files
#' written in the following formats: Medusa (DAT), RedeR (DYN), Pajek (NET)
#' and STRING/STITCH
#' @param networkCoordinatesInput It allows the user to load the coordinate
#' of the nodes the network.
#' @param networkInteractionsInput Parameter available only to
#' STRING/STITCH data format.
#' It allows the user to load the interaction data file of the network.
#' @param geneSymbolnput Column name from gene expression data containing the
#' identifier (gene Symbol, Entrez ID, EMSEMBL, etc).
#' @param readExpColumn Variable from readExpColumn function containing the
#' comparisons of the experiments
#' @param contrastValueInput Numeric value for image contrast. The variable
#' range is 0 to 100. The default value is 50
#' @param zoomValueInput  Numeric value for image zoom. The variable range is 0
#' to 100. The default value is 50.
#' @param resolutionValueInput Numeric value for image resolution. The variable
#' range is 0 to 100. The default value is 50.
#' @param smoothValueInput Numeric value for image smoothness. The variable
#' range is 0 to 100. The default is 50.
#' @param  expressionLog Logical variable to indicate Log2 normalization in the
#' expression levels. The default is FALSE
#' @param contourLevi Logical variable to allow contour lines. The default is
#' TRUE.
#' @param setcolor Select the color palette to build the heatmat. There is
#'two options the **Multicolor** has 20 color levels combined. The
#'**Two colors** has two types of color and the options available are: *purple_pink*,
#'*green_blue*, *blue_yellow*, *pink_green*, *orange_purple*, *green_marine*.
#' @details Integrates the biological network and gene expression levels
#' (or other type of data)
#' @author Isabelle Mira da Silva (isabelle.silva@unesp.br), José Rafael Pilan (rafael.pilan@unesp.br)
#' @return Return a ggplot object and print a image (heatmat).
#' @examples
#'template_network <- file.path(system.file(package="levi"),"extdata",
#'    "medusa.dat", fsep = .Platform$file.sep)
#'
#'template_expression <- file.path(system.file(package="levi"),
#'    "extdata","expression.dat", fsep = .Platform$file.sep)
#'
#'multicolor <- levi(networkCoordinatesInput = template_network,
#'    expressionInput = template_expression, fileTypeInput = "dat",
#'    geneSymbolnput = "ID",
#'    readExpColumn=readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
#'    contrastValueInput = 50, resolutionValueInput  = 50, zoomValueInput = 50,
#'    smoothValueInput = 50, expressionLog = FALSE, contourLevi = TRUE)
#'
#'twocolors <- levi(networkCoordinatesInput = template_network,
#'    expressionInput = template_expression, fileTypeInput = "dat",
#'    geneSymbolnput = "ID",
#'    readExpColumn = readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
#'    setcolor = "pink_green", contourLevi = FALSE)
#'@export
levi <- function(expressionInput, fileTypeInput, networkCoordinatesInput,
    networkInteractionsInput = NA, geneSymbolnput, readExpColumn,
    contrastValueInput = 50, zoomValueInput = 50, resolutionValueInput = 50,
    smoothValueInput = 50, expressionLog = FALSE,
    contourLevi  = FALSE, setcolor = "default"){
        levi_function(expressionInput, fileTypeInput, networkCoordinatesInput,
            networkInteractionsInput, geneSymbolnput, readExpColumn,
            contrastValueInput, zoomValueInput, resolutionValueInput,
            smoothValueInput, expressionLog, contourLevi, setcolor)

}

