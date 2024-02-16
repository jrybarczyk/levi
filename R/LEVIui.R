#' @title LEVIui
#' @description Launch the Levi Graphical User Interface (GUI) in local machine.
#' @details This function launch the LEVI Graphical User Interface. The
#' interface provides the same tools available in the script mode. There are
#' two tools only available in the user interface: 1) Selection of area from
#' heatmap to calculate the gene expression levels in the area selected;
#' 2) Selection of the genes in some specific area from the image.
#' @usage LEVIui(browser)
#' @param browser This argument is necessary to launch Levi GUI. To launch Levi
#' in the web browser the argument required "TRUE". To launch Levi in the R
#' environment the argument required "FALSE". The default is "FALSE"
#' @return return a GUI
#' @export
#' @author José Rafael Pilan <rafael.pilan@unesp.br> &
#' Isabelle Mira da Silva (isabelle.silva@unesp.br)
#' @examples
#' LEVIui(browser)
#' #LEVIui(browser=TRUE)  #Launch Levi to Browser.
#' #LEVIui(borwser=FALSE) #Launch Levi to R environment.
#'

LEVIui <-function(browser){
tryCatch({
    if (!is(browser, "logical"))
        stop("'browser' must be a TRUE or FALSE")
    if (browser) {
        shiny::runApp(system.file("shiny",package="levi"),
        launch.browser = TRUE)
    } else {
        shiny::runApp(system.file("shiny",package="levi"))
    }
},
error = function(e) {print("Parameter must be TRUE or FALSE")
}
)}


