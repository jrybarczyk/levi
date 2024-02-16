suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(igraph))
suppressMessages(require(colorspace))
suppressMessages(require(grid))
suppressMessages(require(Rcpp))
suppressMessages(require(shiny))
suppressMessages(require(shinyjs))
suppressMessages(require(shinydashboard))
suppressMessages(require(dplyr))
suppressMessages(require(DT))
suppressMessages(require(xml2))



if(getRversion() >= "3.4.0") utils::globalVariables(c("V1", "X", "Y",
"Expression", "Gene", "SigCoordPiso", "matrix_entrada", "matrix_saida",
"renderText", "matrixOutFun"))


levi_shiny = function(expression, fileType, networkCoord,
networkInterac, geneSymbol,baseTest,baseControl, alphaValue, betaValue,
backValue, smoothValue, expressionLog){

    #Configuration of contrast, resolution, smothing and zoom
    #contrast
    if (missing(alphaValue)) {alphaValue <- 50}
    if (alphaValue < 0) {alphaValue <- 0}
    if (alphaValue > 100) {alphaValue <- 99}
    alphaValue<-(alphaValue/100)
    alphaValue<-0.1-(0.1*alphaValue)

    #resolution
    if (missing(backValue)) {backValue <- 50}
    backValue<-as.integer((backValue/100)*210+30)

    #zoom
    if (missing(betaValue)) {betaValue <- 50}
    betaValue<-(betaValue/100)
    betaValue<-(0.2*betaValue)-0.2

    #smothing
    if (missing(smoothValue)) {smoothValue <- 50}
    smoothValue=(smoothValue/100)*18
    smoothValue=as.integer(smoothValue)
    if (smoothValue <= 0) {smoothValue <- 1}

    a<-sqrt(betaValue*betaValue)
    b<-1+a-betaValue
    gamaValue<-sqrt(b^2 + b^2)
    increase<-b/(backValue-1)

    #Remove"NA" and "-" from expression file
    expression <- filter(expression,expression[,paste(geneSymbol)] != "NA")
    expression <- filter(expression,expression[,paste(geneSymbol)]   != "-")
    expression <- unique(expression)
    headExpress = as.list(names(expression))
    if (baseControl == " ") {
        baseControl <- baseTest
    }
    arguments <- list(geneSymbol, baseTest, baseControl)
    for (i in seq(arguments)){
        if (!is.element(arguments[i], headExpress) ) {
            stop(paste0("This argument do not exist in this dataframe: ",
                        arguments[i]))}
    }
    #if (expressionLog) {
    #    expressSelect =expression[,c(geneSymbol, baseTest, baseControl)]
    #    expressSelect[,2:3] <- 2^expressSelect[,2:3]
    #} else {
    expressSelect =expression[,c(geneSymbol, baseTest, baseControl)]
    #}
    newExpression <- aggregate(x = expressSelect[c(baseControl,baseTest)],
                               by = expressSelect[c(geneSymbol)],
                               FUN = function(media_valor){mean(media_valor)
                               })



    if (fileType == "dat"){edges <-dplyr::select(networkCoord, 1,2)}
    if (fileType == "net"){edges <-networkCoord}
    if (fileType == "dyn"){edges <-networkCoord}
    if (fileType == "stg"){edges <-networkCoord}

    if (fileType == "dat"){
        nodes <- as.data.frame(networkInterac)
        nodes[, c(2)] <- sapply(nodes[, c(2)], as.double)
        nodes <- arrange(nodes, V1)}

    if (fileType == "net") {
        nodes <- as.data.frame(networkInterac)
        nodes <- arrange(nodes, V1)
    }

    if (fileType == "dyn") {
        nodes <- as.data.frame(networkInterac)
        nodes <- arrange(nodes, V1)
    }

    if (fileType == "stg"){
        nodes <- as.data.frame(networkInterac)
        nodes <- arrange(nodes, V1)
    }


    #signalCoordMerge have values of controle and test
    signalCoordMerge <- merge(nodes, newExpression, by.x = "V1",
                              by.y = geneSymbol,
                              all.x = TRUE)
    #signalCoordMerge[is.na(signalCoordMerge)] <- 0
    listLink<- unique(edges[,c(1,2)])

    graph_edge <- graph.edgelist(as.matrix(listLink), directed = FALSE)
    edgesGraph <- as_long_data_frame(graph = graph_edge)

    colnames(edgesGraph) <- NULL
    colnames(edgesGraph) <- c("a", "b", "c", "V1")
    nodes$V1 = as.character(nodes$V1)
    nodes$V2 = as.numeric(nodes$V2)
    nodes$V3 = as.numeric(nodes$V3)
    nodesCoord <- aggregate(x = nodes[,c(2:3)], by = nodes[1], FUN = mean)


    ###########################################################################
    #merge edgesGraph and nodesCoord
    edgesNodesMerge <- merge(edgesGraph, nodesCoord, by.x = "c", by.y = "V1",
                             all.x = TRUE)
    edgesNodesMerge <- merge(edgesNodesMerge, nodesCoord, by.x = "V1",
                             by.y = "V1",
                             all.x = TRUE)
    edgesNodesMerge <- edgesNodesMerge[,c(3,4,2,1,5,6,7,8)]

    edgesNodesMerge$V2 <- (edgesNodesMerge[,c(5)] + edgesNodesMerge[,c(7)])/2
    edgesNodesMerge$V3 <- (edgesNodesMerge[,c(6)] + edgesNodesMerge[,c(8)])/2

    edgesSignalMerge <- merge(edgesGraph, signalCoordMerge, by.x = "c",
                              by.y = "V1", all.x = TRUE)
    edgesSignalMerge <- merge(edgesSignalMerge, signalCoordMerge, by.x = "V1",
                              by.y = "V1", all.x = TRUE)
    naColumnA <- as.matrix(
        edgesSignalMerge[!complete.cases(edgesSignalMerge[,5]),2],
        stringsAsFactors = FALSE)
    naColumnB <- as.matrix(
        edgesSignalMerge[!complete.cases(edgesSignalMerge[,9]),1],
        stringsAsFactors = FALSE)
    naTotal <- unique(rbind(naColumnA, naColumnB))

    #Create title for chart
    if (baseTest == baseControl) {
        titleChart <- baseTest
    } else {
        titleChart <- paste(baseTest,baseControl, sep = '-')
    }

    #Creates log if exists nodes without expression value
    if (length(naTotal) > 0) {
        showNotification(paste0("There are ",nrow(naTotal),
                                " nodes without expression value, see log
                                in path: ", file.path(tempdir(),titleChart,
                                "levi.log")), type = "warning",
                         closeButton = FALSE, duration = 10)

        if (!file.exists(file.path(tempdir(), titleChart))){
            dir.create(file.path(tempdir(), titleChart))
        }

        file.path(tempdir(),titleChart, "levi.log")
        levi_log <- file(file.path(tempdir(),titleChart, "levi.log"),
                         open = "wt")
        sink(levi_log)
        sink(levi_log, type = "message")
        print(as.vector(naTotal))
        sink(type = "message")
        sink()
    }

    #edgesSignalMerge <- edgesSignalMerge[complete.cases(edgesSignalMerge),]

    edgesSignalMerge <- edgesSignalMerge[,c(3,4,2,1,5,6,7,8,9,10,11,12)]

    edgesSignalMerge$V2 <- (edgesSignalMerge[,c(8)] +
                                edgesSignalMerge[,c(12)])/2
    edgesSignalMerge$V3 <- (edgesSignalMerge[,c(7)] +
                                edgesSignalMerge[,c(11)])/2
    ###########################################################################

    nnodes <- nrow(nodes)
    nedges <- nrow(edges)
    numberAll<-nnodes+nedges
    coordAll <- rbind(nodesCoord[,c(2,3)], edgesNodesMerge[,c(9,10)])

    posCoordAll <- rbind(nodesCoord[1], edgesNodesMerge[4])
    coordAll <- rbind(nodesCoord[,c(2,3)], edgesNodesMerge[,c(9,10)])

    signalExpAll <- data.frame(V1 = c(signalCoordMerge[,c(5)],
                                      edgesSignalMerge[,c(13)]))
    signalCtrlAll <- data.frame(V1 = c(signalCoordMerge[,c(4)],
                                       edgesSignalMerge[,c(14)]))


    numberCoord <- numberAll
    posCoord <- posCoordAll
    coord <- as.matrix(coordAll)
    signalExp <- as.matrix(signalExpAll)

    if (baseTest == baseControl){
        signalCtrl <- matrix(data = 1, ncol = 1, nrow = nrow(signalCtrlAll))
    } else {
        signalCtrl <- as.matrix(signalCtrlAll)
    }
    signalExp[is.na(signalExp)] <- 0.5
    signalCtrl[is.na(signalCtrl)] <- 0.5

    if (expressionLog) {
        SignalOut<-scale(signalExp/(signalExp+signalCtrl))
        SignalOut <- abs(min(SignalOut))+SignalOut
        SignalOut <- (SignalOut-min(SignalOut))/(max(SignalOut)-min(SignalOut))
    } else {
        SignalOut<-signalExp/(signalExp+signalCtrl)
        SignalOut <- (SignalOut-min(SignalOut))/(max(SignalOut)-min(SignalOut))
    }


    a <- max(signalExp)
    b <- max(signalCtrl)
    signalExp<-((signalExp/a)*0.95)+0.05
    signalCtrl<-((signalCtrl/b)*0.95)+0.05

    # normalization and centralization
    a<-min(coord[,c(1)])
    b<-max(coord[,c(1)])
    centroX<-(a+b)/2
    c<-min(coord[,c(2)])
    d<-max(coord[,c(2)])
    centroY<-(c+d)/2


    if (b >= d) {e <- b}
    if (d > b) {e <- d}
    centroX <- centroX/e
    centroY <- centroY/e

    coord[,c(1)]<- (coord[,c(1)]/e)+(0.5-centroX)
    coord[,c(2)]<-(coord[,c(2)]/e)+(0.5-centroY)


    #Applies the calculation and takes the smallest value for coordinates coord


    cppFunction('NumericMatrix SigCoordPiso(NumericMatrix coord,
                int backValue, double gamaValue, double increase,
                double alphaValue, double betaValue, int numberCoord) {

                double r1 = gamaValue, a=0, b=0, c=0, p=betaValue,
                q= betaValue, t = 0;

                NumericMatrix coordPiso(backValue,backValue);

                for(int i=0; i<backValue; i++){
                q = betaValue;
                for(int j=0; j<backValue; j++){
                for(int m=0; m<numberCoord; m++){
                a = (double) (p-coord(m,0))*(p-coord(m,0));
                b = (double) (q-coord(m,1))*(q-coord(m,1));
                t = (double) a + (double) b;
                c = (double) sqrt((double) t);

                if ((double) r1 > (double) c) {r1 = c;}
                }
                if((double) r1 > (double) alphaValue) {coordPiso(i,j) = 0;}
                else {coordPiso(i,j) = 10;}
                r1 = gamaValue;
                q = q+increase;
                }
                p = p+increase;
                }
                return coordPiso;
                }')

    coordPiso <- SigCoordPiso(coord= coord,  backValue = backValue,
                              gamaValue= gamaValue,  increase = increase,
                              alphaValue = alphaValue, betaValue=betaValue,
                              numberCoord=numberCoord)


    cppFunction('List matrix_entrada(NumericMatrix coordPiso,
                NumericMatrix SignalOut,NumericMatrix signalExp,
                NumericMatrix signalCtrl, NumericMatrix coord, int backValue,
                CharacterVector posCoord, double increase, double betaValue,
                int numberCoord) {

                List resultado;
                double p=betaValue, q= betaValue;
                int h =0;

                NumericMatrix matrixIn(numberCoord+(backValue*backValue),5);
                CharacterVector matrixInPos(numberCoord+(backValue*backValue));
                NumericMatrix bandcoord (numberCoord, 1);

                for(int i=0; i<backValue; i++){
                for(int j=0; j<backValue; j++){
                for(int m=0; m<numberCoord; m++){
                if ((coord(m,0) < q) & (bandcoord(m,0) == 0)) {
                matrixIn(h,0)=coord(m,0);
                matrixIn(h,1)=coord(m,1);
                matrixIn(h,2)=SignalOut(m,0);
                matrixIn(h,3)=signalExp(m,0);
                matrixIn(h,4)=signalCtrl(m,0);
                matrixInPos(h)=posCoord(m);
                h=h+1;
                bandcoord(m,0) = 10;
                } else {
                if ((coord(m,0) == q) & (coord(m,1) < p) &
                (bandcoord(m,0) == 0)) {
                matrixIn(h,0)=coord(m,0);
                matrixIn(h,1)=coord(m,1);
                matrixIn(h,2)=SignalOut(m,0);
                matrixIn(h,3)=signalExp(m,0);
                matrixIn(h,4)=signalCtrl(m,0);
                matrixInPos(h)=posCoord(m);
                h=h+1;
                bandcoord(m,0) = 10;
                }
                }
                }
                if (coordPiso(i,j) == 0) {
                matrixIn(h,0)=q;
                matrixIn(h,1)=p;
                matrixIn(h,2)=0;
                matrixIn(h,3)=0;
                matrixIn(h,4)=0;
                h=h+1;
                }
                p=p+increase;
                }
                p=betaValue;
                q=q+increase;
                }
                resultado["m1"] = matrixIn;
                resultado["m2"] = matrixInPos;
                resultado["m3"] = h;

                return resultado;
                }')

    matrix_resultado <- matrix_entrada(coordPiso= coordPiso,
                        SignalOut= SignalOut, coord= coord, backValue=
                        backValue,
                        signalExp = signalExp, signalCtrl = signalCtrl,
                        posCoord= as.matrix(posCoord), increase= increase,
                        betaValue= betaValue, numberCoord= numberCoord)

    matrixIn <- as.matrix(matrix_resultado$m1)
    matrixInPos <- as.matrix(matrix_resultado$m2)
    h <- matrix_resultado$m3

    h<-h-1 #

    cppFunction('List matrixOutFun(NumericMatrix matrixIn, int backValue,
                double smoothValue, double gamaValue, double increase,
                double betaValue, int h, CharacterMatrix matrixInPos) {

                List resultado;
                double a=0, b=0, c=0, p=betaValue, q= betaValue,
                t = 0, s1 = 0, s2 = 0, s3 = 0;

                NumericMatrix vetInPos(h+1,1);
                NumericMatrix matrixOut(backValue,backValue);
                NumericMatrix matrixOutExp(backValue,backValue);
                NumericMatrix matrixOutCtrl(backValue,backValue);
                CharacterMatrix matrixOutPos(backValue,backValue);
                NumericMatrix vetOut(smoothValue+1,1);

                for(int i=0; i<backValue; i++){
                q = betaValue;
                for(int j=0; j<backValue; j++){
                for(int m=0; m<h; m++){
                a = (p-matrixIn(m,0))*(p-matrixIn(m,0));
                b = (q-matrixIn(m,1))*(q-matrixIn(m,1));
                c = sqrt(a+b);
                vetInPos(m,0) = c;
                }


                for(int m=0; m<smoothValue; m++){
                vetOut(m,0)=smoothValue+1;
                vetInPos(vetOut(m,0),0) = gamaValue;
                for(int n=0; n<h; n++){
                for(int l=0; l<m; l++){
                t = 0;
                if (n == vetOut(l,0)){
                t = 1;
                }
                }
                if (((vetInPos(n,0) < vetInPos(vetOut(m,0),0)) & (t!=1))){
                vetOut(m,0) = n;
                matrixOutPos(i,j) = matrixInPos(vetOut(m,0),0);}
                }
                }
                s1 = 0;
                s2 = 0;
                s3 = 0;
                for(int m=0; m<smoothValue; m++){
                s1=s1+matrixIn(vetOut(m,0),2);
                s2=s2+matrixIn(vetOut(m,0),3);
                s3=s3+matrixIn(vetOut(m,0),4);

                }

                //Performs smoothing by joining values in vetOut
                //divided by the smoothing value

                matrixOut(i,j) = s1/smoothValue;
                matrixOutExp(i,j) = s2/smoothValue;
                matrixOutCtrl(i,j) = s3/smoothValue;

                q = q+increase;
                }
                p = p+increase;
                }
                resultado["m1"] = matrixOut;
                resultado["m2"] = matrixOutExp;
                resultado["m3"] = matrixOutCtrl;
                resultado["m4"] = matrixOutPos;

                return resultado;
                }')


    matrixFinal <- matrixOutFun(matrixIn =  matrixIn,backValue = backValue,
                   gamaValue =  gamaValue, increase = increase,
                   betaValue = betaValue, h = h,
                   matrixInPos = matrixInPos, smoothValue = smoothValue)


    matrixOut <- matrixFinal$m1
    matrixOutExp <- matrixFinal$m2
    matrixOutCtrl <- matrixFinal$m3
    matrixOutPos <- matrixFinal$m4


    n<-backValue
    ExpCtrlPos <- matrix(data = 0, ncol = n, nrow = n)

    i <- seq_len(n)
    j <- seq_len(n-1)
    b <- max(matrixOutExp[i, j+1])
    c <- max(matrixOutCtrl[i, j+1])


    matrixOutExp[i, i] <- matrixOutExp[i, i]/b
    matrixOutCtrl[i, i] <- matrixOutCtrl[i, i]/c

    ExpCtrl <- matrixOut[i, rev(i)]
    exp <- matrixOutExp[i, rev(i)]
    ctrl <- matrixOutCtrl[i, rev(i)]
    ExpCtrlPos <- matrixOutPos[i, rev(i)]


    if (baseTest == baseControl){
        landgraph <- melt(exp, value.name = "z")
    } else {
        landgraph <- melt(ExpCtrl, value.name = "z")
    }
    landgraphPos <- melt(ExpCtrlPos, value.name = "z1")
    landgraphFinal <- as.data.frame(cbind(landgraph[,c(1,2,3)],
                                          landgraphPos[,c(3)]))
    landgraphFinal$z<-round(landgraphFinal$z, 2)
    colnames(landgraphFinal) <- c("X", "Y", "Expression", "Gene")

    return(landgraphFinal)
}

# UI as responsive application
ui <- fluidPage(
    shinyjs::useShinyjs(),
    br(),

    # Sidebar layout
    sidebarLayout(
        # Inputs
        sidebarPanel(

            #adds ID to be used as Shiny input to know which tab will be selected.
            sidebarMenu(id = "tab"
            ),

            tabsetPanel(id = "tabset_id", selected = "t1",
                    tabPanel("File", value = "t1",
                            br(),
                            selectInput("fileType", "Network input type:",
                                        c("Medusa (DAT)" = "dat",
                                          "RedeR (DYN)" = "dyn","Pajek (NET)" =
                                           "net", "STRING / STITCH" = "stg")),
                            br(),
                            uiOutput("out2"),
                            fileInput("file","Upload the expression file:"),
                            checkboxInput("log",
                                "Expression values in log scale", FALSE),
                            radioButtons("fields", "Selected fields:",
                                         c("Two Samples" = "twofields",
                                           "One Sample" = "onefield"),
                                         selected = "twofields",
                                         inline = TRUE),
                                 uiOutput("out1"),
                                 helpText("Choose only one"),
                                 checkboxInput(inputId = "multi",
                                               label = "Multicolor"),
                                 checkboxInput(inputId = "twoc",
                                               label = "Two colors",
                                               value = FALSE),
                                 conditionalPanel("input.twoc",
                                selectInput(inputId = "setcolor",
                                            label = "Representation options",
                                            choices =
                                    c("Purple & Pink" = "pp",
                                    "Green & Blue" = "gb",
                                    "Blue & Yellow" = "by",
                                    "Pink & Green" = "pg",
                                    "Orange & Purple" = "op",
                                    "Blue-green & Marine" = "bm"))),
                                 checkboxInput("contour",
                                               "Chart with contour", FALSE),
                                 tags$hr(),
                                 actionButton('action','Run',
                                              style="float:right"),
                                 br()
                        ),
                        tabPanel("Settings", value = "t2", br(),
                                 helpText("Re-run to apply changes"),
                                 sliderInput("contrast", "Contrast:", min = 0,
                                             max = 100,
                                             value = 50),
                                 sliderInput("size", "Resolution:", min = 1,
                                             max = 100, value = 50),
                                 sliderInput("smooth", "Smoothing:",
                                             min = 0, max = 100, value = 50),
                                 sliderInput("zoom", "Zoom:",
                                             min = 0, max = 100, value = 50)))
        ),

        # Output:
        mainPanel(
            tags$head(tags$style(type="text/css", "#loadmessage
            { position: fixed;top: 0px;left: 0px;width: 100%; padding: 5px 0px
            5px 0px;text-align: center;font-weight: bold;font-size:
            100%;color: #000000;background-color: #81efea;z-index: 105;}")),
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("Loading...",id="loadmessage")),



            plotOutput(outputId = "graph", brush = brushOpts(id = "plotBrush",
            resetOnNew = TRUE)),
            h4(textOutput("expArea")),
            uiOutput("out3"),
            br(),br(),br(),

            # Show data table
            dataTableOutput(outputId = "landdatatable")
            )
            )
        )

# Defines the graph function required to generate the scatterplot
server <- function(input, output, session) {
    options(shiny.maxRequestSize=30*1024^2)
    set.seed(123)


    baseSelect <- reactive({
        if(is.null(input$file)){return()}
        tryCatch({
            read.table(file=input$file$datapath, header =TRUE, sep="\t")
        },
        warning=function(w) {showNotification("Incorrect file format",
        type = "warning", closeButton = FALSE, duration = 5)
            return(NA)
        },
        error=function(e) {showNotification("Incorrect file format",
        type = "error", closeButton = FALSE, duration = 5)
            return(NULL)
        })})

    output$out1 <- renderUI({
        dynUi <- tabPanel("", value = "t1",
                          selectInput("geneSymbol", "Select ID field",
                          c(names(baseSelect())),
                          selected = NULL),
                          selectInput("baseTest", "Select test field",
                                      c(names(baseSelect())),
                                      selected = NULL),
                          if(input$fields == "onefield") {
                              selectInput("baseControl",
                                          "Select control field", c(" "),
                              selected = NULL)
                          } else {
                              selectInput("baseControl", "Select control field",
                              c(names(baseSelect()) ),
                              selected = NULL)},
                          helpText("Chart colors can be changed without
                                   having to re-run")
                          )

        return(dynUi)
    })


    inputString <- observeEvent(input$fileType, {
        if(input$fileType == "stg"){
            output$out2 <- renderUI({
                dynUi2 <- tabPanel("", value = "t1",
                 fileInput("file2","Upload the coordinates file:"),
                 fileInput("file3","Upload the interactions file:"))
                return(dynUi2)
            })}
        else {
            output$out2 <- renderUI({
                dynUi2 <- tabPanel("", value = "t1",
                fileInput("file2","Upload the network file:"))
                return(dynUi2)
            })}
    })

    #Allows the graph to be displayed only after executing the "run"
    v <- reactiveValues()
    observeEvent(input$action, {



        v$data0 <- 1
        v$fileType <- input$fileType
        if (v$fileType == "dat") {
            tryCatch({


                data1 <- input$file
                if(is.null(data1)){return()}
                v$expression <- read.delim(file = data1$datapath, header = TRUE,
                sep = "\t",
                quote = "")


                data2 <- input$file2
                if(is.null(data2)){return()}
                nodeslist <- read.delim(file = data2$datapath, header = FALSE,
                sep = "\t",stringsAsFactors=FALSE, fill = TRUE, col.names =
                paste0("V",seq_len(max(count.fields(data2$datapath,
                sep = '\t')))))


                delimiter <- which(nodeslist == "*nodes")

                networkCoord <- slice(nodeslist, 3:delimiter-1)
                v$networkCoord <- networkCoord[,c(1,2)]
                networkInterac <- slice(nodeslist, delimiter+1:nrow(nodeslist))
                v$networkInterac <- networkInterac[,c(1,2,3)]
            },
            warning=function(w) {showNotification("Incorrect file format",
            type = "warning",
            closeButton = FALSE, duration = 5)
                return(NA)
            },
            error=function(e) {showNotification("Incorrect file format",
            type = "error",
            closeButton = FALSE, duration = 5)
                return(NULL)
            }
            )
        }

        if (v$fileType == "stg") {
            tryCatch({
                data1 <- input$file
                if(is.null(data1)){return()}


                v$expression <- read.delim(file = data1$datapath, header = TRUE,
                sep = "\t", quote = "")

                data2 <- input$file2
                if(is.null(data2)){return()}
                data3 <- input$file3
                if(is.null(data3)){return()}

                networkInterac <- read.delim(file = data2$datapath,
                header = TRUE,
                sep = "\t", stringsAsFactors=FALSE, fill = TRUE)
                networkCoord <- read.delim(file = data3$datapath, header = TRUE,
                sep = "\t", stringsAsFactors=FALSE, fill = TRUE)

                networkCoord <- networkCoord[,c(1,2)]
                networkInterac <- networkInterac[,c(1,2,3)]
                colnames(networkCoord) <- c("V1", "V2")
                colnames(networkInterac) <- c("V1", "V2", "V3")
                v$networkCoord <- networkCoord
                v$networkInterac <- networkInterac

            },
            warning=function(w) {showNotification("Incorrect file format
            (STRING/STITCH) or missing file", type = "warning",
            closeButton = FALSE, duration = 5)
                return(NA)
            },
            error=function(e) {showNotification("Incorrect file format
            (STRING/STITCH) or missing file", type = "error",
            closeButton = FALSE, duration = 5)
                return(NULL)
            }
            )
        }

        if (v$fileType == "net") {
            tryCatch({


                data1 <- input$file
                if(is.null(data1)){return()}


                v$expression <- read.delim(file = data1$datapath, header = TRUE,
                sep = "\t", quote = "")


                data2 <- input$file2
                if(is.null(data2)){return()}
                netRead <- read.delim(file = data2$datapath, header = FALSE,
                stringsAsFactors=FALSE)

                delimiter_edge <- which(netRead == "*Edges")
                networkCoord <- data_frame()

                edgesSL <- as.data.frame(slice(netRead,
                delimiter_edge+1:nrow(netRead)))

                delimiter_nodes_end <- which(netRead == "*Edges")
                nodesSl <- as.data.frame(
                    slice(netRead, 3:delimiter_nodes_end-1))
                networkInterac <- data_frame()
                for (i in seq_len(nrow(nodesSl))) {
                    nodesRt <- read.table(text = as.character(nodesSl[i,1]),
                                          sep = " ")

                    nodesFt <- Filter(function(x)!all(is.na(x)), nodesRt)
                    nodesFt <sss- nodesFt[,c(1,2,3,4)]
                    colnames(nodesFt) <- c("V1", "V2","V3", "V4")
                    nodesFt <- data.frame(lapply(nodesFt, function(x)
                    {gsub("FALSE", "F", x)}), stringsAsFactors = FALSE)
                    nodesFt <- data.frame(lapply(nv$matrixSizeodesFt,
                    function(x)
                    {gsub("TRUE", "T", x)}), stringsAsFactors = FALSE)
                    networkInterac <-rbind(networkInterac, nodesFt)
                }


                for (i in seq_len(nrow(edgesSL))) {
                    edgesRt <- read.table(text = as.character(edgesSL[i,1]),
                                          sep = " ")
                    edgesFt <- Filter(function(x)!all(is.na(x)), edgesRt)
                    edgesFt <- edgesFt[,c(1,2)]
                    colnames(edgesFt) <- c("V1", "V2")
                    networkCoord <-rbind(networkCoord, edgesFt)

                }

                netMerge<- merge(networkCoord, networkInterac, by.x = "V1",
                                 by.y = "V1",
                                 all.x = FALSE)
                colnames(netMerge) <- c("a", "b", "c", "d", "e")
                netMerge<- merge(netMerge, networkInterac, by.x = "b",
                                 by.y = "V1", all.x = FALSE)
                networkCoord <- netMerge[,c(3,6)]
                colnames(networkCoord) <- c("V1", "V2")

                networkInterac <- networkInterac[,c(2,3,4)]
                colnames(networkInterac) <- c("V1", "V2", "V3")
                v$networkInterac <-networkInterac
                v$networkCoord <- networkCoord


            },
            warning=function(w) {showNotification("Incorrect file format",
            type = "warning", closeButton = FALSE, duration = 5)
                return(NA)
            },
            error=function(e) {showNotification("Incorrect file format",
            type = "error", closeButton = FALSE, duration = 5)
                return(NULL)
            }
            )}


        if (v$fileType == "dyn") {
            tryCatch({


                data1 <- input$file
                if(is.null(data1)){return()}


                v$expression <- read.delim(file = data1$datapath, header = TRUE,
                                           sep = "\t", quote = "")


                data2 <- input$file2
                if(is.null(data2)){return()}
                tf <- tempfile(tmpdir = tdir <- tempdir())
                dynFiles <- unzip(data2$datapath, exdir = tdir)
                dynRead <- read_xml(dynFiles , stringsAsFactors=FALSE)


                dynLabel <- xml_find_all(dynRead, xpath = "//*/*/@label")
                vals <- trimws(xml_text(dynLabel))
                dynDf = as.data.frame(vals, stringsAsFactors = FALSE)

                dynId <- xml_find_all(dynRead, xpath = "//*/*/@id")
                valsId <- trimws(xml_text(dynId))
                networkInterac <- as.data.frame(slice(dynDf, 1:length(valsId)))
                networkInterac$V1 <- seq(0,nrow(networkInterac)-1)

                dynSource <- xml_find_all(dynRead, xpath = "//*/*/@source")
                dynTarget <- xml_find_all(dynRead, xpath = "//*/*/@target")
                dynX <- xml_find_all(dynRead, xpath = "//*/*/@x")
                dynY <- xml_find_all(dynRead, xpath = "//*/*/@y")

                datasource_tmp <- as.data.frame(lapply(dynSource, gsub,
                pattern = "source=", replacement = "", fixed = TRUE))
                datasource <- as.data.frame(lapply(datasource_tmp, gsub,
                pattern = "\"", replacement = "", fixed = TRUE),
                stringsAsFactors = FALSE)
                colnames(datasource) <- NULL

                datatarget_tmp <- as.data.frame(lapply(dynTarget, gsub,
                pattern = "target=", replacement = "", fixed = TRUE))
                datatarget <- as.data.frame(lapply(datatarget_tmp, gsub,
                pattern = "\"",
                replacement = "", fixed = TRUE), stringsAsFactors = FALSE)
                colnames(datatarget) <- NULL
                datateste <- as.data.frame(cbind(t(datasource), t(datatarget)),
                stringsAsFactors = FALSE)

                dataXTmp <- as.data.frame(lapply(dynX, gsub, pattern = "x=",
                replacement = "", fixed = TRUE))
                datax <- as.data.frame(lapply(dataXTmp, gsub, pattern = "\"",
                replacement = "", fixed = TRUE), stringsAsFactors = FALSE)
                colnames(datax) <- NULL

                dataYTmp <- as.data.frame(lapply(dynY, gsub, pattern = "y=",
                replacement = "", fixed = TRUE))
                datay <- as.data.frame(lapply(dataYTmp, gsub, pattern = "\"",
                replacement = "", fixed = TRUE), stringsAsFactors = FALSE)
                colnames(datay) <- NULL

                networkCoord <- datateste
                networkCoord$V1 <- as.numeric(networkCoord$V1)
                networkCoord$V2 <- as.numeric(networkCoord$V2)
                t1<- merge(networkCoord, networkInterac, by.x = "V1",
                           by.y = "V1")
                colnames(t1) <- c("a", "b", "c")
                t1<- merge(t1, networkInterac, by.x = "b", by.y = "V1",
                           all.x = FALSE)
                networkCoord <- as.matrix(t1[,c(3,4)])
                colnames(networkCoord) <- c("V1", "V2")

                v$networkCoord <- networkCoord

                networkInterac <- as.data.frame(cbind(networkInterac[,1],
                t(datax),t(datay)), stringsAsFactors = FALSE)

                networkInterac$V1 = as.character(networkInterac$V1)
                networkInterac$V2 = as.numeric(networkInterac$V2)
                networkInterac$V3 = as.numeric(networkInterac$V3)
                v$networkInterac <- networkInterac

            },
            warning=function(w) {showNotification("Incorrect file format",
            type = "warning", closeButton = FALSE, duration = 5)
                return(NA)
            },
            error=function(e) {showNotification("Incorrect file format",
            type = "error", closeButton = FALSE, duration = 5)
                return(NULL)
            }
            )}
        })

    # Creates the object of the graph
    observeEvent(input$action, {
        if(is.null(v$data0)){
            return()
        }
        tryCatch({

            input$log
            v$func_ne_return <-levi_shiny(v$expression, v$fileType,
            v$networkCoord, v$networkInterac, input$geneSymbol,
            input$baseTest,input$baseControl , alphaValue = input$contrast,
            betaValue = input$zoom,backValue= input$size,
            smoothValue = input$smooth, expressionLog = input$log)

            output$out3 <- renderUI({

                dynUi <- tabPanel("", value = "t1",
                selectInput("plotType", "", c("TIFF" = "tiff","BMP" = "bmp",
                "JPEG" = "jpeg", "PNG" = "png"), width = '135px'),
                downloadButton('downloadPlot', 'Download Plot',
                style="float:left"))

                return(dynUi)
            })
        },
        warning=function(w) {showNotification("Error", type = "warning",
            closeButton = FALSE, duration = 5)
            return(NA)
        },
        error=function(e) {showNotification("Error", type = "error",
            closeButton = FALSE, duration = 5)
            return(NULL)
        }
        )
    })

    v$graphVC <- renderPlot({
        if(is.null(v$func_ne_return)){
            return()
        }
        v$matrixSize <- sqrt(NROW(v$func_ne_return))
        if (input$multi) {v$multicolor <- c("#180052", "#0c0083",
        "#0000b4","#0000e4","#0010ff", "#0041ff", "#0072ff", "#00A3FF",
        "#00D4FF", "#00FF49","#5AFF00", "#FFE400", "#FFC400", "#FFA300",
        "#FF8300", "#FF6200", "#FF4100","#FF2100", "#FF0000", "#E40000")}

        v$matrixSize <- sqrt(NROW(v$func_ne_return))
        if (input$setcolor == "pp"){v$colorSet <- c("#4e5052", "#b387e6",
            "#ff0000")}
        if (input$setcolor == "gb"){v$colorSet <- c("#4e5052", "#a4db56",
            "#1d02c9")}
        if (input$setcolor == "by"){v$colorSet <- c("#4e5052", "#27b0cf",
            "#ffec2b")}
        if (input$setcolor == "pg"){v$colorSet <- c("#4e5052", "#f757ca",
            "#04ff00")}
        if (input$setcolor == "op"){v$colorSet <- c("#4e5052", "#fcbb63",
            "#7300c4")}
        if (input$setcolor == "bm"){v$colorSet <- c("#4e5052", "#5df0b0",
            "#360d94")}
        v$graficof <-ggplot(data = v$func_ne_return, aes(x = X, y = Y))+
            geom_raster(aes(fill = Expression), interpolate = TRUE, hjust = 0.5,
                        vjust = 0.5) +
            theme_void() +
            theme(legend.margin=margin(0,0,0,-20)) +
            annotate("text", x = c(v$matrixSize*1.045,v$matrixSize*1.045),
                y = c(v$matrixSize*0.42,v$matrixSize*0.6),
                     label = c("decrease", "increase") , size=3 , angle=90) +
            annotate("segment", x = v$matrixSize*1.07, xend = v$matrixSize*1.07,
                y = v$matrixSize*0.49, yend = v$matrixSize*0.31,
                colour = "black", size=0.2, alpha=0.6,
                arrow=arrow(type = "closed",length = unit(x = c(0.2),
                units = "cm"))) +
            annotate("segment", x = v$matrixSize*1.07, xend = v$matrixSize*1.07,
                y = v$matrixSize*0.51, yend = v$matrixSize*0.71,
                colour = "black", size=0.2, alpha=0.6,
                arrow=arrow(type =
                "closed",length = unit(x = c(0.2),units = "cm"))) +
            coord_fixed(ratio = 1) +
            if (input$contour == TRUE) {geom_contour(aes(z = Expression))}

        v$landDataTableTmp <- v$func_ne_return
        if (input$multi) {
            v$graficof + scale_fill_gradientn(colours=v$multicolor,
            values=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
            0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
            breaks=seq(0,1,0.2), limits=c(0,1),
            guide = guide_colorbar(title="Expression Level",
            title.position = "right", title.hjust = 0.5,
            title.theme = element_text(angle = 270, size = 9),
            barwidth=1, barheight = 10))
        }else{
            v$graficof + scale_fill_gradientn(colours=v$colorSet,
            values=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
            0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
            breaks=seq(0,1,0.2), limits=c(0,1),
            guide = guide_colorbar(title="Expression Level",
            title.position = "right", title.hjust = 0.5,
            title.theme = element_text(angle = 270, size = 9),
            barwidth=1, barheight = 10))
        }

    }, height="auto")


    # Displays the data of the selected dataframe in the graph
    v$landDataTableTmp <- setNames(data.frame(matrix(ncol = 1, nrow = 0)),
                                   c("Gene"))
    isolate(output$landdatatable <- DT::renderDataTable({

        tryCatch({
            brushedPoints(v$landDataTableTmp, brush = input$plotBrush) %>%
                select(Gene) %>% filter(Gene != "") %>% unique()

        },
        warning=function(w) {showNotification("Error while selecting chart",
            type = "warning", closeButton = FALSE, duration = 5)
            return(NA)
        },
        error=function(e) {showNotification("Error while selecting chart",
            type = "error", closeButton = FALSE, duration = 5)
            return(NULL)
        })
    }, rownames = FALSE))


    observeEvent(input$action, {
        output$graph <- v$graphVC
    })

    output$expArea <- renderText({
        tryCatch({
            expArea <- sum(brushedPoints(v$func_ne_return,
            brush = input$plotBrush) %>% select(Expression),na.rm=TRUE)
            return({paste0("Expression area: ", expArea)})
        },
        error=function(e) {return({paste0("Expression area: ")})

        })})

    #Saves the selected dataframe data in the graph
    downloadNameFun <- reactive({

        if(input$plotType=="png") filename <- paste0("chart",".png",sep="")
        if(input$plotType=="tiff") filename <- paste0("chart",".tif",sep="")
        if(input$plotType=="jpeg") filename <- paste0("chart",".jpg",sep="")
        if(input$plotType=="bmp") filename <- paste0("chart",".bmp",sep="")
        return(filename)

    })

    fn_download <- function(){
        tryCatch({
            df <- v$graficof + scale_fill_gradientn(colours=v$colorSet,
            values=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
            0.55, 0.6,0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
            breaks=seq(0,1,0.2), limits=c(0,1),
            guide = guide_colorbar(title="Expression Level",
            title.position = "right", title.hjust = 0.5,
            title.theme = element_text(angle = 270, size = 9),
            barwidth=1, barheight = 10))

            if(input$plotType=="png") png(downloadNameFun(), width = 480,
                                          height = 480, units = "px")
            if(input$plotType=="tiff") tiff(downloadNameFun(), width = 480,
                                            height = 480,units = "px")
            if(input$plotType=="jpeg") jpeg(downloadNameFun(), height=480,
                                            width=480,units="px", quality=100)
            if(input$plotType=="bmp") bmp(downloadNameFun(), width = 480,
                                          height = 480,units = "px")
            plot(df)
            dev.off()
        },
        error=function(e) {showNotification("Error while selecting chart",
                                            type = "error", closeButton = FALSE,
                                            duration = 5)})}

    #Save plot
    output$downloadPlot <- downloadHandler(
        filename = downloadNameFun,
        content = function(file) {
            try({
                fn_download()
                file.copy(downloadNameFun(), file, overwrite=TRUE)
            })})

}

shinyApp(ui, server)
