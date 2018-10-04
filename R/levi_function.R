levi_function <- function(expressionInput, fileTypeInput, networkNodesInput,
                          networkEdgesInput, geneSymbolnput, readExpColumn,
                          contrastValueInput,zoomValueInput,
                          resolutionValueInput, smoothValueInput, expressionLog,
                          contourLevi, setcolor){

    if (missing(setcolor)) {setcolor <- "default"}

    if (setcolor == "default"){colorSet <- c("#180052", "#0c0083","#0000b4",
        "#0000e4","#0010ff", "#0041ff", "#0072ff", "#00A3FF", "#00D4FF",
        "#00FF49","#5AFF00", "#FFE400", "#FFC400", "#FFA300", "#FF8300",
        "#FF6200", "#FF4100", "#FF2100", "#FF0000", "#E40000")}
    else
        if (setcolor == "terrain"){colorSet <- terrain.colors(n = 20)}
        else
            if (setcolor == "rainbow"){colorSet <- rainbow(n = 20)}
            else
                if (setcolor == "heat"){colorSet <- heat.colors(n = 20)}
                else
                    if (setcolor == "topo"){colorSet <-
                        topo.colors(n = 20)}
                    else
                        if (setcolor == "cm"){colorSet <-
                            cm.colors(n = 20)}
                        else
                            {print("Incorrect color name.")
                                stop()}




    for (i in seq(2,length(readExpColumn))) {

        columnComb<- do.call('rbind',
            strsplit(as.character(readExpColumn[i]),'-',
            fixed=TRUE))

            if (missing(networkEdgesInput)) {networkEdgesInput <- NA}

            if (missing(expressionLog)) {expressionLog <- FALSE}

            #Configuration of contrast, resolution, smothing and zoom
            #contrast
            if (missing(contrastValueInput)) {contrastValue <- 50}
            else
                {contrastValue <- contrastValueInput}
            if (contrastValue < 0) {contrastValue <- 0}
            if (contrastValue > 100) {contrastValue <- 99}
            contrastValue<-(contrastValue/100)
            contrastValue<-0.1-(0.1*contrastValue)

            #resolution
            if (missing(resolutionValueInput)) {resolutionValue <- 50}
            else
                {resolutionValue <- resolutionValueInput}
            if (resolutionValue > 100) {resolutionValue <- 100}
            if (resolutionValue < 1) {resolutionValue <- 1}
            resolutionValue<-as.integer((resolutionValue/100)*210+30)

            #zoom
            if (missing(zoomValueInput)) {zoomValue <- 50}
            else
                {zoomValue <- zoomValueInput}
            if (zoomValue < 0) {zoomValue <- 0}
            if (zoomValue > 100) {zoomValue <- 100}
            zoomValue<-(zoomValue/100)
            zoomValue<-(0.2*zoomValue)-0.2

            #smothing
            if (missing(smoothValueInput)) {smoothValue <- 50}
            else
                {smoothValue <- smoothValueInput}
            smoothValue=(smoothValue/100)*18
            smoothValue=as.integer(smoothValue)
            if (smoothValue <= 0) {smoothValue <- 1}


            a<-sqrt(zoomValue*zoomValue)
            b<-1+a-zoomValue
            gamaValue<-sqrt(b^2 + b^2)
            increase<-b/(resolutionValue-1)

            nameBase <- expressionInput
            networkNodes<- networkNodesInput
            networkEdges <- networkEdgesInput
            geneSymbol<- geneSymbolnput
            baseTest<- columnComb[,1]
            baseControl<- columnComb[,2]
            fileType <- fileTypeInput

        switch(fileType,
                dat={
                    tryCatch({

                    networkNodes <- read.delim(file = networkNodes,
                    header = FALSE, sep = "\t",
                    stringsAsFactors=FALSE, fill = TRUE, col.names =
                    paste0("V",seq_len(max(count.fields(networkNodes,
                    sep = '\t')))))


                    delimiter <- which(networkNodes == "*nodes")

                    edges <- slice(networkNodes, 3:delimiter-1)
                    edges <- edges[,c(1,2)]
                    nodes <- slice(networkNodes,
                        delimiter+1:nrow(networkNodes))
                    nodes <- nodes[,c(1,2,3)]
                    },
                    warning=function(w) {print("Incorrect file format (dat)")
                        return(NA)
                    },
                    error=function(e) {print("Incorrect file format (dat)")
                        return(NULL)
                    }
                    )},

                stg={
                    tryCatch({

                    nodes <- read.delim(file = networkNodes, header = TRUE,
                        sep = "\t", stringsAsFactors=FALSE, fill = TRUE)
                    edges <- read.delim(file = networkEdges, header = TRUE,
                        sep = "\t", stringsAsFactors=FALSE, fill = TRUE)

                    edges <- edges[,c(1,2)]
                    nodes <- nodes[,c(1,2,3)]
                    colnames(edges) <- c("V1", "V2")
                    colnames(nodes) <- c("V1", "V2", "V3")
                    },
                    warning=function(w) {print("Incorrect file format (String)
                    or missing file")
                        return(NA)
                    },
                    error=function(e) {print("Incorrect file format (String)
                    or missing file")
                        return(NULL)
                    }
                    )},

                net={
                    tryCatch({

                    net_read <- read.delim(file = networkNodes, header = FALSE,
                        stringsAsFactors=FALSE)

                    delimiter_edge <- which(net_read == "*Edges")
                    edges <- data_frame()

                    edges_sl <- as.data.frame(slice(net_read,
                        delimiter_edge+1:nrow(net_read)))

                    delimiter_nodes_end <- which(net_read == "*Edges")
                    nodes_sl <- as.data.frame(slice(net_read,
                        3:delimiter_nodes_end-1))
                    nodes <- data_frame()
                        for (i in seq_len(nrow(nodes_sl))) {
                            nodes_rt <- read.table(text =
                            as.character(nodes_sl[i,1]), sep = " ")

                            nodes_ft <- Filter(function(x)!all(is.na(x)),
                            nodes_rt)
                            nodes_ft <- nodes_ft[,c(1,2,3,4)]
                            colnames(nodes_ft) <- c("V1", "V2","V3", "V4")
                            nodes_ft <- data.frame(lapply(nodes_ft, function(x)
                                {gsub("FALSE", "F", x)}),
                                stringsAsFactors = FALSE)
                            nodes_ft <- data.frame(lapply(nodes_ft, function(x)
                                {gsub("TRUE", "T", x)}),
                                stringsAsFactors = FALSE)
                           nodes <-rbind(nodes, nodes_ft)
                        }


                        for (i in seq_len(nrow(edges_sl))) {
                            edges_rt <- read.table(text =
                            as.character(edges_sl[i,1]), sep = " ")
                            edges_ft <- Filter(function(x)!all(is.na(x)),
                            edges_rt)
                            edges_ft <- edges_ft[,c(1,2)]
                            colnames(edges_ft) <- c("V1", "V2")
                            edges <-rbind(edges, edges_ft)

                       }



                        net_mg<- merge(edges, nodes, by.x = "V1", by.y = "V1",
                        all.x = FALSE)
                        colnames(net_mg) <- c("a", "b", "c", "d", "e")
                        net_mg<- merge(net_mg, nodes, by.x = "b", by.y = "V1",
                        all.x = FALSE)
                        edges <- net_mg[,c(3,6)]
                        colnames(edges) <- c("V1", "V2")

                        nodes <- nodes[,c(2,3,4)]
                        colnames(nodes) <- c("V1", "V2", "V3")
                        nodes <-nodes
                        edges <- edges


                    },
                    warning=function(w) {print("Incorrect file format  (net)")
                        return(NA)
                    },
                    error=function(e) {print("Incorrect file format (net)")
                        return(NULL)
                    }
                    )},


                dyn={
                    tryCatch({

                        tf <- tempfile(tmpdir = tdir <- tempdir())
                        dyn_files <- unzip(networkNodes, exdir = tdir)
                        dyn_read <- read_xml(dyn_files , stringsAsFactors=FALSE)


                        dyn_label <- xml_find_all(dyn_read,
                        xpath = "//*/*/@label")
                        vals <- trimws(xml_text(dyn_label))
                        dyn_df = as.data.frame(vals, stringsAsFactors = FALSE)

                        dyn_id <- xml_find_all(dyn_read, xpath = "//*/*/@id")
                        vals_id <- trimws(xml_text(dyn_id))
                        nodes <- as.data.frame(slice(dyn_df, 1:length(vals_id)))
                        nodes$V1 <- seq(0,nrow(nodes)-1)

                        dyn_source <- xml_find_all(dyn_read,
                        xpath = "//*/*/@source")
                        dyn_target <- xml_find_all(dyn_read,
                        xpath = "//*/*/@target")
                        dyn_x <- xml_find_all(dyn_read, xpath = "//*/*/@x")
                        dyn_y <- xml_find_all(dyn_read, xpath = "//*/*/@y")

                        datasource_tmp <- as.data.frame(lapply(dyn_source, gsub,
                        pattern = "source=",
                        replacement = "", fixed = TRUE))
                        datasource <- as.data.frame(lapply(datasource_tmp, gsub,
                        pattern = "\"",
                        replacement = "", fixed = TRUE),
                        stringsAsFactors = FALSE)
                        colnames(datasource) <- NULL

                        datatarget_tmp <- as.data.frame(lapply(dyn_target, gsub,
                        pattern = "target=",
                        replacement = "", fixed = TRUE))

                        datatarget <- as.data.frame(lapply(datatarget_tmp, gsub,
                        pattern = "\"",
                        replacement = "", fixed = TRUE),
                        stringsAsFactors = FALSE)
                        colnames(datatarget) <- NULL
                        datateste <- as.data.frame(cbind(t(datasource),
                        t(datatarget)),
                        stringsAsFactors = FALSE)

                        datax_tmp <- as.data.frame(lapply(dyn_x, gsub,
                        pattern = "x=",
                        replacement = "", fixed = TRUE))
                        datax <- as.data.frame(lapply(datax_tmp, gsub,
                        pattern = "\"",
                        replacement = "", fixed = TRUE),
                        stringsAsFactors = FALSE)
                        colnames(datax) <- NULL

                        datay_tmp <- as.data.frame(lapply(dyn_y, gsub,
                        pattern = "y=",
                        replacement = "", fixed = TRUE))
                        datay <- as.data.frame(lapply(datay_tmp, gsub,
                        pattern = "\"",
                        replacement = "", fixed = TRUE),
                        stringsAsFactors = FALSE)
                        colnames(datay) <- NULL

                        edges <- datateste
                        edges$V1 = as.numeric(edges$V1)
                        edges$V2 = as.numeric(edges$V2)

                        t1<- merge(edges, nodes, by.x = "V1", by.y = "V1")
                        colnames(t1) <- c("a", "b", "c")
                        t1<- merge(t1, nodes, by.x = "b", by.y = "V1",
                        all.x = FALSE)
                        edges <- as.matrix(t1[,c(3,4)])
                        colnames(edges) <- c("V1", "V2")

                        edges <- edges

                        nodes <- as.data.frame(cbind(nodes[,1], t(datax),
                        t(datay)), stringsAsFactors = FALSE)

                        nodes$V1 = as.character(nodes$V1)
                        nodes$V2 = as.numeric(nodes$V2)
                        nodes$V3 = as.numeric(nodes$V3)
                        nodes <- nodes

                    },
                    warning=function(w) {print("Incorrect file format (dyn)")
                        return(NA)
                    },
                    error=function(e) {print("Incorrect file format (dyn)")
                        return(NULL)
                    }
                    )},
                {
                    print('Input a valid file format')
                }
         )


         tryCatch({

             #Remove"NA" and "-" from expression file
            expression <- read.delim(file = nameBase, header = TRUE,
            sep = "\t", quote = "")

            expression <- subset(expression,expression[,paste(geneSymbol)] !=
            "NA")
            expression <- subset(expression,expression[,paste(geneSymbol)] !=
            "-")
            expression <- unique(expression)
            head_express = as.list(names(expression))
            if (baseControl == " ") {
                baseControl <- baseTest
            }
            arguments <- list(geneSymbol, baseTest, baseControl)
            for (i in seq(arguments)){
                if (!is.element(arguments[i], head_express) ) {
                    stop(paste0("This argument do not exist in this dataframe:
                    ", arguments[i]))}
            }
            if (expressionLog) {
                expressSelect =expression[,c(geneSymbol, baseTest, baseControl)]
                expressSelect[,2:3] <- 2^expressSelect[,2:3]
            } else {
                expressSelect =expression[,c(geneSymbol, baseTest, baseControl)]
            }
            newExpression <- aggregate(x = expressSelect[c
                 (baseControl,baseTest)],
                 by = expressSelect[c(geneSymbol)],
                 FUN = function(media_valor){
                     mean(media_valor)
                     })
        },
        warning=function(w) {print("Incorrect Expression file or Gene Symbol
            column")
            return(NA)
        },
        error=function(e) {print("Incorrect Expression file or Gene Symbol
            column")
            return(NULL)
        }
        )


        if (fileType == "dat"){
            edges <-select(edges, -3)
            nodes <- as.data.frame(nodes)
            nodes[, c(2)] <- sapply(nodes[, c(2)], as.double)
        }

        if (fileType == "net"){
            nodes <- as.data.frame(nodes)
        }

        if (fileType == "dyn"){
            nodes <- as.data.frame(nodes)
        }

        if (fileType == "stg"){
            nodes <- as.data.frame(nodes)
        }

        nodes <- arrange(nodes, nodes[,c(1)])


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


        #######################################################################
        #merge edgesGraph and nodesCoord
        edgesNodesMerge <- merge(edgesGraph, nodesCoord, by.x = "c",
            by.y = "V1",
            all.x = TRUE)
        edgesNodesMerge <- merge(edgesNodesMerge, nodesCoord, by.x = "V1",
            by.y = "V1",
            all.x = TRUE)
        edgesNodesMerge <- edgesNodesMerge[,c(3,4,2,1,5,6,7,8)]

        edgesNodesMerge$V2 <- (edgesNodesMerge[,c(5)] +
            edgesNodesMerge[,c(7)])/2
        edgesNodesMerge$V3 <- (edgesNodesMerge[,c(6)] +
            edgesNodesMerge[,c(8)])/2

        edgesSignalMerge <- merge(edgesGraph, signalCoordMerge, by.x = "c",
            by.y = "V1", all.x = TRUE)
        edgesSignalMerge <- merge(edgesSignalMerge, signalCoordMerge,
            by.x = "V1", by.y = "V1", all.x = TRUE)
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
            print(paste0(
                "There are ",nrow(naTotal)," nodes without expression value,
                see log in path: ",
                file.path(tempdir(),titleChart, "levi.log")))
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

        edgesSignalMerge <- edgesSignalMerge[,c(3,4,2,1,5,6,7,8,9,10,11,12)]

        edgesSignalMerge$V2 <- (edgesSignalMerge[,c(8)] +
        edgesSignalMerge[,c(12)])/2
        edgesSignalMerge$V3 <- (edgesSignalMerge[,c(7)] +
        edgesSignalMerge[,c(11)])/2
        ##############################################################################

        nnodes <- nrow(nodes)
        nedges <- nrow(edges)
        numberAll<-nnodes+nedges
        coordAll <- rbind(nodesCoord[,c(2,3)], edgesNodesMerge[,c(9,10)])

        signalExpAll <- data.frame(V1 = c(signalCoordMerge[,c(5)],
        edgesSignalMerge[,c(13)]))
        signalCtrlAll <- data.frame(V1 = c(signalCoordMerge[,c(4)],
        edgesSignalMerge[,c(14)]))


        numberCoord <- numberAll
        coord <- as.matrix(coordAll)
        signalExp <- as.matrix(signalExpAll)

        if (baseTest == baseControl){
            signalCtrl <- matrix(data = 1, ncol = 1, nrow = nrow(signalCtrlAll))
        } else {
            signalCtrl <- as.matrix(signalCtrlAll)
        }
        signalExp[is.na(signalExp)] <- 0.5
        signalCtrl[is.na(signalCtrl)] <- 0.5
        SignalOut<-signalExp/(signalExp+signalCtrl)


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

        #Applies the calculation and takes the smallest value for coordinates
        #coord


        coordPiso <- SigCoordPiso(coord= coord,
        resolutionValue = resolutionValue,gamaValue= gamaValue,
        increase = increase,  contrastValue = contrastValue,
        zoomValue=zoomValue, numberCoord=numberCoord)


        matrix_resultado <- matrix_entrada(coordPiso= coordPiso,
        SignalOut= SignalOut,
        coord= coord, resolutionValue= resolutionValue,
        signalExp = signalExp, signalCtrl = signalCtrl, increase= increase,
        zoomValue= zoomValue, numberCoord= numberCoord)

        matrixIn <- as.matrix(matrix_resultado$m1)
        h <- matrix_resultado$m3

        h<-h-1


        matrixFinal <- matrix_saida(matrixIn =  matrixIn,
        resolutionValue = resolutionValue, gamaValue =  gamaValue,
        increase = increase, zoomValue = zoomValue, h = h,
        smoothValue = smoothValue)


        matrixOut <- matrixFinal$m1
        matrixOutExp <- matrixFinal$m2
        matrixOutCtrl <- matrixFinal$m3


        n<-resolutionValue
        ExpCtrl <- matrix(data = 0, ncol = n, nrow = n)
        exp <- matrix(data = 0, ncol = n, nrow = n)
        ctrl <- matrix(data = 0, ncol = n, nrow = n)

        b<-matrixOutExp[1,1]
        c<-matrixOutCtrl[1,1]

        for (i in 1:resolutionValue){
            for (j in 1:(resolutionValue-1)){
                if(b < matrixOutExp[i,(j+1)]){b<-matrixOutExp[i,(j+1)]}
                if(c < matrixOutCtrl[i,(j+1)]){c<-matrixOutCtrl[i,(j+1)]}
            }
        }

        for (i in seq_len(resolutionValue)) {
            k<-resolutionValue
            for (j in seq_len((resolutionValue))) {
                matrixOutExp[i,j]<-matrixOutExp[i,j]/b
                matrixOutCtrl[i,j]<-matrixOutCtrl[i,j]/c


                ExpCtrl[i,k]<-matrixOut[i,j]
                exp[i,k]<-matrixOutExp[i,j]
                ctrl[i,k]<-matrixOutCtrl[i,j]
                k<-k-1
            }}

        if (baseTest == baseControl){
            landgraph <- melt(exp, value.name = "z")
        } else {
            landgraph <- melt(ExpCtrl, value.name = "z")
        }

        landgraphFinal <- as.data.frame(landgraph[,c(1,2,3)])


        matrixSize <- sqrt(NROW(landgraphFinal))


        if (missing(contourLevi) || contourLevi == TRUE) {
            landgraphChart <-ggplot(data = landgraphFinal,
                aes(x = Var1, y = Var2))+
                geom_raster(aes(fill = z), interpolate = TRUE, hjust = 0.5,
                vjust = 0.5) +
                geom_contour(aes(z = z)) +
                scale_fill_gradientn(colours=colorSet,
                values=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
                breaks=seq(0,1,0.2), limits=c(0,1),
                guide = guide_colorbar(title="Expression Level",
                title.position = "right", title.hjust = 0.5,
                title.theme = element_text(angle = 270, size = 9),
                barwidth= 1,barheight = 10)) +
                theme_void() +
                ggtitle(titleChart) +
                theme(plot.title = element_text(
                margin = margin(t = 10, b = -10), hjust = 0.5,lineheight=.8,
                face="bold"),legend.margin=margin(0,0,0,-20)) +
                annotate("text", x = c(matrixSize+2,matrixSize+2),
                y = c(matrixSize*0.42,matrixSize*0.6),
                label = c("decrease", "increase"),size=3 , angle=90) +
                annotate("segment", x = matrixSize*1.07, xend = matrixSize*1.07,
                y = matrixSize*0.49, yend = matrixSize*0.35, colour = "black",
                size=0.2, alpha=0.6, arrow=arrow(type = "closed",
                length = unit(x = c(0.2), units = "cm"))) +
                annotate("segment", x = matrixSize*1.07, xend = matrixSize*1.07,
                y = matrixSize*0.51, yend = matrixSize*0.67, colour = "black",
                size=0.2,alpha=0.6, arrow=arrow(type = "closed",
                length = unit(x = c(0.2), units = "cm"))) +
                coord_fixed(ratio = 1)
        } else {
            landgraphChart <-ggplot(data = landgraphFinal,
                aes(x = Var1, y = Var2))+
                geom_raster(aes(fill = z), interpolate = TRUE,
                hjust = 0.5, vjust = 0.5) +
                scale_fill_gradientn(colours=colorSet,
                values=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
                breaks=seq(0,1,0.2), limits=c(0,1),
                guide = guide_colorbar(title="Expression Level",
                title.position = "right", title.hjust = 0.5,
                title.theme = element_text(angle = 270, size = 9), barwidth= 1,
                barheight = 10)) +
                theme_void() +
                ggtitle(titleChart) +
                theme(plot.title =
                element_text(margin = margin(t = 10, b = -10), hjust = 0.5,
                lineheight=.8, face="bold"), legend.margin=margin(0,0,0,-20)) +
                annotate("text", x = c(matrixSize+2,matrixSize+2),
                y = c(matrixSize*0.42,matrixSize*0.6),
                label = c("decrease", "increase"), size=3 , angle=90) +
                annotate("segment", x = matrixSize*1.07, xend = matrixSize*1.07,
                y = matrixSize*0.49, yend = matrixSize*0.35, colour = "black",
                size=0.2, alpha=0.6, arrow=arrow(type = "closed",
                length = unit(x = c(0.2), units = "cm"))) +
                annotate("segment", x = matrixSize*1.07, xend = matrixSize*1.07,
                y = matrixSize*0.51, yend = matrixSize*0.67, colour = "black",
                size=0.2, alpha=0.6, arrow=arrow(type = "closed",
                length = unit(x = c(0.2), units = "cm"))) +
                coord_fixed(ratio = 1)
        }

        landgraphChart + coord_fixed(ratio = 1)

        print(landgraphChart)
    }
}

