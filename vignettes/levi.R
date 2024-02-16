## ----eval=FALSE, fig.height=6, fig.width=6------------------------------------
#  library(levi)
#  LEVIui(browser=TRUE)  #Launch Levi to Browser.
#  LEVIui(browser=FALSE) #Launch Levi to R environment.
#  

## ----eval=TRUE, fig.height=6, fig.width=6-------------------------------------
library(levi)

template_network <- file.path(system.file(package="levi"),"extdata",
                                "medusa.dat", fsep = .Platform$file.sep)

template_expression <- file.path(system.file(package="levi"),
                                "extdata","expression.dat", 
                                fsep = .Platform$file.sep)

multicolor <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolnput = "ID", 
                readExpColumn=
                readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"), 
                contrastValueInput = 50, resolutionValueInput  = 50, 
                zoomValueInput = 50, smoothValueInput = 50, contourLevi = TRUE)

twocolors <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolnput = "ID", 
                readExpColumn=
                readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
                setcolor = "pink_green", contourLevi = FALSE)


## ----eval=TRUE,fig.height=6, fig.width=6--------------------------------------
library(levi)
base <- readExpColumn("TumorFormerSmoker-NormalFormerSmoker", 
                        "TumorNeverSmoker-TumorNeverSmoker")

template_network <- file.path(system.file(package="levi"),"extdata",
                                "medusa.dat", fsep = .Platform$file.sep)

template_expression <- file.path(system.file(package="levi"),
                                "extdata","expression.dat", 
                                fsep = .Platform$file.sep)

multicolor <- levi(networkCoordinatesInput = template_network,
                    expressionInput = template_expression, 
                    fileTypeInput = "dat",
                    geneSymbolnput = "ID", readExpColumn= base, 
                    contrastValueInput = 50, resolutionValueInput  = 50, 
                    zoomValueInput = 50, smoothValueInput = 50, 
                    contourLevi = FALSE)

twocolors <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolnput = "ID", 
                readExpColumn= base,
                setcolor = "pink_green", contourLevi = FALSE)


## ---- eval=TRUE, label='Session information', , echo=FALSE--------------------
sessionInfo()

