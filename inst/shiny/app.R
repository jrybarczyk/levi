suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(igraph))
suppressMessages(require(colorspace))
suppressMessages(require(grid))
suppressMessages(require(shiny))
suppressMessages(require(shinyjs))
suppressMessages(require(shinydashboard))
suppressMessages(require(dplyr))
suppressMessages(require(DT))
suppressMessages(require(xml2))

# plotly is a suggested dependency: without it the interface still works, just
# without the 3D surface. Both the control and the output are left out rather
# than failing when the app starts.
has_plotly <- requireNamespace("plotly", quietly = TRUE)

if(getRversion() >= "3.4.0") utils::globalVariables(c(
    "V1", "X", "Y", "Expression", "Gene",
    "renderText", "scoreTable", "LandscapeScore", "Rank",
    "pval_over", "pval_under", "perm_side", "sig_level",
    "regions", "inference_unit", "PeakRow", "PeakCol", "Label", "xend", "yend",
    "Var1", "Var2", "pval", "peakTable",
    "x", "y", "label", "Type",
    "signal_mode", "logfc_k"))

# The landscape core lives in src/landscape_gauss.cpp and is compiled with
# the package. The interface calls the very same implementation script mode
# does, instead of recompiling private copies with cppFunction() each session.

levi_shiny <- function(expression, fileType, networkCoord,
    networkInterac, geneSymbol, baseTest, baseControl,
    alphaValue, betaValue, backValue, smoothValue, expressionLog,
    n_perm = 0L, perm_side = "both", sig_level = 0.05,
    signal_mode = "ratio", logfc_k = 1,
    inference_unit = c("region", "cell")) {
    inference_unit <- match.arg(inference_unit)

    if (missing(alphaValue)) alphaValue <- 50
    if (missing(betaValue)) betaValue <- 50
    if (missing(backValue)) backValue <- 50
    if (missing(smoothValue)) smoothValue <- 50
    if (baseControl == " ") baseControl <- baseTest
    has_progress <- !is.null(shiny::getDefaultReactiveDomain())
    res <- levi:::levi_function(
        expression, fileType, networkInterac, networkCoord, geneSymbol,
        levi::readExpColumn(paste(baseTest, baseControl, sep = "-")),
        alphaValue, betaValue, backValue, smoothValue, expressionLog,
        FALSE, "default", n_perm = n_perm, sig_level = sig_level,
        perm_side = perm_side, signal_mode = signal_mode, logfc_k = logfc_k,
        # "region" (default) outlines significant areas; "cell" (legacy)
        # draws contours from cell-wise p-values.
        inference_unit = inference_unit,
        .parsed_network = list(nodes = networkInterac, edges = networkCoord),
        .draw = FALSE,
        .progress = if (has_progress) function(i, n)
            shiny::incProgress(1 / n, detail = sprintf("Permutation %d / %d", i, n)))
    grid <- res$metadata$grid
    n <- grid$resolution
    node_idx <- levi:::nearest_node_grid(res$metadata$node_coordinates,
                                         n, grid$zoom, grid$increase)
    labels <- matrix(res$metadata$nodes$V1[node_idx], nrow = n)
    surface <- matrix(res$landscape$z, nrow = n)
    landscape <- data.frame(X = res$landscape$Var1, Y = res$landscape$Var2,
                            Expression = round(res$landscape$z, 2),
                            Gene = as.vector(labels[, rev(seq_len(n))]))
    list(landscape = landscape, surface = surface, scores = res$scores,
         peaks = res$peaks, pval_over = res$pvalues$over,
         pval_under = res$pvalues$under, regions = res$regions,
         inference_unit = inference_unit, perm_side = perm_side,
         sig_level = sig_level, metadata = res$metadata)
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
    shinyjs::useShinyjs(),
    br(),
    sidebarLayout(
        sidebarPanel(
            sidebarMenu(id = "tab"),
            tabsetPanel(id = "tabset_id", selected = "t1",
                tabPanel("File", value = "t1",
                    br(),
                    selectInput("fileType", "Network input type:",
                                c("Medusa (DAT)" = "dat",
                                  "RedeR (DYN)"  = "dyn",
                                  "Pajek (NET)"  = "net",
                                  "STRING / STITCH" = "stg")),
                    br(),
                    uiOutput("out2"),
                    tags$hr(),
                    radioButtons("exprSource",
                                 "Expression data source:",
                                 c("Upload file (.tsv/.dat)" = "file",
                                   "Bioconductor object (.rds)" = "bioc"),
                                 selected = "file", inline = TRUE),
                    conditionalPanel("input.exprSource == 'file'",
                        fileInput("file",
                                  "Upload the expression file:"),
                        checkboxInput("log",
                                      "Expression values in log scale",
                                      FALSE),
                        radioButtons("fields", "Selected fields:",
                                     c("Two Samples" = "twofields",
                                       "One Sample"  = "onefield"),
                                     selected = "twofields", inline = TRUE),
                        uiOutput("out1")
                    ),
                    conditionalPanel("input.exprSource == 'bioc'",
                        fileInput("bioc_rds",
                                  paste0("Upload .rds file ",
                                         "(SummarizedExperiment, ",
                                         "SingleCellExperiment, ",
                                         "ExpressionSet):")),
                        uiOutput("bioc_ui")
                    ),
                    radioButtons(inputId = "palette", label = "Colour scale",
                                 choices = c("Multicolor" = "multi",
                                             "Two colors" = "twoc"),
                                 selected = "multi"),
                    conditionalPanel("input.palette == 'twoc'",
                        selectInput(inputId = "setcolor",
                                    label   = "Representation options",
                                    choices = c(
                                        "Purple & Pink"       = "pp",
                                        "Green & Blue"        = "gb",
                                        "Blue & Yellow"       = "by",
                                        "Pink & Green"        = "pg",
                                        "Orange & Purple"     = "op",
                                        "Blue-green & Marine" = "bm"))),
                    checkboxInput("contour", "Chart with contour", FALSE),
                    if (has_plotly)
                        checkboxInput("plot3d", "3D surface", FALSE),
                    tags$hr(),
                    textInput("geneSearch", "Highlight genes (names, comma-separated):",
                              value = "", placeholder = "e.g. OGDHL, IDH1"),
                    selectInput("highlightColor", "Highlight colour",
                                choices = c("Yellow" = "yellow", "White" = "white",
                                            "Black" = "black", "Red" = "red",
                                            "Cyan" = "cyan", "Magenta" = "magenta",
                                            "Orange" = "orange"),
                                selected = "yellow"),
                    checkboxInput("showPeakLabels", "Label peaks on map", value = TRUE),
                    tags$hr(),
                    actionButton("action", "Run", style = "float:right"),
                    br()
                ),
                tabPanel("Settings", value = "t2",
                    br(),
                    helpText("Re-run to apply changes"),
                    sliderInput("contrast", "Contrast:",
                                min = 0, max = 100, value = 50),
                    sliderInput("size",     "Resolution:",
                                min = 1, max = 100, value = 50),
                    sliderInput("smooth",   "Smoothing:",
                                min = 0, max = 100, value = 50),
                    sliderInput("zoom",     "Zoom:",
                                min = 0, max = 100, value = 50),
                    tags$hr(),
                    helpText("Signal transformation"),
                    selectInput("signal_mode", "Signal mode:",
                                choices = c(
                                    "Ratio Test/(Test+Control) — counts, TPM, FPKM" = "ratio",
                                    "Sigmoid logFC — RMA, VST, log2-proteomics, scRNA-seq FC" = "logfc",
                                    "Z-score (pnorm) — multi-condition / heterogeneous scale" = "zscore"
                                )),
                    numericInput("logfc_k",
                                 "logFC steepness k (logfc mode only):",
                                 value = 1, min = 0.1, max = 5, step = 0.1),
                    helpText(
                        "k < 1: smoother gradient (large FC datasets, e.g. scRNA-seq).",
                        br(),
                        "k > 1: sharper gradient (tight FC, e.g. microarray ±1).",
                        br(),
                        "For single-column logFC: use readExpColumn('FC-FC')."
                    ),
                    tags$hr(),
                    helpText("Significance test (increases run time)"),
                    sliderInput("nperm",
                                "Permutations (0 = disabled):",
                                min = 0, max = 500, value = 0, step = 50),
                    selectInput("inference_unit", "Test unit:",
                                choices = c(
                                    "Regions (maximum regional mass)" = "region",
                                    "Grid cells (legacy, BY-adjusted)" = "cell")),
                    selectInput("perm_side", "Test direction:",
                                choices = c(
                                    "Both (over & under)"  = "both",
                                    "Over-expression only" = "over",
                                    "Under-expression only" = "under")),
                    numericInput("sig_level_ui", "Significance level:",
                                 value = 0.05, min = 0.001,
                                 max = 0.2, step = 0.005),
                    helpText(
                        "Regions: white outline and label = area whose mass ",
                        "beats the null maximum (p <= level).", br(),
                        "Cells: dashed contour = over-expressed, ",
                        "dotted = under-expressed (p < level).")
                )
            )
        ),
        mainPanel(
            tags$head(tags$style(type = "text/css",
                "#loadmessage { position:fixed; top:0px; left:0px; width:100%;
                 padding:5px 0px 5px 0px; text-align:center;
                 font-weight:bold; font-size:100%; color:#000000;
                 background-color:#81efea; z-index:105; }")),
            conditionalPanel(
                condition = "$('html').hasClass('shiny-busy')",
                tags$div("Loading...", id = "loadmessage")),

            # The 2D map and the 3D surface live in two tabs so either can
            # take the full width; everything below (selected area, tables)
            # is shared and unchanged. The brush belongs to the 2D map: the
            # surface is for reading the relief, the map is what you select
            # regions on.
            tabsetPanel(id = "viewTabs", selected = "view2d",
                tabPanel("2D landscape", value = "view2d",
                    br(),
                    plotOutput(outputId = "graph",
                               brush = brushOpts(id = "plotBrush",
                                                 resetOnNew = TRUE)),
                    # Format selector and "Download Plot" for the 2D map;
                    # the 3D tab has its own HTML download. The button
                    # floats left, so the expression area (sum over the
                    # brushed cells of this map) is cleared below it.
                    uiOutput("out3"),
                    div(style = "clear:both"),
                    h4(textOutput("expArea"))),
                tabPanel("3D surface", value = "view3d",
                    br(),
                    if (has_plotly) conditionalPanel(
                        condition = "!input.plot3d",
                        helpText("Tick \"3D surface\" in the File tab to ",
                                 "build the surface.")),
                    if (has_plotly) conditionalPanel(
                        condition = "input.plot3d",
                        plotly::plotlyOutput(outputId = "graph3d",
                                             height = "520px"),
                        helpText("Rotate the surface, then use the camera ",
                                 "icon on the plot toolbar to save a PNG of ",
                                 "the current view, or download it as an ",
                                 "interactive HTML file that keeps this ",
                                 "view."),
                        downloadButton("download3D", "Download 3D (HTML)"),
                        helpText("Camera of the current view, ready to ",
                                 "paste into leviSave3D(camera = ...):"),
                        verbatimTextOutput("camera3dCode")),
                    if (!has_plotly) helpText(
                        "Install the plotly package to enable the 3D ",
                        "surface."))
            ),

            br(),

            # The tables in tabs of their own. "Genes" lists what was brushed
            # on the 2D map, so it is hidden while the 3D tab is in front;
            # the other three describe the landscape and stay available in
            # both views.
            tabsetPanel(id = "tableTabs", selected = "tabGenes",
                tabPanel("Genes", value = "tabGenes",
                    br(),
                    helpText("Genes under the area selected on the 2D map."),
                    dataTableOutput(outputId = "landdatatable")),
                tabPanel("Node scores", value = "tabScores",
                    br(),
                    uiOutput("scoreTableHeader"),
                    uiOutput("scoreDownloadUI"),
                    DT::dataTableOutput(outputId = "scoreTableOutput")),
                tabPanel("Peaks and valleys", value = "tabPeaks",
                    br(),
                    uiOutput("peakTableHeader"),
                    uiOutput("peakDownloadUI"),
                    DT::dataTableOutput(outputId = "peakTableOutput")),
                tabPanel("Regions", value = "tabRegions",
                    br(),
                    uiOutput("regionTableHeader"),
                    uiOutput("regionDownloadUI"),
                    DT::dataTableOutput(outputId = "regionTableOutput"))
            )
        )
    )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
    options(shiny.maxRequestSize = 30 * 1024^2)
    set.seed(123)

    # ── Reactive state (all fields initialised here) ──────────────────────
    v <- reactiveValues(
        expression     = NULL,
        surface        = NULL,
        networkCoord   = NULL,
        networkInterac = NULL,
        fileType       = NULL,
        func_ne_return = NULL,
        graficof       = NULL,
        colorSet       = c("#4e5052", "#0000b4", "#FF0000"),
        multicolor     = c("#180052","#0c0083","#0000b4","#0000e4","#0010ff",
                           "#0041ff","#0072ff","#00A3FF","#00D4FF","#00FF49",
                           "#5AFF00","#FFE400","#FFC400","#FFA300","#FF8300",
                           "#FF6200","#FF4100","#FF2100","#FF0000","#E40000"),
        matrixSize     = NULL,
        landDataTableTmp = setNames(
            data.frame(matrix(ncol = 4, nrow = 0)),
            c("X", "Y", "Expression", "Gene")),
        scoreTable     = data.frame(Gene = character(0),
                                    LandscapeScore = numeric(0),
                                    Rank = integer(0)),
        peakTable      = data.frame(),
        pval_over      = NULL,
        pval_under     = NULL,
        regions        = NULL,
        inference_unit = "region",
        perm_side_res  = NULL,
        sig_level_res  = NULL,
        currentPlot    = NULL
    )

    # ── Reactive: read expression file for column names ───────────────────
    baseSelect <- reactive({
        req(input$file)
        # A warning here does not mean the file is unreadable. The commonest
        # one, "incomplete final line found", comes from a file whose last
        # line has no newline -- which plenty of editors and pipelines produce,
        # and which read.table parses correctly anyway. Catching it with
        # tryCatch discarded the data and reported "Incorrect file format",
        # about the one thing that was fine.
        withCallingHandlers(
            tryCatch(
                read.table(file = input$file$datapath, header = TRUE,
                           sep = "\t"),
                error = function(e) {
                    showNotification("Incorrect file format",
                                     type = "error", closeButton = FALSE,
                                     duration = 5)
                    NULL
                }),
            warning = function(w) {
                showNotification(conditionMessage(w), type = "warning",
                                 closeButton = FALSE, duration = 5)
                invokeRestart("muffleWarning")
            })
    })

    # ── Reactive: load Bioconductor object from uploaded RDS ──────────────
    bioc_obj <- reactive({
        req(input$bioc_rds, input$exprSource == "bioc")
        tryCatch(
            readRDS(input$bioc_rds$datapath),
            error = function(e) {
                showNotification(
                    paste("Error loading RDS:", conditionMessage(e)),
                    type = "error", duration = 8)
                NULL
            }
        )
    })

    # ── Dynamic UI: assay / condition selectors for Bioc object ──────────
    output$bioc_ui <- renderUI({
        req(bioc_obj())
        obj <- bioc_obj()
        cls <- class(obj)[1]

        if (methods::is(obj, "SummarizedExperiment")) {
            cd           <- as.data.frame(
                SummarizedExperiment::colData(obj))
            assays_avail <- SummarizedExperiment::assayNames(obj)
            tagList(
                tags$p(tags$b("Detected: "), cls),
                selectInput("bioc_assay", "Assay:",
                            choices  = assays_avail,
                            selected = assays_avail[1]),
                selectInput("bioc_cond_col",
                            "Condition column (colData):",
                            choices = colnames(cd)),
                uiOutput("bioc_levels_ui"),
                checkboxInput("bioc_log",
                              "Log2-transform counts (log2(x + 1))",
                              value = FALSE)
            )
        } else if (methods::is(obj, "ExpressionSet")) {
            if (!requireNamespace("Biobase", quietly = TRUE)) {
                return(tags$p("Biobase package required for ExpressionSet.",
                              "Install: BiocManager::install('Biobase')"))
            }
            pd <- Biobase::pData(obj)
            tagList(
                tags$p(tags$b("Detected: "), "ExpressionSet"),
                selectInput("bioc_cond_col",
                            "Condition column (pData):",
                            choices = colnames(pd)),
                uiOutput("bioc_levels_ui")
            )
        } else {
            tags$p(paste0(
                "Unsupported class: ", cls, ". ",
                "Supported: SummarizedExperiment, ",
                "SingleCellExperiment, ExpressionSet."))
        }
    })

    output$bioc_levels_ui <- renderUI({
        req(bioc_obj(), input$bioc_cond_col)
        obj <- bioc_obj()
        col <- input$bioc_cond_col

        levs <- if (methods::is(obj, "SummarizedExperiment")) {
            cd <- as.data.frame(SummarizedExperiment::colData(obj))
            sort(unique(as.character(cd[[col]])))
        } else if (methods::is(obj, "ExpressionSet")) {
            pd <- Biobase::pData(obj)
            sort(unique(as.character(pd[[col]])))
        } else {
            return(NULL)
        }

        tagList(
            selectInput("bioc_test_level", "Test condition:",
                        choices  = levs,
                        selected = levs[1]),
            selectInput("bioc_ctrl_level", "Control condition:",
                        choices  = levs,
                        selected = if (length(levs) > 1) levs[2] else levs[1])
        )
    })

    output$out1 <- renderUI({
        tabPanel("", value = "t1",
            selectInput("geneSymbol", "Select ID field",
                        c(names(baseSelect())), selected = NULL),
            selectInput("baseTest", "Select test field",
                        c(names(baseSelect())), selected = NULL),
            if (input$fields == "onefield") {
                selectInput("baseControl", "Select control field",
                            c(" "), selected = NULL)
            } else {
                selectInput("baseControl", "Select control field",
                            c(names(baseSelect())), selected = NULL)
            },
            helpText("Chart colors can be changed without having to re-run"))
    })

    observeEvent(input$fileType, {
        if (input$fileType == "stg") {
            output$out2 <- renderUI({
                tabPanel("", value = "t1",
                    fileInput("file2", "Upload the coordinates file:"),
                    fileInput("file3", "Upload the interactions file:"))
            })
        } else {
            output$out2 <- renderUI({
                tabPanel("", value = "t1",
                    fileInput("file2", "Upload the network file:"))
            })
        }
    })

    # ── Single observeEvent: parse files then compute ─────────────────────
    observeEvent(input$action, {
        ft        <- input$fileType
        bioc_mode <- isTRUE(input$exprSource == "bioc")

        # Column symbols differ between modes
        gs_sym <- if (bioc_mode) "GeneID"  else input$geneSymbol
        bt_sym <- if (bioc_mode) "Test"    else input$baseTest
        bc_sym <- if (bioc_mode) "Control" else input$baseControl

        # — Bioconductor object → expression data.frame —
        if (bioc_mode) {
            obj <- bioc_obj()
            if (is.null(obj)) {
                showNotification("No Bioconductor object loaded.",
                                 type = "error")
                return()
            }
            bioc_expr <- tryCatch({
                if (methods::is(obj, "SummarizedExperiment")) {
                    leviFromSE(obj,
                        assay_name    = input$bioc_assay,
                        condition_col = input$bioc_cond_col,
                        test_level    = input$bioc_test_level,
                        ctrl_level    = input$bioc_ctrl_level,
                        gene_col      = "GeneID",
                        log_transform = isTRUE(input$bioc_log))
                } else if (methods::is(obj, "ExpressionSet")) {
                    leviFromExpressionSet(obj,
                        condition_col = input$bioc_cond_col,
                        test_level    = input$bioc_test_level,
                        ctrl_level    = input$bioc_ctrl_level,
                        gene_col      = "GeneID")
                } else {
                    stop("Unsupported class: ", class(obj)[1])
                }
            }, error = function(e) {
                showNotification(conditionMessage(e), type = "error",
                                 duration = 10)
                NULL
            })
            if (is.null(bioc_expr)) return()
        }

        # — file parsing (expression file OR bioc result; network always file) —
        parsed <- tryCatch({
            if (bioc_mode) {
                expr <- bioc_expr
            } else {
                data1 <- input$file
                req(data1)
                expr <- read.delim(file = data1$datapath, header = TRUE,
                                   sep = "\t", quote = "")
            }

            # One parser per format, shared with script mode. The interface
            # used to carry its own copy of all four; sharing them is what
            # keeps a fix from landing on one path only.
            #
            # .parseNetwork() names its outputs after what they are, while the
            # interface keeps the older names: netCoord holds the edges and
            # netInterac the nodes.
            data2 <- input$file2; req(data2)
            edgeFile <- if (ft == "stg") {
                data3 <- input$file3; req(data3); data3$datapath
            } else NA

            parsedNet  <- levi:::.parseNetwork(data2$datapath, edgeFile, ft)
            netCoord   <- parsedNet$edges
            netInterac <- parsedNet$nodes

            list(expr = expr, netCoord = netCoord, netInterac = netInterac)
        },
        error = function(e) {
            showNotification("Incorrect file format", type = "error",
                             closeButton = FALSE, duration = 5)
            NULL
        })

        if (is.null(parsed)) return()

        v$fileType     <- ft
        v$expression   <- parsed$expr
        v$networkCoord <- parsed$netCoord
        v$networkInterac <- parsed$netInterac

        # — computation —
        #
        # Warnings are shown and the run continues. They used to be caught by
        # tryCatch, which unwinds the stack: any warning killed the whole
        # computation. That is how a deprecation notice from igraph 2.0 left
        # the interface unable to produce anything at all, reported to the user
        # as "Incorrect file format".
        withCallingHandlers(
        tryCatch({
          withProgress(message = "Computing landscape...", value = 0, {
            shiny::incProgress(0.05, detail = "Initialising")
            leviResult <- levi_shiny(
                expression   = v$expression,
                fileType     = v$fileType,
                networkCoord = v$networkCoord,
                networkInterac = v$networkInterac,
                geneSymbol   = gs_sym,
                baseTest     = bt_sym,
                baseControl  = bc_sym,
                alphaValue   = input$contrast,
                betaValue    = input$zoom,
                backValue    = input$size,
                smoothValue  = input$smooth,
                expressionLog = if (bioc_mode) FALSE else input$log,
                n_perm      = as.integer(input$nperm),
                perm_side   = input$perm_side,
                sig_level   = input$sig_level_ui,
                signal_mode = input$signal_mode,
                logfc_k     = input$logfc_k,
                inference_unit = input$inference_unit)

            v$func_ne_return  <- leviResult$landscape
            v$surface         <- leviResult$surface
            v$signal_meaning  <- leviResult$metadata$meaning
            v$landDataTableTmp <- leviResult$landscape
            v$scoreTable      <- leviResult$scores
            v$peakTable       <- leviResult$peaks
            v$pval_over       <- leviResult$pval_over
            v$pval_under      <- leviResult$pval_under
            v$regions         <- leviResult$regions
            v$inference_unit  <- leviResult$inference_unit
            v$perm_side_res   <- leviResult$perm_side
            v$sig_level_res   <- leviResult$sig_level
            v$metadata        <- leviResult$metadata
            shiny::incProgress(0.05, detail = "Done")
          }) # end withProgress

            output$out3 <- renderUI({
                tabPanel("", value = "t1",
                    selectInput("plotType", "",
                                c("TIFF" = "tiff", "BMP" = "bmp",
                                  "JPEG" = "jpeg", "PNG" = "png"),
                                width = "135px"),
                    downloadButton("downloadPlot", "Download Plot",
                                   style = "float:left"))
            })
        },
        error = function(e) {
            showNotification(conditionMessage(e), type = "error",
                             closeButton = FALSE, duration = 5)
        }),
        warning = function(w) {
            showNotification(conditionMessage(w), type = "warning",
                             closeButton = FALSE, duration = 8)
            invokeRestart("muffleWarning")
        })
    })

    # ── 3D surface ────────────────────────────────────────────────────────
    # Built by the same routine script mode uses, so the two views agree --
    # including the significance boundary traced onto the relief when a
    # permutation test ran.
    surface3d <- reactive({
        req(v$surface)
        colours <- if (identical(input$palette, "multi")) v$multicolor else v$colorSet
        pvals <- if (!is.null(v$pval_over))
            list(over = v$pval_over, under = v$pval_under) else NULL

        levi:::.buildSurface3D(
            zmat       = v$surface,
            colours    = colours,
            titleChart = NULL,
            pvals      = pvals,
            i          = if (!is.null(pvals)) seq_len(nrow(v$surface)),
            sig_level  = v$sig_level_res %||% 0.05,
            perm_side  = v$perm_side_res %||% "both")
    })

    if (has_plotly) output$graph3d <- plotly::renderPlotly({
        req(isTRUE(input$plot3d))
        # The relayout event carries the camera every time the user rotates
        # or zooms the surface; the download below reuses it so the saved
        # file opens at the view the user chose.
        plotly::event_register(surface3d(), "plotly_relayout")
    })

    # Ticking "3D surface" brings the 3D tab to the front; unticking goes
    # back to the map.
    observeEvent(input$plot3d, {
        updateTabsetPanel(session, "viewTabs",
                          selected = if (isTRUE(input$plot3d)) "view3d"
                                     else "view2d")
    }, ignoreInit = TRUE)

    # The "Genes" table only makes sense next to the map it is brushed on.
    observeEvent(input$viewTabs, {
        if (identical(input$viewTabs, "view3d")) {
            if (identical(input$tableTabs, "tabGenes"))
                updateTabsetPanel(session, "tableTabs", selected = "tabScores")
            hideTab("tableTabs", "tabGenes")
        } else {
            showTab("tableTabs", "tabGenes")
        }
    })

    if (has_plotly) observeEvent(plotly::event_data("plotly_relayout"), {
        ev <- plotly::event_data("plotly_relayout")
        if (!is.null(ev[["scene.camera"]])) v$camera3d <- ev[["scene.camera"]]
    }, ignoreNULL = TRUE)

    # The camera as R code. plotly's default view is shown until the user
    # rotates the surface, so the box is never empty.
    output$camera3dCode <- renderText({
        cam <- v$camera3d %||% list(eye = list(x = 1.25, y = 1.25, z = 1.25))
        fmt <- function(v) paste0("list(x = ", signif(v$x, 3), ", y = ",
                                  signif(v$y, 3), ", z = ", signif(v$z, 3), ")")
        parts <- c(eye = "eye", center = "center", up = "up")
        parts <- parts[names(parts) %in% names(cam)]
        paste0("camera <- list(\n",
               paste0("  ", parts, " = ", vapply(cam[parts], fmt, ""),
                      collapse = ",\n"),
               "\n)")
    })

    output$download3D <- downloadHandler(
        filename = function() paste0("levi_surface3D_", Sys.Date(), ".html"),
        content  = function(file) {
            fig <- surface3d()
            if (!is.null(v$camera3d))
                fig <- plotly::layout(fig, scene = list(camera = v$camera3d))
            htmlwidgets::saveWidget(fig, file,
                selfcontained = requireNamespace("rmarkdown", quietly = TRUE) &&
                    rmarkdown::pandoc_available())
        })

    # ── Plot (declared once; reacts to v$* changes) ───────────────────────
    output$graph <- renderPlot({
        req(v$func_ne_return)

        # colour palettes
        # Palettes come from the package rather than being written out again
        # here. They had already drifted: "Pink & Green" was #f757ca/#04ff00 in
        # the interface against #e854d9/#90db56 in script mode, so the same
        # option produced different figures depending on how levi was called.
        palette_of <- c(pp = "purple_pink", gb = "green_blue",
                        by = "blue_yellow", pg = "pink_green",
                        op = "orange_purple", bm = "green_marine")
        if (!is.null(input$setcolor) && input$setcolor %in% names(palette_of))
            v$colorSet <- levi:::.colorSet(palette_of[[input$setcolor]])

        n <- v$matrixSize <- sqrt(NROW(v$func_ne_return))

        # Same drawing script mode produces. The interface passes its own
        # column names and its own colour vector -- "Multicolor" is not one of
        # the named palettes -- and no title, since the comparison is shown in
        # the page. Everything else, including the scale arrows, now comes from
        # one place: they used to sit at different heights on each side.
        base_plot <- levi:::.buildLandscapeChart(
            landgraphFinal = v$func_ne_return,
            setcolor       = "default",
            titleChart     = NULL,
            contourLevi    = isTRUE(input$contour),
            colours        = if (identical(input$palette, "multi")) v$multicolor else v$colorSet,
            cols           = c("X", "Y", "Expression")) +
            ggplot2::labs(caption = v$signal_meaning)

        v$graficof <- base_plot

        # Regional inference: outline the areas that pass sig_level and
        # label every area with its p-value, exactly as script mode does.
        regional <- identical(v$inference_unit, "region") &&
            !is.null(v$regions) && !is.null(v$regions$summary$PSpatial)
        if (regional) {
            sl   <- v$sig_level_res %||% 0.05
            side <- v$perm_side_res %||% "both"
            summ <- v$regions$summary
            selected <- summ$Region[summ$PSpatial <= sl &
                (side == "both" | summ$Direction == side)]
            boundary <- levi:::.regionBoundaries(v$regions, n, selected)
            if (nrow(boundary)) base_plot <- base_plot +
                ggplot2::geom_segment(data = boundary,
                    aes(x = x, y = y, xend = xend, yend = yend),
                    colour = "white", linewidth = 0.8, inherit.aes = FALSE)
            summ <- summ[summ$Region %in% selected, , drop = FALSE]
            if (nrow(summ)) {
                summ$Label <- sprintf("%s\np = %.3f", summ$Region, summ$PSpatial)
                base_plot <- base_plot + ggplot2::geom_label(data = summ,
                    aes(x = PeakRow, y = n + 1 - PeakCol, label = Label),
                    size = 3, inherit.aes = FALSE)
            }
        }

        # Significance contours (legacy cell mode)
        if (!regional && (!is.null(v$pval_over) || !is.null(v$pval_under))) {
            i_seq <- seq_len(n)
            sl    <- v$sig_level_res %||% 0.05
            side  <- v$perm_side_res %||% "both"

            if (!is.null(v$pval_over) && side %in% c("both", "over"))
                base_plot <- levi:::.significanceContour(base_plot,
                    v$pval_over, i_seq, sl, "dashed", "higher score")
            if (!is.null(v$pval_under) && side %in% c("both", "under"))
                base_plot <- levi:::.significanceContour(base_plot,
                    v$pval_under, i_seq, sl, "dotted", "lower score")
        }

        # ── Peak labels ───────────────────────────────────────────────────
        if (isTRUE(input$showPeakLabels) && !regional &&
            !is.null(v$peakTable) && nrow(v$peakTable) > 0) {
            label_df <- data.frame(
                x     = v$peakTable$MatrixRow,
                y     = n + 1 - v$peakTable$MatrixCol,
                label = v$peakTable$NearestGene,
                Type  = v$peakTable$Type,
                stringsAsFactors = FALSE)
            if (requireNamespace("ggrepel", quietly = TRUE)) {
                base_plot <- base_plot +
                    ggrepel::geom_text_repel(
                        data = label_df,
                        aes(x = x, y = y, label = label),
                        colour = "white", size = 2.5, fontface = "bold",
                        box.padding = 0.3, max.overlaps = 20,
                        inherit.aes = FALSE)
            } else {
                base_plot <- base_plot +
                    geom_text(data = label_df,
                              aes(x = x, y = y, label = label),
                              colour = "white", size = 2.5, fontface = "bold",
                              inherit.aes = FALSE)
            }
        }

        # ── Gene search highlight ─────────────────────────────────────────
        # One circle at the node's own grid cell. The landscape table maps
        # every grid cell to its nearest node, so matching on it painted the
        # whole Voronoi region of the gene instead of its position.
        # Several genes may be given, separated by comma, semicolon or space.
        gene_query <- unique(trimws(strsplit(input$geneSearch, "[,;[:space:]]+")[[1]]))
        gene_query <- gene_query[nzchar(gene_query)]
        if (length(gene_query) > 0 && !is.null(v$metadata)) {
            md    <- v$metadata
            names <- as.character(md$nodes[[1]])
            idx   <- which(toupper(names) %in% toupper(gene_query))
            missing <- gene_query[!toupper(gene_query) %in% toupper(names)]
            if (length(missing) > 0)
                showNotification(paste0("Gene(s) not found: ",
                                        paste(missing, collapse = ", ")),
                                 type = "warning", duration = 3)
            if (length(idx) > 0) {
                coords <- md$node_coordinates[idx, , drop = FALSE]
                gx <- pmin(pmax(round((coords[, 1] - md$grid$zoom) /
                                          md$grid$increase) + 1L, 1L), n)
                gy <- pmin(pmax(round((coords[, 2] - md$grid$zoom) /
                                          md$grid$increase) + 1L, 1L), n)
                hit <- data.frame(x = gx, y = n + 1L - gy,
                                  label = names[idx])
                hl_col <- input$highlightColor %||% "yellow"
                base_plot <- base_plot +
                    geom_point(data = hit, aes(x = x, y = y),
                               colour = hl_col, fill = NA,
                               size = 6, shape = 21, stroke = 1.5,
                               inherit.aes = FALSE) +
                    geom_text(data = hit, aes(x = x, y = y, label = label),
                              colour = hl_col, size = 3, fontface = "bold",
                              vjust = -1.2, inherit.aes = FALSE)
            }
        }

        v$currentPlot <- base_plot
        base_plot

    }, height = "auto")

    # ── Brush selection table ─────────────────────────────────────────────
    output$landdatatable <- DT::renderDataTable({
        # Same reasoning as the file read above: a warning from the selection
        # is worth showing, but it does not mean there is nothing selected.
        withCallingHandlers(
            tryCatch(
                brushedPoints(v$landDataTableTmp, brush = input$plotBrush) %>%
                    select(Gene) %>% filter(Gene != "") %>% unique(),
                error = function(e) {
                    showNotification("Error while selecting chart",
                                     type = "error", closeButton = FALSE,
                                     duration = 5)
                    NULL
                }),
            warning = function(w) {
                showNotification(conditionMessage(w), type = "warning",
                                 closeButton = FALSE, duration = 5)
                invokeRestart("muffleWarning")
            })
    }, rownames = FALSE)

    # ── Score table ───────────────────────────────────────────────────────
    output$scoreTableHeader <- renderUI({
        req(nrow(v$scoreTable) > 0)
        h4("Node Landscape Scores (ranked)")
    })

    output$scoreDownloadUI <- renderUI({
        req(nrow(v$scoreTable) > 0)
        downloadButton("downloadScores", "Download CSV",
                       style = "margin-bottom:6px")
    })

    output$scoreTableOutput <- DT::renderDataTable({
        req(nrow(v$scoreTable) > 0)
        DT::datatable(
            v$scoreTable,
            rownames = FALSE,
            options  = list(pageLength = 10, order = list(list(2, "asc"))),
            caption  = "Score 0-1: 1 = highly over-expressed, 0 = highly under-expressed")
    }, server = TRUE)

    output$downloadScores <- downloadHandler(
        filename = function() paste0("levi_scores_", Sys.Date(), ".csv"),
        content  = function(file) write.csv(v$scoreTable, file, row.names = FALSE))

    # ── Peak / valley table ───────────────────────────────────────────────
    output$peakTableHeader <- renderUI({
        req(!is.null(v$peakTable) && nrow(v$peakTable) > 0)
        h4("Detected Peaks and Valleys")
    })

    output$peakDownloadUI <- renderUI({
        req(!is.null(v$peakTable) && nrow(v$peakTable) > 0)
        downloadButton("downloadPeaks", "Download CSV",
                       style = "margin-bottom:6px")
    })

    output$peakTableOutput <- DT::renderDataTable({
        req(!is.null(v$peakTable) && nrow(v$peakTable) > 0)
        DT::datatable(
            v$peakTable,
            rownames = FALSE,
            options  = list(pageLength = 10),
            caption  = "Peaks = local over-expression maxima; Valleys = local minima")
    }, server = TRUE)

    output$downloadPeaks <- downloadHandler(
        filename = function() paste0("levi_peaks_", Sys.Date(), ".csv"),
        content  = function(file) write.csv(v$peakTable, file, row.names = FALSE))

    # ── Regions table ─────────────────────────────────────────────────────
    region_summary <- reactive({
        req(!is.null(v$regions), nrow(v$regions$summary) > 0)
        s <- v$regions$summary
        keep <- intersect(c("Region", "Direction", "Cells", "Area", "Mass",
                            "PeakScore", "PSpatial", "Significant"), names(s))
        s <- s[, keep, drop = FALSE]
        for (col in intersect(c("Area", "Mass", "PeakScore", "PSpatial"), keep))
            s[[col]] <- round(s[[col]], 4)
        s
    })

    output$regionTableHeader <- renderUI({
        req(!is.null(v$regions), nrow(v$regions$summary) > 0)
        h4(if (!is.null(v$regions$summary$PSpatial))
            "Landscape Regions (regional permutation test)" else
            "Landscape Regions (descriptive)")
    })

    output$regionDownloadUI <- renderUI({
        req(!is.null(v$regions), nrow(v$regions$summary) > 0)
        downloadButton("downloadRegions", "Download CSV",
                       style = "margin-bottom:6px")
    })

    output$regionTableOutput <- DT::renderDataTable({
        DT::datatable(
            region_summary(),
            rownames = FALSE,
            options  = list(pageLength = 10),
            caption  = paste("Eight-connected areas beyond the neutral score.",
                "PSpatial compares each area's mass with the largest mass",
                "seen under node-label permutation, over both directions."))
    }, server = TRUE)

    output$downloadRegions <- downloadHandler(
        filename = function() paste0("levi_regions_", Sys.Date(), ".csv"),
        content  = function(file) write.csv(v$regions$summary, file,
                                            row.names = FALSE))

    # ── Expression area (brush) ───────────────────────────────────────────
    output$expArea <- renderText({
        tryCatch({
            expArea <- sum(
                brushedPoints(v$func_ne_return, brush = input$plotBrush) %>%
                    select(Expression), na.rm = TRUE)
            paste0("Expression area: ", expArea)
        },
        error = function(e) "Expression area: ")
    })

    # ── Download ──────────────────────────────────────────────────────────
    downloadNameFun <- reactive({
        switch(input$plotType,
               png  = "chart.png",
               tiff = "chart.tif",
               jpeg = "chart.jpg",
               bmp  = "chart.bmp",
               "chart.png")
    })

    fn_download <- function() {
        req(v$currentPlot)
        tryCatch({
            fname <- downloadNameFun()
            if (input$plotType == "png")  png( fname, width=480, height=480, units="px")
            if (input$plotType == "tiff") tiff(fname, width=480, height=480, units="px")
            if (input$plotType == "jpeg") jpeg(fname, width=480, height=480, units="px",
                                               quality=100)
            if (input$plotType == "bmp")  bmp( fname, width=480, height=480, units="px")
            plot(v$currentPlot)
            dev.off()
        },
        error = function(e) showNotification("Error saving plot",
            type = "error", closeButton = FALSE, duration = 5))
    }

    output$downloadPlot <- downloadHandler(
        filename = downloadNameFun,
        content  = function(file) {
            try({
                fn_download()
                file.copy(downloadNameFun(), file, overwrite = TRUE)
            })
        })
}

# ── Null-coalescing helper (backport for R < 4.4) ─────────────────────────────
`%||%` <- function(x, y) if (!is.null(x)) x else y

shinyApp(ui, server)
