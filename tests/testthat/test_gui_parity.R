library(levi)

# The Shiny interface used to carry its own copy of the pipeline: parsing,
# permutation and plotting were written twice. Sharing the compiled core was
# not enough to keep the two in step -- fixes landed on one side and not the
# other. These tests call both paths on the same input and demand the same
# numbers, so a future divergence fails here instead of going unnoticed.

# levi_shiny() lives in the app file, which also builds the UI. Source only the
# function so no server is started.
gui_fun <- local({
    cached <- NULL
    function() {
        if (!is.null(cached)) return(cached)
        app <- system.file("shiny", "app.R", package = "levi")
        if (!nzchar(app)) return(NULL)

        src <- readLines(app, warn = FALSE)
        start <- grep("^levi_shiny <- function", src)
        if (length(start) != 1L) return(NULL)
        # the closing brace of the function, at column 0
        ends <- grep("^\\}$", src)
        end <- ends[ends > start][1]
        if (is.na(end)) return(NULL)

        # The app is sourced by Shiny with dplyr and igraph attached; give the
        # function the same unqualified names it expects.
        env <- new.env(parent = globalenv())
        env$filter  <- dplyr::filter
        env$arrange <- dplyr::arrange
        env$select  <- dplyr::select
        env$slice   <- dplyr::slice
        env$melt    <- reshape2::melt
        env$graph_from_edgelist <- igraph::graph_from_edgelist
        env$as_long_data_frame <- igraph::as_long_data_frame

        eval(parse(text = paste(src[seq(start, end)], collapse = "\n")),
             envir = env)
        cached <<- env$levi_shiny
        cached
    }
})

skip_without_gui <- function() {
    skip_if_not_installed("dplyr")
    skip_if_not_installed("igraph")
    if (is.null(gui_fun())) skip("levi_shiny() could not be sourced")
}

# The interface hands levi_shiny() a network that is already parsed, and its
# two network arguments are named the other way round: networkCoord carries the
# edges and networkInterac the nodes.
parsed_hub <- function() {
    levi:::.parseNetwork(
        system.file("extdata", "hub_network.dat", package = "levi"),
        NA, "dat")
}

run_script <- function(...) {
    levi(networkCoordinatesInput = system.file("extdata", "hub_network.dat",
             package = "levi"),
         expressionInput = system.file("extdata", "hub_expression.dat",
             package = "levi"),
         fileTypeInput   = "dat",
         geneSymbolInput = "ID",
         readExpColumn   = readExpColumn("Test-Control"),
         resolutionValueInput = 40,
         smoothValueInput     = 50,
         contrastValueInput   = 50,
         zoomValueInput       = 50, inference_unit = "cell", ...)
}

run_gui <- function(...) {
    net <- parsed_hub()
    gui_fun()(
        expression   = read.delim(system.file("extdata", "hub_expression.dat",
                                              package = "levi")),
        fileType     = "dat",
        networkCoord = net$edges,
        networkInterac = net$nodes,
        geneSymbol   = "ID",
        baseTest     = "Test",
        baseControl  = "Control",
        alphaValue   = 50,
        betaValue    = 50,
        backValue    = 40,
        smoothValue  = 50,
        expressionLog = FALSE, ...)
}


test_that("GUI and script mode agree on the node scores", {
    skip_without_gui()
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    s <- run_script()
    g <- run_gui()

    a <- setNames(s$scores$LandscapeScore, s$scores$Gene)
    b <- setNames(g$scores$LandscapeScore, g$scores$Gene)

    expect_setequal(names(a), names(b))
    expect_equal(a[names(a)], b[names(a)], tolerance = 1e-8)
})

test_that("GUI and script mode agree on the landscape surface", {
    skip_without_gui()
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    s <- run_script()
    g <- run_gui()

    # Both keep the surface in long form, under different column names: the
    # script calls the value z, the interface calls it Expression. The script
    # keeps full precision while the interface rounds to two decimals for the
    # tooltip, so the comparison is made at that precision.
    zs <- round(s$landscape$z, 2)
    zg <- g$landscape$Expression

    expect_equal(length(zs), length(zg))
    expect_equal(sort(zs, na.last = TRUE), sort(zg, na.last = TRUE),
                 tolerance = 1e-8)
})

test_that("GUI and script mode agree on the permutation p-values", {
    skip_without_gui()
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    # Both sides draw from the same generator, so the same seed has to give the
    # same null. This is what the shared .permutationPvalues() buys.
    set.seed(99)
    s <- run_script(n_perm = 39)
    set.seed(99)
    g <- run_gui(n_perm = 39, inference_unit = "cell")

    expect_false(is.null(s$pvalues))
    expect_false(is.null(g$pval_over))

    expect_equal(s$pvalues$over,  g$pval_over,  tolerance = 1e-12)
    expect_equal(s$pvalues$under, g$pval_under, tolerance = 1e-12)
})

test_that("the interface reuses the package helpers rather than copying them", {
    app <- system.file("shiny", "app.R", package = "levi")
    skip_if(!nzchar(app), "app.R not installed")
    src <- paste(readLines(app, warn = FALSE), collapse = "\n")

    # Shared: a change to any of these has to reach both paths at once.
    for (helper in c("levi_function", ".colorSet", ".parseNetwork",
                     ".buildLandscapeChart"))
        expect_true(grepl(helper, src, fixed = TRUE),
                    info = paste("the interface no longer calls", helper))

    # The old private copy of the permutation stack must not come back.
    expect_false(grepl("permLandscapes", src, fixed = TRUE))
    expect_false(grepl("obs3d", src, fixed = TRUE))
})


test_that("the two charts are built from the same recipe", {
    skip_without_gui()
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    # The interface keeps the surface under different column names and passes
    # an explicit colour vector; everything else has to come out the same.
    # The scale arrows in particular had drifted, sitting at different heights
    # on each side.
    s <- run_script()
    g <- run_gui()

    gui_df <- g$landscape
    script_chart <- levi:::.buildLandscapeChart(
        s$landscape, "default", NULL, FALSE)
    gui_chart <- levi:::.buildLandscapeChart(
        gui_df, "default", NULL, FALSE,
        colours = levi:::.colorSet("default"),
        cols = c("X", "Y", "Expression"))

    expect_equal(length(script_chart$layers), length(gui_chart$layers))

    bs <- ggplot2::ggplot_build(script_chart)$data
    bg <- ggplot2::ggplot_build(gui_chart)$data

    # The annotation layers -- the decrease/increase text and the two arrows --
    # carry no data of their own, so they must match exactly.
    for (k in seq_along(bs)[-1])
        expect_equal(bs[[k]], bg[[k]],
                     info = paste("annotation layer", k, "differs"))
})

test_that("the interface can build the same 3D surface as script mode", {
    skip_if_not_installed("plotly")
    skip_without_gui()
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    # The GUI used to have no 3D at all, while LEVIui's help claimed it exposed
    # everything script mode offers. Both now go through .buildSurface3D().
    set.seed(7)
    s <- run_script(n_perm = 49, plot3d = TRUE)
    set.seed(7)
    g <- run_gui(n_perm = 49, inference_unit = "cell")

    expect_false(is.null(g$surface))

    gui3d <- levi:::.buildSurface3D(
        zmat    = g$surface,
        colours = levi:::.colorSet("default"),
        pvals   = list(over = g$pval_over, under = g$pval_under),
        i       = seq_len(nrow(g$surface)),
        sig_level = g$sig_level, perm_side = g$perm_side)

    bs <- plotly::plotly_build(s$plot3d)$x$data
    bg <- plotly::plotly_build(gui3d)$x$data

    expect_equal(length(bs), length(bg))
    expect_equal(vapply(bs, function(d) d$type, character(1)),
                 vapply(bg, function(d) d$type, character(1)))

    surf_s <- bs[[which(vapply(bs, function(d) d$type, character(1)) ==
                        "surface")[1]]]
    surf_g <- bg[[which(vapply(bg, function(d) d$type, character(1)) ==
                        "surface")[1]]]
    expect_equal(surf_s$z, surf_g$z, tolerance = 1e-10)
})

test_that("the interface leaves out the 3D control when plotly is missing", {
    app <- system.file("shiny", "app.R", package = "levi")
    skip_if(!nzchar(app), "app.R not installed")
    src <- paste(readLines(app, warn = FALSE), collapse = "\n")

    # plotly is suggested, so the app has to start without it. Every use is
    # behind has_plotly; a bare plotly:: call in the UI would stop the app from
    # launching for anyone who has not installed it.
    expect_true(grepl("has_plotly <- requireNamespace", src, fixed = TRUE))
    expect_true(grepl("if (has_plotly) output$graph3d", src, fixed = TRUE))
    expect_true(grepl("if (has_plotly) conditionalPanel", src, fixed = TRUE))
})


# --- the app itself, driven through shiny::testServer -------------------------

test_that("GUI and script mode agree on the regional permutation test", {
    skip_without_gui()
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    # The interface defaults to the regional test, like levi(). Same seed,
    # same node-label draws, same regions and the same p-values.
    set.seed(11)
    s <- levi(networkCoordinatesInput = system.file("extdata", "hub_network.dat",
                  package = "levi"),
              expressionInput = system.file("extdata", "hub_expression.dat",
                  package = "levi"),
              fileTypeInput = "dat", geneSymbolInput = "ID",
              readExpColumn = readExpColumn("Test-Control"),
              resolutionValueInput = 40, smoothValueInput = 50,
              contrastValueInput = 50, zoomValueInput = 50, n_perm = 19)
    set.seed(11)
    g <- run_gui(n_perm = 19)

    expect_equal(g$inference_unit, "region")
    expect_null(g$pval_over)
    expect_equal(g$regions$summary$Region, s$regions$summary$Region)
    expect_equal(g$regions$summary$PSpatial, s$regions$summary$PSpatial)
    expect_equal(g$regions$null_max_mass, s$regions$null_max_mass)
})

# Everything above exercises levi_shiny() in isolation. These drive the whole
# server: inputs in, reactive values out. Without them the interface had no
# end-to-end coverage at all, which is how it came to be unable to compute
# anything at all under igraph >= 2.0 without that being noticed.

app_file <- function() system.file("shiny", "app.R", package = "levi")

test_that("the app computes a landscape from uploaded files", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("dplyr")
    skip_if(!nzchar(app_file()), "app.R not installed")

    net  <- system.file("extdata", "hub_network.dat",    package = "levi")
    expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(
            exprSource = "file", fileType = "dat",
            file  = list(datapath = expr, name = "hub_expression.dat"),
            file2 = list(datapath = net,  name = "hub_network.dat"),
            geneSymbol = "ID", baseTest = "Test", baseControl = "Control",
            contrast = 50, zoom = 50, size = 40, smooth = 50, log = FALSE,
            nperm = 0, sig_level_ui = 0.05, perm_side = "both",
            signal_mode = "ratio", logfc_k = 1, setcolor = "pp",
            palette = "multi", contour = FALSE, plot3d = FALSE,
            geneSearch = "", showPeakLabels = FALSE)

        session$setInputs(action = 1)

        # A deprecation warning used to abort the whole run here, so the most
        # valuable assertion is simply that something came out.
        expect_equal(nrow(v$scoreTable), 9)
        expect_false(is.null(v$surface))
        # size = 40 maps to a grid of (40/100)*210 + 30 = 114 cells a side.
        expect_equal(nrow(v$func_ne_return), 114L * 114L)
        expect_equal(dim(v$surface), c(114L, 114L))
        expect_true(all(v$scoreTable$LandscapeScore >= 0 &
                        v$scoreTable$LandscapeScore <= 1))
    })
})

test_that("the app agrees with script mode end to end", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("dplyr")
    skip_if(!nzchar(app_file()), "app.R not installed")
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    net  <- system.file("extdata", "hub_network.dat",    package = "levi")
    expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    from_script <- run_script()
    a <- setNames(from_script$scores$LandscapeScore, from_script$scores$Gene)

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(
            exprSource = "file", fileType = "dat",
            file  = list(datapath = expr, name = "hub_expression.dat"),
            file2 = list(datapath = net,  name = "hub_network.dat"),
            geneSymbol = "ID", baseTest = "Test", baseControl = "Control",
            contrast = 50, zoom = 50, size = 40, smooth = 50, log = FALSE,
            nperm = 0, sig_level_ui = 0.05, perm_side = "both",
            signal_mode = "ratio", logfc_k = 1, setcolor = "pp",
            palette = "multi", contour = FALSE, plot3d = FALSE,
            geneSearch = "", showPeakLabels = FALSE)
        session$setInputs(action = 1)

        b <- setNames(v$scoreTable$LandscapeScore, v$scoreTable$Gene)
        expect_equal(a[names(a)], b[names(a)], tolerance = 1e-8)
    })
})

test_that("a warning does not abort the computation", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("dplyr")
    skip_if(!nzchar(app_file()), "app.R not installed")

    # ratio mode on signed data warns. The interface used to catch warnings
    # with tryCatch, which unwinds: the warning killed the run and the user was
    # told "Incorrect file format". It must now notify and carry on.
    net <- system.file("extdata", "hub_network.dat", package = "levi")
    signed <- file.path(tempdir(), "signed_expression.dat")
    write.table(
        data.frame(ID = c("HUB", paste0("N", 1:8)),
                   Test = c(4.3, 4.1, 4.4, 4.2, 4.3, -5.3, -5.1, -5.4, -5.2),
                   Control = rep(0, 9)),
        signed, sep = "\t", row.names = FALSE, quote = FALSE)

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(
            exprSource = "file", fileType = "dat",
            file  = list(datapath = signed, name = "signed_expression.dat"),
            file2 = list(datapath = net, name = "hub_network.dat"),
            geneSymbol = "ID", baseTest = "Test", baseControl = "Control",
            contrast = 50, zoom = 50, size = 30, smooth = 50, log = FALSE,
            nperm = 0, sig_level_ui = 0.05, perm_side = "both",
            signal_mode = "ratio", logfc_k = 1, setcolor = "pp",
            palette = "multi", contour = FALSE, plot3d = FALSE,
            geneSearch = "", showPeakLabels = FALSE)
        session$setInputs(action = 1)

        expect_equal(nrow(v$scoreTable), 9)
    })
})

test_that("the interface delegates network calculations to the shared pipeline", {
    skip_if(!nzchar(app_file()), "app.R not installed")
    src <- paste(readLines(app_file(), warn = FALSE), collapse = "\n")

    # graph.edgelist() has been deprecated since igraph 2.0.0. Combined with
    # the tryCatch above, its deprecation notice was fatal.
    expect_false(grepl("graph.edgelist(", src, fixed = TRUE))
    expect_true(grepl("levi:::levi_function(", src, fixed = TRUE))
})


# --- the rest of the interface -----------------------------------------------

test_that("no tryCatch in the app treats a warning as fatal", {
    skip_if(!nzchar(app_file()), "app.R not installed")

    # tryCatch unwinds the stack when it catches, so a warning handler there
    # aborts whatever was running. That is how a deprecation notice from igraph
    # left the interface unable to compute, and how a file without a trailing
    # newline was reported as "Incorrect file format". Warnings belong in
    # withCallingHandlers, which notifies and carries on.
    findBad <- function(e, acc = list()) {
        if (is.call(e)) {
            if (identical(as.character(e[[1]])[1], "tryCatch") &&
                "warning" %in% names(as.list(e)))
                acc[[length(acc) + 1L]] <- deparse(e)[1]
            for (part in as.list(e)[-1])
                if (!missing(part)) acc <- findBad(part, acc)
        }
        acc
    }

    bad <- list()
    for (e in parse(app_file())) bad <- c(bad, findBad(e))
    expect_equal(length(bad), 0L)
})

test_that("a file with no trailing newline is still read", {
    skip_if_not_installed("shiny")
    skip_if(!nzchar(app_file()), "app.R not installed")

    # read.table warns "incomplete final line found" and parses the file
    # correctly. The interface used to discard the data and blame the format.
    f <- file.path(tempdir(), "no_newline.dat")
    con <- file(f, "wb")
    writeBin(charToRaw("ID\tTest\tControl\nA\t1\t2\nB\t3\t4"), con)
    close(con)

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(file = list(datapath = f, name = "no_newline.dat"))
        got <- baseSelect()

        expect_s3_class(got, "data.frame")
        expect_equal(nrow(got), 2L)
        expect_identical(names(got), c("ID", "Test", "Control"))
    })
})

test_that("brushing the map returns the genes under the selection", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("dplyr")
    skip_if(!nzchar(app_file()), "app.R not installed")

    net  <- system.file("extdata", "hub_network.dat",    package = "levi")
    expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(
            exprSource = "file", fileType = "dat",
            file  = list(datapath = expr, name = "hub_expression.dat"),
            file2 = list(datapath = net,  name = "hub_network.dat"),
            geneSymbol = "ID", baseTest = "Test", baseControl = "Control",
            contrast = 50, zoom = 50, size = 30, smooth = 50, log = FALSE,
            nperm = 0, sig_level_ui = 0.05, perm_side = "both",
            signal_mode = "ratio", logfc_k = 1, setcolor = "pp",
            palette = "multi", contour = FALSE, plot3d = FALSE,
            geneSearch = "", showPeakLabels = FALSE)
        session$setInputs(action = 1)

        # A brush covering the whole map has to reach every named node.
        n <- sqrt(nrow(v$func_ne_return))
        session$setInputs(plotBrush = list(
            xmin = 0, xmax = n + 1, ymin = 0, ymax = n + 1,
            direction = "xy", mapping = list(x = "X", y = "Y"),
            domain = list(left = 0, right = n + 1, bottom = 0, top = n + 1),
            range = list(left = 0, right = n + 1, bottom = 0, top = n + 1),
            log = list(x = NULL, y = NULL), outputId = "graph"))

        expect_match(output$expArea, "Expression area")
    })
})

test_that("the 3D output is produced when the box is ticked", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("plotly")
    skip_if(!nzchar(app_file()), "app.R not installed")

    net  <- system.file("extdata", "hub_network.dat",    package = "levi")
    expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(
            exprSource = "file", fileType = "dat",
            file  = list(datapath = expr, name = "hub_expression.dat"),
            file2 = list(datapath = net,  name = "hub_network.dat"),
            geneSymbol = "ID", baseTest = "Test", baseControl = "Control",
            contrast = 50, zoom = 50, size = 30, smooth = 50, log = FALSE,
            nperm = 0, sig_level_ui = 0.05, perm_side = "both",
            signal_mode = "ratio", logfc_k = 1, setcolor = "pp",
            palette = "multi", contour = FALSE, plot3d = TRUE,
            geneSearch = "", showPeakLabels = FALSE)
        session$setInputs(action = 1)

        expect_false(is.null(output$graph3d))
    })
})


# --- Bioconductor mode and downloads -----------------------------------------

# The two surfaces the audits left uncovered. The interface reads a serialised
# Bioconductor object instead of a tab-delimited file, and hands the score and
# peak tables back as CSV.

gui_inputs <- function(net, ...) {
    c(list(fileType = "dat",
           file2 = list(datapath = net, name = "hub_network.dat"),
           contrast = 50, zoom = 50, size = 30, smooth = 50, log = FALSE,
           nperm = 0, sig_level_ui = 0.05, perm_side = "both",
           signal_mode = "ratio", logfc_k = 1, setcolor = "pp",
           palette = "multi", contour = FALSE, plot3d = FALSE,
           geneSearch = "", showPeakLabels = FALSE),
      list(...))
}

test_that("the app computes from a SummarizedExperiment loaded as RDS", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("SummarizedExperiment")
    skip_if(!nzchar(app_file()), "app.R not installed")

    genes <- c("HUB", paste0("N", 1:8))
    counts <- matrix(
        c(200, 200, 200, 200, 200,   5,   5,   5,   5,
          190, 210, 195, 205, 200,   6,   4,   5,   5,
           10,  10,  10,  10,  10, 200, 200, 200, 200),
        nrow = 9,
        dimnames = list(genes, c("tumor1", "tumor2", "normal1")))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(counts = counts),
        colData = data.frame(condition = c("Tumor", "Tumor", "Normal"),
                             row.names = colnames(counts)))

    rds <- file.path(tempdir(), "se_for_gui.rds")
    saveRDS(se, rds)

    net <- system.file("extdata", "hub_network.dat", package = "levi")

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        do.call(session$setInputs, gui_inputs(net,
            exprSource = "bioc",
            bioc_rds = list(datapath = rds, name = "se_for_gui.rds"),
            bioc_assay = "counts", bioc_cond_col = "condition",
            bioc_test_level = "Tumor", bioc_ctrl_level = "Normal",
            bioc_log = FALSE))

        # The object has to come back before anything is computed from it.
        expect_s4_class(bioc_obj(), "SummarizedExperiment")

        session$setInputs(action = 1)

        expect_equal(nrow(v$scoreTable), 9)
        expect_true(all(v$scoreTable$LandscapeScore >= 0 &
                        v$scoreTable$LandscapeScore <= 1))

        # The core is over-expressed and the corners repressed in this data,
        # the same conclusion script mode reaches through leviFromSE().
        sc <- setNames(v$scoreTable$LandscapeScore, v$scoreTable$Gene)
        expect_gt(sc[["HUB"]], 0.5)
        expect_lt(sc[["N5"]], 0.5)
    })
})

test_that("the app refuses an RDS that is not a supported class", {
    skip_if_not_installed("shiny")
    skip_if(!nzchar(app_file()), "app.R not installed")

    rds <- file.path(tempdir(), "not_bioc.rds")
    saveRDS(data.frame(a = 1:3), rds)
    net <- system.file("extdata", "hub_network.dat", package = "levi")

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        do.call(session$setInputs, gui_inputs(net,
            exprSource = "bioc",
            bioc_rds = list(datapath = rds, name = "not_bioc.rds")))

        session$setInputs(action = 1)

        # Nothing is computed, and the app is still standing. The reactive
        # starts life as an empty data.frame, not NULL.
        expect_equal(nrow(v$scoreTable), 0L)
    })
})

test_that("the download handlers write the tables as CSV", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("dplyr")
    skip_if(!nzchar(app_file()), "app.R not installed")

    net  <- system.file("extdata", "hub_network.dat",    package = "levi")
    expr <- system.file("extdata", "hub_expression.dat", package = "levi")

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        do.call(session$setInputs, gui_inputs(net,
            exprSource = "file",
            file = list(datapath = expr, name = "hub_expression.dat"),
            geneSymbol = "ID", baseTest = "Test", baseControl = "Control"))
        session$setInputs(action = 1)

        # testServer runs the handler and hands back the path of the file it
        # wrote, so the download can be read exactly as a user would get it.
        got <- read.csv(output$downloadScores)

        expect_equal(nrow(got), 9)
        expect_true(all(c("Gene", "LandscapeScore", "Rank") %in% names(got)))
        expect_equal(got$Gene, v$scoreTable$Gene)
        expect_equal(got$LandscapeScore, v$scoreTable$LandscapeScore,
                     tolerance = 1e-8)

        peaks <- read.csv(output$downloadPeaks)
        expect_equal(nrow(peaks), nrow(v$peakTable))
        expect_true("Type" %in% names(peaks))
    })
})

test_that("the download file names carry the date", {
    skip_if_not_installed("shiny")
    skip_if(!nzchar(app_file()), "app.R not installed")

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        expect_match(basename(output$downloadScores),
                     "^levi_scores_\\d{4}-\\d{2}-\\d{2}\\.csv$")
        expect_match(basename(output$downloadPeaks),
                     "^levi_peaks_\\d{4}-\\d{2}-\\d{2}\\.csv$")
    })
})
