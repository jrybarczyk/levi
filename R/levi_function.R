utils::globalVariables(c("x", "y", "xend", "yend", "PeakRow", "PeakCol", "Label"))

levi_function <- function(expressionInput, fileTypeInput, networkCoordinatesInput,
    networkInteractionsInput, geneSymbolInput, readExpColumn,
    contrastValueInput, zoomValueInput, resolutionValueInput,
    smoothValueInput, expressionLog, contourLevi, setcolor,
    plot3d = FALSE, n_perm = 0, sig_level = 0.05,
    perm_side = c("both", "over", "under"),
    signal_mode = c("ratio", "logfc", "zscore"), logfc_k = 1,
    p_adjust_method = "BY", region_threshold = 0.1, region_min_cells = 3L,
    .parsed_network = NULL, .draw = TRUE,
    .progress = NULL, inference_unit = c("region", "cell"), perm_strata = NULL,
    edge_weighting = c("midpoint", "degree", "none")) {
    inference_unit <- match.arg(inference_unit)
    edge_weighting <- match.arg(edge_weighting)
    perm_side   <- match.arg(perm_side)
    signal_mode <- match.arg(signal_mode)
    p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
    .validateScalar(n_perm, "n_perm", 0, integer = TRUE)
    .validateScalar(sig_level, "sig_level", 0, 1)
    .validateScalar(logfc_k, "logfc_k", .Machine$double.eps)
    .validateScalar(region_threshold, "region_threshold", .Machine$double.eps,
                    0.5 - .Machine$double.eps)
    .validateScalar(region_min_cells, "region_min_cells", 1, integer = TRUE)




    fileType <- match.arg(fileTypeInput, c("dat", "dyn", "stg", "net"))


    leviResults <- vector("list", length(readExpColumn) - 1)

    # ---------------------------------------------------------------------
    # Invariant setup: image parameters, network parsing, graph and
    # coordinates do not depend on the comparison, so they are done only once
    # instead of on every iteration of the loop over readExpColumn.
    # ---------------------------------------------------------------------
    #Configuration of contrast, resolution, smothing and zoom
    #contrast -> silhouette occupancy threshold.
    # More contrast = higher threshold = tighter silhouette around the network.
    {contrastValue <- contrastValueInput}
    if (contrastValue < 0) {contrastValue <- 0}
    if (contrastValue > 100) {contrastValue <- 100}
    occFrac <- 0.002 + 0.096 * (contrastValue/100)

    #resolution
    {resolutionValue <- resolutionValueInput}
    if (resolutionValue > 100) {resolutionValue <- 100}
    if (resolutionValue < 1) {resolutionValue <- 1}
    resolutionValue<-as.integer((resolutionValue/100)*210+30)

    #smothing -> sigma of the Gaussian kernel, in grid cells.
    # Proportional to resolutionValue so that the apparent smoothing does not
    # change when the user alters the resolution.
    {smoothValue <- smoothValueInput}
    if (smoothValue < 0) {smoothValue <- 0}
    if (smoothValue > 100) {smoothValue <- 100}
    sigmaCells <- (0.005 + 0.08 * (smoothValue/100)) * resolutionValue
    if (sigmaCells < 0.35) {sigmaCells <- 0.35}

    #zoom -> margin of the grid around the network, in coordinate units
    # (the network spans [0, 1] on its longer axis).
    # The silhouette reaches beyond the outermost node as far as the kernel
    # carries occupancy above occFrac, about sigma * sqrt(-2 log occFrac)
    # cells, plus one cell for the bilinear deposit. The grid always leaves
    # at least that much room, so the landscape is never clipped; the zoom
    # adds up to 20% of the network extent on top of it (zoom 0 = widest
    # frame, zoom 100 = frame that just fits the silhouette). Before 2.0.0
    # the margin was the zoom alone and the default cut the silhouette
    # whenever a node sat near the border.
    {zoomValue <- zoomValueInput}
    if (zoomValue < 0) {zoomValue <- 0}
    if (zoomValue > 100) {zoomValue <- 100}
    zoomFraction <- zoomValue/100
    # The margin itself is computed once the support points are known, below.

    nameBase <- expressionInput
    networkNodes <- networkCoordinatesInput
    networkEdges <- networkInteractionsInput
    geneSymbol<- geneSymbolInput
    parsed <- .parsed_network %||% .parseNetwork(networkNodes, networkEdges, fileType)
    nodes <- parsed$nodes
    edges <- parsed$edges


    #Remove"NA" and "-" from expression file
        if (is.character(nameBase)) {
            expression <- read.delim(file = nameBase, header = TRUE,
                sep = "\t", quote = "")
        } else if (is.data.frame(nameBase) || is.matrix(nameBase)) {
            expression <- as.data.frame(nameBase)
        } else {
            stop("expressionInput must be a file path, data.frame, or matrix")
        }

        # Validated here, before the first indexing below: the later check in
        # the comparison loop only runs after these subsets, which failed with
        # the unhelpful "undefined columns selected".
        if (!is.character(geneSymbol) || length(geneSymbol) != 1L ||
            is.na(geneSymbol) || !(geneSymbol %in% names(expression))) {
            stop("'geneSymbolInput' must name one column of the expression ",
                "data. Got \"", paste(geneSymbol, collapse = ", "),
                "\"; the columns available are: ",
                paste(names(expression), collapse = ", "), call. = FALSE)
        }

        expression <- subset(expression,expression[,paste(geneSymbol)] !=
        "NA")
        expression <- subset(expression,expression[,paste(geneSymbol)] !=
        "-")
        expression <- unique(expression)
        head_express <- as.list(names(expression))
    if (fileType == "dat"){
        edges <- edges[, c(1, 2)]
        nodes <- as.data.frame(nodes)
        nodes[, c(2)] <- vapply(nodes[, c(2)], as.double, numeric(1))
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

    nodes <- nodes[order(nodes[, 1]), ]


    nodes$V1 <- as.character(nodes$V1)
    nodes$V2 <- as.numeric(nodes$V2)
    nodes$V3 <- as.numeric(nodes$V3)
    nodesCoord <- aggregate(nodes[, 2:3], nodes[1], mean)
    nodesForMerge <- nodesCoord
    edges <- unique(edges[, 1:2])
    edge_index <- cbind(match(as.character(edges[, 1]), nodesCoord$V1),
                        match(as.character(edges[, 2]), nodesCoord$V1))
    if (anyNA(edge_index))
        stop("Every edge endpoint must have network coordinates.")
    nnodes <- nrow(nodesCoord)
    coord <- rbind(as.matrix(nodesCoord[, 2:3]),
        (as.matrix(nodesCoord[edge_index[, 1], 2:3]) +
         as.matrix(nodesCoord[edge_index[, 2], 2:3])) / 2)
    if (!all(is.finite(coord))) stop("Network coordinates must be finite.")
    support_weights <- .supportWeights(nnodes, edge_index, edge_weighting)

    # normalization and centralization
    minCoordX <- min(coord[,c(1)])
    maxCoordX <- max(coord[,c(1)])
    centroX <- (minCoordX+maxCoordX)/2
    minCoordY <- min(coord[,c(2)])
    maxCoordY <- max(coord[,c(2)])
    centroY <- (minCoordY+maxCoordY)/2

    coordRange <- max(maxCoordX - minCoordX, maxCoordY - minCoordY)
    if (coordRange == 0) coordRange <- 1
    centroX <- centroX/coordRange
    centroY <- centroY/coordRange

    coord[,c(1)] <- (coord[,c(1)]/coordRange)+(0.5-centroX)
    coord[,c(2)] <- (coord[,c(2)]/coordRange)+(0.5-centroY)

    # Grid margin. The silhouette keeps the cells whose occupancy exceeds
    # occFrac times that of an isolated point; W points stacked near the
    # border push it out to about sigma * sqrt(2 log(W / occFrac)) cells,
    # plus one cell for the bilinear deposit. Solving
    # m >= reachCells * (1 + 2m) / (resolutionValue - 1) for the margin m.
    # W is the largest kernel-weighted stack of support points around any
    # point (the occupancy at that point relative to an isolated one). For
    # very large networks the total weight is used as an upper bound instead
    # of the n x n distance matrix.
    sigmaCoord <- sigmaCells * 1.4 / (resolutionValue - 1)
    stackWeight <- if (nrow(coord) <= 3000L) {
        d2 <- as.matrix(stats::dist(coord))^2
        max(colSums(support_weights * exp(-d2 / (2 * sigmaCoord^2))))
    } else sum(support_weights)
    stackWeight <- max(stackWeight, 1)
    reachCells <- sigmaCells * sqrt(2 * log(stackWeight / occFrac)) + 1
    denominator <- resolutionValue - 1 - 2 * reachCells
    reachMargin <- if (denominator > 0) reachCells / denominator else 1
    reachMargin <- min(reachMargin, 1)
    gridFor <- function(margin) {
        zoomValue <- -margin
        list(zoom = zoomValue, increase = (1 + 2 * margin) / (resolutionValue - 1))
    }
    # That bound is generous, so the silhouette is measured once on the
    # generous grid (occupancy depends on coordinates and weights only, not
    # on the signal) and the frame is tightened to what it actually needs,
    # plus two cells of safety. zoom 100 then frames the silhouette exactly
    # and zoom 0 adds 20% of the network extent around it.
    probeGrid <- gridFor(reachMargin)
    dummy <- matrix(0.5, nrow(coord), 1L)
    probe <- landscape_gauss(coord = coord, SignalOut = dummy, signalExp = dummy,
        signalCtrl = dummy, resolutionValue = resolutionValue,
        zoomValue = probeGrid$zoom, increase = probeGrid$increase,
        sigma = sigmaCells, occFrac = occFrac, weights = support_weights)$m1
    inside <- which(!is.na(probe), arr.ind = TRUE)
    if (nrow(inside)) {
        lo <- probeGrid$zoom + (apply(inside, 2, min) - 1) * probeGrid$increase
        hi <- probeGrid$zoom + (apply(inside, 2, max) - 1) * probeGrid$increase
        needed <- max(0, -lo, hi - 1) + 2 * probeGrid$increase
        reachMargin <- min(needed, reachMargin)
    }
    finalGrid <- gridFor(reachMargin + 0.2 * (1 - zoomFraction))
    zoomValue <- finalGrid$zoom
    increase <- finalGrid$increase


    for (k in seq(2,length(readExpColumn))) {
        rng_state <- if (exists(".Random.seed", envir = .GlobalEnv))
            get(".Random.seed", envir = .GlobalEnv) else NULL
        columnComb<- do.call('rbind',
            strsplit(as.character(readExpColumn[k]),'-',
            fixed=TRUE))

            baseTest<- columnComb[,1]
            baseControl<- columnComb[,2]

            if (baseControl == " ") {
                baseControl <- baseTest
            }
            arguments <- list(geneSymbol, baseTest, baseControl)
            for (i in seq(arguments)){
                if (!is.element(arguments[i], head_express) ) {
                    stop("Column not found in the expression data: ",
                         arguments[i])}
            }
            expressSelect <- expression[, c(geneSymbol, baseTest, baseControl)]
            # Back-transform from log2 only for ratio mode.
            # logfc/zscore modes compute logFC = Test - Control directly,
            # so log2 values must NOT be exponentiated first.
            if (expressionLog && signal_mode == "ratio") {
                expressSelect[, 2:3] <- 2^expressSelect[, 2:3]
            } else if (expressionLog && signal_mode != "ratio") {
                message("Note: expressionLog = TRUE is ignored when ",
                        "signal_mode = '", signal_mode, "'. ",
                        "Log2 values are used directly to compute logFC.")
            }
            newExpression <- aggregate(x = expressSelect[c
                 (baseControl,baseTest)],
                 by = expressSelect[c(geneSymbol)],
                 FUN = function(media_valor){
                     mean(media_valor)
                     })

            # Repeated identifiers are averaged, which is what mapping probes
            # onto symbols calls for. Say so: two contradictory measurements
            # average to the neutral point, and a silently neutral gene reads
            # exactly like a gene that genuinely did not change.
            collapsed <- nrow(expressSelect) - nrow(newExpression)
            if (collapsed > 0)
                message(collapsed, " duplicated identifier(s) in the ",
                    "expression data were averaged into ",
                    nrow(newExpression), " unique entries.")



        #signalCoordMerge have values of controle and test
        signalCoordMerge <- merge(nodesForMerge, newExpression, by.x = "V1",
            by.y = geneSymbol,
            all.x = TRUE)
        #signalCoordMerge[is.na(signalCoordMerge)] <- 0
        # Nodes whose identifier found no match in the expression table.
        # The check has to run on signalCoordMerge, which holds one row per
        # node: edgesSignalMerge repeats a node once per edge, and its column
        # positions shift with every merge, so fixed indices into it silently
        # stopped pointing at the expression columns.
        naTotal <- as.matrix(unique(
            signalCoordMerge[!is.finite(signalCoordMerge[[baseTest]]) |
                             !is.finite(signalCoordMerge[[baseControl]]), 1]))

        #Create title for chart
        if (baseTest == baseControl) {
            titleChart <- baseTest
        } else {
            titleChart <- paste(baseTest,baseControl, sep = '-')
        }

        #Creates log if exists nodes without expression value
        if (length(naTotal) > 0) {
            logDir <- file.path(tempdir(), titleChart)
            if (!dir.exists(logDir)) dir.create(logDir, recursive = TRUE)
            logPath <- file.path(logDir, "levi.log")
            writeLines(as.vector(naTotal), logPath)

            if (nrow(naTotal) == nrow(signalCoordMerge)) {
                # No identifier matched at all. The landscape comes out
                # uniformly neutral, which is indistinguishable from a dataset
                # with no variation, so this has to be a warning and not a note.
                warning("None of the ", nrow(signalCoordMerge),
                    " network nodes matched an identifier in the expression ",
                    "data, so the landscape is uniformly neutral and carries ",
                    "no information. Check that 'geneSymbolInput' names the ",
                    "right column and that the network and the expression ",
                    "data use the same identifier type. Node names: ",
                    logPath, call. = FALSE)
            } else {
                message("There are ", nrow(naTotal), " nodes without ",
                    "expression value, see log in path: ", logPath)
            }
        }

        single_col <- (baseTest == baseControl)
        node_values <- cbind(signalCoordMerge[[baseTest]],
                             signalCoordMerge[[baseControl]])
        resolved_strata <- NULL
        if (!is.null(perm_strata)) {
            if (!is.atomic(perm_strata))
                stop("'perm_strata' must be an atomic vector.", call. = FALSE)
            if (!is.null(names(perm_strata)))
                resolved_strata <- perm_strata[as.character(signalCoordMerge$V1)]
            else if (length(perm_strata) == nrow(node_values))
                resolved_strata <- perm_strata
            else stop("Unnamed 'perm_strata' must have one value per network node.",
                      call. = FALSE)
            measured_here <- rowSums(is.finite(node_values)) == ncol(node_values)
            if (anyNA(resolved_strata[measured_here]))
                stop("'perm_strata' is missing a stratum for a measured node.",
                     call. = FALSE)
            resolved_strata <- as.character(resolved_strata)
        }
        signals <- .networkSignals(node_values, edge_index, single_col,
                                    signal_mode, logfc_k)
        SignalOut <- signals$signal
        signalExp <- signals$test
        signalCtrl <- signals$control
        numberCoord <- nrow(SignalOut)

        # Landscape by normalised convolution. Each cell value is the average
        # of the signals weighted by a Gaussian kernel, rather than the mean
        # over the k nearest points. Because the sum of the weights appears
        # in the denominator, the result is always a convex average of the
        # signals:
        # the scale no longer depends on the smoothing or the network density.
        matrixFinal <- landscape_gauss(
            coord           = coord[seq_len(numberCoord), , drop = FALSE],
            SignalOut       = SignalOut,
            signalExp       = signalExp,
            signalCtrl      = signalCtrl,
            resolutionValue = resolutionValue,
            zoomValue       = zoomValue,
            increase        = increase,
            sigma           = sigmaCells,
            occFrac         = occFrac,
            weights         = support_weights)

        matrixOut <- matrixFinal$m1
        n <- resolutionValue
        i <- seq_len(n)
        ExpCtrl <- matrixOut[i, rev(i)]

        landgraph <- melt(ExpCtrl, value.name = "z")

        landgraphFinal <- as.data.frame(landgraph[,c(1,2,3)])


        landgraphChart <- .buildLandscapeChart(landgraphFinal, setcolor,
            titleChart, contourLevi) +
            ggplot2::labs(caption = .signalMeaning(signal_mode, single_col))

        # -- 1. Node landscape scores --
        nodeCoordNorm <- coord[seq_len(nnodes), , drop = FALSE]
        xIdx <- pmin(pmax(
            round((nodeCoordNorm[, 1] - zoomValue) / increase) + 1L, 1L),
            resolutionValue)
        yIdx <- pmin(pmax(
            round((nodeCoordNorm[, 2] - zoomValue) / increase) + 1L, 1L),
            resolutionValue)
        # A node whose cell falls outside the silhouette returns NA; in that
        # case the node's own signal is used, which is the value the convolution
        # would tend to give.
        nodeScores <- mapply(function(xi, yi) matrixOut[xi, yi], xIdx, yIdx)
        naScore <- is.na(nodeScores)
        if (any(naScore))
            nodeScores[naScore] <- SignalOut[seq_len(nnodes), 1][naScore]

        scoreTable <- data.frame(
            Gene           = as.character(nodesCoord[, 1]),
            X              = as.numeric(nodesCoord[, 2]),
            Y              = as.numeric(nodesCoord[, 3]),
            LandscapeScore = round(nodeScores, 4),
            stringsAsFactors = FALSE
        )
        scoreTable <- scoreTable[
            order(scoreTable$LandscapeScore, decreasing = TRUE), ]
        scoreTable$Rank <- seq_len(nrow(scoreTable))
        rownames(scoreTable) <- NULL

        # -- 2. Automatic peak / valley detection --
        peakTable <- .detectPeaks2D(
            matrixOut, coord, nnodes, nodesCoord,
            zoomValue, increase, resolutionValue
        )
        regions <- .extractLandscapeRegions(
            matrixOut, zoomValue, increase,
            threshold = region_threshold, min_cells = region_min_cells)

        # -- 3. Permutation significance test --
        pvalMatrix <- rawPvalues <- NULL
        if (n_perm > 0L) {
            message("Permutation test: ", n_perm, " iterations for '",
                    titleChart, "'...")

            pvalMatrix <- .permutationPvalues(
                coord           = coord[seq_len(numberCoord), , drop = FALSE],
                SignalOut       = SignalOut,
                signalExp       = signalExp,
                signalCtrl      = signalCtrl,
                matrixOut       = matrixOut,
                resolutionValue = resolutionValue,
                zoomValue       = zoomValue,
                increase        = increase,
                sigma           = sigmaCells,
                occFrac         = occFrac,
                n_perm          = n_perm,
                progress        = .progress,
                node_values     = node_values,
                edge_index      = edge_index,
                single_col      = single_col,
                signal_mode     = signal_mode,
                logfc_k         = logfc_k,
                regions         = if (inference_unit == "region") regions else NULL,
                perm_strata     = resolved_strata,
                weights         = support_weights)
            if (inference_unit == "region") {
                regions <- pvalMatrix
                pvalMatrix <- NULL
                regions$summary$Significant <- regions$summary$PSpatial <= sig_level
                selected <- regions$summary$Region[regions$summary$Significant &
                    (perm_side == "both" | regions$summary$Direction == perm_side)]
                boundary <- .regionBoundaries(regions, n, selected)
                if (nrow(boundary)) landgraphChart <- landgraphChart +
                    ggplot2::geom_segment(data = boundary,
                        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                        colour = "white", linewidth = 0.8, inherit.aes = FALSE)
            } else {
            rawPvalues <- pvalMatrix
            pvalMatrix <- .adjustLandscapePvalues(rawPvalues, p_adjust_method)

            # Dashed contour = significantly over-expressed region,
            # dotted = significantly under-expressed.
            if (perm_side %in% c("both", "over"))
                landgraphChart <- .significanceContour(landgraphChart,
                    pvalMatrix$over, i, sig_level, "dashed",
                    "over-expression")
            if (perm_side %in% c("both", "under"))
                landgraphChart <- .significanceContour(landgraphChart,
                    pvalMatrix$under, i, sig_level, "dotted",
                    "under-expression")
            }
        }

        # -- Peak labels on the plot --
        if (inference_unit == "region" && nrow(regions$summary)) {
            area_labels <- regions$summary
            # With a permutation test, label only the regions that pass
            # sig_level (the ones that are outlined); every region carrying
            # "p = 1.000" made the figure unreadable. Without a test the
            # regions are descriptive and keep their names.
            if ("PSpatial" %in% names(area_labels)) {
                area_labels <- area_labels[area_labels$PSpatial <= sig_level &
                    (perm_side == "both" | area_labels$Direction == perm_side), ,
                    drop = FALSE]
                area_labels$Label <- sprintf("%s\np = %.3f", area_labels$Region,
                                             area_labels$PSpatial)
            } else {
                area_labels$Label <- area_labels$Region
            }
            if (nrow(area_labels)) landgraphChart <- landgraphChart + ggplot2::geom_label(
                data = area_labels,
                ggplot2::aes(x = PeakRow, y = n + 1L - PeakCol, label = Label),
                size = 3, inherit.aes = FALSE)
        }
        if (inference_unit != "region" && !is.null(peakTable) && nrow(peakTable) > 0) {
            label_df <- data.frame(
                x     = peakTable$MatrixRow,
                y     = n + 1L - peakTable$MatrixCol,
                label = peakTable$NearestGene,
                stringsAsFactors = FALSE)
            if (requireNamespace("ggrepel", quietly = TRUE)) {
                landgraphChart <- landgraphChart +
                    ggrepel::geom_text_repel(
                        data = label_df,
                        aes(x = x, y = y, label = label),
                        colour = "white", size = 2.5, fontface = "bold",
                        box.padding = 0.3, max.overlaps = 20,
                        inherit.aes = FALSE)
            } else {
                landgraphChart <- landgraphChart +
                    geom_text(data = label_df,
                              aes(x = x, y = y, label = label),
                              colour = "white", size = 2.5, fontface = "bold",
                              inherit.aes = FALSE)
            }
        }

        if (.draw) methods::show(landgraphChart)

        # -- 3D surface (optional) --
        fig3d <- NULL
        if (isTRUE(plot3d)) {
            z_matrix <- ExpCtrl
            fig3d <- .buildSurface3D(z_matrix, .colorSet(setcolor), titleChart,
                pvalMatrix, i, sig_level, perm_side)
            if (!is.null(fig3d)) methods::show(fig3d)
        }

        leviResults[[k - 1]] <- structure(
            list(
                comparison = titleChart,
                landscape  = landgraphFinal,
                scores     = scoreTable,
                peaks      = peakTable,
                regions    = regions,
                pvalues    = pvalMatrix,
                plot       = landgraphChart,
                # The plotly surface used to be shown and then thrown away,
                # so there was no way to save or restyle it the way $plot
                # allows for the 2D figure. NULL when plot3d = FALSE.
                plot3d     = fig3d,
                raw_pvalues = rawPvalues,
                metadata = list(
                    signal_mode = signal_mode, single_col = single_col,
                    logfc_k = logfc_k, expressionLog = expressionLog && signal_mode == "ratio",
                    meaning = .signalMeaning(signal_mode, single_col),
                    nodes = nodesCoord, edges = edge_index,
                    node_coordinates = nodeCoordNorm,
                    node_signal = as.numeric(SignalOut[seq_len(nnodes), 1]),
                    edge_weighting = edge_weighting,
                    support_weights = support_weights,
                    grid = list(resolution = resolutionValue, zoom = zoomValue,
                                increase = increase, sigma = sigmaCells, occupancy = occFrac),
                    region_threshold = region_threshold,
                    region_min_cells = region_min_cells,
                    region_definition = "8-connected grid cells beyond neutral score",
                    missing_genes = as.character(naTotal),
                    n_perm = n_perm, sig_level = sig_level,
                    p_adjust_method = p_adjust_method,
                    inference_unit = inference_unit,
                    permutation_strata = resolved_strata,
                    inference = if (inference_unit == "region")
                        "node-label randomisation; maximum regional mass over both directions" else
                        "node-label randomisation; both tails over occupied grid cells",
                    rng_kind = RNGkind(), rng_state = rng_state,
                    versions = c(R = as.character(getRversion()),
                        levi = as.character(utils::packageVersion("levi")),
                        Rcpp = as.character(utils::packageVersion("Rcpp"))))
            ),
            class = "levi_result"
        )
    }

    if (length(leviResults) == 1L) {
        invisible(leviResults[[1]])
    } else {
        invisible(leviResults)
    }
}
