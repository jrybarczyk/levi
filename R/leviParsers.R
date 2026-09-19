# Network parsers, one per accepted format. Kept apart from levi_function()
# because they share nothing with the rest of the pipeline: each takes the two
# input paths and gives back the node and edge tables.
#
# Every parser returns a list with:
#   nodes  data.frame V1 (name), V2, V3 (coordinates)
#   edges  two columns naming the endpoints of each interaction

.parseNetwork <- function(networkNodes, networkEdges, fileType) {
    parser <- switch(fileType,
        dat = .parseDat,
        stg = .parseStg,
        net = .parseNet,
        dyn = .parseDyn)
    parser(networkNodes, networkEdges)
}

# Medusa (DAT)
#   One file: an *edges section followed by a *nodes section carrying
#   "name  x  y".
.parseDat <- function(networkNodes, networkEdges) {
    nodes <- NULL
    edges <- NULL



    networkNodes <- read.delim(file = networkNodes,
    header = FALSE, sep = "\t",
    stringsAsFactors=FALSE, fill = TRUE, col.names =
    paste0("V",seq_len(max(count.fields(networkNodes,
    sep = '\t')))))


    delimiter <- which(networkNodes == "*nodes")
    if (length(delimiter) != 1L)
        stop("The network file has no '*nodes' section, so it is ",
            "not in Medusa (DAT) format. Check ",
            "'fileTypeInput': it is set to \"dat\".",
            call. = FALSE)

    # DAT format: row 1 is the "*edges" header; the edge list
    # runs from row 2 up to the row before "*nodes".
    n_edges <- delimiter - 2L
    if (n_edges < 1L)
        message("DAT parser: no edge found between *edges and ",
                "*nodes. Check the network file format.")
    edges <- dplyr::slice(networkNodes, seq.int(2L, delimiter - 1L))
    edges <- edges[,c(1,2)]
    nodes <- dplyr::slice(networkNodes,
        seq(delimiter+1, nrow(networkNodes)))
    nodes <- nodes[,c(1,2,3)]

    list(nodes = nodes, edges = edges)
}

# STRING / STITCH (STG)
#   Two inputs: a node table with coordinates and a separate interaction
#   table. Both accept a file path or a data.frame.
.parseStg <- function(networkNodes, networkEdges) {
    nodes <- NULL
    edges <- NULL


    if (is.data.frame(networkNodes) || is.matrix(networkNodes)) {
        nodes <- as.data.frame(networkNodes)
    } else {
        nodes <- read.delim(file = networkNodes, header = TRUE,
            sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    }
    if (is.data.frame(networkEdges) || is.matrix(networkEdges)) {
        edges <- as.data.frame(networkEdges)
    } else {
        edges <- read.delim(file = networkEdges, header = TRUE,
            sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    }
    edges <- edges[, c(1, 2)]
    nodes <- nodes[, c(1, 2, 3)]
    colnames(edges) <- c("V1", "V2")
    colnames(nodes) <- c("V1", "V2", "V3")

    list(nodes = nodes, edges = edges)
}

# Pajek (NET)
#   A *Vertices block with "id label x y" followed by an *Edges block.
.parseNet <- function(networkNodes, networkEdges) {
    nodes <- NULL
    edges <- NULL


    net_read <- read.delim(file = networkNodes, header = FALSE,
        stringsAsFactors=FALSE)

    delimiter_edge <- which(net_read == "*Edges")
    if (length(delimiter_edge) != 1L)
        stop("The network file has no '*Edges' section, so it is ",
            "not in Pajek (NET) format. Check 'fileTypeInput': ",
            "it is set to \"net\".", call. = FALSE)
    edges <- data.frame()

    edges_sl <- as.data.frame(dplyr::slice(net_read,
        seq(delimiter_edge+1, nrow(net_read))))

    delimiter_nodes_end <- which(net_read == "*Edges")
    nodes_sl <- as.data.frame(dplyr::slice(net_read,
        seq.int(3L, delimiter_nodes_end - 1L)))
    nodes <- data.frame()
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

    list(nodes = nodes, edges = edges)
}

# RedeR (DYN)
#   A zip holding XML; node labels and coordinates come from attributes.
.parseDyn <- function(networkNodes, networkEdges) {
    nodes <- NULL
    edges <- NULL


    tf <- tempfile(tmpdir = tdir <- tempdir())
    dyn_files <- unzip(networkNodes, exdir = tdir)
    dyn_read <- read_xml(dyn_files , stringsAsFactors=FALSE)


    dyn_label <- xml_find_all(dyn_read,
    xpath = "//*/*/@label")
    vals <- trimws(xml_text(dyn_label))
    dyn_df <- as.data.frame(vals, stringsAsFactors = FALSE)

    dyn_id <- xml_find_all(dyn_read, xpath = "//*/*/@id")
    vals_id <- trimws(xml_text(dyn_id))
    nodes <- as.data.frame(dplyr::slice(dyn_df, seq_len(length(vals_id))))
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
    edges$V1 <- as.numeric(edges$V1)
    edges$V2 <- as.numeric(edges$V2)

    t1<- merge(edges, nodes, by.x = "V1", by.y = "V1")
    colnames(t1) <- c("a", "b", "c")
    t1<- merge(t1, nodes, by.x = "b", by.y = "V1",
    all.x = FALSE)
    edges <- as.matrix(t1[,c(3,4)])
    colnames(edges) <- c("V1", "V2")

    nodes <- as.data.frame(cbind(nodes[,1], t(datax),
    t(datay)), stringsAsFactors = FALSE)

    nodes$V1 <- as.character(nodes$V1)
    nodes$V2 <- as.numeric(nodes$V2)
    nodes$V3 <- as.numeric(nodes$V3)

    list(nodes = nodes, edges = edges)
}
