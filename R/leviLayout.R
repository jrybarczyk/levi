# Layout computation, shared by leviFromSTRING() and leviFromEdges().
#
# levi draws a landscape over the plane the network occupies, so it needs a
# coordinate per node. When the network arrives as a bare list of interactions
# there are none, and this is what supplies them.

# Positions the vertices of a graph and scales the result to [1, 100], the
# positive range levi expects.
#
# Returns a data.frame with name, x and y, in vertex order.
.layoutNodes <- function(g, layout = c("fr", "kk", "lgl", "dh", "circle"),
                         names = NULL) {
    layout <- match.arg(layout)

    lay_fn <- switch(layout,
        fr     = igraph::layout_with_fr,
        kk     = igraph::layout_with_kk,
        lgl    = igraph::layout_with_lgl,
        dh     = igraph::layout_with_dh,
        circle = igraph::layout_in_circle
    )
    coords <- lay_fn(g)

    # Scale to [1, 100] - levi expects positive coordinates. A graph whose
    # layout collapses on one axis (a single node, or a perfect line) would
    # divide by zero, so those land in the middle.
    scale_01 <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (diff(rng) == 0) return(rep(50, length(x)))
        (x - rng[1]) / diff(rng) * 99 + 1
    }

    # A graph built from an edge list always names its vertices, but one built
    # from an adjacency matrix or with make_empty_graph() may not, and the
    # missing attribute comes back as NULL. Fall back to the vertex index so
    # every node still gets a label.
    vnames <- if (!is.null(names)) {
        names
    } else if (!is.null(igraph::V(g)$name)) {
        igraph::V(g)$name
    } else {
        as.character(seq_len(igraph::vcount(g)))
    }

    data.frame(
        name = vnames,
        x    = scale_01(coords[, 1]),
        y    = scale_01(coords[, 2]),
        stringsAsFactors = FALSE
    )
}


#' @title leviFromEdges
#'
#' @description Build a levi-ready network from a plain list of interactions.
#' \code{levi()} draws its landscape over the plane the network occupies, so it
#' needs a coordinate for every node. When all you have is a list of pairs --
#' which is what most interaction databases and most collaborators hand over --
#' this computes a layout with \pkg{igraph} and returns the node and edge
#' tables ready to pass straight to \code{levi()}.
#'
#' @param edges A \code{data.frame} or \code{matrix} whose first two columns
#'   name the endpoints of each interaction, or an \pkg{igraph} graph. Further
#'   columns are ignored.
#' @param layout Character. Layout algorithm: \code{"fr"} (Fruchterman-Reingold,
#'   the default), \code{"kk"} (Kamada-Kawai), \code{"lgl"}, \code{"dh"} or
#'   \code{"circle"}.
#' @param directed Logical. Whether to treat the pairs as directed when
#'   building the graph. The layout ignores direction either way. Default
#'   \code{FALSE}.
#' @param cols Which two columns of \code{edges} hold the endpoints, by
#'   position or by name. Default is the first two.
#'
#' @return Invisibly, a list with:
#' \describe{
#'   \item{nodes}{data.frame with \code{name}, \code{x} and \code{y}, the
#'     coordinates scaled to \code{[1, 100]}.}
#'   \item{edges}{data.frame with \code{V1} and \code{V2}, the endpoints.}
#'   \item{graph}{the \pkg{igraph} object the layout was computed on.}
#' }
#' Pass \code{nodes} to \code{networkCoordinatesInput}, \code{edges} to
#' \code{networkInteractionsInput} and \code{"stg"} to \code{fileTypeInput}.
#'
#' @details Force-directed layouts are stochastic: \code{"fr"}, \code{"kk"},
#' \code{"lgl"} and \code{"dh"} give a different arrangement on every call.
#' Set \code{set.seed()} beforehand to reproduce a figure. \code{"circle"} is
#' deterministic.
#'
#' The layout is an analytical choice, not only an aesthetic one: the same
#' expression values arranged differently produce a different landscape,
#' because levi reads neighbourhoods. Where a curated layout exists, prefer it
#' over a computed one.
#'
#' @seealso \code{\link{levi}}, \code{\link{leviFromSTRING}}
#'
#' @examples
#' # A small network given only as interactions, with no coordinates.
#' interactions <- data.frame(
#'     from = c("HUB", "HUB", "HUB", "HUB", "N1", "N2", "N3", "N4"),
#'     to   = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"),
#'     stringsAsFactors = FALSE)
#'
#' set.seed(42)
#' net <- leviFromEdges(interactions, layout = "kk")
#' head(net$nodes)
#'
#' expression <- data.frame(
#'     ID      = c("HUB", paste0("N", 1:8)),
#'     Test    = c(200, 200, 200, 200, 200, 5, 5, 5, 5),
#'     Control = c(10, 10, 10, 10, 10, 200, 200, 200, 200))
#'
#' res <- levi(
#'     networkCoordinatesInput  = net$nodes,
#'     networkInteractionsInput = net$edges,
#'     fileTypeInput            = "stg",
#'     expressionInput          = expression,
#'     geneSymbolInput          = "ID",
#'     readExpColumn            = readExpColumn("Test-Control"),
#'     resolutionValueInput     = 20,
#'     smoothValueInput         = 50)
#'
#' head(res$scores)
#'
#' @author Jose Rybarczyk Filho (jose.luiz@@unesp.br)
#'
#' @export
leviFromEdges <- function(edges,
                          layout   = c("fr", "kk", "lgl", "dh", "circle"),
                          directed = FALSE,
                          cols     = seq_len(2L)) {
    layout <- match.arg(layout)

    if (!requireNamespace("igraph", quietly = TRUE))
        stop("Package 'igraph' is required. ",
             "Install with: install.packages('igraph')", call. = FALSE)

    if (inherits(edges, "igraph")) {
        g <- edges
    } else {
        if (!(is.data.frame(edges) || is.matrix(edges)))
            stop("'edges' must be a data.frame, a matrix or an igraph graph.",
                 call. = FALSE)

        edges <- as.data.frame(edges, stringsAsFactors = FALSE)
        if (ncol(edges) < 2L)
            stop("'edges' needs at least two columns naming the endpoints ",
                 "of each interaction; got ", ncol(edges), ".", call. = FALSE)
        if (length(cols) != 2L)
            stop("'cols' must name exactly two columns.", call. = FALSE)

        el <- cbind(as.character(edges[[cols[1]]]),
                    as.character(edges[[cols[2]]]))

        keep <- stats::complete.cases(el) & nzchar(el[, 1]) & nzchar(el[, 2])
        if (!all(keep))
            message(sum(!keep), " interaction(s) with a missing endpoint ",
                    "were dropped.")
        el <- el[keep, , drop = FALSE]

        if (nrow(el) == 0L)
            stop("No usable interaction left in 'edges'.", call. = FALSE)

        g <- igraph::graph_from_edgelist(el, directed = directed)
    }

    if (igraph::vcount(g) == 0L)
        stop("The network has no nodes.", call. = FALSE)

    nodes <- .layoutNodes(g, layout)

    el_out <- igraph::as_edgelist(g)
    out_edges <- data.frame(V1 = el_out[, 1], V2 = el_out[, 2],
                            stringsAsFactors = FALSE)

    message("Network laid out with '", layout, "': ", igraph::vcount(g),
            " nodes, ", igraph::ecount(g), " edges.")

    invisible(list(nodes = nodes, edges = out_edges, graph = g))
}
