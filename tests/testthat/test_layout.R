library(levi)

# levi draws over the plane the network occupies, so it needs a coordinate per
# node. Until now the only way to get one without supplying it by hand was
# leviFromSTRING(), which needs the STRING database and network access.
# leviFromEdges() takes the interactions a user already has.

star_edges <- function() {
    data.frame(
        from = c("HUB", "HUB", "HUB", "HUB", "N1", "N2", "N3", "N4"),
        to   = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"),
        stringsAsFactors = FALSE)
}

star_expression <- function() {
    data.frame(
        ID      = c("HUB", paste0("N", 1:8)),
        Test    = c(200, 200, 200, 200, 200, 5, 5, 5, 5),
        Control = c(10, 10, 10, 10, 10, 200, 200, 200, 200))
}


test_that("leviFromEdges: returns nodes and edges levi can consume", {
    skip_if_not_installed("igraph")

    set.seed(1)
    net <- suppressMessages(leviFromEdges(star_edges()))

    expect_named(net, c("nodes", "edges", "graph"))
    expect_identical(names(net$nodes), c("name", "x", "y"))
    expect_identical(names(net$edges), c("V1", "V2"))

    expect_equal(nrow(net$nodes), 9)
    expect_equal(nrow(net$edges), 8)
    expect_setequal(net$nodes$name, c("HUB", paste0("N", 1:8)))
})

test_that("leviFromEdges: coordinates land in the range levi expects", {
    skip_if_not_installed("igraph")

    for (lay in c("fr", "kk", "circle")) {
        set.seed(3)
        net <- suppressMessages(leviFromEdges(star_edges(), layout = lay))

        expect_true(all(net$nodes$x >= 1 & net$nodes$x <= 100),
                    info = paste("x out of range for layout", lay))
        expect_true(all(net$nodes$y >= 1 & net$nodes$y <= 100),
                    info = paste("y out of range for layout", lay))
        expect_true(all(is.finite(net$nodes$x)), info = lay)
    }
})

test_that("leviFromEdges: the layout feeds a landscape that reads correctly", {
    skip_if_not_installed("igraph")
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    set.seed(42)
    net <- suppressMessages(leviFromEdges(star_edges(), layout = "kk"))

    res <- levi(networkCoordinatesInput  = net$nodes,
                networkInteractionsInput = net$edges,
                fileTypeInput            = "stg",
                expressionInput          = star_expression(),
                geneSymbolInput          = "ID",
                readExpColumn            = readExpColumn("Test-Control"),
                resolutionValueInput     = 20,
                smoothValueInput         = 50)

    sc <- setNames(res$scores$LandscapeScore, res$scores$Gene)

    expect_equal(nrow(res$scores), 9)
    # The core is over-expressed and the outer ring repressed. Whatever the
    # layout does with the positions, that has to survive.
    expect_true(all(sc[c("HUB", "N1", "N2", "N3", "N4")] > 0.5))
    expect_true(all(sc[c("N5", "N6", "N7", "N8")] < 0.5))
})

test_that("leviFromEdges: accepts a matrix and an igraph graph", {
    skip_if_not_installed("igraph")

    m <- as.matrix(star_edges())
    set.seed(5)
    from_matrix <- suppressMessages(leviFromEdges(m, layout = "circle"))

    g <- igraph::graph_from_edgelist(m, directed = FALSE)
    set.seed(5)
    from_graph <- suppressMessages(leviFromEdges(g, layout = "circle"))

    # circle is deterministic, so both routes have to agree exactly.
    expect_equal(from_matrix$nodes, from_graph$nodes)
    expect_equal(nrow(from_matrix$edges), nrow(from_graph$edges))
})

test_that("leviFromEdges: uses the columns it is told to", {
    skip_if_not_installed("igraph")

    wide <- data.frame(score = 1:8, from = star_edges()$from,
                       to = star_edges()$to, stringsAsFactors = FALSE)

    set.seed(2)
    net <- suppressMessages(leviFromEdges(wide, cols = c("from", "to"),
                                          layout = "circle"))
    expect_setequal(net$nodes$name, c("HUB", paste0("N", 1:8)))

    # Pointing at the wrong columns builds a different graph, not a silent
    # wrong one: the score column becomes nine distinct nodes.
    set.seed(2)
    wrong <- suppressMessages(leviFromEdges(wide, cols = c("score", "from"),
                                            layout = "circle"))
    expect_false(setequal(wrong$nodes$name, net$nodes$name))
})

test_that("leviFromEdges: reports interactions with a missing endpoint", {
    skip_if_not_installed("igraph")

    holed <- star_edges()
    holed$to[2] <- NA

    expect_message(net <- leviFromEdges(holed, layout = "circle"),
                   "1 interaction\\(s\\) with a missing endpoint")
    expect_equal(nrow(net$edges), 7)
})

test_that("leviFromEdges: rejects input it cannot make a network from", {
    skip_if_not_installed("igraph")

    expect_error(leviFromEdges("not a network"), "data.frame")
    expect_error(leviFromEdges(data.frame(only_one = 1:3)),
                 "at least two columns")
    expect_error(leviFromEdges(star_edges(), cols = 1),
                 "exactly two columns")
    expect_error(
        suppressMessages(leviFromEdges(data.frame(a = NA, b = NA))),
        "No usable interaction")
    expect_error(leviFromEdges(star_edges(), layout = "spiral"), "'arg'")
})

test_that("leviFromEdges: a force-directed layout is reproducible under a seed", {
    skip_if_not_installed("igraph")

    set.seed(11)
    a <- suppressMessages(leviFromEdges(star_edges(), layout = "fr"))
    set.seed(11)
    b <- suppressMessages(leviFromEdges(star_edges(), layout = "fr"))

    expect_equal(a$nodes, b$nodes)
})

test_that("leviFromEdges: handles graphs whose vertices carry no names", {
    skip_if_not_installed("igraph")

    # graph_from_edgelist() always names its vertices, but a graph built from
    # an adjacency matrix or with make_empty_graph() may not, and igraph
    # returns NULL for the missing attribute. That used to fail with
    # "arguments imply differing number of rows".
    g <- igraph::make_empty_graph(3, directed = FALSE)
    out <- suppressMessages(leviFromEdges(g))

    expect_equal(nrow(out$nodes), 3)
    expect_identical(out$nodes$name, c("1", "2", "3"))
    expect_true(all(is.finite(c(out$nodes$x, out$nodes$y))))

    m <- matrix(c(0, 1, 0,
                  1, 0, 1,
                  0, 1, 0), nrow = 3, byrow = TRUE)
    from_adj <- suppressMessages(
        leviFromEdges(igraph::graph_from_adjacency_matrix(m,
                                                          mode = "undirected")))
    expect_equal(nrow(from_adj$nodes), 3)
    expect_equal(nrow(from_adj$edges), 2)
})

test_that("leviFromEdges: survives the awkward shapes of real networks", {
    skip_if_not_installed("igraph")

    shapes <- list(
        disconnected = data.frame(a = c("A", "C"), b = c("B", "D")),
        self_loop    = data.frame(a = c("A", "B"), b = c("A", "B")),
        duplicated   = data.frame(a = rep("A", 3), b = rep("B", 3)),
        single_edge  = data.frame(a = "A", b = "B"),
        factors      = data.frame(a = factor("A"), b = factor("B")))

    for (nm in names(shapes)) {
        out <- suppressMessages(leviFromEdges(shapes[[nm]], layout = "circle"))
        expect_true(all(is.finite(c(out$nodes$x, out$nodes$y))), info = nm)
        expect_true(all(out$nodes$x >= 1 & out$nodes$x <= 100), info = nm)
        expect_gt(nrow(out$nodes), 0)
    }
})
