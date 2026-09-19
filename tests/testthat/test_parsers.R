library(levi)

# The four network parsers. Until 2.0.0 only the dat path was reachable from
# the test suite, because extdata ships no Pajek or RedeR file; net and dyn were
# refactored with nothing watching them. The fixtures below are built here so
# every format is exercised.

fixture_dir <- function() {
    d <- file.path(tempdir(), "levi_parser_fixtures")
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    d
}

pajek_file <- function() {
    f <- file.path(fixture_dir(), "toy.net")
    writeLines(c("*Network toy", "*Vertices 3",
                 "1 A 0.10 0.20", "2 B 0.50 0.60", "3 C 0.90 0.30",
                 "*Edges", "1 2 1", "2 3 1"), f)
    f
}

reder_file <- function() {
    d <- fixture_dir()
    writeLines(c("<graph>",
                 "<node id=\"0\" label=\"A\" x=\"0.1\" y=\"0.2\"/>",
                 "<node id=\"1\" label=\"B\" x=\"0.5\" y=\"0.6\"/>",
                 "<node id=\"2\" label=\"C\" x=\"0.9\" y=\"0.3\"/>",
                 "<edge source=\"0\" target=\"1\"/>",
                 "<edge source=\"1\" target=\"2\"/>",
                 "</graph>"), file.path(d, "toy.xml"))

    f <- file.path(d, "toy.dyn")
    if (file.exists(f)) unlink(f)
    old <- setwd(d); on.exit(setwd(old), add = TRUE)
    utils::zip(f, "toy.xml", flags = "-q")
    f
}


test_that("dat parser reads the Medusa sections", {
    out <- levi:::.parseNetwork(
        system.file("extdata", "hub_network.dat", package = "levi"),
        NA, "dat")

    expect_named(out, c("nodes", "edges"))
    expect_equal(nrow(out$nodes), 9)
    expect_equal(nrow(out$edges), 8)
    expect_true("HUB" %in% out$nodes[[1]])
})

test_that("stg parser accepts data.frames of nodes and interactions", {
    nodes <- data.frame(name = c("A", "B", "C"),
                        x = c(0.1, 0.5, 0.9), y = c(0.2, 0.6, 0.3),
                        stringsAsFactors = FALSE)
    edges <- data.frame(from = c("A", "B"), to = c("B", "C"),
                        stringsAsFactors = FALSE)

    out <- levi:::.parseNetwork(nodes, edges, "stg")

    expect_equal(nrow(out$nodes), 3)
    expect_equal(nrow(out$edges), 2)
})

test_that("net parser reads a Pajek file", {
    skip_if_not_installed("dplyr")

    out <- levi:::.parseNetwork(pajek_file(), NA, "net")

    expect_equal(nrow(out$nodes), 3)
    expect_equal(nrow(out$edges), 2)
    expect_true(all(c("A", "B", "C") %in% unlist(out$nodes)))
})

test_that("dyn parser reads a RedeR archive", {
    skip_if_not_installed("xml2")
    skip_if(Sys.which("zip") == "", "the zip binary is needed to build the fixture")

    out <- levi:::.parseNetwork(reder_file(), NA, "dyn")

    expect_equal(nrow(out$nodes), 3)
    expect_equal(nrow(out$edges), 2)
})

test_that("parsers reject a file that is not in the declared format", {
    dat <- system.file("extdata", "hub_network.dat", package = "levi")

    expect_error(levi:::.parseNetwork(dat, NA, "net"), "no '\\*Edges' section")
    expect_error(levi:::.parseNetwork(pajek_file(), NA, "dat"),
                 "no '\\*nodes' section")
})
