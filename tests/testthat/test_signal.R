library(levi)

# The landscape score is documented as living in [0, 1] with 0.5 meaning "no
# change". Until 2.0.0 nothing enforced that: ratio mode on signed data
# inverted the ranking without a word, and Test + Control == 0 returned Inf.

genes <- c("HUB", paste0("N", 1:8))

run_hub <- function(expr, ...) {
    levi(networkCoordinatesInput = system.file("extdata", "hub_network.dat",
             package = "levi"),
         expressionInput = expr,
         fileTypeInput   = "dat",
         geneSymbolInput = "ID",
         readExpColumn   = readExpColumn("Test-Control"),
         resolutionValueInput = 10,
         smoothValueInput     = 5, ...)
}

# A differential-expression table: up in the core, down in the corners.
logfc_table <- function() {
    data.frame(ID   = genes,
               Test = c(4.3, 4.1, 4.4, 4.2, 4.3, -5.3, -5.1, -5.4, -5.2),
               Control = rep(0, 9))
}


test_that("ratio mode warns when the data carry negative values", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expect_warning(run_hub(logfc_table()), "expects non-negative values")
    expect_warning(run_hub(logfc_table()), "signal_mode = \"logfc\"")
})

test_that("logfc mode on the same data raises nothing", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expect_no_warning(run_hub(logfc_table(), signal_mode = "logfc"))
})

test_that("logfc mode separates up from down, which ratio mode cannot", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    res <- run_hub(logfc_table(), signal_mode = "logfc")
    sc <- setNames(res$scores$LandscapeScore, res$scores$Gene)

    # Genes with a positive logFC belong above the neutral point and the
    # negative ones below. Getting this backwards is the failure the warning
    # above exists to prevent.
    expect_true(all(sc[c("HUB", "N1", "N2", "N3", "N4")] > 0.5))
    expect_true(all(sc[c("N5", "N6", "N7", "N8")] < 0.5))
})

test_that("Test + Control == 0 does not escape the [0, 1] range", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expr <- data.frame(ID = genes,
                       Test    = c(rep(5, 5), rep(-5, 4)),
                       Control = c(rep(-5, 5), rep(5, 4)))

    expect_warning(res <- run_hub(expr), "Test \\+ Control == 0")

    s <- res$scores$LandscapeScore
    expect_true(all(is.finite(s)))
    expect_true(all(s >= 0 & s <= 1))
    expect_true(all(s == 0.5))     # every value undefined -> all neutral
})

test_that("every signal mode keeps the score inside [0, 1]", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    for (mode in c("ratio", "logfc", "zscore")) {
        res <- suppressWarnings(run_hub(logfc_table(), signal_mode = mode))
        s <- res$scores$LandscapeScore

        expect_true(all(is.finite(s)),
                    info = paste("non-finite score in mode", mode))
        expect_true(all(s >= 0 & s <= 1),
                    info = paste("score outside [0, 1] in mode", mode))

        z <- res$landscape$z
        expect_true(all(z >= 0 & z <= 1, na.rm = TRUE),
                    info = paste("landscape outside [0, 1] in mode", mode))
    }
})

test_that("duplicated identifiers are reported, not just averaged", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expr <- data.frame(ID = c(genes, "HUB"),
                       Test    = c(rep(200, 5), rep(5, 4), 5),
                       Control = c(rep(10, 5), rep(200, 4), 200))

    expect_message(res <- run_hub(expr), "1 duplicated identifier")

    # The two contradictory rows average out to the neutral point.
    hub <- res$scores$LandscapeScore[res$scores$Gene == "HUB"]
    expect_true(abs(hub - 0.5) < 0.1)
})

test_that("clean data raise neither warning nor message", {
    pdf(NULL); on.exit(dev.off(), add = TRUE)

    expect_no_warning(
        expect_no_message(
            run_hub(system.file("extdata", "hub_expression.dat",
                                package = "levi"))))
})
