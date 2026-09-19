library(levi)

test_that("LEVIui: browser defaults to FALSE", {
    expect_identical(formals(LEVIui)$browser, FALSE)
})

test_that("LEVIui: a non-logical browser is rejected with an error", {
    expect_error(LEVIui("yes"), "single TRUE or FALSE")
    expect_error(LEVIui(1), "single TRUE or FALSE")
    expect_error(LEVIui(NULL), "single TRUE or FALSE")
})

test_that("LEVIui: NA and vectors of length != 1 are rejected", {
    expect_error(LEVIui(NA), "single TRUE or FALSE")
    expect_error(LEVIui(c(TRUE, FALSE)), "single TRUE or FALSE")
    expect_error(LEVIui(logical(0)), "single TRUE or FALSE")
})

test_that("LEVIui: a failure inside the app is not reported as a bad argument", {
    # Previously a blanket tryCatch turned every error -- a missing
    # dependency, a busy port, a fault in app.R -- into the message
    # "Parameter must be TRUE or FALSE", sending the user to the wrong place.
    msg <- testthat::with_mocked_bindings(
        tryCatch(LEVIui(TRUE), error = function(e) conditionMessage(e)),
        runApp = function(...) stop("Port 8100 is already in use"),
        .package = "shiny")

    expect_match(msg, "Port 8100 is already in use", fixed = TRUE)
    expect_false(grepl("TRUE or FALSE", msg, fixed = TRUE))
})

test_that("LEVIui: browser selects whether launch.browser is set", {
    seen <- NULL
    capture_args <- function(...) {
        seen <<- list(...)
        invisible(NULL)
    }

    testthat::with_mocked_bindings(LEVIui(browser = TRUE),
        runApp = capture_args, .package = "shiny")
    expect_true(isTRUE(seen$launch.browser))

    # With browser = FALSE the argument is left unset on purpose, so that
    # RStudio can open the app in its Viewer pane.
    testthat::with_mocked_bindings(LEVIui(browser = FALSE),
        runApp = capture_args, .package = "shiny")
    expect_false("launch.browser" %in% names(seen))
})
