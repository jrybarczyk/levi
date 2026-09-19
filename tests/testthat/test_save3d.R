library(levi)

# leviSave3D(): the HTML path is exercised here; static export needs kaleido,
# which the test machines do not have, so its guard is what gets tested.

make_result <- function() {
    skip_if_not_installed("plotly")
    levi(expressionInput = system.file("extdata", "expression.dat",
                                       package = "levi"),
         fileTypeInput = "dat",
         networkCoordinatesInput = system.file("extdata", "medusa.dat",
                                               package = "levi"),
         geneSymbolInput = "ID",
         readExpColumn = readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
         plot3d = TRUE, n_perm = 0)
}

test_that("leviSave3D writes HTML and stores the camera", {
    skip_if_not_installed("htmlwidgets")
    pdf(NULL); on.exit(dev.off())
    res <- suppressMessages(make_result())
    f <- tempfile(fileext = ".html")
    cam <- list(eye = list(x = 1.6, y = -1.6, z = 0.9))
    out <- suppressMessages(leviSave3D(res, f, camera = cam,
                                       selfcontained = FALSE))
    expect_identical(out, f)
    expect_true(file.exists(f))
    fig <- plotly::layout(res$plot3d, scene = list(camera = cam))
    built <- plotly::plotly_build(fig)
    expect_equal(built$x$layout$scene$camera$eye$x, 1.6)
})

test_that("leviSave3D validates its inputs", {
    pdf(NULL); on.exit(dev.off())
    res <- suppressMessages(make_result())
    expect_error(leviSave3D(res, tempfile(fileext = ".bmp")),
                 "Unsupported extension")
    expect_error(leviSave3D(res, tempfile(fileext = ".html"),
                            camera = list(zoom = 2)), "'camera'")
    res$plot3d <- NULL
    expect_error(leviSave3D(res, tempfile(fileext = ".html")),
                 "plot3d = TRUE")
    expect_error(leviSave3D(list(), tempfile(fileext = ".html")),
                 "plot3d = TRUE|plotly object")
})

test_that("leviSave3D explains the kaleido requirement when it is missing", {
    pdf(NULL); on.exit(dev.off())
    res <- suppressMessages(make_result())
    has_kaleido <- requireNamespace("reticulate", quietly = TRUE) &&
        tryCatch(reticulate::py_module_available("kaleido"),
                 error = function(e) FALSE)
    skip_if(has_kaleido, "kaleido installed; static export would succeed")
    expect_error(leviSave3D(res, tempfile(fileext = ".png")), "kaleido")
})
