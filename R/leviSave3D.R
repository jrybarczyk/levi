#' Save the 3D landscape surface from a chosen viewpoint
#'
#' Writes the interactive surface returned by \code{\link{levi}} (field
#' \code{plot3d}) to a file, optionally after setting the camera position, so
#' that a figure can be reproduced from the same angle in every run. The
#' format is taken from the file extension.
#'
#' @param x A \code{levi_result} produced with \code{plot3d = TRUE}, or a
#'   plotly object such as \code{result$plot3d}.
#' @param file Output path. Supported extensions: \code{html} (interactive,
#'   keeps the view and can still be rotated), \code{png}, \code{jpeg},
#'   \code{jpg}, \code{webp}, \code{svg}, \code{pdf} (static, written by
#'   \code{plotly::save_image}, which requires the Python packages
#'   \code{kaleido} and \code{plotly}) and \code{tiff}/\code{tif} (rendered as PNG and converted
#'   with the \code{magick} package).
#' @param camera Optional camera specification, a list with any of
#'   \code{eye}, \code{center} and \code{up}, each a list with \code{x},
#'   \code{y} and \code{z}. \code{eye} is the camera position relative to the
#'   centre of the surface; larger values move it away. \code{NULL} keeps the
#'   camera already stored in the plot (plotly's default view otherwise).
#' @param width,height Image size in pixels for static formats.
#' @param scale Resolution multiplier for static raster formats.
#' @param selfcontained Logical; for \code{html}, embed all dependencies in a
#'   single file. Requires pandoc, and is turned off automatically when pandoc
#'   is not available.
#' @param compression TIFF compression passed to \code{magick::image_write};
#'   \code{"LZW"} by default.
#'
#' @return The path of the written file, invisibly.
#'
#' @details Static export goes through \code{plotly::save_image}, which drives
#' a headless browser via the \code{kaleido} and \code{plotly} Python
#' packages. Install them once with
#' \code{reticulate::py_install(c("kaleido", "plotly"))}. On this route a
#' 1600 x 1200 export at scale 2 takes about ten seconds. Without them, save the HTML
#' and use the camera icon of the plotly toolbar in a browser, or convert a
#' PNG exported that way with \code{magick::image_write}.
#'
#' @examples
#' res <- levi(expressionInput = system.file("extdata", "expression.dat",
#'                                           package = "levi"),
#'             fileTypeInput = "dat",
#'             networkCoordinatesInput = system.file("extdata", "medusa.dat",
#'                                                   package = "levi"),
#'             geneSymbolInput = "ID",
#'             readExpColumn = readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
#'             plot3d = TRUE)
#' cam <- list(eye = list(x = 1.6, y = -1.6, z = 0.9))
#' html <- tempfile(fileext = ".html")
#' leviSave3D(res, html, camera = cam, selfcontained = FALSE)
#' file.exists(html)
#' \donttest{
#' # Static formats need the kaleido and plotly Python packages
#' leviSave3D(res, tempfile(fileext = ".png"), camera = cam)
#' leviSave3D(res, tempfile(fileext = ".tiff"), camera = cam)  # also magick
#' }
#'
#' @seealso \code{\link{levi}}
#' @export
leviSave3D <- function(x, file, camera = NULL, width = 1600, height = 1200,
                       scale = 2, selfcontained = TRUE, compression = "LZW") {
    if (!requireNamespace("plotly", quietly = TRUE))
        stop("The 'plotly' package is required. Install with: ",
             "install.packages('plotly')", call. = FALSE)
    if (!is.character(file) || length(file) != 1L || !nzchar(file))
        stop("'file' must be a single file path.", call. = FALSE)

    fig <- if (inherits(x, "levi_result")) x$plot3d else x
    if (is.null(fig))
        stop("No 3D surface found. Run levi() with plot3d = TRUE.",
             call. = FALSE)
    if (!inherits(fig, "plotly"))
        stop("'x' must be a levi_result or a plotly object.", call. = FALSE)

    if (!is.null(camera)) {
        if (!is.list(camera) ||
            !all(names(camera) %in% c("eye", "center", "up")))
            stop("'camera' must be a list with elements among ",
                 "'eye', 'center' and 'up'.", call. = FALSE)
        fig <- plotly::layout(fig, scene = list(camera = camera))
    }

    ext <- tolower(tools::file_ext(file))
    static <- c("png", "jpeg", "jpg", "webp", "svg", "pdf")

    if (ext == "html") {
        if (!requireNamespace("htmlwidgets", quietly = TRUE))
            stop("The 'htmlwidgets' package is required for HTML output.",
                 call. = FALSE)
        sc <- isTRUE(selfcontained) &&
            requireNamespace("rmarkdown", quietly = TRUE) &&
            rmarkdown::pandoc_available()
        if (isTRUE(selfcontained) && !sc)
            message("pandoc not found: writing HTML with a dependency ",
                    "folder next to it.")
        htmlwidgets::saveWidget(fig, file, selfcontained = sc)
    } else if (ext %in% static) {
        .save3DStatic(fig, file, width, height, scale)
    } else if (ext %in% c("tiff", "tif")) {
        if (!requireNamespace("magick", quietly = TRUE))
            stop("The 'magick' package is required for TIFF output. ",
                 "Install with: install.packages('magick')", call. = FALSE)
        tmp <- tempfile(fileext = ".png")
        on.exit(unlink(tmp), add = TRUE)
        .save3DStatic(fig, tmp, width, height, scale)
        img <- magick::image_read(tmp)
        magick::image_write(img, file, format = "tiff",
                            compression = compression)
    } else {
        stop("Unsupported extension '", ext, "'. Use html, png, jpeg, ",
             "webp, svg, pdf or tiff.", call. = FALSE)
    }
    invisible(file)
}

# plotly::save_image needs the kaleido Python package; say so in plain words
# instead of letting the reticulate error surface.
.save3DStatic <- function(fig, file, width, height, scale) {
    py_ok <- function(m) requireNamespace("reticulate", quietly = TRUE) &&
        tryCatch(reticulate::py_module_available(m), error = function(e) FALSE)
    missing <- c("kaleido", "plotly")[!c(py_ok("kaleido"), py_ok("plotly"))]
    if (length(missing))
        stop("Static export needs the Python package(s) ",
             paste(sQuote(missing, FALSE), collapse = " and "),
             ". Install once with reticulate::py_install(c('kaleido', ",
             "'plotly')), or save as .html and use the camera icon of the ",
             "plotly toolbar.", call. = FALSE)
    plotly::save_image(fig, file, width = width, height = height,
                       scale = scale)
}
