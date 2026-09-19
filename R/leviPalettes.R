# Colour palettes offered through levi()'s setcolor argument.
#
# The multicolour scales return 20 steps and the two-tone ones return 3, which
# is what scale_fill_gradientn() interpolates over.
#
# Split out of levi_function() during the 2.0.0 refactor. The old signature
# took a first argument that no function in the family ever read, and that was
# never defined in the caller: it survived only because R never forced it.

.colorSet <- function(colorType = c("default", "terrain", "rainbow","heat",
    "topo", "cm", "purple_pink", "green_blue", "blue_yellow", "pink_green",
    "orange_purple", "green_marine")) {
    colorType <- match.arg(colorType)
    defaultColors <- function(n) {
        c(
            "#180052", "#0c0083", "#0000b4", "#0000e4", "#0010ff", "#0041ff",
            "#0072ff", "#00A3FF", "#00D4FF", "#00FF49", "#5AFF00", "#FFE400",
            "#FFC400", "#FFA300", "#FF8300", "#FF6200", "#FF4100", "#FF2100",
            "#FF0000", "#E40000"
        )
    }
    purple_pink <- function(n) {
        c("#4e5052", "#b387e6", "#ff0000")
    }
    green_blue <- function(n) {
        c("#4e5052", "#a4db56", "#1d02c9")
    }
    blue_yellow <- function(n) {
        c("#4e5052", "#27b0cf", "#ffec2b")
    }
    orange_purple <- function(n){
        c("#4e5052", "#fcbb63", "#7300c4")
    }
    green_marine <- function(n){
        c("#4e5052", "#5df0b0", "#360d94")
    }
    pink_green <- function(n){
        c("#4e5052", "#e854d9", "#90db56")
    }
    color_list <- list(
        default = defaultColors,
        terrain = terrain.colors,
        rainbow = rainbow,
        heat = heat.colors,
        topo = topo.colors,
        cm = cm.colors,
        purple_pink = purple_pink,
        green_blue = green_blue,
        blue_yellow = blue_yellow,
        orange_purple = orange_purple,
        green_marine = green_marine,
        pink_green = pink_green
    )

    if ((colorType == "default") || (colorType == "terrain") ||
        (colorType == "rainbow") || (colorType == "heat") ||
        (colorType == "topo") || (colorType == "cm")) {
        colorSetRange <- color_list[[colorType]](20)
    } else {
        colorSetRange <- color_list[[colorType]](3)
    }
    return(colorSetRange)
}
