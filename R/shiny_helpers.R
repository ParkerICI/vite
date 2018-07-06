
#' @export
scaffold_GUI <- function(...) {
    shiny::runApp(appDir = file.path(system.file(package = "scgraphs"), "shinyGUI"), ...)
}