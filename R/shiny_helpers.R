

#' Start an analysis from the GUI
#'
#' When the GUI starts you will be prompted to select a working directory. Select any file
#' in that directory, and the directory that contains the file will be selected.
#'
#' @param ... Additional arguments passed to shiny::runApp
#' @export
#'
scgraphs_GUI <- function(...) {
    shiny::runApp(appDir = file.path(system.file(package = "scgraphs"), "shinyGUI"), ...)
}