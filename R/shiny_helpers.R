

#' Start an analysis from the GUI
#'
#' When the GUI starts you will be prompted to select a working directory. 
#' This directory must contain all the files that you want to include in the analysis.
#' Select any file in that directory, and the directory that contains the file will be selected as 
#' working directory 
#'
#' @param ... Additional arguments passed to shiny::runApp
#' @export
#'
vite_GUI <- function(...) {
    shiny::runApp(appDir = file.path(system.file(package = "vite"), "shinyGUI"), ...)
}