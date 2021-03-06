% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_helpers.R
\name{vite_GUI}
\alias{vite_GUI}
\title{Start an analysis from the GUI}
\usage{
vite_GUI(...)
}
\arguments{
\item{...}{Additional arguments passed to shiny::runApp}
}
\description{
When the GUI starts you will be prompted to select a working directory. 
This directory must contain all the files that you want to include in the analysis.
Select any file in that directory, and the directory that contains the file will be selected as 
working directory
}
