
#' Get the columns that are common to a set of input tabular files
#'
#' @param files.list A vector of input file names. Each file should be a tab-separated table, with the first row
#'   representing column headers
#'
#' @return Returns a vector of column names that are present in all the files in \code{files.list}

get_common_columns <- function(files.list) {
    l <- list()
    for(f in files.list) {
        temp <- read.table(f, header = T, sep = "\t", check.names = F, quote = "", nrows = 1)
        l <- c(l, list(names(temp)))
    }
    return(Reduce(intersect, l))
}



