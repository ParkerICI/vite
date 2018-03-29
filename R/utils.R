get_common_columns <- function(files.list) {
    l <- list()
    for(f in files.list) {
        temp <- read.table(f, header = T, sep = "\t", check.names = F, quote = "", nrows = 1)
        l <- c(l, list(names(temp)))
    }
    return(Reduce(intersect, l))
}



