


# Add explanation of what happens in the different cases with or without the saample column, and the pooled issue
#' Write clustering output
#'
#' @param clustered.data A \code{data.frame} containing the clustered data. Each row corresponds to cell. The \code{data.frame}
#'   must include a column called \code{cellType}, indicating cluster membership
#' @param base.name The base name for naming output files. This is only used in two cases:
#'   \itemize{
#'     \item If \code{clustered_data} does not contain a column called \code{sample}
#'     \item If \code{clustered_data} represents a pooled clustering result (e.g. derived from \code{grappolo::cluster_fcs_files_groups}),
#'       the pooled data will be written in a directory called pooled/\code{base.name}
#'   }
#' @param output.dir The output directory
#'
#' @export
write_clusters_data <- function(clustered.data, base.name, output.dir = "./") {

    cluster.data.dir <- file.path(output.dir, "clusters_data")
    write.pooled <- FALSE


    if(!is.null(clustered.data$sample)) {
        samples <- unique(clustered.data$sample)
        sapply(samples, function(x) {dir.create(file.path(cluster.data.dir, x), recursive = TRUE, showWarnings = FALSE)})
        if(length(samples) > 1) {
            dir.create(file.path(cluster.data.dir, "pooled", base.name), recursive = TRUE, showWarnings = FALSE)
            write.pooled <- TRUE
        }
    } else {
        dir.create(file.path(cluster.data.dir, base.name), recursive = TRUE, showWarnings = FALSE)
    }

    plyr::d_ply(clustered.data, ~cellType, function(x) {
        if(is.null(x$sample))
            saveRDS(x, file = file.path(cluster.data.dir, base.name, sprintf("c%d.rds", x$cellType[1])))
        else {
            if(write.pooled)
                saveRDS(x, file = file.path(cluster.data.dir, "pooled", base.name, sprintf("c%d.rds", x$cellType[1])))

            plyr::d_ply(x, ~sample, function(df) {
                saveRDS(df, file = file.path(cluster.data.dir, df$sample[1], sprintf("c%d.rds", df$cellType[1])))
            })


        }
    })

    return(invisible(NULL))
}


#' Write a graph as GraphML file
#'
#' This function writes a graph as a GraphML file. Please use this function when
#' writing graphs created using this package
#'
#' @param G The \code{igraph} object representing the graph
#' @param out.name The name of the output file
#' @param ... Additional arguments passed to \code{igraph::write.graph}
#'
#' @export
#'
write_graph <- function(G, out.name, ...) {
    G <- igraph::set.graph.attribute(G, "fname", basename(out.name))
    igraph::write.graph(G, out.name, format = "graphml", ...)
}



