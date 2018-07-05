#' Calculate cosine similarity between a vector and the rows of a matrix
#'
#' @param x A numeric vector of length \code{P}
#' @param m An \code{NxP} matrix
#'
#' @return Returns a vector of length \code{N} containing the cosine similarity between
#'   the vector \code{x} and all the rows of \code{m}
#'
cosine_similarity_from_matrix <- function(x, m) {
    x <- x / as.vector(sqrt(crossprod(x)))
    return(as.vector((m %*% x) / sqrt(rowSums(m^2))))
}


#' Calculate cosine similarity between the rows of a matrix
#'
#' This function calculates the cosine similarity between the rows of an input matrix, according to the values
#' of the variables in the columns
#'
#' @param m An \code{N x P} matrix
#'
#' @return Returns an \code{N x N} matrix with the cosine similarity between the corresponding rows in \code{m}
#'
cosine_similarity_matrix <- function(m){
    ret <- m %*% t(m) / (sqrt(rowSums(m^2) %*% t(rowSums(m^2))))
    return(ret)
}


#' Builds a graph from a data.frame of nodes
#'
#' @param tab A \code{data.frame} of graph nodes. All the columns of \code{tab} will become vertex properties of the output graph
#' @param col.names A character vector indicating which columns of \code{tab} should be used to calculate distances
#' @param filtering_T The threshold used to filter edges in the graph. If this value is \code{< 1} the edges are filtered based on the
#'   cosine similarity value. If this is an integer \code{>= 1} they are filtered based on rank for each node (i.e. for each node only
#'   the edges with rank less than the threshold are retained)
#'
#' @return Returns an \code{igraph} object
#' @export
build_graph <- function(tab, col.names, filtering_T = 0.8) {
    m <- as.matrix(tab[, col.names])
    row.names(m) <- tab$cellType
    dd <- cosine_similarity_matrix(m)
    diag(dd) <- 0
    dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest


    # Both these filtering functions modify the dd matrix directly, so we don't need to reassign the return value

    if(filtering_T >= 1)
        filter_matrix_by_rank(dd, filtering_T)
    else
        filter_matrix(dd, filtering_T)

    G <- igraph::graph.adjacency(dd, mode = "undirected", weighted = T)

    for(i in names(tab))
        G <- igraph::set.vertex.attribute(G, name = i, value = tab[, i])

    return(G)
}





#' Builds an unsupervised force-directed graph
#'
#' @inheritParams build_graph
#' @param filtering.threshold The threshold used to filter edges in the graph
#'
#' @return Returns an \code{igraph} graph object. This function also performs graph-based clustering. The resulting cluster membership
#'   is contained in the \code{community_id} vertex attribute of the resulting graph
#'
#' @export
get_unsupervised_graph <- function(tab, col.names, filtering.threshold) {
    message("Building graph...")
    flush.console()
    G <- build_graph(tab, col.names, filtering_T = filtering.threshold)

    for(i in names(tab))
        G <- igraph::set.vertex.attribute(G, name = i, value = tab[, i])

    cc <- igraph::multilevel.community(G)
    V(G)$community_id <- as.character(cc$membership)
    V(G)$name <- seq_along(V(G)$name)
    V(G)$type <- "cluster"
    V(G)$Label <- paste("c", V(G)$cellType, sep = "")

    message("Running ForceAtlas2...")
    flush.console()
    G <- complete_forceatlas2(G, first.iter = 50000, overlap.method = NULL, ew.influence = 5)
    message("ForceAtlas2 done")
    flush.console()

    return(G)
}

#' Builds an unsupervised graph starting from a list of input files
#'
#' This function is similar to \code{get_unsupervised_graph}, except that the input in this case is a list of file names containing
#' the vertices of the graph. Each file should be a tab-separated table of nodes, similar to the \code{tab} input of
#' \code{get_unsupervised_graph}. This function will \code{rbind} all the files in \code{files.list}, and a build a single graph
#' with the union of all the nodes. Optionally, file-level metadata can also be incorporated in the vertex properties of the resulting
#' graph, using the \code{metadata.tab} parameter
#'
#' @param files.list The list of files to process. The function will first determine the set of columns that are common to all
#'   the files. Only the common vertex properties will feature in the output graph
#' @inheritParams get_unsupervised_graph
#' @param metadata.tab Optional. If specified, a table of file-level metadata, to be added as vertex properties in the graph. Each row
#'   should specify metadata for a single file, with the columns of \code{metadata.tab} representing metadata values. All the vertices
#'   derived from that file will have the corresponding metadatata value
#' @param metadata.filename.col The name of the column in \code{metadata.tab} that contains the file name to be matched to the files
#'   in \code{files.list}
#' @param use.basename The resulting graph will contain an additional vertex property called \code{sample} identifying which file
#'   the vertex was derived from. If \code{use.basename} is \code{TRUE} the \code{basename} will be used, otherwise the full path
#'   as specified in \code{files.list}. Moreover if \code{use.basename} is \code{TRUE} the matching with the file names contained
#'   in \code{metadata.tab} will be based on the \code{basename} only
#' @param process.clusters.data If this is \code{TRUE} this function will look for a file with extension \code{.all_events.rds} for
#'   each file in \code{files.list} (see the Documentation of \code{scfeatures::cluster_fcs_files}). This file contains single-cell
#'   data (i.e. each row represent a cell, as opposed to the files in \code{files.list} where each row represents a cluster). Each file
#'   will be processed using the \code{\link{write_clusters_data}} function. This processing is used for downstream data visualization
#'   but it is not strictly necessary to create the graph. If \code{use.basename} is \code{TRUE} the \code{basename} of the files in
#'   \code{files.list} will be used for processing.
#' @param downsample.to The target number of events for downsampling. Only used if \code{process.clusters.data == TRUE}. This is only
#'   used for downstream data visualization and does not affect the construction of the graph
#'
#' @return See the return value of \code{get_unsupervised_graph}
#'
#' @export
get_unsupervised_graph_from_files <- function(files.list, col.names, filtering.threshold,
                                              metadata.tab = NULL, metadata.filename.col = NULL, use.basename = TRUE, process.clusters.data = TRUE, downsample.to = 1000) {

    tab <- NULL
    ret <- list(graphs = list())

    common.cols <- scfeatures::get_common_columns(files.list, file.type = "txt")

    if(use.basename && !is.null(metadata.tab) && !is.null(metadata.filename.col))
        metadata.tab[, metadata.filename.col] <- basename(metadata.tab[, metadata.filename.col])

    message("Loading data...")
    flush.console()

    for(f in files.list) {
        temp <- read.table(f, header = T, sep = "\t", check.names = F, quote = "", stringsAsFactors = F)
        temp <- temp[, common.cols]
        temp$sample <- ifelse(use.basename, basename(f), f)

        if(!is.null(metadata.tab) && !is.null(metadata.filename.col))
            temp <- merge(temp, metadata.tab, by.x = "sample", by.y = metadata.filename.col, all.x = T)
        tab <- rbind(tab, temp)
    }

    G <- get_unsupervised_graph(tab, col.names, filtering.threshold)

    #Add downsampling option

    if(process.clusters.data) {
        message("Processing clusters data...")
        for(f in files.list) {
            rds.file <- gsub("txt$", "all_events.rds", f)
            tab <- readRDS(rds.file)
            tab <- downsample_by(tab, "cellType", downsample.to)
            if(use.basename)
                f <- basename(f)
            write_clusters_data(tab, f)
        }
    }

    return(G)
}


# Add explanation of what happens in the different cases with or without the saample column, and the pooled issue
#' Write clustering output
#'
#' @param clustered.data A \code{data.frame} containing the clustered data. Each row corresponds to cell. The \code{data.frame}
#'   must include a column called \code{cellType}, indicating cluster membership
#' @param base.name The base name for naming output files. This is only used if \code{clustered_data} does not contain
#'   a column called \code{sample}. Otherwise the \code{sample} names are used for this purpose
#' @param output.dir The output directory
#'
#' @export
write_clusters_data <- function(clustered.data, base.name, output.dir = "./") {

    cluster.data.dir <- file.path(output.dir, "clusters_data")
    write.pooled <- FALSE


    if(!is.null(clustered.data$sample)) {
        sapply(unique(clustered.data$sample), function(x) {dir.create(file.path(cluster.data.dir, x), recursive = TRUE, showWarnings = FALSE)})
        if(length(unique(clustered.data$sample)) > 1) {
            dir.create(file.path(cluster.data.dir, "pooled"), recursive = TRUE, showWarnings = FALSE)
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
                saveRDS(x, file = file.path(cluster.data.dir, "pooled", sprintf("c%d.rds", x$cellType[1])))

            plyr::d_ply(x, ~sample, function(df) {
                saveRDS(df, file = file.path(cluster.data.dir, df$sample[1], sprintf("c%d.rds", df$cellType[1])))
            })


        }
    })

    return(invisible(NULL))
}

