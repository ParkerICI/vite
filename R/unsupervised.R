



cosine_similarity_from_matrix <- function(x, m) {
    x <- x / sqrt(crossprod(x))
    return(as.vector((m %*% x) / sqrt(rowSums(m^2))))
}


#' Calculate cosine similarity between the rows of a matrix
#'
#' This function calculates the cosine similarity between the rows of an input matrix, according to the values
#' of the variables in the columns
#'
#' @param m An \code{N x P} matrix
#'
#' @return Returns a \code{N x N} matrix with the cosine similarity between the corresponding rows in \code{m}
#'
cosine_similarity_matrix <- function(m){
    ret <- m %*% t(m) / (sqrt(rowSums(m^2) %*% t(rowSums(m^2))))
    return(ret)
}


filter_similarity_matrix <- function(m, T) {
    ret <- t(apply(m, 1, function(x) {
        if(max(x) <= T)
            x[x < max(x)] <- 0
        else
            x[x < T] <- 0
        return(x)
    }))
    return(ret)
}

filter_similarity_matrix_by_rank <- function(m, T) {
    ret <- t(apply(m, 1, function(x) {
        r <- rank(x, ties.method = "first")
        r <- max(r) - r + 1
        x[r > T] <- 0
        return(x)
    }))
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

    if(filtering_T >= 1)
        dd <- filter_similarity_matrix_by_rank(dd, filtering_T)
    else
        dd <- filter_similarity_matrix(dd, filtering_T)

    G <- igraph::graph.adjacency(dd, mode = "undirected", weighted = T)
    n.vertices <- length(V(G))
    lay <- igraph::layout.kamada.kawai(G)
    colnames(lay) <- c("x", "y")
    G <- igraph::set.vertex.attribute(G, name = "x", value = lay[, "x"])
    G <- igraph::set.vertex.attribute(G, name = "y", value = lay[, "y"])
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
get_unsupervised_graph <- function(tab, col.names, filtering.threshold, output.name) {
    message("Building graph...")
    flush.console()
    G <- build_graph(tab, col.names, filtering_T = filtering.threshold)

    for(i in names(tab))
        G <- igraph::set.vertex.attribute(G, name = i, value = tab[, i])

    cc <- igraph::multilevel.community(G)
    V(G)$community_id <- as.character(cc$membership)
    V(G)$name <- seq_along(V(G)$name)

    message("Running ForceAtlas2...")
    flush.console()
    G <- complete_forceatlas2(G, first.iter = 50000, overlap.iter = 1, overlap_method = NULL, ew_influence = 5)
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
#' @param metadata.tab Optional. If specified, a table of file-level metadata, to be added as vertex properties in the graph. Each row
#'   should specify metadata for a single file, with the columns of \code{metadata.tab} representing metadata values. All the vertices
#'   derived from that file will have the corresponding metadatata value
#' @param metadata.filename.col The name of the column in \code{metadata.tab} that contains the file name to be matched to the files
#'   in \code{files.list}
#' @param use.basename The resulting graph will contain an additional vertex property called \code{file_name} identifying which file
#'   the vertex was derived from. If \code{use.basename} is \code{TRUE} the \code{basename} will be used, otherwise the full path
#'   as specified in \code{files.list}. Moreover if \code{use.basename} is \code{TRUE} the matching with the file names contained
#'   in \code{metadata.tab} will be based on the \code{basename} only
#' @inheritDotParams get_unsupervised_graph -tab
#'
#' @return See the return value of \code{get_unsupervised_graph}
#'
#' @export
get_unsupervised_graph_from_files <- function(files.list, metadata.tab = NULL, metadata.filename.col = NULL, use.basename = T, ...) {

    tab <- NULL
    ret <- list(graphs = list())

    common.cols <- get_common_columns(files.list)

    if(use.basename && !is.null(metadata.tab) && !is.null(metadata.filename.col))
        metadata.tab[, metadata.filename.col] <- basename(metadata.tab[, metadata.filename.col])

    message("Loading data...")
    flush.console()

    for(f in files.list) {
        temp <- read.table(f, header = T, sep = "\t", check.names = F, quote = "", stringsAsFactors = F)
        temp <- temp[, common.cols]
        temp$file_name <- ifelse(use.basename, basename(f), f)

        if(!is.null(metadata.tab) && !is.null(metadata.filename.col))
            temp <- merge(temp, metadata.tab, by.x = "file_name", by.y = metadata.filename.col, all.x = T)
        tab <- rbind(tab, temp)
    }

    G <- get_unsupervised_graph(tab, ...)
    return(G)

}
