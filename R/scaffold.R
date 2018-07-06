#' Downsample an input data.frame
#'
#' This function downsamples an input \code{data.frame}, so that for each value of \code{col.name}, there are
#' at most \code{size} rows.
#'
#' @param tab Input \code{data.frame}
#' @param col.name The name of a column in tab (typically the data in the column should either be characters or factors).
#'   For each value of \code{tab$col.name}, there will be at most \code{size} row in the output table
#' @param size The number of rows to select for each value of \code{tab$col.name}. If size is larger than the total number of rows
#'   available for a particular value, all rows for that value will be returned
#' @return Returns a \code{data.frame}
#'
downsample_by <- function(tab, col.name, size) {
    return(plyr::ddply(tab, col.name, function(x, size) {
        if(nrow(x) <= size)
            return(x)
        else
            return(x[sample(1:nrow(x), size),])
    }, size = size))
}




#' Load landmark populations from a list of FCS files
#'
#' @param files.list The list of files to load. The population names will be derived by splitting the file name
#'   using \code{"_"} as separator, and taking the last field
#' @param asinh.cofactor The cofactor to use for \code{asinh} transformation. If this is \code{NULL} no transformation
#'   is applied
#' @param ... Additional arguments passed to \code{flowCore::read.FCS}
#' @return Returns a list with the following elements:
#'   \itemize{
#'     \item{\code{landmarks.data}}: a \code{data.frame} containing the entirety of the data
#'     \item{\code{tab.landmarks}}: a \code{data.frame} with the data for the landmark nodes. Each row is a separate landmark
#'       population, and the columns correspond to the median expression values of each marker. The \code{data.frame} also
#'       contains a column called \code{cellType} which contains the name of the corresponding landmark. In order to derive
#'       the name of the landmark. the name of the corresponding FCS file is split on the \code{"_"} character, and the last
#'       field is used for the landmark name.
#'   }
#'
#'
#' @export
load_landmarks <- function(files.list, asinh.cofactor, transform.data, ...) {
    res <- NULL
    for(f in files.list) {
        population <- tail(strsplit(f, "_")[[1]], n = 1)
        population <- gsub(".fcs", "", population)
        fcs <- flowCore::read.FCS(f, ...)
        tab <- scfeatures::convert_fcs(fcs, asinh.cofactor, transform.data, clip.at.zero = T)

        tab <- data.frame(tab, cellType = population, stringsAsFactors = F, check.names = F)
        res <- rbind(res, tab)
    }

    tab.landmarks <- plyr::ddply(res, ~cellType, plyr::colwise(median))
    return(list(landmarks.data = res, tab.landmarks = tab.landmarks))
}


#' Load landmark populations from a directory
#'
#' @param dir A directory name. All the FCS files in the directory will be loaded
#' @inheritDotParams load_landmarks -files.list
#' @inherit load_landmarks return
#'
#' @export
load_landmarks_from_dir <- function(dir, ...) {
    files.list <- list.files(dir, full.names = T)
    load_landmarks(files.list, ...)
}


#' Remove connections to landmark nodes based on an expression threshold
#'
#' This function removes connections to a landmark node, in cases where the expression of all markers
#' for that landmark is less than a threshold. This is useful when mapping data that contains markers which
#' are all expected to be absent in certain landmarks
#'
#' @param dd A matrix of cosine similarity values between cluster nodes (rows) and landmark nodes (columns)
#' @param tab A \code{data.frame} containing the data (median expression values) for the landmark nodes
#' @param col.names A vector of colum names to be used for the computation
#' @param thresh The threshold to use. Connections to landmarks for which the expression of all the markers
#'   is \code{< thresh}, will be removed
#' @return Returns the same dd matrix, but with the values corresponding to the connections that have been
#'   removed set to 0
#'
distance_from_attractor_hard_filter <- function(dd, tab, col.names, thresh = 0.5) {
    tab <- tab[, col.names]
    w <- apply(tab[,col.names], 1, function(x, thresh) {all(x < thresh)}, thresh = thresh)
    if(any(w)) {
        message(sprintf("Hard removing some connections to landmarks using threshold: %f", thresh))
        flush.console()
        dd[, w] <- 0
    }
    return(dd)
}

#' Calculate (and filter) disatances between cluster and landmark ondes
#'
#' This function calculates the cosine similarity between the cluster and landmark nodes, and then filters
#' out the resulting matrix to include only the top-scoring edges. The result of this function can be used
#' as an adjacency matrix for graph construction
#'
#' @param m A \code{data.frame} containing expression values for the cluster nodes
#' @param tab A \code{data.frame} containing expression values for the landmark nodes
#' @param col.names A vector of column names to be used for the computation
#' @param q.thresh The quantile probability to be used for filtering edges. The algorithm
#'   will calculate a weight threshold based on this quantile of the weight distribution
#' @param min.similarity The minimum similarity value to be used for thresholding edges. Irrespective of
#'   the results of the quantile computation, the actual threshold used will never be less than this value
#' @return Returns a matrix with each row corresponding to a cluster node, and each column
#'   corresponding to the similarity between the cluster node and the respective landmark node. This matrix
#'   can be directly used as an adjacency matrix for graph construction
#'
#'
get_distances_from_landmarks <- function(m, tab, col.names, q.thresh, min.similarity = 0.5) {
    att <- as.matrix(tab[, col.names])
    #row.names(att) <- as.character(1:nrow(tab))
    m <- as.matrix(m[, col.names])
    dd <- t(apply(m, 1, function(x, att) {cosine_similarity_from_matrix(x, att)}, att))
    dd <- distance_from_attractor_hard_filter(dd, tab, col.names, thresh = 1)

    dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
    dist.thresh <- max(c(dist.thresh, min.similarity))



    dd[is.na(dd)] <- 0 #This can happen if one of the landmarks has all 0's for the markers of interest

    # This function modifies the dd matrix directly
    filter_matrix(dd, dist.thresh)
    colnames(dd) <- tab$cellType
    return(dd)
}

add_landmarks_labels <- function(G, v) {
    w <- V(G)$type == "landmark"
    V(G)$name[w] <- V(G)$Label[w] <- V(G)$cellType[w]
    return(G)
}

set_visual_attributes <- function(G) {
    att <- V(G)$type == "landmark"
    V(G)$r <- 79
    V(G)$g <- 147
    V(G)$b <- 222
    V(G)$size <- 10

    V(G)[att]$r <- 255
    V(G)[att]$g <- 117
    V(G)[att]$b <- 128
    V(G)[att]$size <- 20

    E(G)$r <- 180
    E(G)$g <- 180
    E(G)$b <- 180

    return(G)
}

#' Add cluster vertices to a landmark graph
#'
#' @param G An \code{igraph} object representing the landmark graph
#' @inheritParams get_scaffold_map
#'
#' @return Returns a new \code{igraph} object with the added cluster vertices
#'
add_vertices_to_landmarks_graph <- function(G, tab.clustered, tab.landmarks, col.names, min.similarity = 0.5) {
    dd <- get_distances_from_landmarks(tab.clustered, tab.landmarks, col.names, min.similarity)
    n <- nrow(dd)
    num.vertices <- length(V(G))
    G <- igraph::add.vertices(G, n)
    v.seq <- (num.vertices + 1):length(V(G))
    V(G)[v.seq]$name <- as.character(v.seq)
    V(G)[v.seq]$Label <- paste("c", tab.clustered$cellType, sep = "")
    row.names(dd) <- as.character(v.seq)

    for(i in 1:nrow(dd)) {
        v <- dd[i,]
        v <- v[v > 0]
        if(length(v) > 0) {
            e.list <- c(rbind(as.character(num.vertices + i), names(v)))
            G <- G + igraph::edges(e.list, weight = v)
        }
    }

    v.count <- igraph::vcount(G)

    V(G)[1:num.vertices]$type <- "landmark"
    V(G)[(num.vertices + 1):v.count]$type <- "cluster"

    for(i in names(tab.clustered))
        G <- igraph::set.vertex.attribute(G, name = i, index = (num.vertices + 1):v.count, value = tab.clustered[, i])

    G <- set_visual_attributes(G)
    return(G)
}


#' Add connections between cluster nodes
#'
#' This function takes an existing Scaffold map and adds connections between the cluster nodes
#'
#' @param G An \code{igraph} object, representing a Scaffold map, as returned by \code{add_vertices_to_landmarks_graph}
#' @param col.names A vector of column names to be used for the comptation of similarities
#' @param weight.factor Weight factor. The edge weights for the inter-cluster connections will be multiplied
#'   by this weight factor. This is useful in case one wants to weight the connections between the clusters differently
#'   from the connections between the clusters and the landmarks
#' @return Returns an \code{igraph} object, representing a Scaffold map with inter-cluster connections
#'
#'
add_inter_clusters_connections <- function(G, col.names, weight.factor) {
    tab <- igraph::as_data_frame(G, what = "vertices")
    tab <- tab[tab$type == "cluster",]
    m <- as.matrix(tab[, col.names])
    row.names(m) <- tab$name
    dd <- cosine_similarity_matrix(m)
    diag(dd) <- 0
    dd[is.na(dd)] <- 0

    dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
    dist.thresh <- max(c(dist.thresh, 0.5))

    # Both these functions modify the matrix dd directly

    filter_matrix(dd, dist.thresh)
    filter_matrix_by_rank(dd, 3)

    e.list <- NULL

    for(i in 1:nrow(dd)) {
        v <- dd[i,]
        v <- v[v > 0]
        if(length(v) > 0)
            e.list <- rbind(e.list, data.frame(a = tab[i, "name"], b = names(v), weight = v, stringsAsFactors = FALSE))
    }

    temp <- as.matrix(e.list[, c("a", "b")])
    temp <- t(apply(temp, 1, sort))
    e.list <- data.frame(temp, weight = e.list$weight, stringsAsFactors = FALSE)
    names(e.list)[1:2] <- c("a", "b")

    e.list <- e.list[!duplicated(e.list[, c("a", "b")]),]
    e.list.igraph <- c(t(as.matrix(e.list[, c("a", "b")])))


    G <- G + igraph::edges(e.list.igraph, weight = e.list$weight * weight.factor)
    return(G)
}

#' Add information about edge types and highest scoring edges
#'
#' @param G The input \code{igraph} object
#' @return Returns \code{G} with 3 added edge properties (\code{cluster_to_landmark}, \code{highest_scoring}, \code{inter_cluster}),
#'   and an added vertex property (\code{highest_scoring_edge}), with the following meaning:
#'   \itemize{
#'     \item{\code{cluster_to_landmark}}: if this is 1, this edge is between a cluster and a landmark node
#'     \item{\code{highest_scoring}}: if this is 1, this is the highest scoring edge among all the \code{"cluster_to_landmark"}
#'       edges of a single cluster
#'     \item{\code{inter_cluster}}: if this is 1, this is edge is between two clusters
#'     \item{\code{highest_scoring_edge}}: a number representing the index of the \code{"highest_scoring"} edge for this vertex
#'   }
#'
get_highest_scoring_edges <- function(G) {
    E(G)$cluster_to_landmark <- 0
    E(G)$highest_scoring <- 0
    E(G)$inter_cluster <- 0

    e <- igraph::get.edges(G, E(G))
    e <- cbind(V(G)$type[e[,1]], V(G)$type[e[,2]])

    inter.cluster <- e[, 1] == "cluster" & e[, 2] == "cluster"
    to.landmark <-(e[, 1] == "cluster" & e[, 2] == "landmark") |  (e[, 1] == "landmark" & e[, 2] == "cluster")

    E(G)$inter_cluster[inter.cluster] <- 1
    E(G)$cluster_to_landmark[to.landmark] <- 1

    # Remove inter-cluster edges for this calculation
    g.temp <- igraph::delete.edges(G, E(G)[inter.cluster])

    V(g.temp)$highest_scoring_edge <- 0

    for(i in 1:(igraph::vcount(g.temp))) {
        if(V(g.temp)$type[i] == "cluster") {
            sel.edges <- igraph::incident(g.temp, i)
            max.edge <- sel.edges[which.max(E(G)[sel.edges]$weight)]
            V(g.temp)$highest_scoring_edge[i] <- max.edge
            E(G)$highest_scoring[max.edge] <- 1
        }
    }

    V(G)$highest_scoring_edge <- V(g.temp)$highest_scoring_edge
    return(G)
}



#' Creates a Scaffold map
#'
#' @param tab.clustered A \code{data.frame} containing the clusters to be represented in the Scaffold map. Each
#'   row corresponds to a cluster, and each column to a different marker
#' @param col.names A character vector of column names to be used for the analysis (these must correspond to the names
#'   of the columns in \code{tab.clustered})
#' @param tab.landmarks A \code{data.frame} containing the data corresponding to the landmark nodes. Each row corresponds
#'   to a landmark, and each column to a different cluster
#' @param G.landmarks An \code{igraph} object represeting the identity and position of landmark nodes. If this is
#'   \code{NULL} the position of the landmark nodes will be calculated de-novo as part of the computation.
#' @param min.similarity The minimum similarity value to be used as threshold to filter out edges. See the documentation
#'   for \code{\link{get_distances_from_landmarks}}
#' @param inter.cluster.col.names A vector of column names to be used to calculate inter-cluster edges If this
#'   is \code{NULL} no inter-cluster edges will be present in the Scaffold map (not recommended)
#' @param inter.cluster.weight.factor The weight of the inter-cluster edges will be multplied by this number before
#'   calculating the ForceAtlas2 layout. This can be used to fine tune the relative contribution of cluster-to-landmark vs
#'   cluster-to-cluster edges in the laytout
#' @param overlap.method The method to be used for resolving nodes overlap in the final layout. See the documentation for
#'   \code{\link{complete_forceatlas2}}
#'
#' @return Returns a list with the following components
#'  \itemize{
#'     \item{\code{G.landmarks}}: an \code{igraph} object representing the identity and position of landmark nodes, suitable
#'       for succesive invocations of this function
#'     \item{\code{G.complete}}: an \code{igraph} object representing the Scaffold map
#'   }
#'
#' @export
get_scaffold_map <- function(tab.clustered, col.names, tab.landmarks, G.landmarks = NULL, ew.influence = ceiling(length(col.names) / 3),
                             min.similarity = 0.5, inter.cluster.col.names = col.names, inter.cluster.weight.factor = 0.7, overlap.method = "repel") {

    message(sprintf("Running with Edge weight: %f", ew.influence))
    flush.console()

    if(is.null(G.landmarks)) {
        G.landmarks <- build_graph(tab.landmarks, col.names)

        G.complete <- add_vertices_to_landmarks_graph(G.landmarks, tab.clustered, tab.landmarks, col.names, min.similarity)
        G.complete <- complete_forceatlas2(G.complete, first.iter = 50000,
                                           overlap.iter = 20000, ew.influence = ew.influence, overlap.method = overlap.method)
        if(!is.null(inter.cluster.col.names)) {
            message("Adding inter-cluster connections with markers:")
            message(paste(inter.cluster.col.names, collapse = " "))
            message(sprintf("Weight factor:%f", inter.cluster.weight.factor))
            flush.console()
            G.complete <- add_inter_clusters_connections(G.complete, inter.cluster.col.names, weight.factor = inter.cluster.weight.factor)
            G.complete <- complete_forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                               ew.influence = ew.influence, overlap.method = overlap.method)
        }
        V(G.landmarks)$x <- V(G.complete)$x[1:(igraph::vcount(G.landmarks))]
        V(G.landmarks)$y <- V(G.complete)$y[1:(igraph::vcount(G.landmarks))]
    }
    else {
        G.complete <- add_vertices_to_landmarks_graph(G.landmarks, tab.clustered, tab.landmarks, col.names, min.similarity)

        fixed <- rep(FALSE, igraph::vcount(G.complete))
        fixed[1:(igraph::vcount(G.landmarks))] <- TRUE

        G.complete <- complete_forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                           overlap.method = "repel", ew.influence = ew.influence, fixed = fixed)
        if(!is.null(inter.cluster.col.names)) {
            message("Adding inter-cluster connections with markers:")
            message(paste(inter.cluster.col.names, collapse = " "))
            message(sprintf("Weight factor:%f", inter.cluster.weight.factor))
            flush.console()
            G.complete <- add_inter_clusters_connections(G.complete, inter.cluster.col.names, weight.factor = inter.cluster.weight.factor)
            G.complete <- complete_forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                               overlap.method = overlap.method, ew.influence = ew.influence, fixed = fixed)
        }

    }

    G.complete <- get_highest_scoring_edges(G.complete)
    G.complete <- add_landmarks_labels(G.complete)
    V(G.complete)$name <- gsub(".fcs", "", V(G.complete)$name)
    return(list(G.landmarks = G.landmarks, G.complete = G.complete))
}


#' Write the result of a Scaffold analysis
#'
#' This function writes the result of a Scaffold analysis. This includes both the Scaffold map itself as a
#' \code{graphml} file, as well as (optionally) downsampled single-cell data for both the clusters and the landmarks
#' that can be used for visualization
#'
#' @param G An \code{igraph} object representing the Scaffold map
#' @param out.dir The name of the output directory
#' @param out.name The name of the output Scaffold map (The extensions \code{".graphml"} will be added to this name)
#' @param clusters.data A \code{data.frame} containing the single-cell clustered data. Each row of this \code{data.frame}
#'   corresponds to a cell. The \code{data.frame} must contain a column called \code{cellType} indicating cluster membership
#'   for each cell. If the data was clustered with the \code{scfeatures} package, this \code{data.frame} corresponds
#'   to the \code{"all_events.rds"} output file. If this is \code{NULL} no clusters data is written
#' @param landmarks.data The landmarks data, as returned by \code{\link{load_landmarks}} or \code{\link{load_landmarks_from_dir}}
#'   If this is \code{NULL} no landmarks data is written
#' @param downsample.to The number of events to target for downsampling, only used if either \code{clusters.data} or
#'   \code{landmarks.data} are provided
#' @export
#'
write_scaffold_output <- function(G, out.dir, out.name, clusters.data = NULL, landmarks.data = NULL, downsample.to = 1000) {
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    # This is necessary because igraph does not seem to handle NA's correctly when writing graphML

    attrs <- igraph::list.vertex.attributes(G)

    for(x in attrs) {
        v <- igraph::get.vertex.attribute(G, x)
        v[is.na(v)] <- 0
        G <- igraph::set.vertex.attribute(G, name = x, value = v)
    }

    igraph::write.graph(G, file.path(out.dir, sprintf("%s.graphml", out.name)), format = "graphml")

    if(!is.null(clusters.data)) {
        message(sprintf("Downsampling to %d events", downsample.to))
        clusters.data <- downsample_by(clusters.data, "cellType", downsample.to)
        write_clusters_data(clusters.data, out.name, out.dir)
    }

    if(!is.null(landmarks.data))
        write_landmarks_data(landmarks.data, out.dir, downsample.to)

    return(invisible(NULL))
}


#' Write landmarks data
#'
#' @param landmarks.data The landmarks data, as returned by \code{\link{load_landmarks}} or \code{\link{load_landmarks_from_dir}}
#' @param out.dir The output directory
#' @param downsample.to The number of events to be used as downsampling target
#'
#'
write_landmarks_data <- function(landmarks.data, out.dir, downsample.to = 1000) {

    landmark.data.dir <- file.path(out.dir, "landmarks_data")
    dir.create(landmark.data.dir, recursive = TRUE, showWarnings = FALSE)
    downsampled.data <- downsample_by(landmarks.data$landmarks.data, "cellType", downsample.to)

    plyr::d_ply(downsampled.data, ~cellType, function(x) {
        saveRDS(x, file = file.path(landmark.data.dir, sprintf("%s.rds", x$cellType[1])))
    })

}

#' High level wrapper for performing a Scaffold analysis
#'
#' This function represent a high-level entry-point for performing a Scaffold analysis
#'
#' @param files.list Character vector. The list of files to process
#' @param ref.file The name of the file to be used as reference (i.e. the first file to be processed, which will determine
#'   the position of the landmarks in all subsequent maps)
#' @param landmarks.data The landmarks data as returned by \code{\link{load_landmarks}} or \code{\link{load_landmarks_from_dir}}
#' @param col.names A character vector of column (i.e. marker) names to be used for the analysis. These columns have to be present in all the files
#' @param out.dir The name of the output directory
#' @param process.clusters.data If this is \code{TRUE} this function will look for a file with extension \code{.all_events.rds} for
#'   each file in \code{files.list} (see the Documentation of \code{scfeatures::cluster_fcs_files}). This file contains single-cell
#'   data (i.e. each row represent a cell, as opposed to the files in \code{files.list} where each row represents a cluster). Each file
#'   will be processed using the \code{\link{write_clusters_data}} function. This processing is used for downstream data visualization
#'   but it is not strictly necessary to create the graph.
#' @param downsample.to The number of events to target for downsampling when processing the clusters data. This is only used if
#'   \code{process.cluster.data == TRUE} and it does not affect the computation as this data is only used for downstream visualization
#' @param ... Additional parameters passed to \code{\link{get_scaffold_map}}
#'
#' @export
run_scaffold_analysis <- function(files.list, ref.file, landmarks.data, col.names, out.dir = "scaffold_result", process.clusters.data = TRUE, downsample.to = 1000, ...) {
    G.landmarks <- NULL

    for(f in files.list) {
        message(paste("Processing", f, sep = " "))
        flush.console()
        tab <- read.table(f, header = T, sep = "\t", quote = "", check.names = F, comment.char = "", stringsAsFactors = F)

        tab <- tab[!apply(tab[, col.names], 1, function(x) {all(x == 0)}),]

        scaffold.res <- get_scaffold_map(tab, col.names, landmarks.data$tab.landmarks, G.landmarks, ...)

        out.name <- gsub(".clustered.txt", "", basename(f))

        if(process.clusters.data) {
            clusters.data <- readRDS(gsub("txt$", "all_events.rds", f))
            write_scaffold_output(scaffold.res$G.complete, out.dir, out.name, clusters.data = clusters.data, downsample.to = downsample.to)

        }
        else
            write_scaffold_output(scaffold.res$G.complete, out.dir, out.name)

        G.landmarks <- scaffold.res$G.landmarks
    }

    if(process.clusters.data)
        write_landmarks_data(landmarks.data, out.dir, downsample.to)

    return(invisible(NULL))
}


test <- function() {
    setwd("C:/Users/fgherardini/temp/scaffold_demo")

    col.names <- c("CD45", "CD61", "CD7", "CD33", "CD11c", "CD123", "CD14", "CD11b", "CD8",
        "CD4", "CD3", "CD66", "CD16", "CD1c", "BDCA3", "CD45RA", "CD161", "CCR7", "CD19", "IgM", "CD56", "HLA-DR")

    input.files <- list.files(pattern = "*.clustered.txt$")

    landmarks.data <- load_landmarks_from_dir("gated/", asinh.cofactor = 5, transform.data = T)
    run_scaffold_analysis(input.files, input.files[1], landmarks.data, col.names)



}



