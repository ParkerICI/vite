
mark_nearest_neighbour <- function(G) {
    E(G)$nearest_neigh <- 0

    for(i in 1:(igraph::vcount(G))) {
        if(V(G)$type[i] == 2) {
            sel.edges <- igraph::incident(G, i)
            max.edge <- sel.edges[which.max(E(G)[sel.edges]$weight)]
            E(G)[max.edge]$nearest_neigh <- 1
        }
    }
    return(G)
}

cosine_similarity_from_matrix <- function(v, m) {
    m <- as.matrix(m[, names(v), drop = F])
    ret <- apply(m, 1, function(x, v) {return(crossprod(x, v)/sqrt(crossprod(x) * crossprod(v)))}, v)
    return(ret)
}

cosine_similarity_matrix <- function(m) {
    ret <- t(apply(m, 1, function(x, m) {cosine_similarity_from_matrix(x, m)}, m = m))
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



distance_from_attractor_hard_filter <- function(dd, tab, col.names, thresh = 0.5) {
    tab <- tab[, col.names]
    w <- apply(tab[,col.names], 1, function(x, thresh) {all(x < thresh)}, thresh = thresh)
    if(any(w))
        print("Hard removing some connections to unstained landmarks")
    dd[, w] <- 0
    return(dd)
}


get_distances_from_attractors <- function(m, tab, col.names, dist.thresh) {
    att <- as.matrix(tab[, col.names])
    row.names(att) <- as.character(1:nrow(tab))
    m <- as.matrix(m[, col.names])
    dd <- t(apply(m, 1, function(x, att) {cosine_similarity_from_matrix(x, att)}, att))
    dd <- distance_from_attractor_hard_filter(dd, tab, col.names, thresh = 1)

    dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
    dist.thresh <- max(c(dist.thresh, 0.5))



    dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest
    dd <- filter_similarity_matrix(dd, dist.thresh)
    return(dd)
}

add_attractors_labels <- function(G, v) {
    V(G)$name[1:length(v)] <- V(G)$Label[1:length(v)] <- v
    return(G)
}

set_visual_attributes <- function(G) {
    att <- V(G)$type == 1
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



add_vertices_to_attractors_graph <- function(G, tab.clustered, tab.median, col.names, dist.thresh = 0.7) {
    dd <- get_distances_from_attractors(tab.clustered, tab.median, col.names, dist.thresh)
    n <- nrow(dd)
    num.vertices <- length(V(G))
    G <- add.vertices(G, n)
    v.seq <- (num.vertices + 1):length(V(G))
    V(G)[v.seq]$name <- as.character(v.seq)
    V(G)[v.seq]$Label <- paste("c", tab.clustered$groups, sep = "")
    row.names(dd) <- as.character(v.seq)

    for(i in 1:nrow(dd)) {
        v <- dd[i,]
        v <- v[v > 0]
        if(length(v) > 0) {
            e.list <- c(rbind(as.character(num.vertices + i), names(v)))
            G <- G + edges(e.list, weight = v)
        }
    }

    maxx <- maxy <- rep(Inf, vcount(G))
    minx <- miny <- rep(-Inf, vcount(G))

    maxx[1:num.vertices] <- minx[1:num.vertices] <- V(G)$x[1:num.vertices]
    maxy[1:num.vertices] <- miny[1:num.vertices] <- V(G)$y[1:num.vertices]
    lay <- igraph::layout.kamada.kawai(G, minx = minx, maxx = maxx, miny = miny, maxy = maxy)
    colnames(lay) <- c("x", "y")
    G <- igraph::set.vertex.attribute(G, name = "x", value = lay[, "x"])
    G <- igraph::set.vertex.attribute(G, name = "y", value = lay[, "y"])

    V(G)[1:num.vertices]$type <- 1 #attractor
    V(G)[(num.vertices + 1):vcount(G)]$type <- 2 #cell

    for(i in names(tab.clustered))
        G <- igraph::set.vertex.attribute(G, name = i, index = (num.vertices + 1):vcount(G), value = tab.clustered[, i])

    G <- set_visual_attributes(G)
    return(G)
}


get_vertex_table <- function(G) {
    att <- list.vertex.attributes(G)
    ret <- NULL

    for(a in att) {
        d <- data.frame(get.vertex.attribute(G, a), stringsAsFactors = FALSE)
        if(is.null(ret))
            ret <- d
        else
            ret <- cbind(ret, d, stringsAsFactors = FALSE)
    }
    names(ret) <- att
    return(ret)
}

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

add_inter_clusters_connections <- function(G, col.names, weight.factor) {
    tab <- get_vertex_table(G)
    tab <- tab[tab$type == 2,]
    m <- as.matrix(tab[, col.names])
    row.names(m) <- tab$name
    dd <- cosine_similarity_matrix(m)
    diag(dd) <- 0
    dd[is.na(dd)] <- 0

    dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
    dist.thresh <- max(c(dist.thresh, 0.5))
    dd <- filter_similarity_matrix(dd, dist.thresh)
    dd <- filter_similarity_matrix_by_rank(dd, 3)

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



process_data <- function(tab, G.attractors = NULL, tab.attractors = NULL, col.names = NULL, att.labels = NULL, dist.thresh = 0.7,
                         already.clustered = FALSE, inter.cluster.connections = FALSE, col.names.inter_cluster = NULL, inter_cluster.weight_factor = 0.7, ew_influence,
                         overlap_method = NULL)
{

    if(!already.clustered)
    {
        tab <- cluster_data(tab, col.names)
        tab.clustered <- ddply(tab, ~groups, colwise(median))
    }
    else
        tab.clustered <- tab

    if(is.null(col.names.inter_cluster) || col.names.inter_cluster == "")
        col.names.inter_cluster = col.names
    if(is.null(G.attractors))
    {
        G.attractors <- build_graph(tab.attractors, col.names)

        G.complete <- add_vertices_to_attractors_graph(G.attractors, tab.clustered, tab.attractors, col.names, dist.thresh)
        G.complete <- complete.forceatlas2(G.complete, first.iter = 50000,
                                           overlap.iter = 20000, ew_influence = ew_influence, overlap_method = "repel")
        if(inter.cluster.connections)
        {
            print("Adding inter-cluster connections with markers:")
            print(col.names.inter_cluster)
            print(sprintf("Weight factor:%f", inter_cluster.weight_factor))
            G.complete <- add_inter_clusters_connections(G.complete, col.names.inter_cluster, weight.factor = inter_cluster.weight_factor)
            G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                               ew_influence = ew_influence, overlap_method = overlap_method)
        }
        V(G.attractors)$x <- V(G.complete)$x[1:vcount(G.attractors)]
        V(G.attractors)$y <- V(G.complete)$y[1:vcount(G.attractors)]
    }
    else
    {
        G.complete <- add_vertices_to_attractors_graph(G.attractors, tab.clustered, tab.attractors, col.names, dist.thresh)

        fixed <- rep(FALSE, vcount(G.complete))
        fixed[1:vcount(G.attractors)] <- TRUE

        G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                           overlap_method = "repel", ew_influence = ew_influence, fixed = fixed)
        if(inter.cluster.connections)
        {
            print("Adding inter-cluster connections with markers:")
            print(col.names.inter_cluster)
            print(sprintf("Weight factor:%f", inter_cluster.weight_factor))
            G.complete <- add_inter_clusters_connections(G.complete, col.names.inter_cluster, weight.factor = inter_cluster.weight_factor)
            G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                               overlap_method = overlap_method, ew_influence = ew_influence, fixed = fixed)
        }

    }
    G.complete <- add_attractors_labels(G.complete, att.labels)
    V(G.complete)$name <- gsub(".fcs", "", V(G.complete)$name)
    return(list(G.attractors = G.attractors, G.complete = G.complete, tab.attractors = tab.attractors, tab = tab, col.names = col.names))
}





#' Builds an unsupervised force-directed graph
#'
#' @param tab A \code{data.frame} of graph nodes. All the columns of \code{tab} will become vertex properties of the output graph
#' @param col.names A character vector indicating which columns of \code{tab} should be used to calculate distances
#' @param filtering.threshold The threshold used to filter edges in the graph
#'
#' @return Returns an \code{igraph} graph object. This function also performs graph-based clustering. The resulting cluster membership
#'   is contained in the \code{community_id} vertex attribute of the resulting graph
#'
#' @export
get_unsupervised_graph <- function(tab, col.names, filtering.threshold, output.name) {
    G <- build_graph(tab, col.names, filtering_T = filtering.threshold)

    for(i in names(tab))
        G <- igraph::set.vertex.attribute(G, name = i, value = tab[, i])

    cc <- igraph::multilevel.community(G)
    V(G)$community_id <- as.character(cc$membership)
    V(G)$name <- seq_along(V(G)$name)

    message("Running ForceAtlas2...")
    G <- complete_forceatlas2(G, first.iter = 50000, overlap.iter = 1, overlap_method = NULL, ew_influence = 5)
    message("ForceAtlas2 done")

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

