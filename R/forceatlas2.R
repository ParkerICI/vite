adaptive_expand <- function(G, max.iter) {
    message("Starting adaptive expansion")
    flush.console()

    x <- V(G)$x
    y <- V(G)$y
    m <- cbind(x, y)
    ss <- outer(V(G)$size, V(G)$size, "+")

    for(i in 1:max.iter) {
        dd <- as.matrix(dist(m), method = "euclidean")
        dd <- dd - ss
        dd <- dd[upper.tri(dd)]
        if(all(dd >= 0))
            break
        else
            m <- m * 1.2
    }

    message(sprintf("Expansion stopped at iteration: %d", i))
    flush.console()
    V(G)$x <- m[, "x"]
    V(G)$y <- m[, "y"]

    return(G)
}



layout_forceatlas2 <- function(G, ew_influence = 1, kgrav = 1, iter = 1000, prevent.overlap = FALSE, fixed = rep(FALSE, vcount(G)), stopping_tolerance = 0.001, barnes_hut = FALSE) {
    if(vcount(G) >= 2000)
        barnes_hut <- TRUE
    if(vcount(G) > 2000)
        stopping_tolerance <- 0.01
    else if(vcount(G) > 800)
        stopping_tolerance <- 0.005
    else
        stopping_tolerance <- 0.001

    if(is.null(get.vertex.attribute(G, "x")))
        lay <- cbind(x = rnorm(vcount(G)), y = rnorm(vcount(G)))
    else
        lay <- cbind(x = V(G)$x, y = V(G)$y)


    #This is only used with prevent.overlap
    if(is.null(igraph::get.vertex.attribute(G, "size")))
        V(G)$size <- rep(10, igraph::vcount(G))
    mass <- 1 + igraph::degree(G)
    F_att <- (E(G)$weight ^ ew_influence)
    edge_list <- igraph::get.edgelist(G, names = F) - 1 #This is gonna be used in the C code where the indexing is 0-based

    avg_displ <- numeric(iter)
    max_displ <- numeric(iter)

    if(barnes_hut)
        message("Using Barnes-Hut approximation\n")

    message(sprintf("Stopping tolerance: %f\n", stopping_tolerance))
    flush.console()

    layout_forceatlas2Cpp(lay, F_att, mass, V(G)$size, edge_list, avg_displ,
                            kgrav, iter, prevent.overlap, fixed, max_displ, stopping_tolerance, barnes_hut)

    return(list(lay = lay, avg_displ = avg_displ, max_displ = max_displ))
}


#' Performs a Complete cycle of ForceAtlas2
#'
#' This function performs a complete (i.e. possibly including overlap resolution) cycle of the ForceAtlas2 force-directed layout algorithm
#'
#' @param G The input graph
#' @param first.iter The number of iterations in the first cycle, which is performed without overlap resolution
#' @param overlap_method If this is \code{NULL} overlap resolution is not performed. Otherwise this should be a string specifying the
#'   overlap resolution method. Two options are possible
#'   \itemize{
#'     \item{\code{"repel"}}: This is the method used in the original ForceAtlas2 implementation. Using this method, a repulsive force
#'       is applied to nodes that overlap each other. This method can cause problem in cases where the layout is extremely crowded,
#'       as this repulsive force becomes the major determinant of the layout, and the nodes end up being arranged essentially in a grid
#'     \item{\code{"expand"}}: In this method, the graph is linearly expanded, until no two nodes overlap anymore
#'   }
#' @return Returns an \code{igraph} object with two additional vertex attributes \code{x} and \code{y}, containing the x and y coordinates
#'   of the vertices in the final layout
#'
#'
complete_forceatlas2 <- function(G, first.iter = 1000, overlap.iter, overlap_method = NULL, ...) {

    message("First iteration")
    flush.console()

    ret <- layout_forceatlas2(G, prevent.overlap = FALSE, iter = first.iter, ...)
    lay <- ret$lay

    G <- igraph::set.vertex.attribute(G, name = "x", value = lay[, 1])
    G <- igraph::set.vertex.attribute(G, name = "y", value = lay[, 2])
    if(!is.null(overlap_method)) {
        if(overlap_method == "repel") {
            message("Second iteration with prevent overalp")
            flush.console()
            ret <- layout_forceatlas2(G, prevent.overlap = TRUE, iter = overlap.iter, ...)
            lay <- ret$lay

            if(any(is.na(lay)))
                message("Prevent overlap iteration failed")

            else {
                G <- igraph::set.vertex.attribute(G, name = "x", value = lay[, 1])
                G <- igraph::set.vertex.attribute(G, name = "y", value = lay[, 2])
            }
        }
        else if(overlap_method == "expand")
            G <- adaptive_expand(G, overlap.iter)
    }
    return(G)
}











