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


#' ForceAtlas2 force-directed layout
#'
#' @param G The input \code{igraph} object. The graph must have an edge attribute named \code{weight}, representing edge weights
#' @param ew.influence Edge weight influence. The edge weights are set to \code{edge.weight ^ ew.influence} before the
#'   calculation (see original ForceAtlas2 publication)
#' @param kgrav The gravity constant. Higher values will result in more compact graphs (see original ForceAtlas2 publication)
#' @param iter Maximum number of iterations. The algorithm will stop after this
#'   many iterations, or when the average displacement of the nodes between two
#'   iterations is less than the \code{stopping.tolerance} threshold (see below)
#' @param prevent.overlap Set this option to \code{TRUE} to prevent the nodes
#'   from overlapping (see ForceAtlas2 description)
#' @param fixed A boolean vector of length equal to the number of nodes in the
#'   graph which specifies which nodes, need to be held in a fixed
#'   position. If this is \code{NULL} (default), no nodes are held fixed
#' @param stopping.tolerance The algorithm will stop after either \code{iter}
#'   number of iterations, or when the average displacement of the nodes between
#'   two iterations is less than this threshold
#' @param barnes.hut Whether to use the Barnes-Hut approximation for speeding up
#'   the calculations when dealing with large graphs. This option is
#'   automatically set to true when the graph has more than 2000 nodes
#' @return this function returns a list with three elements
#'   \itemize{
#'     \item{\code{lay}}: a numeric matrix with two columns containing the x and y coordinates of each node in the final
#'       layout
#'     \item{\code{avg.displ}}: a numeric vector containing the average displacement of the vertices at each iteration
#'     \item{\code{max.displ}}: a numeric vector containing the maximum displacement between all the vertices after
#'       each iteration
#'   }
#' @references \url{http://gephi.github.io}
#' @references Jacomy M1, Venturini T, Heymann S, Bastian M. ForceAtlas2, a
#'   continuous graph layout algorithm for handy network visualization designed
#'   for the Gephi software. PLoS One. 2014 Jun 10;9(6):e98679. doi:
#'   10.1371/journal.pone.0098679
#' @export
#'
layout_forceatlas2 <- function(G, ew.influence = 1, kgrav = 1, iter = 1000, prevent.overlap = FALSE, fixed = NULL, stopping.tolerance = 0.001, barnes.hut = FALSE) {
    v.count <- igraph::vcount(G)

    if(v.count >= 2000)
        barnes.hut <- TRUE
    if(v.count > 2000)
        stopping.tolerance <- 0.01
    else if(v.count > 800)
        stopping.tolerance <- 0.005
    else
        stopping.tolerance <- 0.001

    if(is.null(fixed))
        fixed <- rep(FALSE, v.count)

    lay <- NULL

    if(is.null(igraph::get.vertex.attribute(G, "x"))) {
        lay <- matrix(ncol = 2, nrow = v.count, data = rnorm(v.count * 2, 10, 2))
        colnames(lay) <- c("x", "y")
    } else {
        lay <- cbind(x = V(G)$x, y = V(G)$y)
        w <- is.na(lay[, "x"])

        if(any(w))
            lay[w,] <- rnorm(sum(w) * 2, 10, 2)

    }





    #This is only used with prevent.overlap
    if(is.null(igraph::get.vertex.attribute(G, "size")))
        V(G)$size <- rep(10, v.count)
    mass <- 1 + igraph::degree(G)
    F_att <- (E(G)$weight ^ ew.influence)
    edge_list <- igraph::get.edgelist(G, names = F) - 1 #This is gonna be used in the C code where the indexing is 0-based

    avg_displ <- numeric(iter)
    max_displ <- numeric(iter)

    if(barnes.hut)
        message("Using Barnes-Hut approximation\n")

    message(sprintf("Stopping tolerance: %f\n", stopping.tolerance))
    flush.console()

    layout_forceatlas2Cpp(lay, F_att, mass, V(G)$size, edge_list, avg_displ,
                            kgrav, iter, prevent.overlap, fixed, max_displ, stopping.tolerance, barnes.hut)

    return(list(lay = lay, avg_displ = avg_displ, max_displ = max_displ))
}


#' Performs a Complete cycle of ForceAtlas2
#'
#' This function performs a complete (i.e. possibly including overlap resolution) cycle of the ForceAtlas2 force-directed layout algorithm
#'
#' @param G The input graph
#' @param first.iter The maximum number of iterations in the first step, which is performed without overlap resolution
#' @param overlap.method If this is \code{NULL} overlap resolution is not performed. Otherwise this should be a string specifying the
#'   overlap resolution method. Two options are possible
#'   \itemize{
#'     \item{\code{"repel"}}: This is the method used in the original ForceAtlas2 implementation. Using this method, a repulsive force
#'       is applied to nodes that overlap each other. This method can cause problem in cases where the layout is extremely crowded,
#'       as this repulsive force becomes the major determinant of the layout, and the nodes end up being arranged essentially in a grid
#'     \item{\code{"expand"}}: In this method, the graph is linearly expanded, until no two nodes overlap anymore
#'   }
#' @param overlap.iter The maximum number of iterations for the overlap resolution step. This is only used if \code{overlap.method} is not
#'   \code{NULL}
#' @return Returns an \code{igraph} object with two additional vertex attributes \code{x} and \code{y}, containing the x and y coordinates
#'   of the vertices in the final layout
#'
#' @export
complete_forceatlas2 <- function(G, first.iter = 1000, overlap.method = NULL, overlap.iter = NULL, ...) {

    message("First iteration")
    flush.console()

    ret <- layout_forceatlas2(G, prevent.overlap = FALSE, iter = first.iter, ...)
    lay <- ret$lay

    G <- igraph::set.vertex.attribute(G, name = "x", value = lay[, 1])
    G <- igraph::set.vertex.attribute(G, name = "y", value = lay[, 2])
    if(!is.null(overlap.method)) {
        if(overlap.method == "repel") {
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
        else if(overlap.method == "expand")
            G <- adaptive_expand(G, overlap.iter)
    }
    return(G)
}











