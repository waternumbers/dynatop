#' Routing by eigen value method
#'
#' @description Functions to setup and then evaluate the eigen vector method of Dummit, 2016
#' @param K matric such that dx/dt = K*s
#' @param v vector of volumes to be redistributed
#' @param dt time step to be taken
#' @param eig an object containing summaries of eigen matrix and vectors
#'
#' @return Either an eig object (eigen_setup) or a revised set of volumes after time step
#'
#' @details The current solution is based on inverting the matrix of eigen vectors in the setup step. This means that on each evaluation requires the sonstant in the solution arecomputed using a matrix multiplication. It may be more robust (and memery efficent) to solve for the initial conditions based on the eigen vector matrix instead.
#' @name eigen_routing
#' @rdname eigen_routing
#' @examples
#' # a real value problem
#' K <- matrix( c(2,0,0,-2,3,-1,4,-2,2),3 )
#' v0 <- c(2,0,2)
#' dt <- 1
#' # solution using eigen_routing
#' tmp <- eigen_routing_setup(K)
#' x <- eigen_routing_step(v0,tmp,dt)
#'
#' # solution should be
#' B <- matrix(c(1, 0, 0, -2,1,1,4,-2,1),3)
#' xhat <- B %*% (solve(B,v0)*exp(-c(2,1,4)*dt))
#'
#' # so this should be small
#' x-xhat
NULL

#' Setup of the eigen value system
#' @rdname eigen_routing
#' @export
eigen_routing_setup <- function(K){
    eig = eigen(K)
    out <- list(U = eig$vectors, #eigen vectors
                lambda = eig$values, # eigen values
                Um1 = solve(eig$vectors) # inverse of eigen matrix
                )
    return(out)
}

#' take one step of the solution
#' @rdname eigen_routing
#' @export
eigen_routing_step <- function(v,eig,dt){
    ## arbritary constants for general solution
    K <- eig$Um1 %*% v
    ## general solution K*exp(L*t)*U
    ## U the eigenmatrix
    ## L the eiegnvalues
    ## t time
    Ext = eig$U %*% (K*exp(eig$lambda*dt))
    return(as.vector(Re(Ext)))
}
