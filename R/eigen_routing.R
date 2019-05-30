#' Routing by eigen value method
#'
#' @description Functions to setup and then evaluate the eigen vector method of Dummit, 2016
#' @description The following notes are taken verbatim from the original implimentation of Peter. They do not appear to match the equations and may be irrelevent but added for completeness...
#'  ## Qin = exS %*% Wsurf
#'  ## Linear storage-discharge relationship
#'  ## Qout = a * vof * ex
#'  ## dExS/dt = Qin - Qbf
#'  ##        = Wv %*% exS
#'  ## ExS is total surface excess storage
#'  ## where Wv = complementary overland routing matrix
#'
#'  ## Solution for surface routing by Eigenvalue method (see Dummit, 2016; )
#'  ## general soln at time t is
#'  ## sigma(exp(lambda_i * t*)ci * ui)
#'  ## where ui are eigenvectors, li eigenvalues and ci constants TBD by
#'  ## boundary conditions
#'
#' @param v vector of volumes to be redistributed
#' @param W matrix of weights for distribution
#' @param vof UNKNOWN
#' @param dt time step to be taken
#' @param eig an object containing summaries of eigen matrix and vectors
#'
#' @return Either an eig object (eigen_setup) or a revised set of volumes after time step
#' @name eigen_routing
#' @rdname eigen_routing
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
