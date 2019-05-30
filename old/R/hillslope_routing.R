#' Routing of surface storage using eigen vector method of (see Dummit, 2016; )
#'
#' @param ex initial surface storage [L]
#' @param area area [L^2]
#' @param U eigenmatrix
#' @param Um1 inverse of eigenmatrix
#' @param L eigenvalues
#' @param t time step
#'
#' @return value of ex at end of time step
route_ex_eigen = function(ex, area, U, Um1, L, t){
    Ex = ex*area;
    ## arbritary constants for general solution
    K = Um1 %*% Ex;
    ## general solution K*exp(L*t)*U
    ## U the eigenmatrix
    ## L the eiegnvalues
    ## t time
    Ext = U %*% (K*exp(L*t))
    return(as.vector(Re(Ext/area)))
}

#' R function to compute gradient of flows relating to downslope distribution in saturated zone
#'
#' @description The function relates the rate of change of flow to unsaturated drainage input, upslope inputs and downslope output
#' @description Th function results from xponential transmissivity assumptions
#' @param t time - unused but required for ode solver
#' @param y numeric vector of discharges from groups
#' @param parms list of:
#' \itemize{
#' \item @param m transmissivity parameters [?]
#' \item @param Wdash complimentary sursurface routing matrix
#' \item @param qbf_max maximum possible flow
#' \item @param uz flow from unsaturated zone [L/T]
#' }
#'
#' @return rate of change of y with time
funR_dqdt <- function(t, y, parms)
{
  # y is q, discharge, and kinematic wave velocity is q/m
  res <-  y/parms$m * (matrix_vect_mult_cpp(parms$Wdash, y) + parms$uz)

  imax = res > parms$qbf_max & res > 0;

  # base flows cannnot increase beyond maximum allowed
  res = res * (1-imax)
  return(list(res))
}

#' R function to compute the recharge rate through unsaturated drainage into saturated zone
#' @description This is based on a time constant proportional to saturated deficit.
#' @description The output is capped so be less then suz/dt
#'
#' @param suz unsaturated zone storage [L]
#' @param sd saturated zone storage deficit [L]
#' @param td constant of proportionality in time constant [T/L]
#' @param dt time step [T]
funR_uz <- function(suz, sd, td, dt = 0.0001){
    uz <- rep(0,length(suz))
    ## location with moisture and space for it to drain to
    idx <- which(sd>0 & suz>0)
    uz[idx] <- pmin( suz[idx] / (sd[idx]*td[idx]) , suz[idx]/dt )
}
