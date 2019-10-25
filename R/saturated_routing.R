#' R function to compute gradient of flows relating to downslope distribution in saturated zone
#'
#' @description The function relates the rate of change of flow to unsaturated drainage input, upslope inputs and downslope output
#' @description The function results from xponential transmissivity assumptions
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
fun_dlex_dt <- function(t, y, parms)
{
    ## y is q, discharge, and kinematic wave velocity is q/m
    res <-  y/parms$m * ( (parms$Wdash %*% y) + parms$quz)
    
    ## res <-  y/parms$m * (matrix_vect_mult_cpp(parms$Wdash, y) + parms$quz)
    
    ##imax = res > parms$lsz_max & res > 0; # original constrains gradient to qbf_max - why??
    imax <- y >= parms$lsz_max & res > 0;
  # base flows cannnot increase beyond maximum allowed
  #res = res * (1-imax) - commented out for testing
  return(list(res))
}
