# determine base flows assuming input output mass balance in steady state
init_base_flows <- function(groups,
                          W,     # flux distribution matrix
                          r,     # specific recharge e.g rainfall and discharge (assummed equal)
                          ichan)
{
  # across each time step
  # input:
  # external recharge e.g. rain, applied equally over each HSU, entering water
  # table as gravity drainage from unsaturated zone
  # base flow distributed out of other HSUs
  # Qint1 <- dt*(r + t(flows$qbf) %*% W)*Ag
  # output at next time:
  # base flow
  # effective evap (ignore)
  # channel discharge - ignore for time being by assumming everything distributed to
  # channel reaches outlet within one step
  # Qout <- dt*flows$qbf*groups$area
  # in steady state these must be equal. Equating
  # Qout = Qin
  # (r + qbf %*% W)*Ag = qbf*Ag

  # t(W)*qbf - qbf = -r
  # (I-t(W))*qbf=r
  # this gives a system of n equations for n unknowns qbf1...qbf which can be
  # easily solved used base::solve. It also gives a quick way to check that the
  # supplied downslope weighting
  # ignore base flow in channel HSUs as these are distributed by the channel
  # delay histogram
  W[ichan,] <- 0 # recharge to river is by redistributed baseflow  only?

  n <- nrow(W)
  ItW <- diag(n)-t(W)

  ra<-groups$area * r    # rep(r, n)*g
  ra[ichan]  <- 0
  # default
  Qbf0 <- r * groups$area
#  ra[ichan] <- 0  # # recharge to river is by redistributed baseflow  only?
  Qbf0 <- tryCatch(Qbf0 <- solve(ItW, ra), silent =TRUE)
  qbf0 <- Qbf0/groups$area
  # if matrix singular somthing is wrong with the way the matrix has been set-up
  # need to check in CheckInput

  return(qbf0)
}
