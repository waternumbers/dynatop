# #############################################################################
# see Beven (2012) eqns B6.1.18 and B6.1.21, pp213-215, for use of areal mean of
# soil-topographic index (Beven, 1986b) to initialise (mean) storage deficts
# -----------------------------------------------------------------------------
init_stores <- function(groups,
                         flows,
                         ichan, qt0)
{
  # storages and deficits for each areal group ([L] per plan area)
  # suz = unsaturated zone storage
  # srz = root zone storage
  # sd = storage deficit (depth to water table)
  # storages all initially zero. Root zone fills from srz0 until reaching
  # "field capacity" (srz_max) - the max RZ storage for the areal group
  stores <- data.frame("id"=groups$id, "tag"=groups$tag, "suz"=0, "srz"=0,
                       "ex"=0,    # excess storage in e.g. overland flow
                       "bf.ex"=0,
                       "sd"=0, "wetting"=TRUE)
  # initialise storage deficit to maximum allowed - allows init base flow to remove some before recharge reduces deficit again
  stores$sd <- groups$sd_max

  # q0  saturated specific areal base flows  estimated from areal average
  # of topographic index and saturated transmissivity
  q0 <- exp(groups$ln_t0-groups$atb.bar)

  # initial base flows estimated from initial recharge assumming steady state
  # see Beven 2012 equation 6.1.21
  sdbar<- pmax(-groups$m * log(flows$qbf/q0), 0, na.rm=TRUE)

  # Assign mean value to areal effective storage deficits (no -ve storages)
  if(any(sdbar<=0))
  {
  #  LogEvent("Calculated initial storage deficit -ve, setting to 0 (saturation)")
  }
  stores$sd <- sdbar
  if(any(stores$sd>groups$sd_max))
  {
    LogEvent("Calculated initial SD below maximum allowed for groups, set to max")
    stores$sd <- pmin(sdbar, groups$sd_max)
  }

  if(any(groups$srz0 > 1)){warning("Initial root zone storage proportion > 1: ignored")}

  # srz0 is proportion 0 - 1 of root zoone initially filled
  stores$srz <- pmin(abs(groups$srz0*groups$srz_max), groups$srz_max)


  # initialise unsat storage now estimates for initial deficit obtained
  # see formulations from Beven and Wood (1983) relating drainage flux (which we take as equal to
  # steady-state recharge) to unsat storage, defict and time delay parameter
  # qv = suz/td.sd
  # what is the maximum unsat storage given a particular sd?
  stores$suz <- qt0*groups$td*stores$sd

  return(stores)
}
