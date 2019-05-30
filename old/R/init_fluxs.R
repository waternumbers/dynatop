#' inialise base flows assumming an initial steady state
init_fluxes <- function(groups,      # hsu definitions
                        W,          # flux distribution matrix (not required anymore)
                        ichan,      # channel identifiers
                        qt0,        # initial (specific) catchment discharge / recharge (m/hr)
                        FLOW_NAMES,
                        dtt=1)      # inner time step

{

  # ex = saturation excess overland flow
  # uz = gravity drainage from unsaturated into saturated zone
  # qof = saturation excess overland flow
  flows <- data.frame("id"=groups$id,
                      "tag"=groups$tag) #

  # add all flow and initialise at zero
  # "qin"=0,    # qin = total input flow into areal groupings (m^3/hr)
  # "qbf"=0,    # qbf = base flow out of saturated zone (m^3/hr/m^2)
  # "ex"=0,     # base flow excess
  # "pex"=0,    # precipitation excess
  # "qof"=0,    #
  # "ae"=0,     # actual evap
  # "rain"=0,
  # "uz"=0)     # uz drainage into water table
  for(nm in FLOW_NAMES)
  {
    flows[,nm]=0;
  }

  # theoretical max for each group
  qbfmax <- exp(groups$ln_t0-groups$atb.bar)

  # solve for base flows assumming initial steady state  recharge / discharge
  flows$qbf <- init_base_flows(groups,  W, r=as.numeric(qt0), ichan)

  # total upslope area for all elements in each group
  # assumme uniform recharge over all upslope area ai for element i
  # then base flow qi for elemnet is rai and total for area is r.sigma(ai)
  # at steady state specific discharge qbf0 is equivalent to rainfall
  # a.bar is specific discharge for entire group per unit recharge
 # flows$qbf <- qt0*groups$sigma.a/15  # aaargh initialisation depends on times steps as these are aggregated over inner loop when distributed to channel
  # initilise to max value and let the loop sort it out

#  flows[ichan,]$qbf <- 0#qt0/15

  # mass balance implies that input is output minus recharge
  # recharge = r + flux distributed from other areas using weighting matrix
 # r <- DistributeDownslopeFluxes(flows$qbf, w)$qin+qt0
#   qin <- dist.flux(groups, flows, W)$qin - qt0*groups$area
#
#   # something awry with init base flows and / or calculated upslope areas
#   if(any(qin<0))
#   {
#     warning("Given initial specific recharge / discharge qt0 and accumulated uplsope areas sigma.a give -ve input fluxes - check")
#   }
  flows$qin <- 0 #pmax(qin, 0)

  # in steady state drainage to unsat zone is equal to recharge (this *is* correct)
  flows$uz <- as.numeric(qt0)

  return(flows)
}
