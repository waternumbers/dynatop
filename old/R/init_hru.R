#' use the specified soil /hydraulic properties if supplied or assign from defaults
init_hru<- function(groups,
                    params=def_hsu_par(),
                    chan.par=def_chan_par(),
							      ichan=1, rain=NULL)
{
  # groups is table of area groupings and params is list of default parameter values
  # Search for columns for saturation transmissivity ln_t0, max root zone storage srMax
  # and exponential transmissivity profile factor m
  # initialise missing parameters with defaults
  params <- def_hsu_par()
  nms <- setdiff(names(params), colnames(groups))
  params<- params[nms]
  groups <- apply_params(groups, params)

  # add missing params using defaults. groups' parameters override defaults
#	nms <- setdiff(colnames(groups), names(params))


  #groups[igroup,] <- as.vector(merge_lists(params, def_hsu_par()))
#    }
	if(length(ichan)>0)
	{
	  # the "river" has zero root zone storage and there is no time delay in drainage reaching baseflow
	  # Storage is limited only by bankfull level - set SD and SDMax to large values to simulate
	  # Could apply physically realistic max storage (i.e average river depth at bankfull level) to simulate overbank events
	  groups[ichan,]$srz0 <- chan.par$srz0
	  groups[ichan,]$srz_max <- chan.par$srz_max   # rainfall goes straight to the unsat zone infiltration
	  groups[ichan,]$sd_max <- chan.par$sd_max   # maximum allowable SD
	  groups[ichan,]$vof <-  groups[ichan,]$vchan   # routing velocity

	#  groups[1:nchan,]$SD <- 3  # this is a huge value reflecting the channel's storage capacity
	  groups[ichan,]$td <- chan.par$td  # all drainage will be routed directly to base flow
	}
	# checj taht groups wetness is in decreasing order
	#if(!all(order(groups[-ichan,$atb.bar, decreasing=TRUE)==groups[-ichan,]$order)){warning("Wetness index not in decreasing order of HRU number")}
  # labels for groups default to ids (shoudl these be sequential?)
  if(is.null(groups$tag)){groups$tag<-groups$id}

	n.gauge <- ifelse(is.null(rain),
										1,
										ncol(rain))
  # ensure gauge reading supplied is withing range
  groups$gauge.id <- pmin(groups$gauge.id, n.gauge)

  return(groups)
}
