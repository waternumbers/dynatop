################################################################################
# initialisation and input checks for Dynamic TOPMODEL
#-------------------------------------------------------------------------------
init_input <- function(groups, dt, ntt,
                      weights, rain, pe,
                      FLOW_NAMES,
                      STORE_NAMES,
                      sim.start, sim.end, # PJS remove graphics disp.par,
                      routing=NULL,
					            ichan=1, i.out=1,
                      qobs=NULL,   # specific or absolute discharge
                      qt0=1e-4,    # initial discharge in m/hr
                      calling.env = NULL)
{
  if(is.null(i.out)){i.out<-1}

  # need this for $ notation
  groups<-data.frame(groups)
  # PJS remove graphics qmax <- disp.par$max.q

  if(length(rain)==0)
  {
    stop("Rainfall input required")
  }

  tz <- ""
  # use the rainfall data to infer the time step and run start and end times,
  # if not otherwise specified
  if(is.null(dt))
  {
    dt <- as.numeric(difftime(index(rain)[2], index(rain)[1], units="hour"))
  }

  # ensure using POSXIct - fails
  index(rain) <- as.POSIXct(index(rain))

  ## PJS remove graphics dispStart <- sim.start + disp.par$"graphics.delay"*3600 # hours -> seconds

  # note that DTM uses an input time step of hours
  # everything is based around the rainfall time series
  tms <- index(rain)
  sim.start <- tms[1]
  sim.end <- tms[length(tms)]

  # check that the river reach distribution is of the correct dimensions -
  # should be ngroup x nreach, where nth row determines the distribution of the
  # nth group into the river network
  ngroups <- nrow(groups)
  # number of river reaches determined from flow length distribution - each row
  # represents a reach length and proportion of channel within that reach
  # nreach <- nrow(routing)

  # inner time step
  dtt <- dt/ntt

  nchan <- length(ichan)
  # channel identifiers
#  ichan <- 1:nchan
  weights <- as.matrix(weights)

  # times for output series, dt in hours
  nmax <- length(tms)

  # create dummy pe and qobs if none supplied, Q) can be used to set absolute discharge max limit
  if(is.null(pe))
  {
  	pe <- rain
  	pe[] <-0
  }

  # ensure everything in same units
  rain <- convert_values(rain)
  qt0 <- convert_values(qt0)
  pe <- convert_values(pe)

  #qobs <- GetTimeSeriesInputRange(qobs, sim.start, sim.end, verbose=FALSE)

  # dummy series if nothing supplied or outside range of simulation
  if(length(qobs)>0)
  {
		qobs <- xts(qobs)

		index(qobs) <- as.POSIXct(index(qobs))
		# want a specific input, observed values often (always?) in cubic m^3/s
		qobs <- convert_values(qobs, msg="Observed flows converted to m/hr: check that specfic discharge in m/hr has been supplied")

		# check (qobs should be in m/hr)
    qin <- sum(rain, na.rm = TRUE)
    qout <- sum(qobs[,1], na.rm=TRUE)

    if((qout-qin)/qout> 0.25)
    {
      #warning("Water balance out by >25%, check rainfall and /or observations")
    }
		if(length(qt0)==0)
		{
			# if the initial value not suppplied then take as first observation
			qt0 <- as.numeric(qobs[1, 1])
			message(paste0("Initialising fluxes using first discharge observation of ", round(qt0*1000, 1), " mm/hr"))
		}
  }

  # initialise storages and deficits for each group add static parameters to
  # areal groups if a soils map provided then fix up the soil parameters from
  # here, otherwise continue to use default
  groups <- init_hru(groups, ichan=ichan, rain=rain)

    # add the actual evap
  evap <- cbind("pe"=pe, "ae"=NA)

  try(check_and_fix_input(environment()))

	# calculated time delay matrix for channel routing. removes initial column with flow lengths
	routing <- get.routing(groups,
											 weights,
											 routing,
											 ichan,	dt)
  # max number of time steps that river flux can be shifted
  # Time series of total output flows from output reach(s). add enough extra elements to hold time delayed flows
  n.shift <- nrow(routing) + 100
  times.ex <-seq(sim.start, sim.end+n.shift*dt*3600, by = 3600*dt)

  nms <- groups[ichan,]$id

  # numbet of outlet reaches to produce flow output
  nout <- ncol(routing)

  # create the output series with enough entries to account for extra flow routed
  Qr <- xts(matrix(0,nrow=n.shift+nmax, ncol=nout), order.by=times.ex)

# total overland flow contribution
  qof <- xts(matrix(0,nrow=nmax,ncol=1), order.by=tms)

  #*****************************************************************************
  # initialisation of flows, storages and deficits for each areal group flows
  # (at time step t, per plan area except where stated)
#  message("... initialising fluxes and storages...\n")

  # can determine a theorectical max for discharge: set all base flow to max ,
  # controlled by saturation transmissivity (z=0)  add in any rainfall that
  # could be transferred from overland flow, distribute to channel by flux routing matrix,
  # then route to outlet by time delay table
  # see discussion of TOPMODEL theory in Beven (2012), Chapter 6, Box 6.1, p210
  # gamma is the average soil-topographic index defined in  Beven (1986a)
  # typical values for ln_t0 are in range (-7,8).
  # Could examine catchment and its observed discharges to calibrate this
  groups$qbf_max <- exp(groups$ln_t0-groups$atb.bar)

  # using the slope - correct!
  if(!is.null(groups$sbar))
  {
    # limiting transmissivty x mean slope
    groups$qbf_max <- exp(groups$ln_t0)*groups$sbar
  }

  # huge value ensures no excess in river (or set to zero to ensure that ay base flow in river produces "surface" run-off)
  groups$qbf_max[ichan]<-1e6

  # estimate for initial base fluxes given a steady state
  flows <- init_fluxes(groups=groups,
                      W=weights,             # flux distribution matrix
                      ichan=ichan,           # channel identifiers
                      FLOW_NAMES=FLOW_NAMES, # names
                      qt0=qt0,               # initial recharge / discharge, assummed equal in steady state, taking into account evapotranspiration
                      dtt=dtt)

  # storage deficits estimated from relationship with base flows estimated above
  stores <- init_stores(groups, flows, ichan,
                         qt0)

  # initial determined from storages vs. deficit

  # index in time series after which to start showing results
  i.start <- which(tms >= sim.start)[1]

  if(is.na(i.start))
  {stop(paste("Error determining start of simulation period "))} ## PJS remove graphics, dispStart))}

  # structure to save the average storage deficit. ae, river input, water balance etc etc
  dat <- data.frame("ae"=rep(0,nmax))
                   # "wb"=rep(0,nmax),
                    #"qr"=rep(0,nmax),
                    #"qriv.in"=rep(0,nmax),   # riveer input fluxes
                    #"ae"=rep(0,nmax))
  data.out <- xts(dat, order.by=tms)

  qriver <- xts(matrix(rep(0,nmax,each=nchan), ncol=nchan), order.by=tms)
  hsus <- groups[-ichan,]
  nhsu <- nrow(hsus)

  # array to hold state of response unit fluxes through time
  fluxes <- array(dim=c(nmax, ngroups, length(FLOW_NAMES)));

  # array that holds state of response unit storages through simulation
  storages <- array(dim=c(nmax, ngroups, length(STORE_NAMES)));
  # graphics interval expressed in hours: convert to time steps

  # start of outer time stepping logic using a time step at each rainfall observation
  time<- sim.start
  it <- 1

  # initialise excess (overland) flow to 0
  q.ex <- rep(0, nrow(groups))

  # set  first discharge estimates to value specified (total, not specific)
  # extend into lag period whilst flux is stil in transit. If this isn't done then it appears that
  # there is a period when no flux enters and we see a initial dip in the hydrograph
  nroute <- length(routing)
  flows1 <- flows
  # steady state - (total) flow into the channel equal those at the outlet
  flows1$qin[ichan] <- qt0*sum(groups$area)

  Qr1 <- Qr
  # distribute initial flow to future flows so that each time step has same
  for(it1 in 1:nroute)
  {
    Qr1 <- route.channel.flows(groups, flows1, stores=stores,
    													 delays=routing,
                                weights, Qr1, it=it1, dt, ichan)
  }
  Qr[1:nroute,]<- qt0*sum(groups$area)

  # remove distributed flows from initial steps of river outlet series
  Qr[1:nroute,] <- Qr[1:nroute,] - Qr1[1:nroute,]

# set channel deficts to a typical value
	stores[ichan,]$sd <- groups[ichan, ]$sd_max/2

	# initialise surface excess storage to zero
	surf.ex <- stores$ex
	pe.dist <- xts(order.by=index(pe), matrix(rep(pe, ngroups), ncol=ngroups))

  # open graphics windows, calculate interval etc
  ## PJS remoove graphics disp.par <- setup.output(disp.par, rain=rain, qobs=qobs)

  rain <- as.xts(rain)
  pe <- as.xts(pe)
#  pe[which(rain>0)]<-0

  # ensure that if it's raining there is nil evapotranspiration
#  tm.rain <- index(rain)[which(rain>0)]
#  pe[tm.rain]<- 0

  # trim the simulation time to the time series
  sim.end <- as.POSIXct(min(max(index(pe)), max(index(rain)), sim.end))
  cur.env <-  environment()
  # copy values to the specified environment, typicaly the  current environment in
  # the calling context
  if(!is.null(calling.env))
  {
    for(n in ls(cur.env, all.names=TRUE))
    {
      # assign name-value in this envir to that specified
      assign(n, get(n, cur.env), calling.env)
    }
  }

  # return as environment that can be loaed into the main routine
  return(cur.env)
}
