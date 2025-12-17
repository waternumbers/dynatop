## ## This contains computation of catchment properties not required by dyantop but perhaps useful for historic reasons

##         ## function to do property calculations on an upwards pass (low to high DEM values)
##         ## if we go up in height order then we are working from near the channel to the heighest point
##         ## could add back in flow distances here
##         apply_upward_pass = function(verbose){

##             rq <- c("filled_dem","channel","channel_fraction")
##             stopifnot(
##                 "Not all required input layers have been generated \n Try running sink_fill first" =
##                     all( rq %in% names( private$brk) )
##             )

##             ## rasterize channel band to start
##             rbnd <- terra::rasterize(private$chn, private$brk[["catchment"]],field = "band",touches=TRUE)
##             rbnd <- terra::mask(rbnd,private$brk$catchment) ## ensure channel bands are within the catchment - else later code fails
##             names(rbnd) <- "band"

##             ## load raster layer
##             d <- terra::as.matrix( private$brk[["filled_dem"]], wide=TRUE )
##             chn <- terra::as.matrix( private$brk[["channel"]], wide=TRUE )
##             chn_frc <- terra::as.matrix( private$brk[["channel_fraction"]], wide=TRUE )
##             bnd <- terra::as.matrix( rbnd,  wide=TRUE )

##             if( verbose ){ print("Computing upward pass") }

##             idx <- order(d,na.last=NA) ## search order
##             nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1) ## neighbours

##             ## set up printing variables
##             if(verbose){
##                 print_step <- c(1,rep(round(length(idx)/20,2)),length(idx)) # current, next print, step, total
##             }

##             ## main loop
##             dz <- rep(NA,8)
##             for(ii in idx){

##                 if( is.finite(chn[ii]) ){ ## then cell is a channel
##                     if(chn_frc[ii]<1){
##                         ## since it is mixed cell - partly landuse , partly channel
##                         bnd[ii] <- bnd[ii] + 1
##                     }
##                 }else{
##                     ## it is not a channel
##                     jdx <- ii+delta
##                     dz[] <- d[ii] - d[jdx]
##                     is_lower <- is.finite(dz) & dz>0
##                     bnd[ii] <- max( bnd[jdx[is_lower]] ) + 1
##                 }

##                 if(verbose){
##                     print_step[1] <- print_step[1] + 1
##                     if( print_step[1] > print_step[2] ){
##                         cat(round(100*print_step[1] / print_step[4],1),
##                             "% complete","\n")
##                         print_step[2] <- print_step[2] + print_step[3]
##                     }
##                 }

##             }

##             ## save
##             terra::values(rbnd) <- bnd
##             private$brk <- c(private$brk,rbnd)
##             private$save_project()


##         },
##         apply_downward_pass = function(min_grad,verbose){

##             if( verbose ){ print("Loading data for downward pass") }
##             rq <- c("filled_dem","channel")
##             stopifnot(
##                 "Not all required input layers have been generated \n Try running sink_fill first" =
##                     all( rq %in% names( private$brk) )
##             )

##             ## load raster layer
##             d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE)
##             ch <- terra::as.matrix( private$brk[["channel"]] , wide=TRUE)
##             ch_frc <- terra::as.matrix( private$brk[["channel_fraction"]] , wide=TRUE)

##             if( verbose ){ print("Setting up computation") }

##             ## distance between cell centres
##             rs <- terra::res( private$brk )
##             dxy <- rep(sqrt(sum(rs^2)),8)
##             dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
##             dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs) ## assumes square
##             nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)

##             ## initialise output
##             gr <- upa <- atb <- d*NA
##             upa <- prod(rs)*(1-ch_frc) ## initialise upslope area from resolution

##             idx <- order(d,na.last=NA,decreasing=TRUE) ## search order

##             ## set up printing variables
##             if(verbose){
##                 print_step <- c(1,rep(round(length(idx)/20,2)),length(idx)) # current, next print, step, total
##             }

##             if( verbose ){ print("Computing hillslope") }

##             uA <- private$chn$area ## upsteam areas for channels

##             ## loop downslope
##             w <- rep(0,8)
##             for(ii in idx){
##                 if( is.finite(ch[ii]) ){
##                     ## pass on upslope area to channel
##                     uA[ ch[ii] ] <- uA[ ch[ii] ] + upa[ii]
##                     if( ch_frc[ii] < 1 ){
##                         ## mixed cells
##                         ## work out gradient from cells flowing in
##                         jdx <- ii+delta
##                         grd <- (d[ii]-d[jdx])/dxy
##                         gcl <- grd*dcl
##                         is_higher <- is.finite(gcl) & gcl<0 #& is.finite(cjdx) & cjdx==ctch[ii]
##                         if( any(is_higher) ){
##                             sum_gcl <- sum( gcl[is_higher] )
##                             sum_dcl <- sum( dcl[is_higher] )
##                             gr[ii] <- max(-sum_gcl / sum_dcl,min_grad)
##                         }else{ ## nothing flows into the cell
##                             gr[ii] <- min_grad
##                         }
##                         atb[ii] <- log(upa[ii]/gr[ii])
##                     }else{ ## pure water cell
##                         upa[ii] <- NA
##                         gr[ii] <- NA
##                         atb[ii] <- NA
##                     }

##                 }else{

##                     ## it is not a channel
##                     w[] <- 0
##                     jdx <- ii+delta
##                     grd <- (d[ii]-d[jdx])/dxy
##                     gcl <- grd*dcl
##                     is_lower <- is.finite(gcl) & gcl>0 #& is.finite(cjdx) & cjdx==ctch[ii]
##                     sum_gcl <- sum( gcl[is_lower] )
##                     sum_dcl <- sum( dcl[is_lower] )
##                     w[is_lower] <- gcl[is_lower] / sum_gcl

##                     if( !any(w>0) ){ stop(paste("Cell",ii,"is a hillslope cell with no outflows")) }

##                     ## gradient
##                     gr[ii] <- max(sum_gcl / sum_dcl,min_grad)
##                     ## topographic index
##                     atb[ii] <- log(upa[ii]/gr[ii]) #log( upa[ii] / sum(gcl) )
##                     ## propogate area downslope
##                     upa[ jdx ]  <- upa[ jdx ] + w*upa[ii]
##                 }

##                 ## verbose output here
##                 if(verbose){
##                     print_step[1] <- print_step[1] + 1
##                     if( print_step[1] > print_step[2] ){
##                         cat(round(100*print_step[1] / print_step[4],1),
##                             "% complete","\n")
##                         print_step[2] <- print_step[2] + print_step[3]
##                     }
##                 }
##             }

##             if( verbose ){ print("Computing channel") }
##             sN <- private$chn$startNode
##             eN <- private$chn$endNode


##             ## merge upslope areas into the channel object
##             ##ch_upa <- tapply(upa,ch,sum)
##             ##ch_upa <- ch_upa[setdiff(names(ch_upa),"NaN")]
##             ##idx <- match(names(ch_upa),paste(private$chn$id)) #,names(ch_upa))
##             ##uA[idx] <- as.numeric(ch_upa)

##             ## remove channel area bit from hillslope upslope area
##             ##upa[ch_frc==1] <- NA

##             ## compute catchment area to each reach
##             for(ii in length(sN):1){
##                 idx <- sN == eN[ii]
##                 uA[idx] <- uA[idx] + uA[ii] / sum(idx) ## TO CHECK not sure how channel routing fractions originally done
##             }
##             stopifnot(
##                 "All channel upstream areas should be finite" = all(is.finite(uA)),
##                 "All channel upstream areas should be non-negative" = all(uA>0)
##             )

##             private$chn$upstream_area <- uA

##             ## save output
##             private$brk <- c(private$brk,
##                              terra::rast( private$brk[["dem"]], names="gradient", vals=gr ),
##                              terra::rast( private$brk[["dem"]], names="upslope_area", vals=upa ),
##                              terra::rast( private$brk[["dem"]], names="atb", vals=atb )
##                              )

##             private$save_project(chn=TRUE)
##         },
##         ## split_to_class
##         apply_classify = function(base_layer,cuts,layer_name){

##             rq <- c("channel",base_layer)
##             stopifnot(
##                 "Missing channel or base layer to classify" = all( rq %in% names(private$brk) ),
##                 "layer_name has zero length" = nchar(layer_name)>0,
##                 "layer_name is already used" = !(layer_name %in% names(private$brk)),
##                 "layer_name is reserved" = !(layer_name %in% names(private$reserved_layers))
##             )

##             ## load base layer and mask out channel
##             ##x <-  terra::mask( private$brk[[base_layer]], private$brk[["channel"]], inverse=TRUE)
##             x <- private$brk[[base_layer]]

##             ## work out breaks
##             brk <- as.numeric(cuts)
##             rng <- as.numeric( terra::global(x, fun="range",na.rm=TRUE) )
##             if( length(brk)==1 ){
##                 ## this defines brks in the same way as cut would otherwise
##                 brk <- seq(rng[1],rng[2],length=brk+1)
##             }else{
##                 brk <- sort(brk)
##                 if( brk[1] > rng[1]){ brk <- c(rng[1],brk) }
##                 if( tail(brk,1) < rng[2]){ brk <- c(brk,rng[2]) }
##             }
##             if( any(is.na(brk)) ){ stop("NA value in brk") }
##             M <- cbind( head(brk,-1), tail(brk,-1), 1:(length(brk)-1) )

##             ## cut the raster and save

##             return( terra::classify(x,rcl=M,include.lowest=TRUE,names=layer_name) )
##         },
##         ## split_to_class
##         apply_combine_classes = function(layer_name,pairs,burns){

##             stopifnot(
##                 "Missing layers in pairs list" = all(pairs %in% names(private$brk)),
##                 "Missing layers in burns list" = all(burns %in% names(private$brk)),
##                 "There must be at lest one entry in pairs" = length(pairs)>0
##             )

##             x <- terra::as.matrix( private$brk[[ pairs ]] )
##             idx <- rowSums(is.finite(x)) == ncol(x) ## thses are the valid cells
##             xstr <- apply(x,1,function(r){paste(r,collapse="_")})

##             ## add burns sequentally
##             if(length(burns)>0){
##                 y <- terra::as.matrix( private$brk[[ burns ]] )
##                 y[y<=0] <- NA
##                 for(ii in 1:ncol(y)){
##                     jdx <- is.finite(y[,ii])
##                     xstr[jdx] <- paste("burn",ii,y[jdx,ii],sep="_")
##                 }
##             }

##             ## make numeric class
##             uxstr <- unique(xstr[idx])
##             ux <- setNames(1:length(uxstr),uxstr)
##             z <- rep(NA,nrow(x))
##             z[idx] <- ux[xstr[idx]]

##             return( terra::rast( private$brk[["dem"]], names=layer_name, vals=z ) )
##         },
