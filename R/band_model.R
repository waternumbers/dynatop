#' Simplifies a dynamic TOPMODEL by merging computational bands
#'
#' @description take a model and tries to merge the computational bands to produce a simplier description
#'
#' @param model a RasterBrick as created by create_model
#' @param cuts number of bands to merge
#'
#' @return a simplified model
#' @export
band_model <- function(model, cuts){

    ## see if the input mode is valid
    ##dynatop::check_model(model)

    ##model <- readRDS("./dynatop/data/Swindale_ordered.rds")
    ##cuts <- seq(0,2760,by=10)
    hs <- model$hillslope
    ##hs$class <- hs$split_id

    ## create new uids
    sz_grp <- as.numeric( cut(hs$sz_band,cuts) )
    sf_grp <- as.numeric( cut(hs$sf_band,cuts) )
    if(any(!is.finite(sz_grp)) | any(!is.finite(sf_grp))){
        stop("Cut does not assign each HRU to a new band")
    }
    # unique new uid by cantor pair
    grp <- 0.5*(hs$class+sz_grp)*(hs$class+sz_grp+1)+sz_grp 
    grp <- 0.5*(grp+sf_grp)*(hs$class+sf_grp+1)+sf_grp 
    ugrp <- unique(grp)
    
    hs[,'new_id'] <- NA
    for(ii in 1:length(ugrp)){
        hs[grp==ugrp[ii],'new_id'] <- ii
    }
    hs[,'new_id']  <-  hs[,'new_id'] + max(model$channel$id)

    ## create the new table
    hs2 <- data.frame(
        id = unique(hs[,'new_id']),
        area = as.numeric(by(hs$area,hs$new_id,sum)),
        atb_bar = as.numeric(by(hs$atb_bar*hs$area,hs$new_id,sum)),
        s_bar = as.numeric(by(hs$s_bar*hs$area,hs$new_id,sum)),
        delta_x = as.numeric(by(hs$delta_x,hs$new_id,sum)),
        class = as.numeric(by(hs$class,hs$new_id,unique)),
        sz_band = as.numeric(by(hs$sz_band,hs$new_id,min)),
        sf_band = as.numeric(by(hs$sf_band,hs$new_id,min)),
        precip="unknown",
        pet="unknown",
        q_sfmax="q_sfmax_default",
        s_rzmax="s_rzmax_default",
        s_rz0="s_rz0_default",
        ln_t0="ln_t0_default",
        m="m_default",
        t_d="t_d_default",
        t_sf="t_sf_default",
        stringsAsFactors=FALSE
    )
    hs2$s_bar <- hs2$s_bar/hs2$area
    hs2$atb_bar <- hs2$atb_bar/hs2$area

    ## adapt the redistribution matrices
    K <- sparseMatrix(i=c(model$channel$id,hs$id),
                      j=c(model$channel$id,hs$new_id),
                      x=1,
                      dims=c(max(hs$id),max(hs2$id)))

    
    A <- Diagonal(max(hs$id),c(model$channel$area,hs$area))
    
    Dsz <- t(K)%*%model$Fsz%*%A%*%K ## fractions combined with area
    Asz <- t(K)%*%A%*%K # area in each
    Dsz <- Dsz %*% Diagonal(ncol(Asz),1/diag(Asz)) # explicit inverse since diagonal
    
    Dsf <- t(K)%*%model$Fsf%*%A%*%K
    Asf <- t(K)%*%A%*%K # area in each
    Dsf <- Dsf %*% Diagonal(ncol(Asf),1/diag(Asf))
    
    ## impose limit can't drain to same or lower band - but with same area draining out
    tDsz <- colSums(Dsz) # current sum of fractions
    tDsf <- colSums(Dsz) # current sum of fractions

    ## remove elements which can't be evlauated since drain to the same band of lower...
    for(ii in unique(hs2$sz_band)){
        idx <- hs2$id[hs2$sz_band<=ii]
        for(jj in idx){
            Dsz[idx,jj] <- 0
        }
    }
    for(ii in hs2$sf_band){
        idx <- hs2$id[hs2$sf_band<=ii]
        for(jj in idx){
            Dsf[idx,jj] <- 0
        }
    }

    ## standardise
    Dsz <- Dsz %*% Diagonal(ncol(Dsz),tDsz/colSums(Dsz))
    Dsf <- Dsf %*% Diagonal(ncol(Dsf),tDsf/colSums(Dsf))
    browser()

    ## handle areas with no drainage - add to the HSU of the same class downslope
    idx <- intersect( which(colSums(Dsf)==0), hs2$id ) # this in an id
    for(ii in idx){
        ## get class and band
        cls <- hs2$class[ hs2$id == ii ]
        bnd <- hs2$sf_band[ hs2$id == ii ]
        
        ## idintify replacement class
        jdx <- hs2$id[ hs2$class ==cls ]
        jdxb <- hs2$sf_band[ hs2$class ==cls ]
        jdx <- jdx[ jdxb > bnd ]
        jdxb <- jdxb[ jdxb > bnd ]
        jdx <- jdx[which.min(jdxb)] # this is the replacement id

        ## get location as logical
        iii <- hs2$id == ii
        jdx <- hs2$id == jdx
        tmp <- hs2$area[jdx] + hs2$area[jdx]
        hs2$atb_bar[jdx] <- ( hs2$atb_bar[jdx]*hs2$area[jdx] +
                              hs2$atb_bar[ii]*hs2$area[iii] )/tmp
        hs2$s_bar[jdx] <- ( hs2$s_bar[jdx]*hs2$area[jdx] +
                            hs2$s_bar[ii]*hs2$area[iii] )/tmp
        hs2$delta_x[jdx] <- hs2$delta_x[jdx] + hs2$delta_x[iii]
        hs2$area[jdx] <- tmp
        Dsz[jdx,] <- Dsz[jdx,] + Dsz[iii,]
    }
    ## drop columns and row
    Dsz
    


    
    ## standardise
    Dsz <- Dsz %*% Diagonal(ncol(Dsz),tDsz/colSums(Dsz))
    Dsf <- Dsf %*% Diagonal(ncol(Dsf),tDsf/colSums(Dsf))

    ## handle areas with no drainage - add to the HSU of the same class downslope
    idx <- intersect( which(colSums(Dsf)==0), hs2$id ) # this in an id
    for(ii in idx){
        ## get class and band
        cls <- hs2$class[ hs2$id == ii ]
        bnd <- hs2$sf_band[ hs2$id == ii ]
        
        ## idintify possible d/s point
    }
    
    ##
    
    ## put back into model
    model$hillslope <- hs2
    model$Fsf <- Dsf
    model$Fsz <- Dsz

    return(model)
}

