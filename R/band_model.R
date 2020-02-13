#' Simplifies a dynamic TOPMODEL by merging computational bands
#'
#' @description take a model and tries to merge the computational bands to produce a simplier description
#'
#' @param model a RasterBrick as created by create_model
#' @param nband number of bands to merge
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
    grp <- as.numeric( cut(hs$band,cuts) )
    if(any(!is.finite(grp))){
        stop("Cut does not assign each HRU to a new band")
    }
    grp <- 0.5*(hs$class+grp)*(hs$class+grp+1)+grp # unique new uid by cantor pair
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
        band = as.numeric(by(hs$band,hs$new_id,min)),
        precip_series="unknown",
        pet_series="unknown",
        qex_max="qex_max_default",
        srz_max="srz_max_default",
        srz_0="srz_0_default",
        ln_t0="ln_t0_default",
        m="m_default",
        td="td_default",
        tex="tex_default",
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
    
    Dsz <- t(K)%*%model$Dsz%*%A%*%K ## fractions combined with area
    Asz <- t(K)%*%A%*%K # area in each
    Dsz <- Dsz %*% Diagonal(ncol(Asz),1/diag(Asz)) # explicit inverse since diagonal
    
    Dex <- t(K)%*%model$Dex%*%A%*%K
    Aex <- t(K)%*%A%*%K # area in each
    Dex <- Dex %*% Diagonal(ncol(Aex),1/diag(Aex))
    
    ## impose limit can;t drain to same band - but with same area draining out
    tDsz <- colSums(Dsz) # current sum of fractions
    tDex <- colSums(Dsz) # current sum of fractions

    ## remove elements not needed
    for(ii in hs2$band){
        idx <- hs2$id[hs2$band==ii]
        for(jj in idx){
            Dsz[idx,jj] <- 0
            Dex[idx,jj] <- 0
        }
    }

    ## standardise
    Dsz <- Dsz %*% Diagonal(ncol(Dsz),tDsz/colSums(Dsz))
    Dex <- Dex %*% Diagonal(ncol(Dex),tDsz/colSums(Dex))

    ## put back into model
    model$hillslope <- hs2
    model$Dex <- Dex
    model$Dsz <- Dsz

    return(model)
}

