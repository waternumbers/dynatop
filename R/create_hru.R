#' Initialises the hru states for a dynatop simulation
#'
#' @description Initialise the states of the different hru types. these can be time varying such as store of fluxes or fixed properties such as parmeter values and areas
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return A list of lists, one list for each HRU type which contains vectors of states and properties for use in computation
#'
#' @name initialise_dynatop
#' @export
create_hru <- function(model){

    ## function to convert tables to lists
    tbl2lst <- function(tbl,type,param){
        print(type)
        hru <- switch(type,
                     "hillslope" = create_hillslope(tbl[['id']]),
                     "channel" = create_channel(tbl[['id']]),
                     )
        h0 <- switch(type,
                     "hillslope"=list(prop=c("area","s_bar","atb_bar","delta_x","band"),
                                      param=c("qex_max","srz_max","srz_0","ln_t0","m","td","tex",
                                              "precip_series","pet_series")),
                     "channel"=list(prop=c("area","band"),
                                    param=c("precip_series","pet_series"))
                     )

        for(ii in setdiff(h0$param,c("precip_series","pet_series"))){
            tbl[,ii] <- param[tbl[,ii]]
        }
        
        #browser()
        ## basic info
        #hru <- rep(list(h0),nrow(tbl))
        #lst <- lapply(tbl[['id']],function(x,y){list(id=x)})
        #hru <- Map(c,hru,lst)

        # hru <- lapply(tbl[['id']],function(x,y){list(id=x,type=y)},y=type)
        ## get properties and parameters
        
        for(jj in c("prop","param")){
            print(jj)
            if(length( h0[[jj]] ) <1 ){
                next
            }
            
            nms <- intersect(h0[[jj]],names(tbl))
            lst <- rep(list(NULL),nrow(tbl))
            for(ii in nms){   
                tmp <- lapply(tbl[[ii]],
                              function(x,nm){setNames(list(x),nm)},
                              nm=ii)
                lst <- Map(c,lst,tmp)
            }
            lst <- lapply(lst,function(x,y){z <- list();z[[y]] <- x;return(z)},y=jj)
            #browser()
            hru <- Map(c,hru,lst)
        }

        return(hru)
    }
    
    ## convert tables to lists
    hru <- c(tbl2lst(model$hillslope,"hillslope",model$param),
             tbl2lst(model$channel,"channel",model$param))

    ## sort hrus by id for applying redistribution
    id <- sapply(hru,function(x){x$id})
    hru <- hru[order(id)]

    ## extract redistribution tables
    redist <- list()
    redist$lex <- as.matrix(summary(model$Dex))
    redist$lsz <- as.matrix(summary(model$Dsz))

    ## work out the reweighting - works since hru sorted by id
    area <- sapply(hru,FUN=function(x){x$prop$area})
    for(jj in names(redist)){
        idx <- split(redist$lex[,1],redist$lex[,2])
        w <- split(redist$lex[,3],redist$lex[,2])
        kdx <- sapply(split(redist$lex[,2],redist$lex[,2]),unique)
        for(ii in 1:length(kdx)){
            k <- kdx[ii]
            hru[[k]]$output[[jj]]$id <- idx[[ii]]
            hru[[k]]$output[[jj]]$w <- w[[ii]]*area[k]/area[idx[[ii]]]
        }
    }

    return(hru)
}

