#' Create an hrus based on the output of dynatopGIS
#'
#' @description This creates and populates the hru table, properties and parameters of the hru
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return A list of lists, one list for each HRU type which contains vectors of states and properties for use in computation
#'
#' @name initialise_dynatop
#' @export
model_to_hru <- function(model){

    ## function to convert tables to lists
    tbl2lst <- function(tbl,type,param){
        print(type)
        ## template of the HRU
        h0 <- template_hru(type)

        ## create vector
        hru <- rep(list(h0),nrow(tbl))

        ## copy id
        for(jj in 1:length(hru)){
            hru[[jj]]$id <- tbl$id[jj]
        }
        
        ## add in parameter values to the table
        for(ii in names(h0$param)){
            #print(ii)
            p0 <- param[tbl[,ii]]
            for(jj in 1:length(hru)){
                hru[[jj]]$param[ii] <- p0[jj]
            }
        }
        ## add in prop values to the hrus
        for(ii in c("id",names(h0$prop))){
            #print(ii)
            p0 <- tbl[,ii]
            for(jj in 1:length(hru)){
                hru[[jj]]$prop[ii] <- p0[jj]
            }
        }

        names(hru) <- paste0(substr(type,1,1),sapply(hru,function(x){x$id}))

        return(hru)
    }
    
    ## convert hillslope and channel tables to lists
    hru <- c(tbl2lst(model$hillslope,"hillslope",model$param),
             tbl2lst(model$channel,"channel",model$param))

    ## sort hrus by id for applying redistribution
    id <- sapply(hru,function(x){x$id})
    hru <- hru[order(id)]
    redist <- model$Dex + model$Dsz
    redist[redist>0] <- 1
    colnames(redist) <- rownames(redist) <- names(hru)

    ## add additional hrus for each input
    inhru <- do.call(rbind,
                     list(
                         model$hillslope[,c("id",'precip_series','pet_series')],
                         model$channel[,c("id",'precip_series','pet_series')]
                     ))
    inhru <- inhru[order(inhru$id),]
    ex_comb <- unique(inhru) # unique combinations
    
    
    

    ## work out levels
    lvl <- setNames(rep(NA,length(hru)),nms)
    nparents <- rowSums(redist)
    to_add <- nms[nparents==0]
    clvl <- 1
    while(length(to_add)>0){
        print(clvl)
        lvl[ to_add ] <- clvl
        nparents[ to_add ]  <- -1
        nparents <- nparents - rowSums(redist[,to_add,drop=FALSE])
        to_add <- nms[nparents==0]
        clvl <- clvl+1
    }
    
    browser()
    ## extract redistribution tables
    redist <- list()
    redist$lex <- as.matrix(summary(model$Dex))
    redist$lsz <- as.matrix(summary(model$Dsz))
    #redist$comb <- do.call(+,redist)
    
    ## work out parents (uptree) nodes
    nms <- names(hru)
    tmp <- by(redist$comb[,2],redist$comb[,1],unique,simplify=FALSE)
    names(tmp) <- nms[unique(redist$comb[,1])]
    prnt <- setNames(rep(list(NULL),length(hru)),nms)
    prnt[names(tmp)] <- tmp

    browser()
    ## convert parent to names and have number
    prnt <- lapply(prnt,function(x,y){list(p=y[x],n=length(x))},y=nms)
    
    ## initialise vector of order
    lvl <- setNames(rep(NA,length(hru)),nms)

    ## work out the first level of hrus
    flg <- TRUE
    clvl <- 1
        
    


    
    ## work out the reweighting - works since hru sorted by id
    browser()

    
    area <- sapply(hru,FUN=function(x){x$prop['area']})

    
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

