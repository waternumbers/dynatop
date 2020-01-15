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
        hru <- rep(h0,nrow(tbl))

        ## add in parameter values to the table
        for(ii in names(h0$param)){
            p0 <- param[tbl[,ii]]
            for(jj in 1:length(hru)){
                hru[[jj]]$param[ii] <- p0[jj]
            }
        }
        ## add in prop values to the hrus
        for(ii in c("id",names(h0$prop))){
            p0 <- tbl[,ii]
            for(jj in 1:length(hru)){
                hru[[jj]]$prop[ii] <- p0[jj]
            }
        }

        names(hru) <- paste0(substr(type,1,1),sapply(hru,function(x){x$id}))
        
        #browser()
        ## basic info
        #hru <- rep(list(h0),nrow(tbl))
        #lst <- lapply(tbl[['id']],function(x,y){list(id=x)})
        #hru <- Map(c,hru,lst)

        # hru <- lapply(tbl[['id']],function(x,y){list(id=x,type=y)},y=type)
        ## get properties and parameters
        
        ## for(jj in c("prop","param")){
        ##     print(jj)
        ##     if(length( h0[[jj]] ) <1 ){
        ##         next
        ##     }
            
        ##     nms <- intersect(h0[[jj]],names(tbl))
        ##     lst <- rep(list(NULL),nrow(tbl))
        ##     for(ii in nms){   
        ##         tmp <- lapply(tbl[[ii]],
        ##                       function(x,nm){setNames(list(x),nm)},
        ##                       nm=ii)
        ##         lst <- Map(c,lst,tmp)
        ##     }
        ##     lst <- lapply(lst,function(x,y){z <- list();z[[y]] <- x;return(z)},y=jj)
        ##     #browser()
        ##     hru <- Map(c,hru,lst)
        ## }

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

