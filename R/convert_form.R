#' Contains functions for converting between the storage form and list form used in simulation
#'
#' @description TODO
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @details TODO
#' 
#' @return TODO
#'
#' @name convert_form
NULL

#' @rdname convert_form
## convert a data frame to a storage list
df_to_store_list <- function(df,type=c("hillslope","channel"),use_states=FALSE){
    
    type <- match.arg(type)
    
    ## get a list of all variables that should exist
    vrbl <- required_vrbl(type)
        
    ## initialise depending upon type
    out <- list(attr=list(),
                input=list(),
                param=list(),
                state=list(),
                flux=list())
    
    for(ii in vrbl$attr){
        out$attr[[ii]] <- unname( df[,ii] )
    }
    for(ii in vrbl$input){
        out$input[[ii]] <-  rep(0,nrow(df))
    }
    for(ii in vrbl$param){
        out$attr[[ii]] <- unname( df[,ii] )
        out$param[[ii]] <- unname( model$param[ df[,ii] ] )
    }
    for(ii in vrbl$state){
        if(use_states){
            out$state[[ii]] <- unname( df[,ii] )
        }else{
            out$state[[ii]] <- rep(0,nrow(df)) ## initialise all stores and fluxes as 0
        }
    }
    for(ii in vrbl$flux){
        out$flux[[ii]] <- rep(0,nrow(df)) ## initialise all stores and fluxes as 0
    }
    
    return(out)
}


#' @rdname convert_form
#' @export
store_to_sim <- function(model, use_states=FALSE){

    for(ii in c("hillslope","channel")){
        model[[ii]] <- df_to_store_list(model[[ii]],ii,use_states)
    }

    ## collapsed matrices as lists
    ## ith element of list contains the parents (and fractions) for the
    ## hsu with id = i
    tmp <- as.matrix(summary(model$Fsf))
    tmp <- lapply(by(tmp[,2:3],tmp[,1],unique,simplify=FALSE),as.list)
    model$sf <- rep(list(NULL),max(as.numeric(names(tmp))))
    model$sf[as.numeric(names(tmp))] <- tmp
    
    tmp <- as.matrix(summary(model$Fsz))
    tmp <- lapply(by(tmp[,2:3],tmp[,1],unique,simplify=FALSE),as.list)
    model$sz <- rep(list(NULL),max(as.numeric(names(tmp))))
    model$sz[as.numeric(names(tmp))] <- tmp

    to_keep <- c("hillslope","channel","sz","sf","param","gauge","point_inflow")

    return(model[to_keep])
}

#' @rdname convert_form
## convert the store back to a table
sim_to_store <- function(model,with_matrix=TRUE){
    
    for(ii in c("hillslope","channel")){
        model[[ii]] <- cbind( as.data.frame(model[[ii]]$attr,stringsAsFactors=FALSE),
                             as.data.frame(model[[ii]]$state,stringsAsFactors=FALSE))
    }

    ## recreate matrices
    tmp <- rep(list(NULL),length(model$sz))
    for(ii in 1:length(model$sz)){
        if(length(model$sz[[ii]]$x)>0){
            tmp[[ii]] <- cbind( i=ii,
                               as.data.frame(model$sz[[ii]],stringsAsFactors=FALSE) )
        }
    }
    tmp <- do.call(rbind,tmp)
    #browser()
    model$Fsz <- Matrix::sparseMatrix(i=tmp$i,j=tmp$j,x=tmp$x)

    tmp <- rep(list(NULL),length(model$sf))
    for(ii in 1:length(model$sf)){
        if(length(model$sf[[ii]]$x)>0){
            tmp[[ii]] <- cbind( i=ii,
                               as.data.frame(model$sf[[ii]],stringsAsFactors=FALSE) )
        }
    }
    tmp <- do.call(rbind,tmp)
    model$Fsf <- Matrix::sparseMatrix(i=tmp$i,j=tmp$j,x=tmp$x)

    to_keep <- c("hillslope","channel","Fsz","Fsf","param","gauge","point_inflow")

    return(model[to_keep])
}
