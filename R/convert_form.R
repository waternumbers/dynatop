#' Functions for converting between the model forms used in simulation
#'
#' @description Function that convert a dynatop model into forms that are
#' more useful for simulation and then extract properties from them.
#'
#' @param model a dynamic TOPMODEL model list object
#' @param used_by which function the transformed model is to be used by
#' @param obj a transformed object
#' @param type of HSU for which states are to be returned (either hillslope or channel)
#'
#' @details These functions are designed to be used internally by other functions within the package. In general \code{convert_form} with an appropriate \code{used_by} should be used although functions are available directly. \code{get_states} extract the states from a converted object.
#'
#' @return Either the data object used in the simulation or a data.frame of states along with the HSU ID.
#'
#' @name convert_form
NULL

#' @rdname convert_form
#' @export
## general function for calling
convert_form <- function(model,used_by="dynatop"){
    used_by <- match.arg(used_by)

    switch(used_by,
           dynatop = convert_dynatop(model),
           stop("Conversion type error"))
}

#' @rdname convert_form
#' @export
## conversion for use in dynatop
convert_dynatop <- function(model){

    check_model(model)

    ## initialise output as a list
    out <- list()

    ## get the descriptions of the variabels to be returned
    desc <- list(hillslope = model_description("hillslope",TRUE,TRUE),
                 channel = model_description("channel",TRUE,TRUE))
    ## TODO - trim desc so only returns what is needed

    ## unpack the matrices into lists
    tmp <- as.matrix(summary(model$Fsf))
    Fsf <- lapply(by(tmp[,2:3],tmp[,1],unique,simplify=FALSE),as.list)
    tmp <- as.matrix(summary(model$Fsz))
    Fsz <- lapply(by(tmp[,2:3],tmp[,1],unique,simplify=FALSE),as.list)

    ## convert into lists
    for(tbl in names(desc)){
        out[[tbl]] <- list()
        for(ii in 1:nrow(desc[[tbl]])){
            nm <- desc[[tbl]]$name[ii]
            out[[tbl]][[ nm ]] <- switch(desc[[tbl]]$role[ii],
                                        attribute = unname( model[[tbl]][[nm]] ),
                                        data_series = unname( model[[tbl]][[nm]] ),
                                        parameter = unname( model$param[ model[[tbl]][[nm]] ]),
                                        state = unname( model[[tbl]][[nm]] ),
                                        tmp = rep(0, nrow(model[[tbl]]))
                                        )
        }
        ## copy the upstream locations to the lists
        out[[tbl]]$Fsz <- Fsz[paste(out[[tbl]]$id)]
        out[[tbl]]$Fsf <- Fsf[paste(out[[tbl]]$id)]
    }

    ## work out the sequences for computing the lateral flux bands these are the index in the vectors NOT the id
    out$sqnc <- list(sf_band=list(),sz_band=list())
    for(ii in names(out$sqnc)){
        tmp <- sort(unique(c(model$hillslope[[ii]],
                             model$channel[[ii]]))) # sorted list of unique bands
        for(jj in 1:length(tmp)){
            out$sqnc[[ii]][[jj]] <- list(
                hillslope = which(model$hillslope[[ii]] == tmp[jj]),
                channel = which(model$channel[[ii]] == tmp[jj]))
        }
    }

    ## storage for lateral fluxes stored as volumes
    out$lateral_flux <- list(sf = rep(0,max(c(model$hillslope$id,model$channel$id))),
                             sz = rep(0,max(c(model$hillslope$id,model$channel$id))))

    return(out)
}

#' @rdname convert_form
#' @export
## convert a data frame to a storage list
get_states <- function(obj,used_by="dynatop",type=c("hillslope","channel")){
    used_by <- match.arg(used_by)

    type <- match.arg(type)

    stt <- model_description(type,TRUE)
    stt <- c("id",stt$name[stt$role=="state"])
    out <- as.data.frame( obj[stt], stringsAsFactors=FALSE )
    return(out)
}

