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
    }
    
    ## work out the sequences for computing the lateral flux bands these are the index in the hillslope vectors NOT the id
    out$sqnc <- list(sf=list(),sz=list())
    for(ii in names(out$sqnc)){
        bnd <- switch(ii,
                      sf = sapply(model$hillslope$sf_dir,function(x){x$band}),
                      sz = sapply(model$hillslope$sz_dir,function(x){x$band})
                      )
        out$sqnc[[ii]] <- by(1:length(model$hillslope$id),bnd,c)
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

