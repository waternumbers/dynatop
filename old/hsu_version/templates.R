#' Functions for generating and checking lists for different HRU types
#'
#' @description Functions for creating, checking, initialising and evolving a hillslope HRU.
#'
#' @param type a string containing the type required
#' @param hru one or a list of dynamic TOPMODEL HRUs
#'
#' @return Either a template list, or if checking will fail with error message
#'
#' @rdname hru
#' @examples
#' tmp <- template_hru()
#' check_hru(tmp)
#' ## add a mock input to combine
#' tmp[['mux']]$new_input <- tmp[['ex']]$output
#' check_hru(tmp)

template_output <- function(){
    c("precip"=0,"pet"=0,"lsz"=0,"lex"=0,"qch"=0)
}

#' @name hru
#' @export 
template_hru <- function(type=NULL){
    ## descriptor has properties:
    ## 1. one of length, min_length, names
    ## 2. class
    ## 3. hasNames

    output <- template_output()

    out <- list(
        "ex"=list("id"=numeric(1),
                  "type"="ex",
                  "series"=setNames(character(length(output)),names(output)),
                  "output"=output),
        "mux"=list("id"=numeric(1),
                   "type"="mux",
                   "output"=output),
        "channel"=list("id"=numeric(1),
                       "type"="channel",
                       "uptree"=character(1),
                       "prop"=c("area"=0),
                       "output"=output),
        "hillslope"=list("id"=numeric(1),
                         "type"="hillslope",
                         "uptree"=character(1),
                         "prop"=c("area"=0,"atb_bar"=0,"s_bar"=0,"delta_x"=0),
                         "param"=c("qex_max"=0,"srz_max"=0,"srz_0"=0,
                                   "ln_t0"=0,"m"=0,"td"=0),
                         "state"=c("ex"=0,"rz"=0,"uz"=0,"sz"=0,
                                   "lsz_in"=0,"lsz"=0,"lsz_max"=0),
                         "output"=output))
    if(length(type)>0){
        if(any(names(out)==type)){
            out <- out[[type]]
        }else{
            stop("Not a valid type")
        }
    }
    return(out)
}

#' @name hru
#' @export
check_hru <- function(hru){
    if( all(c("type","id") %in% names(hru)) ){
        ## probably a single hru so make a list
        hru <- list(hru)
        names(hru) <- "single"
    }
   
    ## check unique names
    if( length(unique(names(hru))) != length(hru) ){
        stop("Not unique HRU names")
    }

    ## check individual HRUS
    templates <- template_hru()
    #browser()
    for(ii in names(hru)){
        
        h <- hru[[ii]]
        if( !any(names(templates)==h$type) ){
   
            stop(paste("HRU",ii,"is not a valid type"))
        }
        
        t <- templates[[h$type]]
        if(!all(names(t) %in% names(h))){
            stop(paste("Not all components presetn in HRU:",ii))
        }
        estr <- NULL
        for(jj in names(t)){
            if( !(class(t[[jj]]) %in% class(h[[jj]])) ){
                estr <- c(estr,paste(jj,"wrong class in HRU",ii))
            }
            if( length(t[[jj]]) != length(h[[jj]]) ){
                c(estr,paste(jj,"wrong length in HRU",ii))
            }
            if( length(names(t[[jj]])) > 0  &&
                !all(names(t[[jj]]) %in% names(h[[jj]])) ){
                c(estr,paste(jj,"should all be named in HRU",ii))
            }    
        }

        ## if mux may have other elements that are vectors for weights
        if( h$type=="mux" ){
            for(jj in setdiff(names(h),names(t))){
                if( !all( class(h[[jj]]) == class(t[['output']]) ,
                         names(h[[jj]])==names(t[['output']]) ) ){
                    c(estr,paste(jj,"is incorrectly specified in in HRU",ii))
                }
            }
        }
        
        if( length(estr) > 0 ){
            stop(paste(estr,collapse="\n"))
        }
    }
}
