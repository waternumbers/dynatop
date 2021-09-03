#' Start profiler
#'
#' @description start gperftools
#'
#' @param fileName name of file to write to
#' @export
startProfiler <- function(fileName){
    start_profiler(fileName)
}

#' Stop profiler
#'
#' @description stop gperftools
#'
#' @export
stopProfiler <- function(){
    stop_profiler()
}
