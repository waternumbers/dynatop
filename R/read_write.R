#' Functions for reading and writing dynatop models
#'
#' @description Helper functions for reading adn writing dynatop models output by dynaGIS.
#'
#' @param model A model object to output
#' @param file path to the model file
#'
#' @return A SpatVect containing polygons of the channel network, with at least the following properties: uid, width, depth, slope, startNode and endNode.
#'
#' @details While the basic dynatop model format is a text file delimited by ";" there are some complications in handling
#'   - the json strings containing the parameters
#'   - keeping the quote approriate for easy import in (Q)GIS
#'
#' @examples
#' mdl <- read_model(system.file("extdata", "mdl", "SwindaleModel.csv"))
#' tmpFile <- tempfile()
#' on.exit(unlink(tmpFile))
#' write_model(mdl,tmpFile)
#' @export
write_model <- function(model,file){
    write.table(model, file,
                sep=";",row.names=FALSE,
                na="",quote=FALSE,
                fileEncoding = "UTF-8")
}
#' @rdname write_model
#' @export
read_model <- function(file){
    read.table(file,
               header=TRUE,
               sep=";",quote="'",
               fileEncoding = "UTF-8")
}
