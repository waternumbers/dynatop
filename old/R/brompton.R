#' Topographic and observation data for running Dynamic TOPMODEL.
#'
#'@description
#' Brompton is a small agricultural catchment in N.Yorkshire, UK. Its eastern
#' edges rise in the North Yorks Moors and it drains southwards, becoming
#' North Beck before joining the Wiske in Northallerton.
#'
#' In the late 19th century the area upstream of Water End was drained and turned
#' over to arable cultivation and has since suffered from infrequent, but severe,
#' flooding due to intense rainfall from weather systems moving in from the
#' North Sea.The last event hat flooded the village was in November 2012; flooding
#' was narrowly avoided in the storms of December 2015 and in a convective event
#' in July 2017.
#'
#' The catchment exhibits high land-channel connectivity due to heavily-modified
#' natural channels and extensive artificial drainage, both surface and
#' subsurface. It has a homogenous land cover, with almost all its area comprising
#' class 1 arable grassland and crops. The terrain is undulating with slightly
#' acid, base-rich loam and clay soils predominating. Distances from the
#' channel appear to exert most influence over the catchment response. It is
#' hypothesised that this is due to the influence of the field drainage.

#'@docType data
#'
#'@usage data(brompton)
#' @references
#' Metcalfe P. (2016). Case study 2. Brompton runoff attenuation modelling. In Hankin, B., Burgess-Gamble, L., Bentley, S., Rose, S. (Eds.). How to model and map catchment processes when flood risk management planning. Science report SC120015/R1, Environment Agency, Bristol, UK.
#'
#' Metcalfe, P., Beven, K., Hankin, B., & Lamb, R. (2017). A modelling framework for evaluation of the hydrological impacts of nature-based approaches to flood risk management, with application to in-channel interventions across a 29 km^2 scale catchment in the United Kingdom. Hydrological Processes, 31(9), 1734-1748.
#'@examples
#' require(dynatop)
#' data(brompton)
#'
#'
#' # Show it
#' # raster::image(brompton$dem)

#'@format An environment comprising the DEM, river network, observed flows from stage data reconstructed at the EA gauge at Water End, and hydrometric (rainfall and pe)
#' necessary to run the model.
#'
#' @seealso \code{\link{run_dtm}}
"brompton"
