#' Time of maximum observation
#'
#' @description Determine the time of the maximum item in the supplied time series.
#' @author Peter Metcalfe
#' @export time_at_peak
#' @param ts Time series
#' @param icol Column index if a multi-column time series
#' @author Peter Metcalfe
#' @examples
#'
#' require(dynatop)
#'
#' data(brompton)
#'
#' with(brompton$storm.run, time_at_peak(qsim))
#'
time_at_peak <- function(ts, icol=1)
{
  tms <- index(ts)
  imax <- which.max(ts)
  return(tms[imax])

  # ttp <- apply(ts, MARGIN=2,
  #   function(x)
  #   {
  #     imax <- which.max(x)
  #     return(tms[imax])
  #   })
  #
  # return(unlist(ttp))
}


#' Nash Sutcliffe Efficiency of a model's output against observations
#'
#' @description Returns the the NSE (NSE, Nash and Sutcliffe, 1970) of the simulated values against the given observations.
#' @param qsim Time series or vector of simulated values
#' @param qobs Time series or vector of observations
#' @param digits No. decimal places in returned value
#' @return A number <= 1 indicating the goodness of fit of the simulated series against observations (1= perfect fit). Values of >0.8 are generally regarded as "behavioural"
#' @author Peter Metcalfe
#' @import zoo
#' @export NSE
#' @references Nash, J., & Sutcliffe, J. V. (1970). River flow forecasting through conceptual models part I-A discussion of principles. Journal of hydrology, 10(3), 282-290.
#' @examples
#'\dontrun{
#' require(dynatop)
#'
#' data(brompton)
#'
#' # Goodness of fit for the storm simulation
#'
#' NSE(brompton$storm.run$qsim, brompton$storm.run$qobs)
#' }
NSE <- function(qsim, qobs, digits=2)
{
  #qobs <- qobs[index(qsim)]

  qobs <- as.vector(qobs)
	qsim <- as.vector(qsim)

	# shrink to the smallest array
	len <- min(length(qobs), length(qsim))
	qobs <- qobs[1:len]
	qsim <- qsim[1:len]

	igood <- which(!is.na(qobs))
	# test only against non-null observations
	qsim <- qsim[igood]
	qobs <- qobs[!is.na(qobs)]

	if (length(qobs) == 0 || length(qsim) == 0)
	{
	  stop("No non-null observations found")
	}

	res <- 1 - (sum((qobs - qsim)^2)/sum((qobs - mean(qobs))^2))

	res <- round(res, digits)
	return(res)
}

# write a log entry console and/or file
LogEvent <- function(msg, tm=NULL, # time can be actual run time or simulation time
                     warn=FALSE)  #, logout=stderr())  #System.getenv("log"))  # destination
{

     # browser()
		if(is.null(tm)){tm=Sys.time()}
    msg <- paste(as.character(tm), ": ", msg)
    if(warn)
    {
      warning(msg, immediate.=TRUE)
    }
    else
    {
      message(msg)
    }

}


