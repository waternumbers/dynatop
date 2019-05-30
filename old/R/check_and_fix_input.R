#' checking the input for potential errors, fix them and report if specified
check_and_fix_input <- function(envir, show.warnings=FALSE, show.errors=FALSE)
{

  errs <- NULL
  warns <- NULL
  if(is.null(envir$qt0) | !is.finite(envir$qt0))
  {
    warns <- c(warns, "Initial specific discharge qt0 required, setting to 1e-5")
    envir$qt0 <- 1e-5
  }
  if(is.null(envir$qmax) | !is.finite(envir$qt0))
  {
    warns <- c(warns, "Max displayed specific discharge required: setting to 5/1000")
    envir$qmax <- 5/100
  }
  if(!(is.data.frame(envir$groups)|is.matrix(envir$groups))){
    errs <- c(errs, "Areal grouping information required as data frame or matrix")
  }
  ngroups <- nrow(envir$groups)

  if(!is.matrix(envir$weights)){
    errs <- c(errs, "Flow weightings - matrix required")
  }
  if(nrow(envir$weights)!= ncol(envir$weights))
  {
    errs <- c(errs, "Square input weighting matrix required")
  }
  if(nrow(envir$weights)!=ngroups | ncol(envir$weights) != nrow(envir$weights))
  {
    errs <- c(errs, "Flow weighting matrix incorrect size - require square matrix of side equal to # areal groupings")
  }
  # ensure we have the expected col names and minimal data requirements
  check_cols(envir$groups, c("id", "area", "atb.bar"))

  if(length(envir$rain)==0)
  {
    errs <- c(errs, "Rainfall input required. Check start / end dates")
  }
  else
  {
    # check that input rain and pe (output)
    if(nrow(envir$rain) != nrow(envir$pe))
    {
      min.len <- min(nrow(envir$rain), nrow(envir$pe))
      warns <- c(warns, "rainfall and pe series differ in length, trimming to shortest")
      envir$rain <- envir$rain[1:min.len,]
      envir$pe <- envir$pe[1:min.len,]
    }
  }
  # each row of weighting matrix is total flow out of the group so should add to 1
  dist.sum <- round(rowSums(envir$weights), 1)
  not.unity <- which(dist.sum !=1)
  if(length(not.unity>0))
  {
    envir$weights <-normalise_rows(envir$weights)
  #  warns <- c(warns, paste("Group ", groups[not.unity]$tag, " weights should add to 1", sep=""))
    # normalise
   # envir$weights[not.unity,] <- normalise.vector(weights[not.unity,])
  }


  return(envir)
}
