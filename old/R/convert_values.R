# if the max value of the given series or value is greater than 1m it implies that the values are in mm/hr.
# ensure consistency by dividing by the given factor. Would have been better do supply inputs in mm/hr. Too late now!
# also replace any NAs with zeroes
convert_values <- function(x,
                           msg="Value(s) converted to m/hr", max=1, fact=1000)
{
  # convert any NAs to zeroes
  ina <- which(is.na(x[]))
  if(length(ina)>0)
  {
    x[ina] <- 0
  }
  if(length(x) > 0 && max(x, na.rm = TRUE) > max)
  {
    #message(msg)
    x <- x/fact
  }
  return(x)
}
