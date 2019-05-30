# check that tabular object has given column names
check_cols <- function(obj, names, stop=TRUE)
{
  msg <- NULL
  for(nm in names)
  {
    if(!(nm %in% colnames(obj)))
    {
      msg <- paste(msg, "Missing column name ", nm, "\n")
    }
  }
  if(length(msg)>0 & stop){stop(msg)}
  return(msg)
}
