# apply a list of named values to the existing data table or list
# if names missing in target then add to target
# if a list of parameters to ignore provided do not overwrite these in the target
apply_params <- function(groups, params, ihru=1:nrow(groups), preserve=NULL,
                         existing.only=FALSE)
{

  nms <- names(params)

  if(existing.only)
  {
    # only apply the parameters in the target
    nms <- intersect(names(params), names(groups))
  }
  #  params <- params[nms]

  # miss out any given named parameters from the target
  if(length(preserve)>0)
  {
    igood <- which(!(nms %in% preserve))
    nms <- nms[igood]
  }

  # igood <- which(!is.na(names(params)))
  params <- params[nms]


  if(length(params)==0)
  {
   # warning("No valid parameters found to apply")
    return(groups)
  }

  for(nm in names(params))
  {
    if(is.list(groups))
    {
      groups[[nm]] <- params[[nm]]
    }
    else if(is.data.frame(groups))
    {
      groups[ihru,nm] <- params[[nm]]
    }
    #  catch_discretisation[[nm]] <- params[[nm]]
  }
  return(groups)

}
