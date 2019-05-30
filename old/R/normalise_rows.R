# ensure sum of all rows in matrix is unity
normalise_rows <- function(w)
{
  if(!is.matrix(w)){warning("matrix input expected to normalise.rows"); return(w)}
  if(any(rowSums(w)!=1))
  {
    cat("Normalising rows...\n")

    #     w <-apply(w, MARGIN=1,
    #               FUN=function(x)
    #               {
    #                 sum <- sum(x, na.rm=TRUE)
    #                 res<- ifelse(sum==0, rep(0, length(x)), x/sum)
    #                 return(res)
    #               }
    #     )


    w <- apply(w, MARGIN=1,
               FUN=function(x)
               {
                 sum <- sum(x, na.rm=TRUE)
                 res<- x/sum
                 return(res)
               })
    w[!is.finite(w)]  <- 0
    #  w[1:nrow(w),]<-w[1:nrow(w),]/rowSums(w)
    # transpose to get in rightg row. col order
    w <- t(w)
  }
  # need to transpose to get in right row-col order
  return(w)
}
