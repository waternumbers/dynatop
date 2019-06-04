## script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- './dynatop_helpers'
Rcpp::compileAttributes(pacPath)
devtools::document(pacPath)
devtools::check(pacPath)
## devtools::build(pacPath)


