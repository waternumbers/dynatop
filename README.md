# dynatop

This R package coontains the core code to perform a dynamic TOPMODEL
simulation.

## TO DO
* Provide R and Rcpp version of all functions
* Check parameterisation and definitions in intialisation
* try to fix notes relating to xts and zoo packages in build
* Alter model structure to match GIS
* Write model structure vignette

## Suggestions

## Coding
* Move to object oriented programming systle making use of S4 or R5 classes
      - May get improved memory usage
	  - Allows for a more managed growth of complexity by imposing a more rigid
        style of programming	
* Further checks on eigen vector routing and stability of solution
* Use RcppArmadillo to speed up matrix operations further & can use sparse
  matrices
* Test speed up of moving further components to Rcpp

## Hydrological
* Improve initialisation of model
* Consider improving calculation of quz
* Consider re-introduction of sd_max
