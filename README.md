# dynatop

This R package coontains the core code to perform a dynamic TOPMODEL
simulation.

## TO DO
* Reformat example to reflect changes in debug example
* Provide R and Rcpp version of all functions
* clean out unused Rcpp functions
* Check function documentation and add missing variables
* What happens to cpp version of route_ex_eigen when eigen values of WV are
imaginary?
* check solution to initialisation when only one chanel
* change from channel to chanel
* Does channnel need vof??
* write checks
* add back in water balance?
* add back in output of states?

## Suggestions
* Move to object oriented?
* Check application of eigen vector routing?
* improve initialisation
* Change how initialisation is passed, currently loaded into workspace is
* allow checks to warn and change rather then fail
* Use RcppArmadillo - spped up matrix operation further + can use sparse matrices
* consider consistancy of sd and qsz
* qsz\_max initialised by atb\_bar rather then mean slope

# dynatop v0.0.0.9000

* Gone through Peter code to find functions actually beign used
* Altered input files to allow for more easy reparameterisation and
  use of multiple inputs (used named columns, parameters not implicit
  by columns order)
* rewritten initialisation and main execution loop of dynamic TOPMODEL
  so that
  * Chanels are handled explicily not as 'special cases'
    (parameterisations) of hillslope units
  * Outputs can be selected

