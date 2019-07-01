# dynatop 0.0.2

## Bug fixes
- Added purl=FALSE where error=TRUE in vignette

## Other changes
- Removed unused Rcpp functions and library dependencies

# dynatop 0.0.1

## Other changes
- Removed vignette data from package to reduce size. Now on git hub
- Minor spelling and format changes to vignette

# dynatop 0.0.0.9000

## Context
This package is the result of an almost complete rewrite of the dynatopmodel package
formally CRAN and the associated development code (not in the public
domain).

This package contains the code for model evaluation and helpers for preparing time
series input data. The package [dynatopGIS](https://waternumbers.github.io/dynatopGIS/) contains the
tools for preparing models from GIS data.

## New Features
- New definition of a dynamic TOPMODEL 'object'
      - seperates out GIS data for more compact size
	  - introduces the concept of Hydrological Response Unit (HRU) types
	  - Altered HRU definition table to allow for more transparent and easy
        reparameterisation
	  - Altered HRU definition table to allow for more transparent use of
        multiple input series (referenced by name, not input column)
- Altered input of time series data to a single xts object containing names
columns
- New output list containing:
      - model object, with final states of the system
	  - xts object of flows to channel HRUs
- Rewritten initialisation and main execution loop of dynamic TOPMODEL so that
      - Chanels are handled explicily not as 'special cases' (parameterisations) of hillslope units
      - Analytical solution of root zone can handle case where both pet and
        precip are positive
	  - Initialisation does not return a new workspace
- Vignettes documenting
      - use of dynamic TOPMODEL
	  - Equations and solution methodology
	  - Performance of code for larger simulations
- 'Hard' checks are now implimented, previously code was performaing various
  'corrections' sometimes without warning.
- Model simulation can be initialised by passed in states allowing for
  chinking of the longer simulations.

## Regressions of note
- Plotting functions and performance calculations have not been replicated
- The time delay historgram method for river routing is not implimented
- Model no longer outputs states and fluxes at every time step
