# dynatopGIS 0.1

## Breaking Changes

- Code base reformulated in an Object orientated form using the `R6`
  package. Except for change below algorithms as for v0.0.4
- Input format of model changed to reflect `dynatopGIS` v0.1
- Removal of banding functions since these are now in dynatopGIS

## New Features

- Reformulation allows for data and code to saved in a single object allowing fuller
  reproducibility
- Additional plotting and data extraction functions

# dynatop 0.0.4

## Breaking changes
- New solutions to the surface and saturated zone means models for previous
  versions of dynatop will not work. See the model description and the
  `dynatopGIS` package (0.0.4) for a method of generating models int he
  revised format

## Other changes
- The surface store, previously labeled $s_{ex}$ has been relabelled $s_{sf}$
  since this is felt to be more logical.
- The `Matrix` package is allowing the use of sparse matrices.
- Computation bands have been introduced for the surface and saturated
  zones. These allow the revised solutions below. See the associated vignette
  on banding HSUs
- A revised approximation to the surface water movement has been implimented
  allowing for larger numbers of HSUs without the performance overhead of
  computing large matrix exponentials.
- A four point kinematic wave solution to the saturated zone is utilised. This
  is both more performant and has better representation of hillslope length
  then the numerical ODE solution used in v0.0.3.

# dynatop 0.0.3

## New features
- Added time delay histogram river routing. See vignettes on theory and use.

## Bug fixes
- Fixed some mass balance issues in dynatop relating to the saturated zone
  solution in dynatop
- Improved handling of case of infinite saturated zone deficit

## Breaking changes
- The model structure has been adapted to allow for the channel connectivity
  and specification of gauges and inputs on the channel network. See the model
  object vignette. If using dynatopGIS then rerun create_model.
  
## Other changes
- Fixed minor formating issues in the vignettes
- Added vignette on checking mass balances

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
