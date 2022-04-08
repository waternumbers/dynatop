# dynatop 0.2.2

- Fixed missing timestep in some Courant number calculations
- Fixed accounting of channel inflow in mass balance computations and added
  matching test. Previously impacted simulations with a substep.
- Added Runoff Attenuation Feature (RAF) representation to surfce of Hillslope HRU

# dynatop 0.2.1

- Release for CRAN

# dynatop 0.2.0.9200

- Fix so using a ingle channel works correctly
- Fix so fluxes passed correctly during initialisation
- Fix bug in get_states() so returns correctly when recrods=FALSE
- Tidying for release to CRAN

# dynatop 0.2.0.9101

## Breaking changes

- New model structure allowing for the transmissivity profile to be specified per HRU
- New transmissivity profile (double exponential) meaning additional parameter columnes (m_2, omega) in the model

## Other changes

- dropping of Boost bisection algorithm for direct implimentation with user
  specified tolerance and maximum iterations
- C++ code refactored for speed improvements, making use of a single hillslope_hru class
- Refinements to the vignettes and documentation to reflect changes 

# dynatop 0.2.0.904*

- changes to improve plotting


# dynatop 0.2.0.9035

- Adapted R and C++ code to pass data frames of hillslope and channel properties
- tidy up of C++ code
- Fix bug in hillslope HRU solution which was indexing incorrect PET value - this will only impact simulations with multiple PET series, if constant then results shoudl be identical to v0.2.0.9030.
- Revert Channel routing to compute average flux over the time step

## New features

- None

## Breaking Changes

- None

# dynatop 0.2.0.9030

- Adapted hillslope solution to Finite Volume to address mass balance issues in v0.2.0.9020
- Improved mass balance checkign and output

## New features

- Addition of two further transmissivity profiles - these implimentations should be considered **experimental**

## Breaking Changes

- Alterations to model sturcture
    - Two new data.frames added for defining Precipitation and PET inputs. These can now be specified as weighted sums of observed series.
    - The parameter vector `$param` is dropped. Numeric parameters are now stored in the `$hillslope` and `$channel` data.frames. This is to allow for cleaner code and future development of different model stroage less dependent of `*.rds` files.


# dynatop 0.2.0.9020

- Adapted to use a contour (cross section) solution. This brings brings the model into line with the original
  dynamic TOPMODEL concept.
- Complete rewrite of the hillslope simulation code making better use of std C++ classes and
  Boost libraries.

## Breaking Changes

- Model structure altered to contain a further data frame containing all the HRU linkages


# dynatop 0.13

## New features

- C++ implimentation of hillslope simulation code can now return states and is
  fully feature compatable with R implimentation.
- R implimentation of hillslope simulations will be depreciated in a future release.

# dynatop 0.12

## New features

- Main hydrological simulation code implimented in C++ usign Rcpp for better performance. 
- The R version of the code can be used instead by setting the use_R input parameter to TRUE. This will be depreciated in a later version.
- Currently C++ code does not return intermediate states, use the R version if these are required


# dynatop 0.11

## New features

- improved calculations within the saturated zone
- improved mass_check
- more complet vignette - model equation and coding notes

# dynatop 0.1

## Breaking Changes

- Code base reformulated in an Object orientated form using the `R6`
  package. Except for change below algorithms as for v0.0.4
- Input format of model changed to reflect `dynatopGIS` v0.1
- Removal of banding functions since these are now in dynatopGIS
- All units now in m and seconds
- Precipitation and PET inputs now expected to be metres accured over the timestep rather then m per hour

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
