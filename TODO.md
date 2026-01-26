## Things to update in v0.5

- memory management in dynaGIS - goes out with UK sized data, could try
  - move some algorithms to C++ with deque to save sorting
  - stop holding the brk read in rasters as required (not sure what
    impact this would make)
  - make better use of writing rasters on the fly rather then at end
    of steps
  - possible gc() calls???
  
- dynaGIS initialisation slow with large projects

- To revert in model build
  - sort out hru numbering - error between channel and hillslope 
  - sort out band numbering - error between channel and hillslope

- Fix percentage complete calcs in all passes in dygis.R
- Get all verbose in dynaGIS to work: 
  - sink_fill upward pass
  - compute_properties
  - aggregate_layer
  
- get initial state paraemters in create_model
  
- merge in and check ability to handle daily PET values in evap_est.R
  from v0.4.999
- finish merging in all functions in convert_channel.R
- do we need processing_helpers.R from v0.4.999
- documentation & vignettes

