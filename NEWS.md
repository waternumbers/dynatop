# dynatop 0.0.0.9000
- Gone through Peter code to find functions actually beign used
- Altered input files to allow for more easy reparameterisation and use of multiple inputs (used named columns, parameters not implicit by columns order)
- rewritten initialisation and main execution loop of dynamic TOPMODEL so that
    - Chanels are handled explicily not as 'special cases' (parameterisations) of hillslope units
    - Outputs can be selected
    - changed analytical solution of root zone to handle case where both pe
    and precip are positive

