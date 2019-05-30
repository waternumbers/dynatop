# default set of hydrological parameters - applied to HSUs where not supplied
# use underscores rather than periods as ESRI shape format doesn'r support this charcter
# in field names of attached data
def_hsu_par <- function()
{
  list("gauge.id"=1,
  	"srz_max"=0.1,
       "ln_t0"=7,     # saturation transmissivity
       "m"=0.01,
       "srz0"=0, "td"=1,
       "vchan"=1000,
       "vof"=100,
  	"cw"=1.7,        # weir coefficient for overflow discharge from surface
       "k0"=1e8,     # surface conductivity, Used in infiltration excess calcs in Buytaert's implementation of TOPMODEL. Set to large value for no infiltration excess
       "cd"=0.1,     # capillary drive (see Morel-Seytoux and Khanji, 1974, Beven, 1984), as above.
       "atb.bar"=0,  # areal average of topographic index
       "sd_max"=0.5, # max storage deficit
  	"pe_fact"=1,
  	"vof_fact"=1,
  	"rain_fact"=1,   # scaling factor for rainfall input
  	"mann.n"=0.01,    # roughness factor for rough grass
  	"s0"=0.1,       # nominal surface slope
  	"ex_max"=1)      # maximum surface excess storage
}

# default set of parameters applied to channel HSU
def_chan_par <- function(vals=list())
{
  res <- list("srz_max"=0, # rainfall goes straight to the unsat zone infiltration
       "ln_t0"=8,
       "m"=0.01,
       "srz0"=0,
       "td"=1e-8 ,
       "vchan"=1500,
       "vof"=10,
       "sd_max"=2)    # depth of rectangular channel?
  merge_lists(res, vals)
}
