# To do list for v0.3

- check initialisation part of HRU in digest_hru
- stop call to initisation failing if add_data has not been called
- handle get_output nicely when there is no output
- check for id in output_defn isn't correct (fails if set to 1 when the only
  id is zero)
- fix digest_obs so that NA values in ob series not used anre allowed
- document flow split on exit since this is confusing (q_sf when s_sz is zero
  at tstart and no inflow to it...)
- document that Dt >= kappa x eta. In kinemtatic case c*Dt >= Dx/2 to produce positive outflows - knock on impact on
  values of maximum downward flux
- look at whether q_sz can go negative
