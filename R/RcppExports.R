# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

dt_init <- function(mdl, vtol, etol, max_it) {
    invisible(.Call(`_dynatop_dt_init`, mdl, vtol, etol, max_it))
}

dt_sim <- function(mdl, out_dfn, keep_states, obs_matrix, mass_balance, out_matrix, state_rec, timestep, n_sub_step, vtol, etol, max_it) {
    invisible(.Call(`_dynatop_dt_sim`, mdl, out_dfn, keep_states, obs_matrix, mass_balance, out_matrix, state_rec, timestep, n_sub_step, vtol, etol, max_it))
}

