// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dt_init
void dt_init(Rcpp::List mdl, double const vtol, double const etol, int const max_it);
RcppExport SEXP _dynatop_dt_init(SEXP mdlSEXP, SEXP vtolSEXP, SEXP etolSEXP, SEXP max_itSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type mdl(mdlSEXP);
    Rcpp::traits::input_parameter< double const >::type vtol(vtolSEXP);
    Rcpp::traits::input_parameter< double const >::type etol(etolSEXP);
    Rcpp::traits::input_parameter< int const >::type max_it(max_itSEXP);
    dt_init(mdl, vtol, etol, max_it);
    return R_NilValue;
END_RCPP
}
// dt_sim
void dt_sim(Rcpp::List mdl, Rcpp::DataFrame out_dfn, std::vector<bool> keep_states, Rcpp::NumericMatrix obs_matrix, Rcpp::NumericMatrix mass_balance, Rcpp::NumericMatrix out_matrix, Rcpp::List state_rec, double const timestep, int const n_sub_step, double const vtol, double const etol, int const max_it);
RcppExport SEXP _dynatop_dt_sim(SEXP mdlSEXP, SEXP out_dfnSEXP, SEXP keep_statesSEXP, SEXP obs_matrixSEXP, SEXP mass_balanceSEXP, SEXP out_matrixSEXP, SEXP state_recSEXP, SEXP timestepSEXP, SEXP n_sub_stepSEXP, SEXP vtolSEXP, SEXP etolSEXP, SEXP max_itSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type mdl(mdlSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type out_dfn(out_dfnSEXP);
    Rcpp::traits::input_parameter< std::vector<bool> >::type keep_states(keep_statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type obs_matrix(obs_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mass_balance(mass_balanceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type out_matrix(out_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type state_rec(state_recSEXP);
    Rcpp::traits::input_parameter< double const >::type timestep(timestepSEXP);
    Rcpp::traits::input_parameter< int const >::type n_sub_step(n_sub_stepSEXP);
    Rcpp::traits::input_parameter< double const >::type vtol(vtolSEXP);
    Rcpp::traits::input_parameter< double const >::type etol(etolSEXP);
    Rcpp::traits::input_parameter< int const >::type max_it(max_itSEXP);
    dt_sim(mdl, out_dfn, keep_states, obs_matrix, mass_balance, out_matrix, state_rec, timestep, n_sub_step, vtol, etol, max_it);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dynatop_dt_init", (DL_FUNC) &_dynatop_dt_init, 4},
    {"_dynatop_dt_sim", (DL_FUNC) &_dynatop_dt_sim, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_dynatop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}