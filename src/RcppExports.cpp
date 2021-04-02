// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dt_exp_init
void dt_exp_init(std::vector<int> id, Rcpp::NumericMatrix states, Rcpp::NumericMatrix attr, Rcpp::NumericMatrix param, std::vector<int> channel_id, std::vector<int> flow_from, std::vector<int> flow_to, std::vector<double> flow_frc, std::vector<double> r_uz_sz_0);
RcppExport SEXP _dynatop_dt_exp_init(SEXP idSEXP, SEXP statesSEXP, SEXP attrSEXP, SEXP paramSEXP, SEXP channel_idSEXP, SEXP flow_fromSEXP, SEXP flow_toSEXP, SEXP flow_frcSEXP, SEXP r_uz_sz_0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type id(idSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type states(statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type attr(attrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type param(paramSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type channel_id(channel_idSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type flow_from(flow_fromSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type flow_to(flow_toSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type flow_frc(flow_frcSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type r_uz_sz_0(r_uz_sz_0SEXP);
    dt_exp_init(id, states, attr, param, channel_id, flow_from, flow_to, flow_frc, r_uz_sz_0);
    return R_NilValue;
END_RCPP
}
// dt_exp_implicit
void dt_exp_implicit(std::vector<int> id, Rcpp::NumericMatrix states, Rcpp::NumericMatrix attr, Rcpp::NumericMatrix param, std::vector<int> channel_id, Rcpp::NumericMatrix channel_attr, std::vector<int> flow_from, std::vector<int> flow_to, std::vector<double> flow_frc, std::vector<int> precip_col, std::vector<int> precip_id, std::vector<double> precip_frc, std::vector<int> pet_col, std::vector<int> pet_id, std::vector<double> pet_frc, Rcpp::NumericMatrix obs, Rcpp::NumericMatrix channel_inflow, Rcpp::NumericMatrix mass_balance, std::vector<bool> keep_states, Rcpp::List state_rec, double timestep, int n_sub_step);
RcppExport SEXP _dynatop_dt_exp_implicit(SEXP idSEXP, SEXP statesSEXP, SEXP attrSEXP, SEXP paramSEXP, SEXP channel_idSEXP, SEXP channel_attrSEXP, SEXP flow_fromSEXP, SEXP flow_toSEXP, SEXP flow_frcSEXP, SEXP precip_colSEXP, SEXP precip_idSEXP, SEXP precip_frcSEXP, SEXP pet_colSEXP, SEXP pet_idSEXP, SEXP pet_frcSEXP, SEXP obsSEXP, SEXP channel_inflowSEXP, SEXP mass_balanceSEXP, SEXP keep_statesSEXP, SEXP state_recSEXP, SEXP timestepSEXP, SEXP n_sub_stepSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type id(idSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type states(statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type attr(attrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type param(paramSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type channel_id(channel_idSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type channel_attr(channel_attrSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type flow_from(flow_fromSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type flow_to(flow_toSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type flow_frc(flow_frcSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type precip_col(precip_colSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type precip_id(precip_idSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type precip_frc(precip_frcSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type pet_col(pet_colSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type pet_id(pet_idSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pet_frc(pet_frcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type channel_inflow(channel_inflowSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mass_balance(mass_balanceSEXP);
    Rcpp::traits::input_parameter< std::vector<bool> >::type keep_states(keep_statesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type state_rec(state_recSEXP);
    Rcpp::traits::input_parameter< double >::type timestep(timestepSEXP);
    Rcpp::traits::input_parameter< int >::type n_sub_step(n_sub_stepSEXP);
    dt_exp_implicit(id, states, attr, param, channel_id, channel_attr, flow_from, flow_to, flow_frc, precip_col, precip_id, precip_frc, pet_col, pet_id, pet_frc, obs, channel_inflow, mass_balance, keep_states, state_rec, timestep, n_sub_step);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dynatop_dt_exp_init", (DL_FUNC) &_dynatop_dt_exp_init, 9},
    {"_dynatop_dt_exp_implicit", (DL_FUNC) &_dynatop_dt_exp_implicit, 22},
    {NULL, NULL, 0}
};

RcppExport void R_init_dynatop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
