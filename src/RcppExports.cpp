// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// inner_gibbs_conditional
IntegerMatrix inner_gibbs_conditional(IntegerMatrix zinit, arma::fcube& cond_weights, IntegerMatrix R, arma::fcube& theta, int ncycles);
RcppExport SEXP _mrf2dbayes_inner_gibbs_conditional(SEXP zinitSEXP, SEXP cond_weightsSEXP, SEXP RSEXP, SEXP thetaSEXP, SEXP ncyclesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type zinit(zinitSEXP);
    Rcpp::traits::input_parameter< arma::fcube& >::type cond_weights(cond_weightsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::fcube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type ncycles(ncyclesSEXP);
    rcpp_result_gen = Rcpp::wrap(inner_gibbs_conditional(zinit, cond_weights, R, theta, ncycles));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mrf2dbayes_inner_gibbs_conditional", (DL_FUNC) &_mrf2dbayes_inner_gibbs_conditional, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_mrf2dbayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}