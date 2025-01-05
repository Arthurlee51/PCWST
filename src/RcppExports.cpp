// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_grad_val_mat
List rcpp_grad_val_mat(NumericVector M_vec, NumericVector n_vec, NumericVector w_vec);
RcppExport SEXP _PCWST_rcpp_grad_val_mat(SEXP M_vecSEXP, SEXP n_vecSEXP, SEXP w_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type M_vec(M_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_vec(n_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_vec(w_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_grad_val_mat(M_vec, n_vec, w_vec));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_val
double rcpp_val(NumericVector M_vec, NumericVector n_vec, NumericVector w_vec);
RcppExport SEXP _PCWST_rcpp_val(SEXP M_vecSEXP, SEXP n_vecSEXP, SEXP w_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type M_vec(M_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_vec(n_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_vec(w_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_val(M_vec, n_vec, w_vec));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_like
double rcpp_like(NumericMatrix W, NumericMatrix Mhat, int n);
RcppExport SEXP _PCWST_rcpp_like(SEXP WSEXP, SEXP MhatSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Mhat(MhatSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_like(W, Mhat, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PCWST_rcpp_grad_val_mat", (DL_FUNC) &_PCWST_rcpp_grad_val_mat, 3},
    {"_PCWST_rcpp_val", (DL_FUNC) &_PCWST_rcpp_val, 3},
    {"_PCWST_rcpp_like", (DL_FUNC) &_PCWST_rcpp_like, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_PCWST(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
