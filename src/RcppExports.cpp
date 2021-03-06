// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// css_statistic
NumericVector css_statistic(NumericVector& y);
RcppExport SEXP _variance_css_statistic(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(css_statistic(y));
    return rcpp_result_gen;
END_RCPP
}
// BS_rcpp
Rcpp::List BS_rcpp(int s, int e, NumericVector& y, float penality);
RcppExport SEXP _variance_BS_rcpp(SEXP sSEXP, SEXP eSEXP, SEXP ySEXP, SEXP penalitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< float >::type penality(penalitySEXP);
    rcpp_result_gen = Rcpp::wrap(BS_rcpp(s, e, y, penality));
    return rcpp_result_gen;
END_RCPP
}
// penalty_fun
float penalty_fun(int n, int params);
RcppExport SEXP _variance_penalty_fun(SEXP nSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(penalty_fun(n, params));
    return rcpp_result_gen;
END_RCPP
}
// variance
float variance(std::vector<float>& v, float& mu);
RcppExport SEXP _variance_variance(SEXP vSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<float>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< float& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(variance(v, mu));
    return rcpp_result_gen;
END_RCPP
}
// cost_function
float cost_function(std::vector<float>& v);
RcppExport SEXP _variance_cost_function(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<float>& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(cost_function(v));
    return rcpp_result_gen;
END_RCPP
}
// slicing
std::vector<float> slicing(std::vector<float> const& v, int X, int Y);
RcppExport SEXP _variance_slicing(SEXP vSEXP, SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<float> const& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(slicing(v, X, Y));
    return rcpp_result_gen;
END_RCPP
}
// argmin
int argmin(std::vector<float>& v);
RcppExport SEXP _variance_argmin(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<float>& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(argmin(v));
    return rcpp_result_gen;
END_RCPP
}
// OP_cpp
std::vector<int> OP_cpp(std::vector<float>& data, int params);
RcppExport SEXP _variance_OP_cpp(SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<float>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(OP_cpp(data, params));
    return rcpp_result_gen;
END_RCPP
}
// PELT_cpp
std::vector<int> PELT_cpp(std::vector<float>& data, int params, int K);
RcppExport SEXP _variance_PELT_cpp(SEXP dataSEXP, SEXP paramsSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<float>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(PELT_cpp(data, params, K));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_variance_css_statistic", (DL_FUNC) &_variance_css_statistic, 1},
    {"_variance_BS_rcpp", (DL_FUNC) &_variance_BS_rcpp, 4},
    {"_variance_penalty_fun", (DL_FUNC) &_variance_penalty_fun, 2},
    {"_variance_variance", (DL_FUNC) &_variance_variance, 2},
    {"_variance_cost_function", (DL_FUNC) &_variance_cost_function, 1},
    {"_variance_slicing", (DL_FUNC) &_variance_slicing, 3},
    {"_variance_argmin", (DL_FUNC) &_variance_argmin, 1},
    {"_variance_OP_cpp", (DL_FUNC) &_variance_OP_cpp, 2},
    {"_variance_PELT_cpp", (DL_FUNC) &_variance_PELT_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_variance(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
