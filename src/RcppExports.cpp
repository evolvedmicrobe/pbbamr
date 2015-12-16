// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// loadpbi
DataFrame loadpbi(std::string filename);
RcppExport SEXP pbbamr_loadpbi(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    __result = Rcpp::wrap(loadpbi(filename));
    return __result;
END_RCPP
}
// timesTwo
CharacterVector timesTwo(NumericVector x);
RcppExport SEXP pbbamr_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = Rcpp::wrap(timesTwo(x));
    return __result;
END_RCPP
}
