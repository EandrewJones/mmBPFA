// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sample_alpha
arma::vec sample_alpha(const arma::mat& x, const arma::mat& lambda, const arma::mat& omega);
RcppExport SEXP _mmBPFA_sample_alpha(SEXP xSEXP, SEXP lambdaSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_alpha(x, lambda, omega));
    return rcpp_result_gen;
END_RCPP
}
// sample_gamma_k
arma::vec sample_gamma_k(const arma::mat& zeros, const arma::mat& lambda, const double c, const double d);
RcppExport SEXP _mmBPFA_sample_gamma_k(SEXP zerosSEXP, SEXP lambdaSEXP, SEXP cSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type zeros(zerosSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_gamma_k(zeros, lambda, c, d));
    return rcpp_result_gen;
END_RCPP
}
// sample_d
double sample_d(const arma::mat& zeros, const arma::vec& gamma_k, const double c, const double c0, const double d0);
RcppExport SEXP _mmBPFA_sample_d(SEXP zerosSEXP, SEXP gamma_kSEXP, SEXP cSEXP, SEXP c0SEXP, SEXP d0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type zeros(zerosSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gamma_k(gamma_kSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const double >::type d0(d0SEXP);
    rcpp_result_gen = Rcpp::wrap(sample_d(zeros, gamma_k, c, c0, d0));
    return rcpp_result_gen;
END_RCPP
}
// sample_IBP_a
double sample_IBP_a(const double ibp_b, const arma::mat& zeros, const double e, const double f);
RcppExport SEXP _mmBPFA_sample_IBP_a(SEXP ibp_bSEXP, SEXP zerosSEXP, SEXP eSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type ibp_b(ibp_bSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type zeros(zerosSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_IBP_a(ibp_b, zeros, e, f));
    return rcpp_result_gen;
END_RCPP
}
// sample_IBP_b
double sample_IBP_b(const double ibp_a, const double ibp_b, const arma::mat& zeros);
RcppExport SEXP _mmBPFA_sample_IBP_b(SEXP ibp_aSEXP, SEXP ibp_bSEXP, SEXP zerosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type ibp_a(ibp_aSEXP);
    Rcpp::traits::input_parameter< const double >::type ibp_b(ibp_bSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type zeros(zerosSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_IBP_b(ibp_a, ibp_b, zeros));
    return rcpp_result_gen;
END_RCPP
}
// sample_lambda_and_zeros
List sample_lambda_and_zeros(arma::mat x, arma::mat lambda, arma::mat zeros, arma::mat omega, arma::vec alpha, arma::vec gamma_k, std::vector<int> dc, const double ibp_a, const double ibp_b, const double tau, const bool sparse, const bool infinite);
RcppExport SEXP _mmBPFA_sample_lambda_and_zeros(SEXP xSEXP, SEXP lambdaSEXP, SEXP zerosSEXP, SEXP omegaSEXP, SEXP alphaSEXP, SEXP gamma_kSEXP, SEXP dcSEXP, SEXP ibp_aSEXP, SEXP ibp_bSEXP, SEXP tauSEXP, SEXP sparseSEXP, SEXP infiniteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zeros(zerosSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_k(gamma_kSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< const double >::type ibp_a(ibp_aSEXP);
    Rcpp::traits::input_parameter< const double >::type ibp_b(ibp_bSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< const bool >::type infinite(infiniteSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_lambda_and_zeros(x, lambda, zeros, omega, alpha, gamma_k, dc, ibp_a, ibp_b, tau, sparse, infinite));
    return rcpp_result_gen;
END_RCPP
}
// sample_omega
arma::mat sample_omega(arma::mat lambda, arma::mat x, arma::vec alpha);
RcppExport SEXP _mmBPFA_sample_omega(SEXP lambdaSEXP, SEXP xSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_omega(lambda, x, alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mmBPFA_sample_alpha", (DL_FUNC) &_mmBPFA_sample_alpha, 3},
    {"_mmBPFA_sample_gamma_k", (DL_FUNC) &_mmBPFA_sample_gamma_k, 4},
    {"_mmBPFA_sample_d", (DL_FUNC) &_mmBPFA_sample_d, 5},
    {"_mmBPFA_sample_IBP_a", (DL_FUNC) &_mmBPFA_sample_IBP_a, 4},
    {"_mmBPFA_sample_IBP_b", (DL_FUNC) &_mmBPFA_sample_IBP_b, 3},
    {"_mmBPFA_sample_lambda_and_zeros", (DL_FUNC) &_mmBPFA_sample_lambda_and_zeros, 12},
    {"_mmBPFA_sample_omega", (DL_FUNC) &_mmBPFA_sample_omega, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_mmBPFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
