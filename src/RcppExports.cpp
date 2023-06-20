// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// docbi_kde
NumericVector docbi_kde(NumericVector data, NumericVector positions);
RcppExport SEXP _BHMOI_docbi_kde(SEXP dataSEXP, SEXP positionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type positions(positionsSEXP);
    rcpp_result_gen = Rcpp::wrap(docbi_kde(data, positions));
    return rcpp_result_gen;
END_RCPP
}
// docbi_clusterKDE
NumericVector docbi_clusterKDE(List data, NumericVector positions);
RcppExport SEXP _BHMOI_docbi_clusterKDE(SEXP dataSEXP, SEXP positionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type positions(positionsSEXP);
    rcpp_result_gen = Rcpp::wrap(docbi_clusterKDE(data, positions));
    return rcpp_result_gen;
END_RCPP
}
// docbi_ovl_continuous
double docbi_ovl_continuous(NumericVector data_1, NumericVector data_2, size_t grid);
RcppExport SEXP _BHMOI_docbi_ovl_continuous(SEXP data_1SEXP, SEXP data_2SEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type data_1(data_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data_2(data_2SEXP);
    Rcpp::traits::input_parameter< size_t >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(docbi_ovl_continuous(data_1, data_2, grid));
    return rcpp_result_gen;
END_RCPP
}
// docbi_ovl_discrete
double docbi_ovl_discrete(NumericVector data_1, NumericVector data_2);
RcppExport SEXP _BHMOI_docbi_ovl_discrete(SEXP data_1SEXP, SEXP data_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type data_1(data_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data_2(data_2SEXP);
    rcpp_result_gen = Rcpp::wrap(docbi_ovl_discrete(data_1, data_2));
    return rcpp_result_gen;
END_RCPP
}
// docbi_kmeansCluster
List docbi_kmeansCluster(size_t K, List dataSets, bool weighted, double b, size_t n_sim, size_t nthreads, size_t nIters, bool showmessage, bool is_pdf, int seed, int grid);
RcppExport SEXP _BHMOI_docbi_kmeansCluster(SEXP KSEXP, SEXP dataSetsSEXP, SEXP weightedSEXP, SEXP bSEXP, SEXP n_simSEXP, SEXP nthreadsSEXP, SEXP nItersSEXP, SEXP showmessageSEXP, SEXP is_pdfSEXP, SEXP seedSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< size_t >::type K(KSEXP);
    Rcpp::traits::input_parameter< List >::type dataSets(dataSetsSEXP);
    Rcpp::traits::input_parameter< bool >::type weighted(weightedSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< size_t >::type n_sim(n_simSEXP);
    Rcpp::traits::input_parameter< size_t >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< size_t >::type nIters(nItersSEXP);
    Rcpp::traits::input_parameter< bool >::type showmessage(showmessageSEXP);
    Rcpp::traits::input_parameter< bool >::type is_pdf(is_pdfSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(docbi_kmeansCluster(K, dataSets, weighted, b, n_sim, nthreads, nIters, showmessage, is_pdf, seed, grid));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4continuous_response_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4tran_bhm_general_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_BHMOI_docbi_kde", (DL_FUNC) &_BHMOI_docbi_kde, 2},
    {"_BHMOI_docbi_clusterKDE", (DL_FUNC) &_BHMOI_docbi_clusterKDE, 2},
    {"_BHMOI_docbi_ovl_continuous", (DL_FUNC) &_BHMOI_docbi_ovl_continuous, 3},
    {"_BHMOI_docbi_ovl_discrete", (DL_FUNC) &_BHMOI_docbi_ovl_discrete, 2},
    {"_BHMOI_docbi_kmeansCluster", (DL_FUNC) &_BHMOI_docbi_kmeansCluster, 11},
    {"_rcpp_module_boot_stan_fit4continuous_response_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4continuous_response_mod, 0},
    {"_rcpp_module_boot_stan_fit4tran_bhm_general_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4tran_bhm_general_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_BHMOI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
