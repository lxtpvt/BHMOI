# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.kde <- function(data, positions) {
    .Call(`_BHMOI_docbi_kde`, data, positions)
}

.clusterKDE <- function(data, positions) {
    .Call(`_BHMOI_docbi_clusterKDE`, data, positions)
}

.ovl_continuous <- function(data_1, data_2, grid) {
    .Call(`_BHMOI_docbi_ovl_continuous`, data_1, data_2, grid)
}

.ovl_discrete <- function(data_1, data_2) {
    .Call(`_BHMOI_docbi_ovl_discrete`, data_1, data_2)
}

.KMeanCluster <- function(K, dataSets, weighted, b, n_sim, nthreads, nIters, showmessage, is_pdf, seed, grid = 1000L) {
    .Call(`_BHMOI_docbi_kmeansCluster`, K, dataSets, weighted, b, n_sim, nthreads, nIters, showmessage, is_pdf, seed, grid)
}

