#' Find the best K clustering
#'
#' Fix the number of clusters K and find the best clustering result through maximizing OCI value.
#'
#' @param K A integer to specify the number of clusters.
#' @param dataSets A list of \eqn{n (n \geq K)} vectors. Each vector stores a dataset that represents a distribution.
#' @param data_type A bool value to indicate the type of distributions, 1 : continuous, 0 : discrete. No mixed type is allowed.
#' @param weighted A bool value to indicate the K-Means algorithm is weighted or not.
#' @param b A real number which is the power parameter of generalized weighted K-Means algorithm.
#' @param n_sim A integer to set the repeat times of weighted K-Means algorithm. K-Means algorithm can stuck at a local optimal. Increase n_sim can increase the probability to find the global optimal.
#' @param nthreads A integer to set the number of threads for computation.
#' @param kmIters A integer to set the max number of iterations in K-Means algorithm.
#' @param seed A integer to set the seed for K-Means algorithm.
#' @param show_progress_bar A bool value to show the progress bar.
#' @param grid A integer the number of discrete points to calculate the probability density.
#'
#' @details There are a set datasets in \code{dataSets}. Each dataset includes discrete data points of a distribution. This function will search the best \code{K} clustering result of the distributions.
#'
#' @return A list to store the results of running the weighted K-Means algorithm. The length of list is \code{n_sim}
#' which is corresponding to the number of running weighted K-Means algorithm. Running multiple times to avoid the
#' initial issue in K-Means algorithm. Each element of the list is a sub-list which includes "clusters","ams"
#' and "centroids" for each simulation. "ams" is a value for simplifying OCI and OBI calculation.
#'
#' @references The weighted K-Means clustering algorithm in the paper
#' \emph{Distribution-free Overlapping Indices for Dynamic Information Borrowing in Bayesian Hierarchical Modeling}.
#'
#' @examples
#'
#'# Discrete distributions
#' set.seed(123)
#' data1 = rbinom(n=50,size=20,prob=0.3)
#' data2 = rbinom(n=30,size=10,prob=0.5)
#' data3 = rbinom(n=20,size=15,prob=0.2)
#' data4 = rbinom(n=28,size=25,prob=0.4)
#' data5 = rpois(35,4)
#' data6 = rpois(45,9)
#' data7 = rpois(15,8)
#' ds = list(data1,data2,data3,data4,data5,data6,data7)
#' res = KMeansCluster(2,ds,0)
#'
#' @export
KMeansCluster <- function(K, dataSets, data_type, weighted=F, b=5, n_sim=5,
                          nthreads=1, kmIters=50, seed=123, show_progress_bar=1, grid=1000){
  res = .KMeanCluster(K=K, dataSets=dataSets, weighted=weighted, b=b, n_sim=n_sim, nthreads=nthreads,
                      nIters=kmIters, showmessage=show_progress_bar,is_pdf=data_type, seed=seed, grid=grid)
  for (i in 1:length(res)) {
    class(res[[i]])<-"clusterResult"
    names(res[[i]])<-c('clusters','ams','centroids')
  }
  return(res)
}
