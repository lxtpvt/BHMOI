#' Find the best clustering (partitioning) result of a set of distributions
#'
#' This function uses the weighted K-Means (see the reference) method to find the best clustering (partitioning) result of a set of distributions.
#'
#' @param dataSets A list of \eqn{n (n \geq K)} vectors. Each vector stores a dataset that represents a distribution.
#' @param data_type A bool value to indicate the type of distributions, 1 : continuous, 0 : discrete. No mixed type is allowed.
#' @param a A real number which is the power parameter of generalized overlapping clustering index. It is used to control the number of clusters. The bigger the \code{a}, the more clusters we get.
#' @param weighted A bool value to indicate the K-Means algorithm is weighted or not.
#' @param b A real number which is the power parameter of generalized weighted K-Means algorithm.
#' @param fast A bool value to select the strategy to calculate the next level of K. 1: stop at the first time K does not increase; 0: stop at OCI_K>OCI_{K+1}.
#' @param levels A list of \eqn{n (n \geq K)} vectors. Each vector stores the values of K that need to calculate OCI.
#' @param K A integer to set the number of cluster.
#' @param n_sim A integer to set the repeat times of weighted K-Means algorithm. K-Means algorithm can stuck at a local optimal. Increase n_sim can increase the probability to find the global optimal.
#' @param nthreads A integer to set the number of threads for computation.
#' @param showmessage A bool value to show messages of process running.
#' @param kmIters A integer to set the max number of iterations in K-Means algorithm.
#' @param seed A integer to set the seed for K-Means algorithm.
#' @param grid A integer the number of discrete points to calculate the probability density.
#'
#'
#' @details There are a set datasets in \code{dataSets}. Each dataset includes discrete data points of a distribution. This function will search the best clustering result of the distributions. The number of clusters can be influenced by parameter \code{a}. The bigger the \code{a}, the more clusters we get.
#'
#' @return A list that stores the best (i.e. the global maximum OCI) clustering results.
#' The list includes the maximum K, corresponding OIC value and the best clustering result(check ?maxOCI_K).
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
#' class(ds)<-"pmdfs"
#' plot(data=ds,distribution_type="pmf")
#' res = clusterMaxOCI(dataSets=ds, data_type=0, a=0.4)
#' res$K
#' res$OCI
#' res$clusterResult$clusters
#' plot(data=ds,distribution_type="pmf",clusterResult=res$clusterResult)
#'
#' @export
clusterMaxOCI <- function(dataSets, data_type, a, weighted=0, b=1, fast=0, levels = 0, K=0, n_sim=5, nthreads=1,
                          showmessage=1, kmIters=50, seed=123, grid=1000){
  if(K>0){
    return(.calMaxOCI(dataSets=dataSets, data_type=data_type, a=a, weighted=weighted, b=b, fast=fast,
                      allKs=K, n_sim=n_sim, nthreads=nthreads, kmIters=kmIters, seed=seed,
                      showmessage=showmessage, grid=grid))
  }else{
    if(levels==0){
      return(.calMaxOCI(dataSets=dataSets, data_type=data_type, a=a, weighted=weighted, b=b, fast=fast,
                        allKs=c(1:length(dataSets)), n_sim=n_sim, nthreads=nthreads, kmIters=kmIters, seed=seed,
                        showmessage=showmessage, grid=grid))
    }else{
      return(.calMaxOCI(dataSets=dataSets, data_type=data_type, a=a, weighted=weighted, b=b, fast=fast,
                        allKs=levels, n_sim=n_sim, nthreads=nthreads, kmIters=kmIters, seed=seed,
                        showmessage=showmessage, grid=grid))
    }
  }
}
