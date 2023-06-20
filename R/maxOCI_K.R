#' Find the best K clustering (partitioning) result
#'
#' Fix the number of clusters K and find the best clustering result through the weighted K-Means algorithm.
#'
#' @param K A integer to specify the number of clusters.
#' @param dataSets A list of \eqn{n (n \geq K)} vectors. Each vector stores a dataset that represents a distribution.
#' @param data_type A bool value to indicate the type of distributions, 1 : continuous, 0 : discrete. No mixed type is allowed.
#' @param a A real number which is the power parameter of generalized overlapping clustering index. It is used to control the number of clusters. The bigger the \code{a}, the more clusters we get.
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
#' @return A list to store the best (i.e. the maximum OCI) clustering results among the \code{n_sim} K-Means clustering results.
#' The list includes the maximum OIC value and the best clustering result which includes the clusters, ams and centroids (check ?KMeansCluster).
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
#' res = maxOCI_K(K=3, dataSets=ds, data_type=0)
#' res$maxOIC_K
#' res$clusterResult$clusters
#' plot(data=ds,distribution_type="pmf",clusterResult=res$clusterResult)
#'
#'# Continuous distributions
#' set.seed(123)
#' data1 = rt(100, df=2, ncp=2)
#' data2 = rnorm(100,3,5)
#' data3 = rnorm(100, 5, 3)
#' data4 = rnorm(89, 4, 2)
#' data5 = rgamma(90,4,1)
#' data6 = rgamma(95,9,2)
#' data7 = rt(100, df=3, ncp=6)
#' ds = list(data1,data2,data3,data4,data5,data6,data7)
#' class(ds)<-"pmdfs"
#' plot(data=ds,distribution_type="pdf")
#' res = maxOCI_K(K=3, dataSets=ds, data_type=1)
#' res$maxOIC_K
#' res$clusterResult$clusters
#' plot(data=ds,distribution_type="pdf",clusterResult=res$clusterResult)
#'
#' @export
maxOCI_K <- function(K, dataSets, data_type, a=1, weighted=F, b=5, n_sim=10,
                     nthreads=1, kmIters=50, seed=123, show_progress_bar=1, grid=1000){
  resKm <- .KMeanCluster(K=K, dataSets=dataSets, weighted=weighted, b=b, n_sim=n_sim, nthreads=nthreads,
                         nIters=kmIters, showmessage=show_progress_bar,is_pdf=data_type, seed=seed, grid=grid)
  n_pmdfs = length(unlist(resKm[[1]][[1]]))
  bestId=-1
  OIC_K = 0
  for (i in 1:length(resKm)) {
    if(weighted){
      p_m = unlist(lapply(resKm[[i]][[1]], length))/n_pmdfs
    }else{
      p_m = rep(1/length(resKm[[i]][[1]]),length(resKm[[i]][[1]]))
    }
    A_m = unlist(resKm[[i]][[2]])
    currentOIC_K = .OCI_K(p_m,A_m,a)
    if(OIC_K<currentOIC_K){
      OIC_K = currentOIC_K
      bestId = i
    }
  }
  clRes = resKm[[bestId]]
  for (i in 1:length(clRes)) {
    names(clRes)<-c('clusters','ams','centroids')
  }
  # ordered the result for plot
  orderedClRes = clRes
  rank(unlist(lapply(clRes$clusters, min)))->ords
  for (i in 1:K) {
    which(ords==i)->id
    orderedClRes$clusters[[i]]<-clRes$clusters[[id]]
    orderedClRes$ams[[i]]<-clRes$ams[[id]]
    orderedClRes$centroids[[i]]<-clRes$centroids[[id]]
  }
  return(list(maxOIC_K=OIC_K, clusterResult = orderedClRes))
}
