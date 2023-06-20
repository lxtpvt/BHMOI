#' Find the best K clustering (partitioning) result from the result of weighted K-Means algorithm
#'
#' Fix the number of clusters K and find the best clustering result from the result of weighted K-Means algorithm.
#'
#' @param kmeansRes A integer to specify the number of clusters.
#' @param a A real number which is the power parameter of generalized overlapping clustering index. It is used to control the number of clusters. The bigger the \code{a}, the more clusters we get.
#' @param weighted A bool value to indicate the K-Means algorithm is weighted or not.

#'
#' @details Fix the number of clusters K and find the best clustering result from the result of weighted K-Means algorithm.
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
#' resKM = KMeansCluster(2,ds,0)
#' res = maxResOCI_K(kmeansRes=resKM, weighted=F, a=0.2)
#' res$maxOIC_K
#' res$clusterResult$clusters
#' plot(data=ds,distribution_type="pmf",clusterResult=res$clusterResult)
#'
#' @export
maxResOCI_K <- function(kmeansRes,weighted,a){
  n_pmdfs = length(unlist(kmeansRes[[1]][[1]]))
  bestId=-1
  OIC_K = 0
  for (i in 1:length(kmeansRes)) {
    if(weighted){
      p_m = unlist(lapply(kmeansRes[[i]][[1]], length))/n_pmdfs
    }else{
      p_m = rep(1/length(kmeansRes[[i]][[1]]),length(kmeansRes[[i]][[1]]))
    }
    A_m = unlist(kmeansRes[[i]][[2]])
    if(OIC_K<.OCI_K(p_m,A_m,a)){
      OIC_K = .OCI_K(p_m,A_m,a)
      bestId = i
    }
  }
  clRes = kmeansRes[[bestId]]
  for (i in 1:length(clRes)) {
    names(clRes)<-c('clusters','ams','centroids')
  }
  class(clRes)<-"clusterResult"
  return(list(maxOIC_K=OIC_K, clusterResult = clRes))
}
