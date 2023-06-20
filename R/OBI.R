#' Calculate OBIs with the result of clustering
#'
#' This function calculate OBIs with the result of clustering.
#' 
#'
#' @param clusterResult A list to store the clustering results (see the return in ?KMeansCluster).
#' 
#' @details This function calculate OBIs with the result of clustering.
#'
#' @return A vector values of OBIs.
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
#' res = maxOCI_K(K=2, dataSets=ds, data_type=0, a=0.2)
#' res$clusterResult$clusters
#' OBI(res$clusterResult)
#'
#' @export
OBI <- function(clusterResult){
  res = NULL
  for (i in 1:length(clusterResult$clusters)) {
    n_m = length(clusterResult$clusters[[i]])
    if(n_m==1){
      res = c(res,1)
    }else{
      A_m = clusterResult$ams[[i]]
      res = c(res,(A_m-1)/(n_m-1))
    }
  }
  return(res)
}
