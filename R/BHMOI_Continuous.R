#' BHMOI for clinical trials with subgroups and Continuous outcomes (CTSC)
#'
#'
#'
#' @param sizeG Number of observations in each subgroup.
#' @param Y A list of vectors. Each vector includes a subgroup of observations.
#' @param nGroups Number of subgroups.
#' @param a A real number in the interval [0,1]. It is used to control the number of clusters. The bigger the \code{a}, the more clusters we get.
#' @param weighted A bool value to indicate the K-Means algorithm is weighted (include parameter b) or not.
#' @param b A real number which is the power parameter of generalized weighted K-Means algorithm. When the weighted=True, it is used to control the number of clusters. The bigger b is, the more extreme clusters (large or small) are generated.
#' @param alphaY A real value denotes the parameter of the gamma prior for the precision of the real data.
#' @param betaY A real value denotes the parameter of the gamma prior for the precision of the real data.
#' @param alpha.min A real value denotes the parameter that correspond approximately to no borrowing in the gamma prior.
#' @param alpha.max A real value denotes the parameter that correspond approximately to full borrowing in the gamma prior.
#' @param beta A real value denotes the parameter in the gamma prior.
#' @param tau0 A real value denotes the precision of the hyper-parameter of the cluster mean.
#' @param mu0 A real value denotes the mean of the hyper-parameter of the cluster mean.
#' @param K A integer to set the number of cluster.
#' @param fast A bool value to select the strategy to calculate the next level of K. 1: stop at the first time K does not increase; 0: stop at OCI_K>OCI_{K+1}.
#' @param n_sim A integer to set the repeat times of weighted K-Means algorithm. K-Means algorithm can stuck at a local optimal. Increase n_sim can increase the probability to find the global optimal.
#' @param seed A integer to set the number of threads for computation.
#' @param nIter A integer to set the number of MCMC iterations.
#' @param kmIters A integer to set the max number of iterations in K-Means algorithm.
#' @param grid A integer the number of discrete points to calculate the probability density.
#'
#' @details There are a set datasets in \code{dataSets}. Each dataset includes discrete data points of a distribution. This function will search the best clustering result of the distributions. The number of clusters can be influenced by parameter \code{a}. The bigger the \code{a}, the more clusters we get.
#'
#' @return A list that stores the results of Bayesian hierarchical model.
#'
#' @references The paper
#' \emph{Distribution-free Overlapping Indices for Efficient and Robust Information Borrowing in Bayesian Hierarchical Modeling}.
#'
#' @examples
#'
#' nGroups = 10
#' W = rep(25,nGroups)
#' J = c(2, 3,  5,  6,  8, 6,  9, 11, 13,7)
#' res = CTSB(J=J, W=W, nGroups=nGroups, a=0.35)
#' res$OBIs
#' plot(data=res$initPost,distribution_type="pdf",ylim = c(0,9))
#' plot(data=res$initPost,distribution_type="pdf",clusterResult=res$maxOCIres$clusterResult,ylim = c(0,9))
#' plot(data=res$borrowedPost,distribution_type="pdf",clusterResult=res$maxOCIres$clusterResult,ylim = c(0,9))
#'
# #' @export
CTSC<-function(Y, sizeG, nGroups, a=0.3, weighted=FALSE, b=5, alphaY=1, betaY=10, alpha.min=1,
               alpha.max=50, beta=10, tau0=0.01, mu0=.logit(0.1), K=0, fast=1, n_sim=3, seed=123, nIter=2000,kmIters=50, grid=1000){
  # Step 1, run initApproxDists()
  init_thetas <- initTranPostContinuous(Y=Y, sizeG=sizeG, nGroups=nGroups, alphaY=alphaY, betaY=betaY, alpha=alpha.min/100,
                              beta=beta, tau0=tau0, mu0=mu0, seed=seed, nIter=nIter)
  # Step 2, get best clustering result by maximizing OCI
  oci.res <- clusterMaxOCI(dataSets=init_thetas, data_type = 1, a=a, weighted=weighted, b=b,
                           fast=fast, levels=0, K=K, n_sim=n_sim, showmessage=1,kmIters=kmIters, seed=seed, grid=grid)

  # get clusters and their ids
  clusters_ids = oci.res$clusterResult$clusters
  # Step 3, calculate OBIs and corresponding borrow strength.
  obis = OBI(oci.res$clusterResult)
  # calculate borrow strength specified by alpha in gamma(alpha, beta) with a fixed beta.
  als = .calGammaAlpha(obis)*(alpha.max-alpha.min)+alpha.min
  #Step 4, run model in (22)
  temp_M = NULL
  for (j in 1:length(clusters_ids)) {
    K_m = length(clusters_ids[[j]])
    if(K_m==1){ # cluster only has one distribution, no information borrowing.
      thetas_m = matrix(init_thetas[[(clusters_ids[[j]])]],ncol = 1)
    }else{ # cluster has more than one distributions
      Y_m = Y[clusters_ids[[j]]]
      sizeG_m = sizeG[clusters_ids[[j]]]
      # cluster centroid mean. How about median and mode?
      c_mu = oci.res$clusterResult$centroids[[j]][[1]][which.max(oci.res$clusterResult$centroids[[j]][[2]])]
      Y_m = unlist(Y_m)
      mydata = list(K=K_m, sizeG=sizeG_m, Y=Y_m, aY=alphaY, bY=betaY, mu0 = c_mu, tau0 = tau0, a = als[j], b = beta)
      if(is.null(seed)){
        thetas_m = rstan::extract(rstan::sampling(stanmodels$continuous_response,mydata, iter = nIter,
                                                  chains = 1, show_messages=F,refresh=F))$theta
      }else{
        thetas_m = rstan::extract(rstan::sampling(stanmodels$continuous_response,mydata,seed=seed,
                                                  iter = nIter, chains = 1, show_messages=F,refresh=F))$theta
      }
    }
    # calculate the posterior mode. We should use mode instead of mean because the skewness.
    temp_M = cbind(temp_M, thetas_m)
  }
  # reorder temp_M
  thetasMat=temp_M
  thetasMat[,unlist(clusters_ids)] = temp_M
  resL = .toList(thetasMat)
  class(resL)<-"pmdfs"
  return(list(initPost=init_thetas, borrowedPost=resL, maxOCIres=oci.res, OBIs=obis))
}

