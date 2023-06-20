#' BHMOI for Meta-analysis
#'
#'
#'
#' @param mus A vector includes the sample means of studies.
#' @param sigmas A vector includes the sample standard devistions of studies.
#' @param nSample A integer denotes the sample size that you want to draw from every study.
#' @param nGroups A integer denotes the number of studies.
#' @param a A real number in the interval [0,1]. It is used to control the number of clusters. The bigger the \code{a}, the more clusters we get.
#' @param weighted A bool value to indicate the K-Means algorithm is weighted (include parameter b) or not.
#' @param b A real number which is the power parameter of generalized weighted K-Means algorithm.
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
#' mus = c(0.7091362, 0.3548641, 1.7911700, 0.1824552, 0.4218509, 0.6300000, 0.7248838, 0.5286638, 0.2840000,
#'          1.2750682, 0.1036082, 0.3883906, 0.5407398, 0.4261593, 0.5153969, 1.4797260, 0.6125782, 0.6000000)
#' nGroups = length(mus)
#' sigmas = c(0.2608202, 0.1963624, 0.3455692, 0.1177874, 0.1448128, 0.1960000, 0.2246641, 0.2104609, 0.1680000,
#'            0.3371997, 0.1947275, 0.2307689, 0.2443133, 0.2579379, 0.3512737, 0.3152817, 0.2266834, 0.2490000)
#' res = Meta(mus, sigmas, nGroups)
#' plot(data=res$initPost,distribution_type="pdf",ylim = c(0,9))
#' plot(data=res$initPost,distribution_type="pdf",clusterResult=res$maxOCIres$clusterResult,ylim = c(0,9))
#' plot(data=res$borrowedPost,distribution_type="pdf",clusterResult=res$maxOCIres$clusterResult,ylim = c(0,9))
#'
# #' @export
Meta<-function(mus, sigmas, nSample, nGroups, a=0.3, weighted=FALSE, b=5, alphaY=1, alpha.min=1, betaY=10, alpha.max=50,
               beta=10, tau0=0.01, mu0=.logit(0.1), K=0, fast=1, n_sim=3, seed=NULL, nIter=2000,kmIters=50, grid=1000){
  # Step 1, generate data
  Y = list()
  sizeG = rep(nSample,length(mus))
  for (i in 1:length(mus)) {
    Y = append(Y, list(rnorm(n=sizeG[i], mean = mus[i], sd = sigmas[i])))
  }

  res <- CTSC(Y, sizeG, nGroups, a=a, weighted=weighted, b=b, alphaY=alphaY, betaY=betaY, alpha.min=alpha.min,
              alpha.max=alpha.max, beta=beta, tau0=tau0, mu0=mu0, K=K, fast=fast, n_sim=n_sim, seed=seed, nIter=nIter,kmIters=kmIters, grid=grid)
  return(res)
}

