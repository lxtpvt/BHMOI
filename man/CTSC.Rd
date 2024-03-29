% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BHMOI_Continuous.R
\name{CTSC}
\alias{CTSC}
\title{BHMOI for clinical trials with subgroups and Continuous outcomes (CTSC)}
\usage{
CTSC(
  Y,
  sizeG,
  nGroups,
  a = 0.3,
  weighted = FALSE,
  b = 5,
  alphaY = 1,
  betaY = 10,
  alpha.min = 1,
  alpha.max = 50,
  beta = 10,
  tau0 = 0.01,
  mu0 = .logit(0.1),
  K = 0,
  fast = 1,
  n_sim = 3,
  seed = 123,
  nIter = 2000,
  kmIters = 50,
  grid = 1000
)
}
\arguments{
\item{Y}{A list of vectors. Each vector includes a subgroup of observations.}

\item{sizeG}{Number of observations in each subgroup.}

\item{nGroups}{Number of subgroups.}

\item{a}{A real number in the interval \link{0,1}. It is used to control the number of clusters. The bigger the \code{a}, the more clusters we get.}

\item{weighted}{A bool value to indicate the K-Means algorithm is weighted (include parameter b) or not.}

\item{b}{A real number which is the power parameter of generalized weighted K-Means algorithm. When the weighted=True, it is used to control the number of clusters. The bigger b is, the more extreme clusters (large or small) are generated.}

\item{alphaY}{A real value denotes the parameter of the gamma prior for the precision of the real data.}

\item{betaY}{A real value denotes the parameter of the gamma prior for the precision of the real data.}

\item{alpha.min}{A real value denotes the parameter that correspond approximately to no borrowing in the gamma prior.}

\item{alpha.max}{A real value denotes the parameter that correspond approximately to full borrowing in the gamma prior.}

\item{beta}{A real value denotes the parameter in the gamma prior.}

\item{tau0}{A real value denotes the precision of the hyper-parameter of the cluster mean.}

\item{mu0}{A real value denotes the mean of the hyper-parameter of the cluster mean.}

\item{K}{A integer to set the number of cluster.}

\item{fast}{A bool value to select the strategy to calculate the next level of K. 1: stop at the first time K does not increase; 0: stop at OCI_K>OCI_{K+1}.}

\item{n_sim}{A integer to set the repeat times of weighted K-Means algorithm. K-Means algorithm can stuck at a local optimal. Increase n_sim can increase the probability to find the global optimal.}

\item{seed}{A integer to set the number of threads for computation.}

\item{nIter}{A integer to set the number of MCMC iterations.}

\item{kmIters}{A integer to set the max number of iterations in K-Means algorithm.}

\item{grid}{A integer the number of discrete points to calculate the probability density.}
}
\value{
A list that stores the results of Bayesian hierarchical model.
}
\description{
BHMOI for clinical trials with subgroups and Continuous outcomes (CTSC)
}
\details{
There are a set datasets in \code{dataSets}. Each dataset includes discrete data points of a distribution. This function will search the best clustering result of the distributions. The number of clusters can be influenced by parameter \code{a}. The bigger the \code{a}, the more clusters we get.
}
\examples{

nGroups = 10
W = rep(25,nGroups)
J = c(2, 3,  5,  6,  8, 6,  9, 11, 13,7)
res = CTSB(J=J, W=W, nGroups=nGroups, a=0.35)
res$OBIs
plot(data=res$initPost,distribution_type="pdf",ylim = c(0,9))
plot(data=res$initPost,distribution_type="pdf",clusterResult=res$maxOCIres$clusterResult,ylim = c(0,9))
plot(data=res$borrowedPost,distribution_type="pdf",clusterResult=res$maxOCIres$clusterResult,ylim = c(0,9))

}
\references{
The paper
\emph{Distribution-free Overlapping Indices for Efficient and Robust Information Borrowing in Bayesian Hierarchical Modeling}.
}
