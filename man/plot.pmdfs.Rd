% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.pmdfs}
\alias{plot.pmdfs}
\title{Plot a set of distributions}
\usage{
\method{plot}{pmdfs}(
  data,
  distribution_type,
  clusterResult = NULL,
  centroid = TRUE,
  centroid.colors = NULL,
  rep = T,
  jump = 1,
  xlim = NULL,
  ylim = NULL,
  title = "",
  lineSize = 0.8,
  lineTypes = NULL,
  centroidLineSize = 1.5,
  int_hue_pal = NULL,
  colors = NULL,
  ...
)
}
\arguments{
\item{data}{A list of vectors. Each vector stores a dataset that represents a distribution.}

\item{distribution_type}{A char variable to indicate the types of distributions, "pdf": continuous or "pmf": discrete.}

\item{clusterResult}{A list stores the best clustering result which includes the clusters, ams and centroids (check ?KMeansCluster).}

\item{centroid}{A bool variable controlling plot the centroid distribution or not.}

\item{xlim}{A vector of two values represent the limit of horizontal axis.}

\item{ylim}{A vector of two values represent the limit of vertical axis.}

\item{lineSize}{A real value to set the size of line.}
}
\description{
Plot distributions that are estimated from the datasets stored in \code{data}.
}
\details{
Plot distributions that are estimated from the datasets stored in \code{data}.
}
\examples{

# Discrete distributions
set.seed(123)
data1 = rbinom(n=50,size=20,prob=0.3)
data2 = rbinom(n=30,size=10,prob=0.5)
data3 = rpois(35,4)
data4 = rpois(45,9)
data5 = rpois(15,8)
ds = list(data1,data2,data3,data4,data5)
plot.pmdfs(data=ds,distribution_type="pmf")

}
