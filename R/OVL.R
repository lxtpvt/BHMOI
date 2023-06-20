#' Calculate the coefficient of overlapping of two distributions
#'
#' The distributions are estimated from \code{data1} and \code{data2}
#' respectively.
#'
#' @param data1 A vector of data for estimating distribution 1.
#' @param data2 A vector of data for estimating distribution 2.
#' @param is.continuous A bool variable to indicate the types of two distributions, 1: continuous or 0: discrete.
#'
#' @details The two distributions must be the same type. The support of two distributions are estimated from their corresponding data.
#'
#' @return A real value between 0 and 1 to measure the coefficient of overlapping between two distributions.
#'
#' @references
#' \emph{https://github.com/duncanmcn/kdepp}.
#'
#' @examples
#' # Continuous distributions
#' # Generate two datasets
#' data1 <- rnorm(100,0.2)
#' data2 <- rbeta(80,2,5)
#'
#' # Calculate the coefficient of two distributions.
#' ovl <- OVL(data1, data2)
#' ovl
#'
#'# Discrete distributions
#' data(iris)
#' data1 <- sample(iris$Species,20)
#' data2 <- sample(iris$Species,30)
#'
#' # Calculate the coefficient of two distributions.
#' ovl <- OVL(data1, data2, F)
#' ovl
#'
#' @export
OVL <- function(data1,data2,is.continuous=TRUE,grids=1000){
  if(is.continuous){
    return(.ovl_continuous(data1,data2,grids))
  }else{
    # Convert the datasets into integers.
    int.datasets = .cates2ints(data1,data2)
    return(.ovl_discrete(int.datasets$int.data1,int.datasets$int.data2))
  }
}
