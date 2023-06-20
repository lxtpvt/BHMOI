#' Estimate empirical probability mass function
#'
#' The pmf is estimated from data \code{x}.
#'
#' @param x A vector of data for estimating the pmf.
#'
#' @return A dataframe with its first column including all categories and the second column including the empirical pmf.
#'
#' @examples
#' # Generate data
#' data(iris)
#' x <- sample(iris$Species, 20)
#'
#' # Estimate the pmf.
#' epmf <- epmf(x)
#' epmf
#'
#' @export
epmf <- function(x){
  epmf=NULL
  allCates = unique(x)
  #if(is.numeric(allCates)){allCates = order(allCates)}
  for (e in allCates) {
    count = 0
    for (ex in x) {
      if(ex==e){
        count=count+1
      }
    }
    epmf = c(epmf, count/length(x))
  }
  res = data.frame(category = allCates,epmf = epmf)
  res = res[order(res$category),]
  return(res)
}

#' Estimate empirical probability density function
#'
#' The pdf is estimated from data \code{x}.
#'
#' @param x A vector of data for estimating the pdf.
#' @param npos A vector of positions to discrete the pdf.
#' @param scale A real value to scale the support of pdf. To keep the accuracy of estimation, it should be less than 0.15.
#'
#' @return A list values (positions, probabilities) of the pdf.
#'
#' @examples
#' # Generate data
#' x <- rnorm(1000)
#'
#' # Estimate the pdf.
#' res <- epdf(x)
#' plot(x=res$positions, y=res$probabilities, type="l")
#'
#' @export
epdf <- function(x, npos=1000, scale=0.1){
  l = min(x)*(1-scale)
  r = max(x)*(1+scale)
  pos = seq(from=l, to=r, by = (r-l)/(npos-1))
  return(list(positions = pos, probabilities = .kde(x,pos)))
}


#' @export
epdf.cluster <- function(cluster_data_list, npos=1000, scale=0.1){
  l = min(unlist(cluster_data_list))*(1-scale)
  r = max(unlist(cluster_data_list))*(1+scale)
  pos = seq(from=l, to=r, by = (r-l)/(npos-1))
  return(list(positions = pos, probabilities = .clusterKDE(cluster_data_list,pos)))
}