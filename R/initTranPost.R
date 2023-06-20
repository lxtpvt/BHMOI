
#' @export
initTranPostBinary <- function(J, W, nGroups, alpha=0.01, beta=10, tau0=0.01, mu0=.logit(0.1), seed=NULL, nIter=2000){

  mydata = list(Y=J, n=W, K=nGroups, a = alpha, b = beta, mu0 = mu0, tau0 = tau0)
  if(is.null(seed)){
      temp_M = rstan::extract(rstan::sampling(stanmodels$tran_bhm_general,mydata,iter = nIter, chains = 1,
                                              show_messages=F,refresh=F))$p_k
  }else{
    temp_M = rstan::extract(rstan::sampling(stanmodels$tran_bhm_general,mydata,iter = nIter, chains = 1,
                                            seed=seed, show_messages=F,refresh=F))$p_k
  }
  resL = .toList(temp_M)
  class(resL)<-"pmdfs"
  return(resL)
}

#' @export
initTranPostContinuous <- function(Y, sizeG, nGroups, alphaY=1, betaY=10, alpha=0.01, beta=10, tau0=0.01,
                         mu0=.logit(0.1), seed=NULL, nIter=2000){
  Y = unlist(Y)
  mydata = list(K=nGroups, sizeG=sizeG, Y=Y, aY=alphaY, bY=betaY, mu0 = mu0, tau0 = tau0, a = alpha, b = beta)
  if(is.null(seed)){
    temp_M = rstan::extract(rstan::sampling(stanmodels$continuous_response,mydata,iter = nIter, chains = 1,
                                            show_messages=F,refresh=F))$theta
  }else{
    temp_M = rstan::extract(rstan::sampling(stanmodels$continuous_response,mydata,iter = nIter, chains = 1,
                                            seed=seed, show_messages=F,refresh=F))$theta
  }
  resL = .toList(temp_M)
  class(resL)<-"pmdfs"
  return(resL)
}
