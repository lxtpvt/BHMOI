# calculate alpha in the gamma(alpha, beta) prior
.calGammaAlpha <- function(obi){return(obi*exp(-5*(1-obi)))}
# logit function
.logit<-function(p){return(log(p/(1-p)))}
# convert matrix to list
.toList<-function(data){
  data.list = list()
  if(is.vector(data)){
    data.list = append(data.list,list(data))
  }else if(is.data.frame(data) | is.matrix(data) | is.array(data)){
    for (i in 1:dim(data)[2]) {
      data.list = append(data.list,list(data[,i]))
    }
  }else{
    print("Acceptable data types include vector, data.frame, matrix or array.")
    data.list = NULL
  }
  return(data.list)
}

# Convert the categories of discrete data into integers
.cates2ints <- function(data1, data2){
  allCates = sort(unique(union(data1,data2))) # it's very important to keep the order for uniqueness
  d1=integer(length = length(data1))
  d2=integer(length = length(data2))
  for (i in 1:length(allCates)) {
    d1[data1==allCates[i]]=i
    d2[data2==allCates[i]]=i
  }
  return(list(int.data1=d1,int.data2=d2, orderCates = allCates))
}

.OCI_K <- function(p_m,A_m,a){
  return(p_m^a %*% A_m)
}

.OBI_m <- function(n_m,A_m){
  return((A_m-1)/(n_m-1))
}

.nextUpperBound <- function(n_pmdfs,K,a,weighted){
  if(weighted){
    return(((n_pmdfs-K)^(a+1)+K)/n_pmdfs^a)
  }else{
    return(n_pmdfs/(K+1)^a)
  }
}

# "fast": stop at the first time that maxOIC_K > maxOIC_K+1
.calMaxOCI <- function(dataSets, data_type, a, weighted, b, fast, allKs, n_sim, nthreads, kmIters,
                       seed, showmessage, grid){
  n_pmdfs = length(dataSets)
  oci.max = 0
  K.max = 1
  res = NULL
  for (K in allKs) {
    # calculate the best clustering results when fix the number of clusters K.
    bestResK = maxOCI_K(K=K, dataSets=dataSets, data_type=data_type, a=a, weighted=weighted, b=b, n_sim=n_sim,
                        nthreads=nthreads, kmIters=kmIters, seed=seed, show_progress_bar=showmessage,grid=grid)
    # when less
    if(oci.max<bestResK$maxOIC_K){
      oci.max = bestResK$maxOIC_K
      res = bestResK$clusterResult
      K.max = K
    }else if(fast){
      # the first time  that maxOIC_K > maxOIC_K+1. If fast is TRUE then return the values at K.
      return(list(K = K.max, OCI = oci.max, clusterResult = res))
    }
    # if fast is FALSE, we first need to check OCI_K > nextUpperBound is TRUE or FALSE.
    nextUB = .nextUpperBound(n_pmdfs,K,a,weighted)
    if(showmessage){
      print(paste0("The max OCI_", toString(K), " = ",toString(round(oci.max,2))," . The upper bound of OCI_",
                   toString(K+1), " = ", toString(round(nextUB,2))," ."))
    }
    if(oci.max>=nextUB | K==allKs[length(allKs)]){
      return(list(K = K.max, OCI = oci.max, clusterResult = res))
    }
  }
}
