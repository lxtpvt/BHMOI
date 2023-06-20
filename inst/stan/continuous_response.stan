data{
  int<lower=0> K;  // number of subgroups
  int sizeG[K]; // the size of each group
  vector[sum(sizeG)] Y;
  real<lower=0> aY;
  real<lower=0> bY;
  real mu0;
  real<lower=0> tau0;
  real<lower=0> a;
  real<lower=0> b;
}
parameters{
  real<lower=0> tauY;
  real theta[K];  // mean of each group
  real mu;
  real<lower=0> tau;
}


model{
  int start = 0;
  // prior
  tauY ~ gamma(aY,bY);
  theta ~ normal(mu,tau^(-0.5));
  mu ~ normal(mu0,tau0^(-0.5));
  tau ~ gamma(a,b);
  // likelihood
  for(i in 1:K){
    if(i==1){
      start = 1;
    }else{
      start = sum(head(sizeG,(i-1))) + 1;
    }
    segment(Y,start,sizeG[i]) ~ normal(theta[i],tauY);
  }
}
