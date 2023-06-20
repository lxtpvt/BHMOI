data{
  int<lower=0> K;  // number of subgroups
  int Y[K];  // number of response
  int n[K];  // number of tials
  real mu0;
  real<lower=0> tau0;
  real<lower=0> a;
  real<lower=0> b;
}

parameters{
  vector<lower=0,upper=1>[K] p_k;
  real mu;
  real<lower=0> tau;
}

transformed parameters{
  vector[K] theta;
  theta = logit(p_k);
}

model{
  // prior
  theta ~ normal(mu,tau^(-0.5));
  mu ~ normal(mu0,tau0^(-0.5));
  tau ~ gamma(a,b);
  // likelihood
  for(j in 1:K){
    Y[j] ~ binomial(n[j],p_k[j]);
  }
}
