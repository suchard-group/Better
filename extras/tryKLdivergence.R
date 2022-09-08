# 09/07/2022
# try out KL divergence computation

library(mclust)

gmmapprox <- function(x, K = 3){
  fit = Mclust(x, G = K, model = 'V')
  params = fit$parameters
  ps = params$pro
  means = params$mean
  sds = sqrt(params$variance$sigmasq)
  
  FUN = function(v){
    ds = numeric(length = length(v))
    for(i in 1:K){
      this.dens = ps[i] * dnorm(v, mean = means[i], sd = sds[i])
      ds = ds + this.dens
    }
    ds
  }
  
  FUN
}

KL_div <- function(x,y, K = 3){
  # x: samples from the base measure
  p = gmmapprox(x, K)
  q = gmmapprox(y, K)
  
  KL = mean(log(p(x)) - log(q(x)))
  
  KL
}