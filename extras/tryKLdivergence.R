# 09/07/2022
# try out KL divergence for consecutive posterior samples

## 1. from stack exchange... ------
## this works pretty bad, numerically??

my.ecdf  <-  function(x)   {
  x   <-   sort(x)
  x.u <-   unique(x)
  n  <-  length(x) 
  x.rle  <-  rle(x)$lengths
  y  <-  (cumsum(x.rle)-0.5) / n
  FUN  <-  approxfun(x.u, y, method="linear", yleft=0, yright=1,
                     rule=2)
  FUN
}  

KL_est  <-  function(x, y)   {
  # x are samples from the base measure...
  dx  <-  diff(sort(unique(x)))
  dy  <-  diff(sort(unique(y)))
  ex  <-  min(dx) ; ey  <-  min(dy)
  e   <-  min(ex, ey)/2
  n   <-  length(x)    
  P  <-   my.ecdf(x) ; Q  <-  my.ecdf(y)
  KL  <-  sum( log( (P(x)-P(x-e))) - log(Q(x)-Q(x-e))) / n - 1
  KL              
}

# ## test a little bit
# x = rnorm(100)
# #y <- rt(100, df=5) # this works...
# y = rnorm(100, mean = 1) # this doesn't!!
# 
# KL_est(x, y)


## 2. use a naive, normal mixture approach
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

# ## test it out....
# x = rnorm(1000)
# y = rnorm(1000)
# KL_div(x,y)
