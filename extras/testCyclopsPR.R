# test Cyclops Poisson regression fit
# 04/14/2022: mystery solved!!
# when ex. (x1,y1) = (1,0); (x2,y2) = c(0, ++)
# log likelihood of PR is monotonously decreasing...
# so will return a very negative coef
# which should actually be -infty

library(Cyclops)

# data directly copied from a specific Historical Comparator anlaysis
data = data.frame(cohortCounts = c(0,402), 
                  personYears = c(1171.9828,12200647.8653), 
                  exposed = c(1,0))

# use Cyclops
cyclopsData = createCyclopsData(cohortCounts ~ exposed, 
                                data = data, 
                                offset = log(personYears),
                                modelType = 'pr')
cyclopsFit <- fitCyclopsModel(cyclopsData)
coef(cyclopsFit)

# use glm
pr = glm(cohortCounts ~ exposed,
         offset =log(personYears),
         data=data, family = "poisson")
pr$coefficients

## -----
# see if zero responses generally create problems for Poisson regression

generateData <- function(n = 1){
  y1 = rep(0,n)
  y0 = sample(1:50, size = n, replace=TRUE)
  
  z1 = sample(50:200, size = n, replace=TRUE)
  z0 = sample(1000:2000, size = n, replace=TRUE)
  
  x = rep(c(1,0),each = n)
  
  data.frame(x=x, y=c(y1,y0), z=c(z1, z0))
}

N = 1000
dat = generateData(N)

## Cyclops fit
cyclopsData = createCyclopsData(y ~ x, 
                                data = dat, 
                                offset = log(z),
                                modelType = 'pr')
cyclopsFit <- fitCyclopsModel(cyclopsData)
coef(cyclopsFit)

## base R glm
pr = glm(y ~ x, 
         data = dat, 
         offset = log(z),
         family = "poisson")
pr$coefficients

# 1. if we have a bunch of zeros (half of the observations): 
#    then will get a really negative coefficient no matter what

# 2. why does base R glm produce different estimates from Cyclops????


# compare Cyclops and glm using regular-looking data  ------
N = 1000
x = rnorm(N, sd = 2)
beta = 0.6
inter = 0.2
gammas = inter + beta * x
y = rpois(N, lambda = exp(gammas))

dat = data.frame(y=y, x=x)

## Cyclops fit
cyclopsData = createCyclopsData(y ~ x, 
                                data = dat, 
                                modelType = 'pr')
cyclopsFit <- fitCyclopsModel(cyclopsData)
coef(cyclopsFit)
## base R glm
pr = glm(y ~ x, 
         data = dat, 
         family = "poisson")
pr$coefficients

# if the data are okay, then Cyclops and glm produce same results!
# WTF?????
