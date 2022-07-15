# 07/14/2022
# simulation to check error rates of MaxSPRT

library(Sequential)
library(dplyr)

# use a Poisson model first
# 1. MaxSPRT -------

## a function to simulate Poisson data and run MaxSPRT
simulatePoissonTest <- function(S, numLooks=12, 
                                Nbatch=100, 
                                expectedCount = 5,
                                effectSize = 1, 
                                alphaSpend = 'Wald',
                                cachePath = './localCache/'){
  
  rejects = 0
  
  for(s in 1:S){
    cat(sprintf('Simulation round %s:...\n', s))
    testName = sprintf('Poisson%s', s)
    
    if(file.exists(file.path(cachepath, sprintf('%s.txt',testName)))){
      file.remove(file.path(cachepath, sprintf('%s.txt',testName)))
    }
    
    AnalyzeSetUp.Poisson(name = testName, 
                         SampleSize = expectedCount * numLooks,
                         M = 4, # recommended by Kulldorff and Silva
                         AlphaSpendType = alphaSpend,
                         address = cachePath)
    
    #AllRes = NULL
    
    for(i in 1:numLooks){
      ## simulate Poisson count data
      obsCounts = sum(rpois(Nbatch, expectedCount/Nbatch * effectSize))
      
      ## test
      res = Analyze.Poisson('Poisson1', test = i, 
                            mu0 = expectedCount, 
                            events = obsCounts)
      
      if(res[1,'Reject H0'] == 'Yes'){
        cat(sprintf('Null rejected at %sth data look.\n', i))
        rejects = rejects + 1
        break
      }
    }
    
    #rm()
  }
  
  return(rejects)
}


## try simulation run

cachepath = '~/Documents/Research/better/localCache/'

# rejects = simulatePoissonTest(1, numLooks = 12, 
#                               Nbatch = 100,
#                               expectedCount = 5,
#                               effectSize = 1, 
#                               alphaSpend = 'Wald',
#                               cachePath = cachepath)
# # somehow there is an internal stuff that doesn't allow me to repeat testing inside the function???


# we'll just run it one by one then...----
nL = 12
Nt = 60
eC = 3
S = 500

# (1) effect = 1

set.seed(41)

effect = 1

allRejects = 0
for(s in 1:S){
  cat(sprintf('\n\n ROUND %s ...\n\n', s))
  reject = simulatePoissonTest(1, numLooks = nL, 
                               Nbatch = Nt,
                               expectedCount = eC,
                               effectSize = effect, 
                               alphaSpend = 'Wald',
                               cachePath = cachepath)
  allRejects = allRejects + reject
}

(Type1error = allRejects/S) # 0.006


# (2) effect = 1.5

effect = 1.5

allRejects = 0
for(s in 1:S){
  cat(sprintf('\n\n ROUND %s ...\n\n', s))
  reject = simulatePoissonTest(1, numLooks = nL, 
                               Nbatch = Nt,
                               expectedCount = eC,
                               effectSize = effect, 
                               alphaSpend = 'Wald',
                               cachePath = cachepath)
  allRejects = allRejects + reject
}

(Type2error15 = 1-allRejects/S)


# (2) effect = 2

effect = 2

allRejects = 0
for(s in 1:S){
  cat(sprintf('\n\n ROUND %s ...\n\n', s))
  reject = simulatePoissonTest(1, numLooks = nL, 
                               Nbatch = Nt,
                               expectedCount = eC,
                               effectSize = effect, 
                               alphaSpend = 'Wald',
                               cachePath = cachepath)
  allRejects = allRejects + reject
}

(Type2error20 = 1-allRejects/S)



# (3) effect = 4

effect = 4

allRejects = 0
for(s in 1:S){
  cat(sprintf('\n\n ROUND %s ...\n\n', s))
  reject = simulatePoissonTest(1, numLooks = nL, 
                               Nbatch = Nt,
                               expectedCount = eC,
                               effectSize = effect, 
                               alphaSpend = 'Wald',
                               cachePath = cachepath)
  allRejects = allRejects + reject
}

(Type2error40 = 1-allRejects/S)


# 2. Bayesian test -----

## a function to simulate data streams
## and also get posterior probs
simulatePoissonBayes <- function(S, numLooks=12, 
                                 Nbatch = 100, 
                                 expectedCount = 5,
                                 effectSize = 1, 
                                 priorAlpha = 0.01,
                                 priorBeta = 0.01,
                                 alpha = 0.05){
  
  # results storage
  P1s = rep(0.01,S*numLooks)
  Rounds = rep(1:S, each = numLooks)
  Period = rep(1:numLooks, S)
  Reject = rep(FALSE, S*numLooks)
  
  # the expected rate
  expectedRate = expectedCount/Nbatch
  
  for(s in 1:S){
    cat(sprintf('\n\nSimulation round %s:...\n', s))
    
    # set up posterior starting point
    postAlpha = priorAlpha
    postBeta = priorBeta

    flag = FALSE
    
    for(i in 1:numLooks){
      ## simulate Poisson count data
      obsCounts = sum(rpois(Nbatch, expectedRate * effectSize))
      
      ## update posterior
      postAlpha = postAlpha + obsCounts
      postBeta = postBeta + Nbatch
      
      ## get P1
      this.P1 = pgamma(expectedRate, 
                       shape = postAlpha, rate = postBeta,
                       lower.tail = FALSE) #evaluate the upper tail prob
      
      ## save the P1 entry
      ind = (s-1)*numLooks + i
      P1s[ind] = this.P1
      
      ## decision making using delta_1 = 1-alpha
      if(this.P1 > 1-alpha){
        Reject[ind] = TRUE
        if(!flag){
          flag = TRUE
          cat(sprintf('\nNull rejected at data look %s!\n',i))
        }
      }

    }
  }
  
  # results saved in a dataframe
  res = data.frame(Simulation = Rounds,
                   Period = Period,
                   P1 = P1s,
                   Reject = Reject)
  
  return(res)
}

## another function to make decision and check rejection rate
getRejectionRate <- function(res, alpha = 0.05){
  
  if(alpha != 0.05){
    res$Reject = res$P1 > 1-alpha
  }
  
  Decisions = res %>% group_by(Simulation) %>%
    summarize(overallReject = any(Reject)) %>%
    ungroup()
  
  ## sanity check
  cat(sprintf('Out of %s total simulations, %s final decision are reached.\n', 
              max(res$Simulation), nrow(Decisions)))
  
  mean(Decisions$overallReject)
}

## try it

nL = 12
Nt = 60
eC = 3

S = 5000

pAlpha = 0.01
pBeta = 0.01
alphaDefault = 0.05

set.seed(42)

# (1) effect = 1

effect = 1

sims1 = simulatePoissonBayes(S, numLooks = nL,
                            Nbatch = Nt, 
                            expectedCount = eC,
                            effectSize = effect,
                            priorAlpha = pAlpha,
                            priorBeta = pBeta,
                            alpha = alpha)

getRejectionRate(sims1, alpha = 0.0155) # calibrated at 0.05; 0.0492


# (2) effect = 1.5

effect = 1.5

sims15 = simulatePoissonBayes(S, numLooks = nL,
                             Nbatch = Nt, 
                             expectedCount = eC,
                             effectSize = effect,
                             priorAlpha = pAlpha,
                             priorBeta = pBeta,
                             alpha = alpha)

1 - getRejectionRate(sims15, alpha = 0.0155) # errorRate = 0.2126


# (3) effect = 2

effect = 2

sims20 = simulatePoissonBayes(S, numLooks = nL,
                              Nbatch = Nt, 
                              expectedCount = eC,
                              effectSize = effect,
                              priorAlpha = pAlpha,
                              priorBeta = pBeta,
                              alpha = alpha)

1 - getRejectionRate(sims20, alpha = 0.0155) # errorRate = 0.0018


# (4) effect = 4

effect = 4

sims40 = simulatePoissonBayes(S, numLooks = nL,
                              Nbatch = Nt, 
                              expectedCount = eC,
                              effectSize = effect,
                              priorAlpha = pAlpha,
                              priorBeta = pBeta,
                              alpha = alpha)

1 - getRejectionRate(sims40, alpha = 0.0155) # errorRate = 0.0




