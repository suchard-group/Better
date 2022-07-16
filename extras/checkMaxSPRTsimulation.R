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

# another function:
## pre-compute critical values (CVs)
## and then simulate to get llr's
maxSprtPoisson <- function(S, numLooks=12, 
                           expectedCount = 5,
                           effectSize = 1,
                           alpha = 0.05){
  
  ## pre-compute CV
  theCV = CV.Poisson(SampleSize = expectedCount * numLooks,
                     M = 4, # recommended by Kulldorff et al.
                     alpha = alpha,
                     GroupSizes = rep(expectedCount, numLooks))
  
  ## results storage
  LLRs = rep(0.01,S*numLooks)
  Rounds = rep(1:S, each = numLooks)
  Period = rep(1:numLooks, S)
  Reject = rep(FALSE, S*numLooks)
  
  for(s in 1:S){
    cat(sprintf('\n\nSimulation round %s:...\n', s))
    
    ## simulate data
    ## across all data looks in one go...
    obsCounts = rpois(numLooks, expectedCount * effectSize)
    
    ## accumulate observations
    cumCounts = cumsum(obsCounts)
    cumExpected = expectedCount * c(1:numLooks)
    
    ## evaluate LLR's
    this.LLRs = dpois(cumCounts, cumCounts, log = TRUE) - 
      dpois(cumCounts, cumExpected, log = TRUE)
    
    ## reject?
    this.reject = this.LLRs > theCV
    
    ## save results
    inds = seq(from = (s-1)*numLooks + 1,
               to = (s-1)*numLooks + numLooks,
               by = 1)
    
    LLRs[inds] = this.LLRs
    Reject[inds] = this.reject
    
    if(any(this.reject)){
      cat(sprintf('Null rejected in simulation round %s!\n', s))
    }
  }
  
  res = data.frame(Simulation = Rounds,
                   Period = Period,
                   LLR = LLRs,
                   Reject = Reject)
  
  return(list(res = res,
              numLooks = numLooks,
              expectedCount = expectedCount))
}

## utils function to get rejection
getMaxSprtRejectRate <- function(resLs, alpha = 0.05){
  
  res = resLs$res
  
  if(alpha != 0.05){
    theCV = CV.Poisson(SampleSize = resLs$expectedCount * resLs$numLooks,
                       M = 4, # recommended by Kulldorff et al.
                       alpha = alpha,
                       GroupSizes = rep(resLs$expectedCount, resLs$numLooks))
    
    res$Reject = res$LLR > theCV
  }
  
  Decisions = res %>% group_by(Simulation) %>%
    summarize(overallReject = any(Reject)) %>%
    ungroup()
  
  ## sanity check
  cat(sprintf('Out of %s total simulations, %s final decisions are reached.\n', 
              max(res$Simulation), nrow(Decisions)))
  
  mean(Decisions$overallReject)
}


## try simulation run----

# cachepath = '~/Documents/Research/better/localCache/'

# rejects = simulatePoissonTest(1, numLooks = 12, 
#                               Nbatch = 100,
#                               expectedCount = 5,
#                               effectSize = 1, 
#                               alphaSpend = 'Wald',
#                               cachePath = cachepath)
# # somehow there is an internal stuff that doesn't allow me to repeat testing inside the function???

# try the updated function --- a lot faster!!
# resLs = maxSprtPoisson(10, numLooks = 12, expectedCount = 3, effectSize = 1, alpha = 0.05)

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



# run simulations using the updated, faster function-----
nL = 12
#Nt = 60
eC = 3

alpha = 0.05

S = 5000

set.seed(41)

# (1) effect = 1

effect = 1

simsMPRT1 = maxSprtPoisson(S, numLooks = nL, 
                           expectedCount = eC,
                           effectSize = effect,
                           alpha = alpha)
getMaxSprtRejectRate(simsMPRT1, alpha = alpha) # >0.14 if alpha = 0.05???
getMaxSprtRejectRate(simsMPRT1, alpha = 0.0245) # 0.0474


# (2) effect = 1.5

effect = 1.5

simsMPRT15 = maxSprtPoisson(S, numLooks = nL, 
                           expectedCount = eC,
                           effectSize = effect,
                           alpha = alpha)
1-getMaxSprtRejectRate(simsMPRT15, alpha = alpha) # 0.229 if alpha = 0.05
1-getMaxSprtRejectRate(simsMPRT15, alpha = 0.0245) # 0.336


# (3) effect = 2

effect = 2

simsMPRT20 = maxSprtPoisson(S, numLooks = nL, 
                            expectedCount = eC,
                            effectSize = effect,
                            alpha = alpha)
1-getMaxSprtRejectRate(simsMPRT20, alpha = alpha) # 0.0022 if alpha = 0.05
1-getMaxSprtRejectRate(simsMPRT20, alpha = 0.0245) # 0.0048


# (4) effect = 4

effect = 4

simsMPRT40 = maxSprtPoisson(S, numLooks = nL, 
                            expectedCount = eC,
                            effectSize = effect,
                            alpha = alpha)
1-getMaxSprtRejectRate(simsMPRT40, alpha = alpha) # 0 alpha = 0.05
1-getMaxSprtRejectRate(simsMPRT40, alpha = 0.0245) # 0





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




