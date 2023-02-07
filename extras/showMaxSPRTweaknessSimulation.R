# 02/06/2023
# a cleaner version of MaxSPRT simulations...

# <i> if things go with the plan
# <ii> if actual sample sizes differ from the plan
# <iii> if data looks change mid-way

# functionalities borrowed from "extras/checkMaxSPRTsimulation.R"
library(Sequential)

# (1) utils function to simulate data and run MaxSPRT----
# 
maxSprtPoisson <- function(S, numLooks=12, 
                           expectedCountByPlan = 10,
                           actualExpectedCount = 15,
                           effectSize = 1,
                           alpha = 0.05){
  
  ## pre-compute CV
  theCV = CV.Poisson(SampleSize = expectedCountByPlan * numLooks,
                     M = 4, # recommended by Kulldorff et al.
                     alpha = alpha,
                     GroupSizes = rep(expectedCountByPlan, numLooks))
  
  ## results storage
  LLRs = rep(0.01,S*numLooks)
  Rounds = rep(1:S, each = numLooks)
  Period = rep(1:numLooks, S)
  Reject = rep(FALSE, S*numLooks)
  
  for(s in 1:S){
    cat(sprintf('\n\nSimulation round %s:...\n', s))
    
    ## simulate data
    ## across all data looks in one go...
    obsCounts = rpois(numLooks, actualExpectedCount * effectSize)
    
    ## accumulate observations
    cumCounts = cumsum(obsCounts)
    cumExpected = actualExpectedCount * c(1:numLooks)
    
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
              expectedCountByPlan = expectedCountByPlan,
              actualExpectedCount = actualExpectedCount))
}

# # try a effect_size = 1 case
# simulated = maxSprtPoisson(S = 500, 
#                            numLooks=12, 
#                            expectedCountByPlan = 10,
#                            actualExpectedCount = 15,
#                            effectSize = 1,
#                            alpha = 0.05)

# (2) utils function to calculate rejection rate ----
getMaxSprtRejectRate <- function(resLs, 
                                 numLooks = 12, # the number of looks planned by MaxSPRT (not actual looks in simulated data)
                                 alpha = 0.05){
  
  res = resLs$res
  
  if(alpha != 0.05 || numLooks != resLs$numLooks){
    # allows different alpha OR different number of looks here....
    theCV = CV.Poisson(SampleSize = resLs$expectedCountByPlan * numLooks,
                       M = 4, # recommended by Kulldorff et al.
                       alpha = alpha,
                       GroupSizes = rep(resLs$expectedCountByPlan, numLooks))
    
    res$Reject = res$LLR > theCV
  }

  # overall rejection rate over all simulations
  Decisions = res %>% group_by(Simulation) %>%
    summarize(overallReject = any(Reject)) %>%
    ungroup()
  
  ## sanity check
  rejectRate = mean(Decisions$overallReject)
  
  cat(sprintf('Out of %s total simulations, %s final decisions are reached. %.1f %% are rejections.\n', 
              max(res$Simulation), nrow(Decisions), round(rejectRate * 100, 1)))
  
  # rejection rates over periods
  DecisionsByPeriod = res %>% 
    group_by(Simulation) %>%
    arrange(Period) %>%
    mutate(cumReject = cumsum(Reject)) %>%
    mutate(rejectSoFar = (cumReject > 0)) %>%
    ungroup() %>% group_by(Period) %>%
    summarise(rejectRate = mean(rejectSoFar)) %>%
    ungroup()
  
  return(list(rejectByPeriod = DecisionsByPeriod, overall = rejectRate))
}



## (3) actually run a very hypothetical experiment... 
set.seed(42)
## (i) everything goes exactly according to the plan
simulated0 = maxSprtPoisson(S = 500, 
                           numLooks=12, 
                           expectedCountByPlan = 10,
                           actualExpectedCount = 10,
                           effectSize = 1,
                           alpha = 0.05)
Rates0 = getMaxSprtRejectRate(simulated0)

## (ii) actual number of data points > anticipated
simulated1 = maxSprtPoisson(S = 500, 
                            numLooks=12, 
                            expectedCountByPlan = 10,
                            actualExpectedCount = 15,
                            effectSize = 1,
                            alpha = 0.05)
Rates1 = getMaxSprtRejectRate(simulated1)

## (iii) actually planned out a 1-year surveillance, but had to do 2-year
simulated2 = maxSprtPoisson(S = 500, 
                            numLooks=24, 
                            expectedCountByPlan = 10,
                            actualExpectedCount = 10,
                            effectSize = 1,
                            alpha = 0.05)
Rates2 = getMaxSprtRejectRate(simulated1, numLooks = 12)


