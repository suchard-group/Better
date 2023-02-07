# 02/06/2023
# a cleaner version of MaxSPRT simulations...

# <i> if things go with the plan
# <ii> if actual sample sizes differ from the plan
# <iii> if data looks change mid-way

# some functionalities borrowed from "extras/checkMaxSPRTsimulation.R"
library(Sequential)
library(dplyr)
library(ggplot2)

# (0) function to compute Bayesian post probs for Poisson data
BayesianPoisson <- function(S, numLooks = 12,
                            expectedCount = 10,
                            effectSize = 1,
                            priorAlpha = 1, 
                            priorBeta = 1,
                            delta1 = 0.95){
  P1s = rep(0.01,S*numLooks)
  Rounds = rep(1:S, each = numLooks)
  Period = rep(1:numLooks, S)
  Reject = rep(FALSE, S*numLooks)
  
  rejectCounter = 0
  
  for(s in 1:S){
    cat(sprintf('\n\nSimulation round %s:...\n', s))
    
    ## simulate data
    ## across all data looks in one go...
    obsCounts = rpois(numLooks, expectedCount * effectSize)
    
    ## accumulate observations
    cumCounts = cumsum(obsCounts)
    cumExpected = expectedCount * c(1:numLooks)
    
    ## evaluate posterior probs of (rate ratio > 1)
    this.P1s = pgamma(1,
                     shape = priorAlpha + cumCounts, 
                     rate = priorBeta + cumExpected,
                     lower.tail = FALSE)
    
    ## reject?
    this.reject = this.P1s > delta1
    
    ## save results
    inds = seq(from = (s-1)*numLooks + 1,
               to = (s-1)*numLooks + numLooks,
               by = 1)
    
    P1s[inds] = this.P1s
    Reject[inds] = this.reject
    
    if(any(this.reject)){
      cat(sprintf('Null rejected in simulation round %s!\n', s))
      rejectCounter = rejectCounter + 1
    }
  }
  
  res = data.frame(Simulation = Rounds,
                   Period = Period,
                   P1 = P1s,
                   Reject = Reject)
  
  cat(sprintf('Overall rejection rate is %.3f\n', rejectCounter/S))
  
  return(list(res = res,
              numLooks = numLooks,
              expectedCount = expectedCount,
              rejectRate = rejectCounter/S,
              delta1 = delta1))
  
}

# # try Bayesian Poisson sequential test
# bayesSim = BayesianPoisson(S = 500,
#                            numLooks = 24,
#                            expectedCount = 10,
#                            delta1 = 0.98)

# (0)-2 utils function to get rejection rates for Bayesian method----
getBayesRejectionRate <- function(resLs, 
                                  numLooks = 12, # the number of looks to include
                                  delta1 = 0.95){
  res = resLs$res
  
  if(numLooks != resLs$numLooks){
    res = res %>% filter(Period <= numLooks)
  }
  
  if(delta1 != resLs$delta1){
    res$Reject = res$P1 > delta1
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

# ## get Bayesian results rate
# RatesB = getBayesRejectionRate(bayesSim,
#                                numLooks = 24,
#                                delta1 = .98)

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



## (3) actually run a very hypothetical experiment...----
# with effect_size = 1 here ----
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
                            actualExpectedCount = 12,
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
Rates2 = getMaxSprtRejectRate(simulated2, numLooks = 12)

## (iv) the Bayesian thing...
bayesSim = BayesianPoisson(S = 500,
                           numLooks = 24,
                           expectedCount = 10,
                           delta1 = 0.98)
RatesB = getBayesRejectionRate(bayesSim,
                               numLooks = 24,
                               delta1 = .98)


## (4) make plots for this very stupid hypothetical experiment----
pdf('~/Documents/Research/betterResults/plots/maxSPRT-weakness-simulation.pdf',
    height = 4, width = 7.5)

type1colors = wes_palette("GrandBudapest1")[c(1,4,3)]

## (i) everything goes according to plan
allRates = Rates0$rejectByPeriod %>% mutate(type = 'A')

ggplot(allRates, aes(x=Period, y=rejectRate)) +
  geom_hline(yintercept = .05, color = 'gray60', 
             size = 1, linetype=2)+
  geom_vline(xintercept = 12, color = 'brown', alpha = 0.6,
             size = 1.2, linetype=2)+
  geom_line(size = 1.5, color = type1colors[1]) +
  labs(x='Analysis period (by months)', y = 'Type 1 error rate') +
  scale_x_continuous(limits = c(0.5,24), breaks = seq(from=3,to=24,by=3)) +
  scale_y_continuous(limits = c(0,0.25), breaks = c(0,0.05, 0.1, 0.15, 0.2)) +
  scale_color_manual(values = type1colors) +
  theme_bw(base_size = 16)

## (ii) if more data points are observed than planned??
allRates = allRates %>% bind_rows(Rates1$rejectByPeriod %>% mutate(type = 'B'))

ggplot(allRates, aes(x=Period, y=rejectRate, color = type)) +
  geom_hline(yintercept = .05, color = 'gray60', 
             size = 1, linetype=2)+
  geom_vline(xintercept = 12, color = 'brown', alpha = 0.6,
             size = 1.2, linetype=2)+
  geom_line(size = 1.5) +
  labs(x='Analysis period (by months)', y = 'Type 1 error rate') +
  scale_x_continuous(limits = c(0.5,24), breaks = seq(from=3,to=24,by=3)) +
  scale_y_continuous(limits = c(0,0.25), breaks = c(0,0.05, 0.1, 0.15, 0.2)) +
  scale_color_manual(values = type1colors) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none')

  
## (iii) if longer surveillance is needed? 
allRates = allRates %>% bind_rows(Rates2$rejectByPeriod %>% mutate(type = 'C'))
ggplot(allRates, aes(x=Period, y=rejectRate, color = type)) +
  geom_hline(yintercept = .05, color = 'gray60', 
             size = 1, linetype=2)+
  geom_vline(xintercept = 12, color = 'brown', alpha = 0.6,
             size = 1.2, linetype=2)+
  geom_line(size = 1.5) +
  labs(x='Analysis period (by months)', y = 'Type 1 error rate') +
  scale_x_continuous(limits = c(0.5,24), breaks = seq(from=3,to=24,by=3)) +
  scale_y_continuous(limits = c(0,0.25), breaks = c(0,0.05, 0.1, 0.15, 0.2)) +
  scale_color_manual(values = type1colors) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none')

## (iv) adding Bayesian method consistent performance curve

allRates = allRates %>% bind_rows(RatesB$rejectByPeriod %>% mutate(type = 'D'))

type1colors = c(type1colors, wes_palette("Darjeeling2")[2])

ggplot(allRates, aes(x=Period, y=rejectRate, color = type)) +
  geom_hline(yintercept = .05, color = 'gray60', 
             size = 1, linetype=2)+
  geom_vline(xintercept = 12, color = 'brown', alpha = 0.6,
             size = 1.2, linetype=2)+
  geom_line(size = 1.5) +
  labs(x='Analysis period (by months)', y = 'Type 1 error rate') +
  scale_x_continuous(limits = c(0.5,24), breaks = seq(from=3,to=24,by=3)) +
  scale_y_continuous(limits = c(0,0.25), breaks = c(0,0.05, 0.1, 0.15, 0.2)) +
  scale_color_manual(values = type1colors) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none')

dev.off()


# (6) run experiment with effectSize > 1?

effect_size = 1.2
set.seed(42)

## (i) everything goes exactly according to the plan
simulated0 = maxSprtPoisson(S = 500, 
                            numLooks=12, 
                            expectedCountByPlan = 10,
                            actualExpectedCount = 10,
                            effectSize = effect_size,
                            alpha = 0.05)
Rates0 = getMaxSprtRejectRate(simulated0)

## (ii) actual number of data points > anticipated
simulated1 = maxSprtPoisson(S = 500, 
                            numLooks=12, 
                            expectedCountByPlan = 10,
                            actualExpectedCount = 12,
                            effectSize = effect_size,
                            alpha = 0.05)
Rates1 = getMaxSprtRejectRate(simulated1)

## (iii) actually planned out a 1-year surveillance, but had to do 2-year
simulated2 = maxSprtPoisson(S = 500, 
                            numLooks=24, 
                            expectedCountByPlan = 10,
                            actualExpectedCount = 10,
                            effectSize = effect_size,
                            alpha = 0.05)
Rates2 = getMaxSprtRejectRate(simulated2, numLooks = 12)

## (iv) the Bayesian thing...
bayesSim = BayesianPoisson(S = 500,
                           numLooks = 24,
                           expectedCount = 10,
                           effectSize = effect_size,
                           delta1 = 0.975)
RatesB = getBayesRejectionRate(bayesSim,
                               numLooks = 24,
                               delta1 = .975)

# (7) plot the stats power thing in one big plot
allRates = bind_rows(Rates0$rejectByPeriod %>% mutate(type = 'A'),
                     Rates1$rejectByPeriod %>% mutate(type = 'B'),
                     Rates2$rejectByPeriod %>% mutate(type = 'C'),
                     RatesB$rejectByPeriod %>% mutate(type = 'D'))

type1colors = c(wes_palette("GrandBudapest1")[c(1,4,3)], 
                wes_palette("Darjeeling2")[2])

ggplot(allRates, aes(x=Period, y=rejectRate, color = type)) +
  geom_hline(yintercept = .05, color = 'gray60', 
             size = 1, linetype=2)+
  geom_vline(xintercept = 12, color = 'brown', alpha = 0.6,
             size = 1.2, linetype=2)+
  geom_line(size = 1.5) +
  labs(x='Analysis period (by months)', y = 'Statistical Power',
       caption = 'Statistical power for RR > 1 signal detection',
       color = 'Power for:') +
  scale_x_continuous(limits = c(0.5,24), breaks = seq(from=3,to=24,by=3)) +
  scale_y_continuous(limits = c(0,1), 
                     breaks = seq(from = 0, to = 1, by = 0.25)) +
  scale_color_manual(values = type1colors,
                     labels = c('MaxSPRT (oracle)',
                                'MaxSPRT (more data)',
                                'MaxSPRT (longer)',
                                'Bayesian')) +
  theme_bw(base_size = 16)+
  theme(legend.position = c(0.8,0.35))
