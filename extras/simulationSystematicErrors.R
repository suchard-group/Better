# 02/15/2023
# run simple simulations to show residual systematic error
# using Historical Comparator, and SCCS

library(foreach)

## 1. function to simulate an individual trajectory ----
## given weekly background rates, and 

simulateOnePerson <- function(weeklyRates,
                              RR = 2,
                              exposureWeek = 1,
                              atRiskWeeks = 6,
                              totalWeeks = 52){
  
  # check legnth of `weekRates`
  if(length(weeklyRates) != totalWeeks){
    stop(sprintf('Need %s total weeks of reference rates, as suggested by totalWeeks!',
                 totalWeeks))
  }
  
  # if there is an exposure...
  if(!is.null(exposureWeek)){
    rates = weeklyRates
    rates[exposureWeek:(exposureWeek + atRiskWeeks - 1)] = 
      weeklyRates[exposureWeek:(exposureWeek + atRiskWeeks - 1)] * RR
  }else{
    rates = weeklyRates
  }
  
  # simulate number of incidence counts for this individual
  counts = rpois(totalWeeks, rates)
  
  return(counts)
}

# ## try it
# refRate = 1/20
# 
# refWeeklyRates = c(rep(refRate, 13),
#                    rep(refRate * 0.6, 13),
#                    rep(refRate * 1.3, 13),
#                    rep(refRate * 0.8, 13))
# 
# (traj = simulateOnePerson(refWeeklyRates, exposureWeek = 5))
# 
# (traj = simulateOnePerson(refWeeklyRates, exposureWeek = NULL))


## 2. function to simulate a lot of individual trajectories (post-vaccination campaign)----
## with a 1-d binary simple covariate X

# Feb 16 2023 update: allow non-uniform exposure week distributions....
simulateManyPersons <- function(N,
                                weeklyRates,
                                exposureRates = NULL,
                                presentTime = TRUE,
                                xEffect = 1.5,
                                totalWeeks = 52,
                                atRiskWeeks = 6,
                                seed = 42, ...){
  
  set.seed(seed)
  
  ## randomly assign binary X
  X = rbinom(N, 1, 0.5)
  Xfactors = sapply(X, function(xi) ifelse(xi == 1, xEffect, 1))
  
  ## engineer exposureRates to match total number of exposure weeks possible
  if(!is.null(exposureRates)){
    exposureRates = exposureRates[1:(totalWeeks-atRiskWeeks)]
  }
  
  ## randomly assign exposure weeks
  if(presentTime){
    exposureWeeks = sample(totalWeeks-atRiskWeeks, N, 
                           replace = TRUE,
                           prob = exposureRates)
  }else{
    exposureWeeks = NULL
  }
  
  
  ## simulate all individual trajectories
  allCounts = foreach(i = 1:N, .combine = cbind) %do% {
    this.weeklyRates = weeklyRates * Xfactors[i]
    if(!is.null(exposureWeeks)){
      this.exposureWeek = exposureWeeks[i]
    }else{
      this.exposureWeek = NULL
    }
    simulateOnePerson(weeklyRates = this.weeklyRates,
                      exposureWeek = this.exposureWeek,
                      totalWeeks = totalWeeks, ...)
  }
  
  return(list(incidenceCounts = allCounts,
              exposureWeeks = exposureWeeks,
              atRiskWeeks = atRiskWeeks,
              X = X,
              currentData = presentTime))
}

## try it-----

refRate = 1/(52*3) # annual background rate ~ 1 in 3 persons

refWeeklyRates = c(rep(refRate, 13),
                   rep(refRate * 0.6, 13),
                   rep(refRate * 1.3, 13),
                   rep(refRate * 0.8, 13))

presentTrajs = simulateManyPersons(N = 1000, 
                                   xEffect = 1,
                                   weeklyRates = refWeeklyRates)
historicTrajs = simulateManyPersons(N = 1000, 
                                    xEffect = 1,
                                    weeklyRates = refWeeklyRates * 0.2,
                                    presentTime = FALSE)


# 3. Data prep functions ----
## (1)a prepare data for HC ----
prepDataHC <- function(trajs){
  exposureWeeks = trajs$exposureWeeks
  atRiskWeeks = trajs$atRiskWeeks
  
  N = length(trajs$X)
  nw = nrow(trajs$incidenceCounts)
  
  if(!is.null(exposureWeeks)){
    exposed = sapply(1:N, function(i){
      expoVec = rep(0, nw); 
      expoVec[exposureWeeks[i]:(exposureWeeks[i]+atRiskWeeks-1)] = 1;
      expoVec
    })
  }else{
    exposed = 0
  }
  
  
  longData = data.frame(week = rep(1:nw, N),
                        x = rep(trajs$X, each = nw),
                        exposure = c(exposed),
                        count = c(trajs$incidenceCounts),
                        indivId = rep(1:N, each = nw))
  
  return(longData)
}

## (1)b further single out cases only for SCCS----
caseOnly <- function(longData){
  caseIds = longData %>% group_by(indivId) %>%
    summarise(totalCases = sum(count)) %>%
    ungroup() %>%
    filter(totalCases > 0) %>%
    select(indivId) %>% pull()
  
  return(longData %>% filter(indivId %in% caseIds))
}

## try this...
longDataHC = prepDataHC(presentTrajs)
longDataSCCS = caseOnly(longDataHC)


# 4. HC estimation, static ----
runHistoricalComparator <- function(presentData, 
                                    historyData,
                                    seasonWeeks = 4){
  presentData = presentData %>% 
    mutate(season = as.factor((week-1) %/% seasonWeeks + 1)) %>%
    mutate(phase = 'present')
  
  historyData = historyData %>% 
    mutate(season = as.factor((week-1) %/% seasonWeeks + 1)) %>%
    mutate(phase = 'history')
  
  allData = bind_rows(presentData, historyData) %>%
    # mutate(chunk = as.factor(week %/% resolutionWeeks + 1)) %>%
    group_by(x, season, exposure) %>%
    summarise(total_count = sum(count),
              total_time = n()) %>%
    ungroup()
  
  mHC = glm(total_count ~ exposure + x + season - 1 + offset(log(total_time)), 
            family = poisson(),
            data = allData)
  
  return(list(model = mHC,
              data = allData))
}

## try this...
longDataHC = prepDataHC(presentTrajs)
historyLongDataHC = prepDataHC(historicTrajs)

HCresults = runHistoricalComparator(longDataHC, 
                              historyLongDataHC)
summary(HCresults$model)

# 5. SCCS estimation, static ----
runSCCS <- function(presentData, filterCases = TRUE){
  if(filterCases){
    presentData = caseOnly(presentData)
  }
  
  allData = presentData %>% 
    mutate(indivId = as.factor(indivId)) %>%
    group_by(indivId, x, exposure) %>%
    summarize(total_count = sum(count),
              total_time = n()) %>%
    ungroup() 
  
  mSCCS = glm(total_count ~ exposure + indivId - 1 + offset(log(total_time)), 
              family = poisson(),
              data = allData)
  
  return(list(model = mSCCS,
              data = allData))
}

## try this ----
SCCSresults = runSCCS(longDataHC, filterCases = TRUE)
summary(SCCSresults$model)

# 6. sequential analysis ----
# (for now: data accrual order by order of vaccination exposure for the individuals)
batchData <- function(longData, numBatches = 13){
  batchInterval = max(longData$week)/numBatches
  
  firstExposureWeeks = longData %>% 
    group_by(indivId) %>%
    filter(exposure==1) %>%
    summarize(firstWeek = min(week)) %>%
    ungroup()
  
  batches = firstExposureWeeks %>%
    mutate(batch = if_else(firstWeek/batchInterval < numBatches,
                           ceiling(firstWeek/batchInterval),
                           numBatches))
  
  longData = longData %>% left_join(batches, by = 'indivId')
  
  return(longData)
}

## try it... 
longDataHC_batches = batchData(longDataHC)

sequentialAnalysis <- function(analysisFunc, 
                               presentData, 
                               numBatches = 13, ...){
  
  presentData = batchData(presentData, numBatches = numBatches)
  
  allEstimates = NULL
  for(b in 1:max(presentData$batch)){
    thisBatchData = presentData %>% filter(batch <= b)
    
    this.model = analysisFunc(thisBatchData, ...)
    
    # save "exposure" estimates and Wald-like 95% CIs
    this.RR.estimate = exp(coef(this.model$model)['exposure'])
    this.RR.confint = tryCatch(
      exp(confint(this.model$model, 'exposure')),
      error = function(e){
        sprintf('Profiling does not work for batch %s! Using Wald-like confidence intervals instead...\n',
                b)
        sd = summary(this.model$model)['exposure','Std. Error']
        return(exp(c(this.RR.estimate + qnorm(.025) * sd, 
                     this.RR.estimate + qnorm(.975) * sd)))
      }
    )
    
    this.res = c(b, this.RR.estimate, this.RR.confint)
    allEstimates = rbind(allEstimates, this.res)
  }
  
  allEstimates = as.data.frame(allEstimates)
  names(allEstimates) = c('batch', 'estimate', 'LB95', 'UB95')
  
  return(list(RRestimates = allEstimates,
              presentData = presentData))
}

## try it 
seqResHC = sequentialAnalysis(analysisFunc = runHistoricalComparator,
                              presentData = longDataHC,
                              historyData = historyLongDataHC,
                              seasonWeeks = 4,
                              numBatches = 13)
seqResHC$RRestimates

longDataSCCS = caseOnly(longDataHC)
seqResSCCS = sequentialAnalysis(analysisFunc = runSCCS,
                                presentData = longDataSCCS,
                                filterCases = FALSE)
seqResSCCS$RRestimates


#### CODE TO RUN ####
#### try everything all at once ----

# annual background rate: ~ 1 in 3 persons
refRate = 1/(52*3) 

# seasonality rate factors
seasonFactors = c(1, 0.6, 1.3, 0.8)

# historical differential rate factor
historyRate = 0.3

# total number of people to simulate
N.pop = 5000

# rate factor for the binary covariate X
xeffect = 1.5

# random seed
theSeed = 42

# construct weekly reference rates with seasonality
refWeeklyRates = c(rep(refRate * seasonFactors[1], 13),
                   rep(refRate * seasonFactors[2], 13),
                   rep(refRate * seasonFactors[3], 13),
                   rep(refRate * seasonFactors[4], 13))

# Feb 16 2023: build non-uniform weekly exposure rates as well
exposureRates = c(rep(seasonFactors[1], 13),
                  rep(seasonFactors[2], 13),
                  rep(seasonFactors[3], 13),
                  rep(seasonFactors[4], 13))

# present-time and historical trajectories... ----
presentTrajs = simulateManyPersons(N = N.pop, 
                                   xEffect = xeffect,
                                   weeklyRates = refWeeklyRates,
                                   exposureRates = exposureRates,
                                   seed = theSeed)
historicTrajs = simulateManyPersons(N = N.pop, 
                                    xEffect = xeffect,
                                    weeklyRates = refWeeklyRates * historyRate,
                                    exposureRates = NULL,
                                    presentTime = FALSE,
                                    seed = theSeed)

# reformat data to long format
# also extract cases only for SCCS
longDataHC = prepDataHC(presentTrajs)
historyLongDataHC = prepDataHC(historicTrajs)

longDataSCCS = caseOnly(longDataHC)

# run sequential analysis----
## HC
seqResHC = sequentialAnalysis(analysisFunc = runHistoricalComparator,
                              presentData = longDataHC,
                              historyData = historyLongDataHC,
                              seasonWeeks = 4,
                              numBatches = 13)
seqResHC$RRestimates

## SCCS 
seqResSCCS = sequentialAnalysis(analysisFunc = runSCCS,
                                presentData = longDataSCCS,
                                filterCases = FALSE)
seqResSCCS$RRestimates
