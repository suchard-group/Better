# functions to reverse engineer critival values
# with a "time extension factor" either between 0 and 1
#                                 or > 1 integer!!!


## 0. overhead main function for SCCS and HC----
computeCvs_reverse<- function(estimates, cvFunction, ex_factor = 1, maxCores = 4) {
  # note: estimates has to be estimates with same databaseId AND method!!!
  cluster <- ParallelLogger::makeCluster(min(10, maxCores))
  ParallelLogger::clusterRequire(cluster, "dplyr")
  on.exit(ParallelLogger::stopCluster(cluster))
  
  subsets <- split(estimates, 
                   paste(estimates$analysisId, estimates$exposureId, estimates$outcomeId))
  cvs <- ParallelLogger::clusterApply(cluster, subsets, cvFunction, ex_factor = ex_factor)
  cvs <- bind_rows(cvs)
  rm(subsets)
  
  estimatesWithCvs <- estimates %>% select(-criticalValue) %>%
    inner_join(cvs, by = c("analysisId", "exposureId", "outcomeId"))
  
  return(estimatesWithCvs)
}


## 1. SCCS ------
computeSccsCv_reverse <- function(subset, ex_factor=1) {
  
  computeTruncatedBinomialCv <- function(n, z, groupSizes) {
    if (n > 250) {
      groupSizes <- round(groupSizes * 250 / n)
      groupSizes <- groupSizes[groupSizes > 0]
      n <- sum(groupSizes)
    }
    # This check is done inside Sequential::CV.Binomial as well, but will throw an error there:
    pst <- 1/(1 + z)
    if (1 - pbinom(n - 1, n, pst) > 0.05) {
      return(NA)
    }
    cv <- try({
      # Setting time-out to escape infinite loops in Sequential:
      if (is.null(getOption("CvTimeout"))) {
        setTimeLimit(60, Inf)
      } else {
        setTimeLimit(getOption("CvTimeout"), Inf)
      }
      Sequential::CV.Binomial(N = n,
                              M = 1,
                              alpha = 0.05,
                              z = z,
                              GroupSizes = groupSizes)$cv
    }, silent = TRUE) 
    if (inherits(cv, "try-error")) {
      cv <- NA
    }
    return(cv)
  }
  
  # subset <- subsets[[1]]
  # reverse engineer back needed columns! 
  subset = subset %>% 
    mutate(outcomeEvents = counterfactualOutcomes + exposureOutcomes,
           daysObserved = counterfactualDays + exposureDays) %>%
    arrange(periodId)
  
  # pretend the plan is at ex_factor * original length
  num_rows = nrow(subset)
  if(num_rows == 0){
    cv = NA
    return(NULL)
  }
  
  if(ex_factor < 1){
    slice_row = ceiling(num_rows * ex_factor)
    subset = subset %>% slice(1:slice_row)
  }else if(ex_factor > 1){
    if(!is.integer(ex_factor)){
      ex_factor = ceiling(ex_factor)
    }
    subset = subset %>% slice(rep(row_number(),ex_factor))
    
    for(i in 2:ex_factor){
      # add on observed days, exposure days, and outcome events
      subset$daysObserved[((i-1)*num_rows+1):(i*num_rows)] = 
        subset$daysObserved[((i-2)*num_rows+1):((i-1)*num_rows)] + 
        subset$daysObserved[1:num_rows]
      
      subset$exposureDays[((i-1)*num_rows+1):(i*num_rows)] = 
        subset$exposureDays[((i-2)*num_rows+1):((i-1)*num_rows)] + 
        subset$exposureDays[1:num_rows]
      
      subset$outcomeEvents[((i-1)*num_rows+1):(i*num_rows)] = 
        subset$outcomeEvents[((i-2)*num_rows+1):((i-1)*num_rows)] + 
        subset$outcomeEvents[1:num_rows]
    }
    
  }
  
  # then compute CVs as usual
  sampleSizeUpperLimit <- max(subset$outcomeEvents , na.rm = TRUE)
  if (sampleSizeUpperLimit == 0) {
    cv <- NA
  } else {
    events <- subset %>%
      pull(outcomeEvents)
    looks <- length(events)
    if (looks > 1) {
      events[2:looks] <- events[2:looks] - events[1:(looks-1)]
      events <- events[events != 0]
    }
    cv <- computeTruncatedBinomialCv(n = sampleSizeUpperLimit,
                                     z = max(subset$daysObserved - subset$exposureDays) / max(subset$exposureDays),
                                     groupSizes = events)
  }
  return(tibble(analysisId = subset$analysisId[1],
                exposureId = subset$exposureId[1],
                outcomeId = subset$outcomeId[1],
                criticalValue = cv))
}



## 2. Historical Comparator ----
computeHistoricalComparatorCv_reverse <- function(subset, ex_factor = 1) {
  
  
  computeTruncatedPoissonCv <- function(n, groupSizes) {
    if (n > 250) {
      groupSizes <- round(groupSizes * 250 / n)
      groupSizes <- groupSizes[groupSizes > 0]
      n <- sum(groupSizes)
    }
    cv <- try({
      # Setting time-out to escape infinite loops in Sequential:
      if (is.null(getOption("CvTimeout"))) {
        setTimeLimit(60, Inf)
      } else {
        setTimeLimit(getOption("CvTimeout"), Inf)
      }
      cv <- Sequential::CV.Poisson(SampleSize = n,
                                   alpha = 0.05,
                                   M = 1,
                                   GroupSizes = groupSizes)
    }, silent = TRUE) 
    if (inherits(cv, "try-error")) {
      cv <- NA
    }
    return(cv)
  }
  
  # subset <- subsets[[2611]]
  
  subset = subset %>% 
    mutate(expectedOutcomes = counterfactualOutcomes/(counterfactualDays) * 
             (exposureDays)) %>%
    arrange(periodId)
  
  # pretend the plan is at ex_factor * original length
  num_rows = nrow(subset)
  if(num_rows == 0){
    cv = NA
    return(NULL)
  }
  
  if(ex_factor < 1){
    slice_row = ceiling(num_rows * ex_factor)
    subset = subset %>% slice(1:slice_row)
  }else if(ex_factor > 1){
    if(!is.integer(ex_factor)){
      ex_factor = ceiling(ex_factor)
    }
    subset = subset %>% slice(rep(row_number(),ex_factor))
    
    for(i in 2:ex_factor){
      # add on expected Outcomes
      subset$expectedOutcomes[((i-1)*num_rows+1):(i*num_rows)] = 
        subset$expectedOutcomes[((i-2)*num_rows+1):((i-1)*num_rows)] + 
        subset$expectedOutcomes[1:num_rows]
    }
    
  }
  
  # compute CV as usual
  expectedOutcomes <- subset %>%
    pull(expectedOutcomes)
  looks <- length(expectedOutcomes)
  if (looks > 1) {
    expectedOutcomes[2:looks] <- expectedOutcomes[2:looks] - expectedOutcomes[1:(looks-1)]
    # Per-look expected counts < 1 can lead to CV.Poisson() getting stuck in infinite loop, so combining smaller looks:
    eos <- c()
    pending <- 0
    for (eo in expectedOutcomes) {
      if (!is.na(eo)) {
        if (eo + pending >= 1) {
          eos <- c(eos, eo + pending)
          pending <- 0
        } else {
          pending <- eo + pending
        }
      }
    }
    expectedOutcomes <- eos
    sampleSizeUpperLimit <- sum(expectedOutcomes)
  }
  if (sampleSizeUpperLimit == 0) {
    cv <- NA
  } else {
    cv <- computeTruncatedPoissonCv(n = sampleSizeUpperLimit,
                                    groupSizes = expectedOutcomes)
  }
  return(tibble(analysisId = subset$analysisId[1],
                exposureId = subset$exposureId[1],
                outcomeId = subset$outcomeId[1],
                criticalValue = cv))
}


# # try SCCS
# (
# subset_small = estimates_CCAE %>%
#   filter(method == 'SCCS',
#          analysisId == 1,
#          exposureId == 211983,
#          outcomeId == 197032)
# )
# 
# computeSccsCv_reverse(subset_small, ex_factor = 3)

# # try HC
# (
# subset_small = estimates_CCAE %>%
#   filter(method == 'HistoricalComparator',
#          analysisId == 2,
#          exposureId == 211983,
#          outcomeId == 197032)
# )
# computeHistoricalComparatorCv_reverse(subset_small, ex_factor = 1)
# # well, more or less okay!!

# # try overhead function
# sections = estimates_CCAE %>%
#   filter(method == 'HistoricalComparator', analysisId == 2, exposureId == 211983)
# 
# test_res = computeCvs_reverse(sections, 
#                               computeHistoricalComparatorCv_reverse,
#                               ex_factor = 2, maxCores = 4)
