# March 25
# code to do Bayesian "calibration"

library(tidyverse)

# 1. with fixed threshold (between H0 and H1), find delta1 threshold to achieve Type I error rate
# 2. with fixed delta1 threshold, find threshold (between H0 and H1) to achieve Type I error rate

# function 1 -----
# using default threshold 0 (testing if beta > 0)
calibrateByDelta1 <- function(database_id,
                              method,
                              analysis_id,
                              exposure_id,
                              prior_id,
                              resPath,
                              summ = NULL,
                              alpha = 0.05,
                              minOutcomes = 5,
                              useAdjusted = FALSE){
  # if summ is provided, directly query from given dataframe
  # otherwise, load it from saved summary file
  if(is.null(summ)){
    fname = sprintf('AllSummary-%s-%s.rds', database_id, method)
    summ = readRDS(file.path(resPath, fname))
  }
  
  # get relevant rows for NCs
  dat  = summ %>% 
    filter(analysis_id == !!analysis_id, 
           exposure_id == !!exposure_id,
           prior_id == !! prior_id,
           negativeControl == TRUE)
  
  # check if total num. of NCs meet `minOutcomes`
  if(nrow(dat) == 0 || length(unique(dat$outcome_id)) < minOutcomes){
    mes = sprintf('Num. of negative controls for analysis %s, exposure %s and prior %s is smaller than minimum %s!\n',
                  analysis_id, exposure_id, prior_id, minOutcomes)
    cat(mes)
    return(NULL)
  }else{
    if(useAdjusted){
      p1s = dat %>% group_by(outcome_id) %>%
        summarize(maxP1 = max(adjustedP1))
    }else{
      p1s = dat %>% group_by(outcome_id) %>%
        summarize(maxP1 = max(P1))
    }
    p1s = p1s %>% ungroup() %>% select(maxP1) %>% pull()
    calibratedThres = quantile(p1s, 1-alpha) %>% as.numeric()
    
    res = data.frame(database_id, method = method, analysis_id = analysis_id,
                     exposure_id = exposure_id, prior_id = prior_id, 
                     calibratedDelta1 = calibratedThres, alpha = alpha,
                     adjusted = useAdjusted)
    return(res)
  }
}

## try it ------
# resultspath = '~/Documents/Research/betterResults/summary'
# calibrateByDelta1(database_id = 'CCAE',
#                   method = 'SCCS',
#                   analysis_id = 1,
#                   exposure_id = 21184,
#                   prior_id = 2,
#                   resPath = resultspath,
#                   useAdjusted = TRUE)


# function 2 -------
## smaller helper function to compute type 1 error rate from a matrix of post samples
computeTypeI <- function(samples, h, delta1 = 0.95){
  p1s = apply(samples, 1, function(x) mean(x > h))
  mean(p1s > delta1)
}
  
## use binary search to determine a proper threshold to achieve Type I error rate
# technically, should always set `useAdjusted = FALSE`!!!
calibrateNull <- function(database_id,
                          method,
                          analysis_id,
                          exposure_id,
                          prior_id,
                          resPath,
                          cachPath,
                          searchRange = c(-2,2),
                          samps = NULL,
                          delta1 = 0.95,
                          alpha = 0.05,
                          tol = 0.002,
                          minOutcomes = 5,
                          useAdjusted = FALSE){
  
  # load samples if not provided as input
  if(is.null(samps)){
    subdir = paste0('samples-',database_id)
    dbname = elseif(database_id %in% c('MDCD','MDCR'),
                    paste0('IBM_',database_id),
                    database_id)
    fnamePattern = sprintf(
      '%s_%s_%s_period[1-9]*_analysis%s_samples.rds',
      dbname,
      method,
      exposure_id,
      analysis_id
    )
    filesToRead = list.files(path = file.path(resPath, subdir), 
                             pattern = fnamePattern)
      
  }
}




## try it
resultspath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/samples-MDCD" 


