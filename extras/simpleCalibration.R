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
## smaller helper function to compute type 1 error rate from 
## a named list of post samples
## name: negative control outcome id
## for a specific negative control outcome
computeTypeI <- function(samps, h, delta1 = 0.95){
  ## yet another smaller helper function for ONE outcome only
  decideOneOutcome <- function(oneSamp, h, delta1 = 0.95){
    postProbs = lapply(oneSamp, function(x) mean(x > h)) %>% unlist()
    max(postProbs) > delta1
  }
  
  decs = lapply(samps, decideOneOutcome, h = h, delta1 = delta1) %>% unlist()
  mean(decs)
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
  
  # load IPC table to obtain NC ids
  NCs = readRDS(file.path(cachePath, 'allIPCs.rds'))$NEGATIVE_CONTROL_ID %>% 
    unique() %>%
    as.character()
  
  # load samples if not provided as input
  if(is.null(samps)){
    subdir = paste0('samples-',database_id)
    dbname = elseif(database_id %in% c('MDCD','MDCR'),
                    paste0('IBM_',database_id),
                    database_id)
    samps = list()
    ## read files in by period
    for(pid in 1:12){
      fname = sprintf(
        '%s_%s_%s_period%s_analysis%s_samples.rds',
        dbname,
        method,
        exposure_id,
        pid,
        analysis_id
      )
      fpath = file.path(resPath, subdir, fname)
      # read in if file exists
      if(file.exists(fpath)){
        this.samps = readRDS(fpath)
        if(useAdjusted){
          this.samps = this.samps[[prior_id]]$adjustedPostSamps
        }else{
          this.samps = this.samps[[prior_id]]$postSamps
        }
        this.NCs = rownames(this.samps)[rownames(this.samps) %in% NCs]
        #this.samps = this.samps[this.NCs,]
        for(nc in this.NCs){
          # add this period samples to the named list
          # name: NC id
          if(!nc %in% names(samps)){
            samps[[nc]] = list()
          }
          samps[[nc]][[as.character(pid)]] = this.samps[nc,]
        }
      }
    }
  }
  
  # process the post. samples
  if(length(samps) < minOutcomes){
    mes = sprintf('Num. of negative controls for analysis %s, exposure %s and prior %s is smaller than minimum %s!\n',
                  analysis_id, exposure_id, prior_id, minOutcomes)
    cat(mes)
    return()
  }else{
    #samps
    
    ## do binary search 
    st = searchRange[1]
    en = searchRange[2]
    type1 = NULL
    while(TRUE){
      stError = computeTypeI(samps, st, delta1) - alpha
      enError = computeTypeI(samps, en, delta1) - alpha
      
      # if range too small, stop...
      if(en - st < tol){
        cat('Search grid gets too narrow before convergence. Interpret with caution!!\n')
        h = elseif(abs(stError) < abs(enError), st, en)
        type1 = elseif(abs(stError) < abs(enError), stError, enError) + alpha
        break
      }

      # check if within error margin
      if(stError > tol & enError < -tol){
        ## take middle point and continue
        mid = (st+en)/2
        midError = computeTypeI(samps, mid, delta1) - alpha
        if(midError > 0){
          st = mid
        }else{
          en = mid
        }
      }else if(abs(stError - alpha) <= tol){
        ## use `st` as answer
        h = st
        type1 = stError + alpha
        break
      }else if(abs(enError - alpha) <= tol){
        ## use `en` as answer
        h = en
        type1 = enError + alpha
        break
      }else{
        ## otherwise, something must have gone wrong
        mes = sprintf('At st=%.3f, en=%.3f, Type I errors are %.4f and %.4f. Running failed!\n\n',
                      st, en, stError+alpha, enError+alpha)
        cat(mes)
        break
      }
    }
    
    # return result
    return(list(h = h, type1 = type1))
  }
}




## try it
resultspath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/samples-MDCD" 


