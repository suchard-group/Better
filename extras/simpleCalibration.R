# March 25
# code to do Bayesian "calibration"

library(tidyverse)

# 1. with fixed threshold (between H0 and H1), find delta1 threshold to achieve Type I error rate
# 2. with fixed delta1 threshold, find threshold (between H0 and H1) to achieve Type I error rate

# function 1 -----
# using default threshold 0 (testing if beta > 0)
# March 30: update with Type 2 error evaluation too
#           while also reporting the actual Type 1 error...


## smaller helper func: compute F1 score 
## (a classification metric -- assumes EQUAL importance of precision and recall)
## (supposedly, higher --> better)
computeF1 <- function(type1, type2){
  if(is.na(type2)){
    f1 = NA
  }else{
    f1 = 2*((1-type1) * (1-type2))/(1-type1 + 1-type2)
  }
  f1
}


# main function for delta1 "calibration"
calibrateByDelta1 <- function(database_id,
                              method,
                              analysis_id,
                              exposure_id,
                              prior_id,
                              resPath,
                              summ = NULL,
                              alpha = 0.05,
                              minOutcomes = 5,
                              useAdjusted = FALSE,
                              evalType2 = TRUE){
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
    
    # evaluate the actual type1 error rate using calibrated threshold
    type1 = mean(p1s > calibratedThres)
    untype1 = mean(p1s > (1-alpha)) # an "uncalibrated" version
    
    # evaluate type 2 error as well if...
    if(evalType2){
      # get relevant rows for PCs first
      pc.dat = summ %>% 
        filter(analysis_id == !!analysis_id, 
               exposure_id == !!exposure_id,
               prior_id == !! prior_id,
               negativeControl == FALSE)
      if(useAdjusted){
        p1s = dat %>% group_by(outcome_id) %>%
          summarize(maxP1 = max(adjustedP1))
      }else{
        p1s = dat %>% group_by(outcome_id) %>%
          summarize(maxP1 = max(P1))
      }
      # evaluate type 2 error rates
      type2 = 1 - mean(p1s > calibratedThres)
      untype2 = 1 - mean(p1s > (1-alpha))
    }else{
      type2 = NA
      untype2 = NA
    }
    
    # compute F1 score (because why not...)
    f1 = computeF1(type1, type2)
    unf1 = computeF1(untype1, untype2)
    
    # return result as a one-row data frame to allow batch run...
    res = data.frame(database_id = database_id, method = method, analysis_id = analysis_id,
                     exposure_id = exposure_id, prior_id = prior_id, 
                     calibratedDelta1 = calibratedThres, alpha = alpha,
                     type1 = type1, type2 = type2,
                     uncalibratedType1 = untype1, uncalibratedType2 = untype2,
                     f1 = f1, uncalibratedF1 = unf1,
                     adjusted = useAdjusted)
    return(res)
  }
}

## try it ------
resultspath = '~/Documents/Research/betterResults/summary'
calibrateByDelta1(database_id = 'CCAE',
                  method = 'SCCS',
                  analysis_id = 1,
                  exposure_id = 21184,
                  prior_id = 2,
                  resPath = resultspath,
                  useAdjusted = TRUE,
                  evalType2 = TRUE)



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

## March 30 update: add function to eval Type 2 error rate from calibrated threshold too!
## a smaller helper function for Type 2 error
computeType2 <- function(samps, h, delta1 = 0.95){
  ## yet another smaller helper function for ONE outcome only
  decideOneOutcome <- function(oneSamp, h, delta1 = 0.95){
    postProbs = lapply(oneSamp, function(x) mean(x > h)) %>% unlist()
    max(postProbs) > delta1
  }
  
  decs = lapply(samps, decideOneOutcome, h = h, delta1 = delta1) %>% unlist()
  # return (1 - rateOfSignal) as Type 2 error rate
  1 - mean(decs)
}

  
## use binary search to determine a proper threshold to achieve Type I error rate
# technically, should always set `useAdjusted = FALSE`!!!
## March 30 update: (optional) eval Type 2 error rate from calibrated threshold too!
calibrateNull <- function(database_id,
                          method,
                          analysis_id,
                          exposure_id,
                          prior_id,
                          resPath,
                          cachePath,
                          searchRange = c(-2,2),
                          samps = NULL,
                          delta1 = 0.95,
                          alpha = 0.05,
                          tol = 0.002,
                          minOutcomes = 5,
                          useAdjusted = FALSE,
                          evalType2 = TRUE){
  
  # load IPC table to obtain NC ids
  NCs = readRDS(file.path(cachePath, 'allIPCs.rds'))$NEGATIVE_CONTROL_ID %>% 
    unique() %>%
    as.character()
  
  # load samples if not provided as input
  if(is.null(samps)){
    subdir = paste0('samples-',database_id)
    dbname = ifelse(database_id %in% c('MDCD','MDCR'),
                    paste0('IBM_',database_id),
                    database_id)
    samps = list()
    if(evalType2){
      pcSamps = list()
    }
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
      #cat(fpath)
      
      # read in if file exists
      if(file.exists(fpath)){
        cat(sprintf('Reading file at path %s.....\n', fpath))
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
        
        ## save positive controls if...
        if(evalType2){
          this.PCs = rownames(this.samps)[!rownames(this.samps) %in% NCs]
          for(pc in this.PCs){
            # add this period samples to the named list for positive controls
            # name: PC id
            if(!pc %in% names(pcSamps)){
              pcSamps[[pc]] = list()
            }
            pcSamps[[pc]][[as.character(pid)]] = this.samps[pc,]
          }
        }
      }
    }
  }
  
  # process the post. samples
  if(length(samps) < minOutcomes){
    mes = sprintf('Num. of negative controls for analysis %s, exposure %s and prior %s is smaller than minimum %s!\n',
                  analysis_id, exposure_id, prior_id, minOutcomes)
    cat(mes)
    return(NULL)
  }else{

    ## do binary search for desired threshold
    st = searchRange[1]
    en = searchRange[2]
    type1 = NULL
    while(TRUE){
      stError = computeTypeI(samps, st, delta1) - alpha
      enError = computeTypeI(samps, en, delta1) - alpha
      
      # if range too small, stop...
      if(en - st < tol){
        mes = sprintf('Search grid range [%.4f, %.4f] gets too narrow before convergence. Interpret with caution!!\n',
                      st, en)
        cat(mes)
        h = ifelse(abs(stError) < abs(enError), st, en)
        type1 = ifelse(abs(stError) < abs(enError), stError, enError) + alpha
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
      }else if(abs(stError) <= tol){
        ## use `st` as answer
        h = st
        type1 = stError + alpha
        break
      }else if(abs(enError) <= tol){
        ## use `en` as answer
        h = en
        type1 = enError + alpha
        break
      }else{
        ## otherwise, something must have gone wrong
        mes = sprintf('At st=%.3f, en=%.3f, Type I errors are %.4f and %.4f. Running failed!\n\n',
                      st, en, stError+alpha, enError+alpha)
        cat(mes)
        h = ifelse(abs(stError) < abs(enError), st, en)
        type1 = ifelse(abs(stError) < abs(enError), stError, enError) + alpha
        break
      }
    }
    
    # report search result
    mes = sprintf('\nSearch finished with h=%.3f and Type I error = %.4f.\n', h, type1)
    cat(mes)
    
    # evaluate type2 error rate with the found threshold if...
    if(evalType2){
      type2 = computeType2(pcSamps, h, delta1)
      mes = sprintf('Type II error after calibration = %.4f\n', type2)
      cat(mes)
    }else{
      # placeholder if type2 not evaluated
      type2 = NA
    }
    
    # obtain uncalibrated error rates as well
    untype1 = computeTypeI(samps, 0, delta1)
    if(evalType2){
      untype2 = computeType2(pcSamps, 0, delta1)
      mes = sprintf('\nWithout calibration, Type I error = %.4f, Type II error = %.4f.\n', 
                    untype1, untype2)
    }else{
      untype2 = NA
      mes = sprintf('\nWithout calibration, Type I error = %.4f.\n', 
                    untype1)
    }
    cat(mes)
    
    # get F1 score too
    if(evalType2){
      f1 = computeF1(type1, type2)
      unf1 = computeF1(untype1, untype2)
      mes = sprintf('With calibration, F1 score = %.4f;\nwithout calibration, F1 score = %.4f.\n',
                    f1, unf1)
      cat(mes)
    }else{
      f1 = NA; unf1 = NA
    }
    
    # return result with a summary data frame
    summ = data.frame(database_id = database_id, 
                      method = method, analysis_id = analysis_id,
                      exposure_id = exposure_id, prior_id = prior_id, 
                      threshold = h, delta1 = delta1, alpha = alpha,
                      type1 = type1, type2 = type2,
                      uncalibratedType1 = untype1, uncalibratedType2 = untype2,
                      f1 = f1, uncalibratedf1 = unf1,
                      adjusted = useAdjusted)
    
    return(list(ncSamps = samps, h = h, type1 = type1, type2 = type2,
                summary = summ))
  }
}




## try it
resultspath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/" 
cachepath = './localCache/'
nullRes = 
  calibrateNull(database_id = 'MDCD',
              method = 'SCCS',
              analysis_id = 15,
              exposure_id = 211983,
              prior_id = 1,
              resPath = resultspath,
              cachePath = cachepath,
              tol = 0.004,
              evalType2 = TRUE)


