# March 25
# code to do Bayesian "calibration"

suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(ggplot2))
library(wesanderson)
library(stringr)

# 1. with fixed threshold (between H0 and H1), find delta1 threshold to achieve Type I error rate
# 2. with fixed delta1 threshold, find threshold (between H0 and H1) to achieve Type I error rate

# function 1 -----
# using default threshold 0 (testing if beta > 0)
# March 30: update with Type 2 error evaluation too
#           while also reporting the actual Type 1 error...


## smaller helper func: compute F1 score 
## (a classification metric -- assumes EQUAL importance of precision and recall)
## (supposedly, higher --> better)
## fix F1 score bug: need number of NCs and IPCs
computeF1 <- function(type1, type2, nn, np){
  if(is.na(type2)){
    f1 = NA
  }else{
    R = 1-type2
    P = (1-type2) * np/((1-type2) * np + type1 * nn)
    f1 = 2*(R*P)/(R + P)
    # f1 = 2*((1-type1) * (1-type2))/(1-type1 + 1-type2)
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
    db_name = ifelse(database_id == 'MDCD', 'IBM_MDCD', database_id)
    fname = sprintf('AllSummary-%s-%s.rds', db_name, method)
    # cat(file.path(resPath, fname))
    # cat('\n')
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
        p1s = pc.dat %>% group_by(outcome_id) %>%
          summarize(maxP1 = max(adjustedP1))
      }else{
        p1s = pc.dat %>% group_by(outcome_id) %>%
          summarize(maxP1 = max(P1))
      }
      p1s = p1s %>% ungroup() %>% select(maxP1) %>% pull()
      # evaluate type 2 error rates
      type2 = 1 - mean(p1s > calibratedThres)
      untype2 = 1 - mean(p1s > (1-alpha))
    }else{
      type2 = NA
      untype2 = NA
    }
    
    # compute F1 score (because why not...)
    if(evalType2){
      np = length(unique(pc.dat$outcome_id))
      nn = length(unique(dat$outcome_id))
      f1 = computeF1(type1, type2, nn, np)
      unf1 = computeF1(untype1, untype2, nn, np)
    }else{
      f1 = NA
      unf1 = NA
    }
    
    # return result as a one-row data frame to allow batch run...
    res = data.frame(database_id = database_id, method = method, analysis_id = analysis_id,
                     exposure_id = exposure_id, prior_id = prior_id, 
                     calibratedDelta1 = calibratedThres, alpha = alpha,
                     type1 = type1, type2 = type2,
                     uncalibratedType1 = untype1, uncalibratedType2 = untype2,
                     f1 = f1, uncalibratedf1 = unf1,
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
#                   useAdjusted = TRUE,
#                   evalType2 = TRUE)


# 04/27/2022: calibrate on delta1 over time periods and see how thresholds and errors change over time
# main function for delta1 "calibration"
tempCalibrateByDelta1 <- function(database_id,
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
    db_name = ifelse(database_id == 'MDCD', 'IBM_MDCD', database_id)
    fname = sprintf('AllSummary-%s-%s.rds', db_name, method)
    # cat(file.path(resPath, fname))
    # cat('\n')
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
  }
  
  # proceed if having enough NCs to work with
  
  ## helper function to process data up to a period
  calibrateUpToPeriod <- function(period){
    # get relevant data up to period
    dat.p = summ %>% 
      filter(analysis_id == !!analysis_id, 
             exposure_id == !!exposure_id,
             prior_id == !! prior_id,
             negativeControl == TRUE,
             period_id <= period)
    # check if having enough data
    if(nrow(dat.p) == 0 || length(unique(dat.p$outcome_id)) < minOutcomes){
      mes = sprintf('Num. of negative controls for period %s, is smaller than minimum %s!\n',
                    period, minOutcomes)
      cat(mes)
      return(NULL)
    }
    
    # proceed if okay
    if(useAdjusted){
      p1s = dat.p %>% group_by(outcome_id) %>%
        summarize(maxP1 = max(adjustedP1))
    }else{
      p1s = dat.p %>% group_by(outcome_id) %>%
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
      pc.dat.p = summ %>% 
        filter(analysis_id == !!analysis_id, 
               exposure_id == !!exposure_id,
               prior_id == !! prior_id,
               negativeControl == FALSE,
               period_id <= period)
      if(useAdjusted){
        p1s = pc.dat.p %>% group_by(outcome_id) %>%
          summarize(maxP1 = max(adjustedP1))
      }else{
        p1s = pc.dat.p %>% group_by(outcome_id) %>%
          summarize(maxP1 = max(P1))
      }
      p1s = p1s %>% ungroup() %>% select(maxP1) %>% pull()
      # evaluate type 2 error rates
      type2 = 1 - mean(p1s > calibratedThres)
      untype2 = 1 - mean(p1s > (1-alpha))
    }else{
      type2 = NA
      untype2 = NA
    }
    
    # no longer returning F1 scores....
    
    # return result as a one-row data frame to allow batch run...
    res = data.frame(database_id = database_id, method = method, analysis_id = analysis_id,
                     exposure_id = exposure_id, prior_id = prior_id, 
                     calibratedDelta1 = calibratedThres, alpha = alpha,
                     type1 = type1, type2 = type2,
                     uncalibratedType1 = untype1, uncalibratedType2 = untype2,
                     adjusted = useAdjusted,
                     period_id = period)
    return(res)
  }
  
  # calibration by periods
  periods = sort(unique(dat$period_id))
  
  allRes = NULL
  for(per in periods){
    res.p = calibrateUpToPeriod(per)
    allRes = rbind(allRes, res.p)
  }
  
  return(allRes)
  
}

# ## try it
# resultspath = '~/Documents/Research/betterResults/summary'
# resByPeriods1 = tempCalibrateByDelta1(database_id = 'CCAE',
#                                      method = 'SCCS',
#                                      analysis_id = 1,
#                                      exposure_id = 21184,
#                                      prior_id = 2,
#                                      resPath = resultspath,
#                                      useAdjusted = TRUE,
#                                      evalType2 = TRUE)
# resByPeriods2 = tempCalibrateByDelta1(database_id = 'CCAE',
#                                       method = 'HistoricalComparator',
#                                       analysis_id = 1,
#                                       exposure_id = 21184,
#                                       prior_id = 2,
#                                       resPath = resultspath,
#                                       useAdjusted = TRUE,
#                                       evalType2 = TRUE)


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
      nn = length(samps); np = length(pcSamps)
      f1 = computeF1(type1, type2, nn, np)
      unf1 = computeF1(untype1, untype2, nn, np)
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
# resultspath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/" 
# cachepath = './localCache/'
# nullRes = 
#   calibrateNull(database_id = 'MDCD',
#               method = 'SCCS',
#               analysis_id = 15,
#               exposure_id = 211983,
#               prior_id = 1,
#               resPath = resultspath,
#               cachePath = cachepath,
#               tol = 0.004,
#               evalType2 = TRUE)



# 3. plot calibration (with one exampple) ------------
# (1) compare type1 and type2 error rates (calibrate hypothesis or decision threshold)
# (2) compare F1 scores
## add option to only return data frames used for plotting, but not show plots
plotCalibration <- function(database_id,
                            method,
                            analysis_id,
                            exposure_id,
                            prior_id,
                            summaryPath,
                            samplePath,
                            cachePath,
                            searchRange = c(-2,2),
                            samps = NULL,
                            delta1 = 0.95,
                            alpha = 0.05,
                            tol = 0.002,
                            minOutcomes = 5,
                            useAdjusted = list(delta1 = TRUE, null=FALSE),
                            showPlots = TRUE){
  # calibrate Delta1
  caliDelta = #invisible(
    calibrateByDelta1(database_id, method, analysis_id, exposure_id, prior_id,
                                resPath = summaryPath,
                                alpha = alpha, minOutcomes = minOutcomes, 
                                useAdjusted = useAdjusted$delta1, evalType2 = TRUE)
  #)
  
  # calibrate null threshold
  caliThres = #invisible(
    calibrateNull(database_id, method, analysis_id, exposure_id, prior_id,
                            resPath = samplePath, cachePath = cachePath, 
                            searchRange = searchRange, 
                            delta1 = delta1,alpha = alpha, tol = tol, 
                            minOutcomes = minOutcomes, useAdjusted = useAdjusted$null, 
                            evalType2 = TRUE)
  #)
  summ = caliThres$summary
  h = caliThres$h
  
  # generate figure caption (with info)
  analysis_name = readRDS(file.path(cachePath,'analyses.rds')) %>%
    filter(method == !!method, analysis_id == !!analysis_id) %>%
    select(description, time_at_risk) %>%
    mutate(description_text = sprintf('%s, time at risk = %s days', description, time_at_risk)) %>%
    select(description_text) %>% pull() %>% as.character()
  
  prior_name = readRDS(file.path(cachePath,'priorTable.rds')) %>%
    filter(prior_id == !!prior_id) %>%
    mutate(priorLabel = sprintf('Mean=%s, SD=%.1f', Mean, Sd)) %>%
    select(priorLabel) %>% pull() %>% as.character()
  
  exposure_name = readRDS(file.path(cachePath,'exposures.rds')) %>%
    filter(exposure_id == !!exposure_id) %>%
    select(exposure_name) %>% pull() %>% as.character()
   
  capt = sprintf('%s\nDatabase: %s\nExposure: %s\nPrior:%s\nalpha=%.2f',
                 analysis_name, database_id, exposure_name, prior_name, alpha)
  
  # 1. Type 1 and Type 2 errors
  # re-arrange data frame
  h.panel = sprintf('Hypothesis threshold: h=%.3f', h)
  delta.panel = sprintf('Decision threshold: delta1=%.3f', caliDelta$calibratedDelta1)
  errorDat = data.frame(error = c(caliDelta$type1, caliDelta$uncalibratedType1,
                                  caliDelta$type2, caliDelta$uncalibratedType2,
                                  summ$type1, summ$uncalibratedType1,
                                  summ$type2, summ$uncalibratedType2),
                        errorType = rep(c('Type 1','Type 1', 
                                          'Type 2', 'Type 2'), 2),
                        calibrate = rep(c('Calibrated','Uncalibrated'), 4),
                        method = rep(c(delta.panel,h.panel), each = 4))
  
  # make plot
  p1 = ggplot(data=errorDat, aes(x=errorType, y=error, fill=calibrate)) +
    geom_bar(stat='identity', position = position_dodge()) +
    geom_hline(yintercept = 0.05, color = 'gray60', 
               size = 1, linetype=2)+
    scale_y_continuous(limits = c(0,1))+
    labs(x='', y='error rate', caption = capt, fill='')+
    facet_grid(.~method) +
    scale_fill_manual(values = wes_palette("Darjeeling2")[c(2,4)]) +
    theme_bw(base_size = 13)
  
  if(showPlots){
    print(p1)
  }
  
  
  # 2. F1 scores
  # re-arrange to data frame
  f1Dat = data.frame(f1 = c(caliDelta$f1, caliDelta$uncalibratedf1,
                               summ$f1, summ$uncalibratedf1),
                     method = rep(c('Decision', 'Hypothesis'), each = 2),
                     calibrate = rep(c('Calibrated','Uncalibrated'), 2))
  
  # plot
  p2 = ggplot(data=f1Dat, aes(x=method, y=f1, fill=calibrate)) +
    geom_bar(stat='identity', position = position_dodge()) +
    # geom_hline(yintercept = 0.05, color = 'gray60', 
    #            size = 1, linetype=2)+
    scale_y_continuous(limits = c(0,1))+
    labs(x='Calibrate on...', y='F1 score', caption = capt, fill='')+
    scale_fill_manual(values = wes_palette("Darjeeling2")[c(2,4)]) +
    theme_bw(base_size = 13)
  
  if(showPlots){
    print(p2)
  }
  
  # return summary results
  return(list(caliDelta = caliDelta,
              caliThresSumm = summ,
              errorDat = errorDat,
              f1Dat = f1Dat,
              caption = capt))
    
  
}

## test it ------
# summarypath = '~/Documents/Research/betterResults/summary'
# samplepath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/"
# cachepath = './localCache/'
# plotCalibration(database_id = 'MDCD',
#                 method = 'HistoricalComparator', # 'SCCS'
#                 analysis_id = 2,
#                 exposure_id = 211983,
#                 prior_id = 1,
#                 summaryPath = summarypath,
#                 samplePath = samplepath,
#                 cachePath = cachepath,
#                 tol = 0.004,
#                 useAdjusted = list(delta1 = TRUE, null=TRUE))




# 4. additional function to present two sets of plots side-by-side (left: un-adjusted; right: adjusted)
plotSideBySide <- function(adjLst, unadjLst, plotToShow = 'both'){
  adjErrorDat = adjLst$errorDat %>% mutate(adjLabel = 'With bias adjustment')
  adjF1Dat = adjLst$f1Dat %>% mutate(adjLabel = 'With bias adjustment')
  
  unadjErrorDat = unadjLst$errorDat %>% mutate(adjLabel = 'Without bias adjustment')
  unadjF1Dat = unadjLst$f1Dat %>% mutate(adjLabel = 'Without bias adjustment')
  
  errorDat = rbind(adjErrorDat, unadjErrorDat)
  f1Dat = rbind(adjF1Dat, unadjF1Dat)
  
  # try to order the categories for the Type 1 and 2 errors
  errorDat$calibrate = factor(errorDat$calibrate, levels = c('Uncalibrated', 'Calibrated'))
  errorDat$adjLabel = factor(errorDat$adjLabel, 
                             levels = c('Without bias adjustment', 'With bias adjustment'))
  
  # add commentary on un-calibrated h and delta1 values
  errorDat = errorDat %>% 
    mutate(methodLabel = ifelse(str_starts(method, 'Decision'), 
                                paste0(method, '\n    with h=0 fixed'),
                                paste0(method, '\n    with delta1=0.95 fixed')))
    
  
  
  capt = unadjLst$caption
  
  p1 = ggplot(data=errorDat, aes(x=errorType, y=error, fill=calibrate)) +
    geom_bar(stat='identity', position = position_dodge()) +
    geom_hline(yintercept = 0.05, color = 'gray60', 
               size = 1, linetype=2)+
    scale_y_continuous(limits = c(0,1))+
    labs(x='', y='error rate', caption = capt, fill='')+
    facet_grid(.~adjLabel + methodLabel) +
    scale_fill_manual(values = wes_palette("Darjeeling2")[c(4,2)]) + # my wesanderson colors
    #scale_fill_manual(values = c("#547BD3", "#D3AD4E")) + # Trevor blue and yellow -- which I don't like
    theme_bw(base_size = 13)
  
  #print(p1)
  
  p2 = ggplot(data=f1Dat, aes(x=method, y=f1, fill=calibrate)) +
    geom_bar(stat='identity', position = position_dodge()) +
    scale_y_continuous(limits = c(0,1))+
    labs(x='Calibrate on...', y='F1 score', caption = capt, fill='')+
    facet_grid(.~adjLabel) +
    scale_fill_manual(values = wes_palette("Darjeeling2")[c(2,4)]) + 
    theme_bw(base_size = 13)
  
  #print(p2)
  
  if(plotToShow == 'both'){
    print(p1)
    print(p2)
  }else if(plotToShow == 'error'){
    print(p1)
  }else if(plotToShow == 'f1'){
    print(p2)
  }else{
    # if "neither", then return the data frames used for plotting
    return(list(errorDat = errorDat, 
                f1Dat = f1Dat))
  }
}


# April 12: plot distribution of biases (as in negative control estimates)
# April 13: for EUMAEUS results, pull se_log_rr and ci bounds too for funnel plots
getBiases <- function(database_id,
                      method,
                      analysis_id,
                      exposure_id,
                      prior_id,
                      resPath,
                      summ = NULL,
                      minOutcomes = 5,
                      estimateType = 'Mean',
                      source = 'Bayesian'){
  
  # load from summary files
  # can take Bayesian results OR EUMAEUS results
  if(source == 'Bayesian'){
    if(is.null(summ)){
      db_name = ifelse(database_id == 'MDCD', 'IBM_MDCD', database_id)
      fname = sprintf('AllSummary-%s-%s.rds', db_name, method)
      # cat(file.path(resPath, fname))
      # cat('\n')
      summ = readRDS(file.path(resPath, fname))
    }
    
    # get relevant rows for NCs
    dat  = summ %>% 
      filter(analysis_id == !!analysis_id, 
             exposure_id == !!exposure_id,
             prior_id == !! prior_id,
             negativeControl == TRUE)
    
    # check if has minOutcomes of results
    if(nrow(dat) == 0 || length(unique(dat$outcome_id)) < minOutcomes){
      mes = sprintf('Num. of negative controls for analysis %s, exposure %s and prior %s is smaller than minimum %s!\n',
                    analysis_id, exposure_id, prior_id, minOutcomes)
      cat(mes)
      return(NULL)
    }
    
    # get latest period results
    maxPeriod = max(dat$period_id)
    dat = dat %>% filter(period_id == maxPeriod)
    
    # select relevant estimates by name
    estName = paste0('post',estimateType)
    estimates = dat %>% select(!!estName) %>% pull()
  }else{
    # if using EUMAEUS results
    # load from local save...
    allNCs = readRDS(file.path(resPath, 'CompNegControls.rds'))
    names(allNCs) = tolower(names(allNCs))
    
    db_name = ifelse(database_id %in% c('MDCD','MDCR'), 
                     sprintf('IBM_%s', database_id), 
                     database_id)
    
    dat = allNCs %>% 
      filter(database_id == db_name,
             method == !!method,
             analysis_id == !!analysis_id, 
             exposure_id == !!exposure_id,
             !is.na(log_rr),
             !is.na(se_log_rr))
    
    # check if has minOutcomes of results
    if(nrow(dat) == 0 || length(unique(dat$outcome_id)) < minOutcomes){
      mes = sprintf('Num. of negative controls for analysis %s, exposure %s and prior %s is smaller than minimum %s!\n',
                    analysis_id, exposure_id, prior_id, minOutcomes)
      cat(mes)
      return(NULL)
    }
    
    # get latest period results
    maxPeriod = max(dat$period_id)
    dat = dat %>% filter(period_id == maxPeriod)
    
    # pull estimates
    estimates = dat %>% select(log_rr) %>% pull()
    ses = dat %>% select(se_log_rr) %>% pull()
    ci_lbs = dat %>% select(ci_95_lb) %>% pull()
    ci_ubs = dat %>% select(ci_95_ub) %>% pull()
  }
  
  # return a result list
  res = 
  list(estimates = estimates, 
       mean = mean(estimates),
       sd = sd(estimates),
       num = length(estimates),
       estimateType = estimateType,
       database_id = database_id,
       method = method,
       analysis_id = analysis_id,
       exposure_id = exposure_id,
       outcome_id = dat$outcome_id,
       prior_id = prior_id,
       period_id = maxPeriod)
  
  if(source != 'Bayesian'){
    res$se = ses
    res$ci_95_lb = ci_lbs
    res$ci_95_ub = ci_ubs
  }
  
  return(res)
}

# ## try it... --------------
# summarypath = '~/Documents/Research/betterResults/summary'
# HCbiases = getBiases(database_id = 'MDCD',
#                      method = 'HistoricalComparator', # 'SCCS'
#                      analysis_id = 2,
#                      exposure_id = 211983,
#                      prior_id = 1,
#                      resPath = summarypath,
#                      minOutcomes = 5,
#                      estimateType = 'Median')

# # another example ----
# cachepath = './localCache/'
# HCbiases = getBiases(database_id = 'MDCD',
#                      method = 'HistoricalComparator', # 'SCCS'
#                      analysis_id = 2,
#                      exposure_id = 211983,
#                      prior_id = 1,
#                      resPath = cachepath,
#                      minOutcomes = 5,
#                      source = 'EUMAEUS')

# April 14: produce funnel plot for systematic errors
plotSystematicErrors <- function(resls, xLabel = 'Rate ratio estimates') {
  d = as.data.frame(resls)
  d$Significant <- d$ci_95_lb > 1 | d$ci_95_ub < 1
  
  oneRow <- data.frame(nLabel = paste0(formatC(nrow(d), big.mark = ","), 
                                       " estimates"),
                       meanLabel = paste0(formatC(100 *
                                                    mean(!d$Significant, na.rm = TRUE), 
                                                  digits = 1, format = "f"), 
                                          "% of CIs includes 1"))
  
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  themeLA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 0)
  
  alpha <- 1 - min(0.95 * (nrow(d)/50000)^0.1, 0.95)
  plot <- ggplot2::ggplot(d, ggplot2::aes(x = estimates, y = se)) +
    ggplot2::geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.5) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1/qnorm(0.025)),
                         colour = rgb(0.8, 0, 0),
                         linetype = "dashed",
                         size = 1,
                         alpha = 0.5) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1/qnorm(0.975)),
                         colour = rgb(0.8, 0, 0),
                         linetype = "dashed",
                         size = 1,
                         alpha = 0.5) +
    ggplot2::geom_point(size = 2, color = rgb(0, 0, 0, alpha = 0.05), alpha = alpha, shape = 16) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_label(x = log(0.3),
                        y = 0.97,
                        alpha = 1,
                        hjust = "left",
                        ggplot2::aes(label = nLabel),
                        size = 5,
                        data = oneRow) +
    ggplot2::geom_label(x = log(0.3),
                        y = 0.87,
                        alpha = 1,
                        hjust = "left",
                        ggplot2::aes(label = meanLabel),
                        size = 5,
                        data = oneRow) +
    ggplot2::scale_x_continuous(xLabel, limits = log(c(0.1,
                                                       10)), breaks = log(breaks), labels = breaks) +
    ggplot2::scale_y_continuous("Standard error", limits = c(0, 1)) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.y = themeRA,
                   axis.text.x = theme,
                   axis.title = theme,
                   legend.key = ggplot2::element_blank(),
                   strip.text.x = theme,
                   strip.background = ggplot2::element_blank(),
                   legend.position = "none")
  return(plot)
}


# 04/27/2022: plot the temporal calibration results (on Delta1 only!)
# (one plot for choosing delta1, one for using default delta1=0.95)
plotTempDelta1 <- function(database_id,
                            method,
                            analysis_id,
                            exposure_id,
                            prior_id,
                            summaryPath,
                            cachePath,
                            alpha = 0.05,
                            minOutcomes = 5,
                            useAdjusted = TRUE,
                            showPlots = TRUE){
  # calibrate Delta1 temporally
  caliDelta = tempCalibrateByDelta1(database_id, method, analysis_id, exposure_id, prior_id,
                                    resPath = summaryPath,
                                    alpha = alpha, minOutcomes = minOutcomes, 
                                    useAdjusted = useAdjusted, evalType2 = TRUE)
  
  
  # generate figure caption (with info)
  analysis_name = readRDS(file.path(cachePath,'analyses.rds')) %>%
    filter(method == !!method, analysis_id == !!analysis_id) %>%
    select(description, time_at_risk) %>%
    mutate(description_text = sprintf('%s, time at risk = %s days', description, time_at_risk)) %>%
    select(description_text) %>% pull() %>% as.character()
  
  prior_name = readRDS(file.path(cachePath,'priorTable.rds')) %>%
    filter(prior_id == !!prior_id) %>%
    mutate(priorLabel = sprintf('Mean=%s, SD=%.1f', Mean, Sd)) %>%
    select(priorLabel) %>% pull() %>% as.character()
  
  exposure_name = readRDS(file.path(cachePath,'exposures.rds')) %>%
    filter(exposure_id == !!exposure_id) %>%
    select(exposure_name) %>% pull() %>% as.character()
  
  capt = sprintf('%s\nExposure: %s\nDatabase: %s\nPrior:%s\nalpha=%.2f; default delta1=%.2f',
                 analysis_name, exposure_name, database_id, prior_name, alpha, 1-alpha)
  
  
  if(useAdjusted){
    capt = paste0(capt,'\nWith bias adjustment')
  }else{
    capt = paste0(capt,'\nWithout bias adjustment')
  }
  
  withChoose = data.frame(period_id = rep(caliDelta$period_id, 3),
                          y = c(caliDelta$calibratedDelta1, 
                                caliDelta$type1,
                                caliDelta$type2),
                          stats = rep(c('delta1', 'type 1', 'type 2'), each = nrow(caliDelta)))
  noChoose = data.frame(period_id = rep(caliDelta$period_id, 3),
                        y = c(rep(1-alpha, nrow(caliDelta)), 
                              caliDelta$uncalibratedType1,
                              caliDelta$uncalibratedType2),
                        stats = rep(c('delta1', 'type 1', 'type 2'), each = nrow(caliDelta)))
  
  res = bind_rows(withChoose, noChoose)
  res$method = rep(c('choosing delta1','not choosing delta1'), 
                   each = nrow(withChoose))
  res$method = factor(res$method, 
                      levels = c('not choosing delta1', 'choosing delta1'))
  
  period_breaks = seq(from = min(caliDelta$period_id),
                      to = max(caliDelta$period_id),
                      by = 2)
  period_labels = as.integer(period_breaks)
  
  p = ggplot(res, aes(x=period_id, y=y, color=stats))+
    geom_line(size = 1.5) +
    geom_point(size=2)+
    geom_hline(data = data.frame(yinter=c(0.05,0.95,0.05), 
                          method = factor(c('choosing delta1',
                                            'choosing delta1',
                                            'not choosing delta1'),
                                     levels = c('not choosing delta1', 'choosing delta1'))),
               mapping = aes(yintercept = yinter), 
               color = 'gray60', 
               size = 1, linetype=2)+
    scale_y_continuous(limits = c(0,1))+
    scale_x_continuous(breaks = period_breaks, labels = period_labels)+
    labs(x='analysis period (months)', y='', caption = capt, color='')+
    scale_color_manual(values = wes_palette("Royal1")[c(1,2,4)]) +
    facet_grid(.~method)+
    theme_bw(base_size = 13)
  
  # show plots if...
  if(showPlots){
    print(p)
  }
  
  return(res)
}

## try it
# summarypath = '~/Documents/Research/betterResults/summary'
# cachepath = './localCache/'
# caliDelta = 
#   calibrateByDelta1(database_id= 'MDCD', 
#                     method = 'HistoricalComparator', 
#                     analysis_id = 2,
#                     exposure_id = 211983, 
#                     prior_id = 1,
#                     resPath = summarypath,
#                     alpha = 0.05, 
#                     minOutcomes = 5, 
#                     useAdjusted = TRUE, evalType2 = TRUE)
# res = plotTempDelta1(database_id = 'MDCD',
#                 method = 'HistoricalComparator', # 'SCCS'
#                 analysis_id = 2,
#                 exposure_id = 211983,
#                 prior_id = 1,
#                 summaryPath = summarypath,
#                 cachePath = cachepath,
#                 useAdjusted = TRUE)
