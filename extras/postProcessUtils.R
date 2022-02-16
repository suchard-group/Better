# functions used to post process cluster run results
# (from using likelihood profiles)

library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(cores = 4)

#####
## function to put together all summaries for each (or a set of) exposure
## and for each method

##also need to produce grouping for negative/positive controls

pullResults <- function(database_id, 
                        method, 
                        exposure_id, 
                        resultsPath, 
                        savePath = NULL,
                        IPCpath='./localCache/'){
  
  # get IPC table
  IPCtable = readRDS(file.path(IPCpath, 'allIPCs.rds'))
  
  res = 
    foreach(expo = exposure_id, .combine = 'bind_rows') %dopar% {
    pullResultsOneExpo(database_id, method, expo, resultsPath, IPCtable)
    } %>%
    mutate(database_id = database_id,
           method = method)
  
  if(!is.null(savePath)){
    fname = paste0(sprintf('all_summary_%s_%s_',
                    database_id,
                    method),
                   paste(exposure_id, collapse = '-'),
                   '.rds')
  }else{
    res
  }
  
}

## function for doing for one method and one exposure only
pullResultsOneExpo <- function(database_id, 
                               method, 
                               exposure_id, 
                               resultsPath,
                               IPCtable,
                               savePath = NULL,
                               verbose = FALSE){
  
  # get a list of negative control outcome ids
  NCs = unique(IPCtable$NEGATIVE_CONTROL_ID)
  
  # produce file name pattern for summaries
  fnamePattern_analyses = sprintf(
    '%s_%s_%s_period[1-9]*_analysis[1-9]*_summary.rds',
    database_id,
    method,
    exposure_id
  )
  
  fnamePattern_periods = sprintf(
    'period_summary_%s_%s_%s_period[1-9]*_analysis.*\\.rds',
    database_id,
    method,
    exposure_id
  )
  
  analyses_files = list.files(path = resultsPath, 
                              pattern = fnamePattern_analyses)
  periods_files = list.files(path = resultsPath, 
                             pattern = fnamePattern_periods)
  
  if(verbose){
    cat(sprintf('For database %s, method %s, exposure %s, found %s separate analysis files and %s period-aggregate files.\n',
                database_id,
                method,
                exposure_id,
                length(analyses_files),
                length(periods_files)))
  }
  
  # read all those files and combine rows into a big data table
  allFiles = c(analyses_files, periods_files)
  #allFiles = analyses_files
  if(length(allFiles) == 0) return()
  
  res = 
    foreach(fname = allFiles, .combine = 'bind_rows') %dopar% {
      if(fname %in% analyses_files){
        this.period = getPeriodID(fname)
        readRDS(file.path(resultsPath, fname)) %>% 
          mutate(period_id = this.period)
      }else{
        readRDS(file.path(resultsPath, fname))
      }
    } %>% 
    distinct() %>%
    mutate(negativeControl = outcome_id %in% NCs,
           exposure_id = exposure_id)
  
  # save it if...
  if(!is.null(savePath)){
    fname = sprintf('exposure_summary_%s_%s_%s.rds',
                    database_id,
                    method,
                    exposure_id)
    saveRDS(res, file.path(savePath, fname))
  }else{
    # return it if not to save
    res
  }
}

## helper function to get the period id from file names
## (apparently for separate analysis files, period ids are not recorded in data frame...)
getPeriodID <- function(fname){
  chunks = stringr::str_split(fname, pattern = '([1-9]*_period)|(_a)') %>% 
    unlist()
  as.numeric(chunks[2])
}

#### examples tried ------------
# ## try it
# IPCs = readRDS('./localCache/allIPCs.rds')
# resPath = '~/Documents/Research/betterResults/betterResults-MDCD/'
# comb_res = pullResultsOneExpo(database_id = 'IBM_MDCD',
#                               method = 'SCCS',
#                               exposure_id = 211981,
#                               resultsPath = resPath,
#                               IPCtable = IPCs, verbose=TRUE)
# 
# all_expos = sort(unique(IPCs$EXPOSURE_ID))
# big_comb_res = pullResults(database_id = 'IBM_MDCD',
#                            method = 'SCCS',
#                            exposure_id = all_expos,
#                            resultsPath = resPath,
#                            IPCpath = './localCache/')



## function to make test decisions using the Bayesian postprob thresholds
## re-arrange the columns a little, and allow saving the dataframe w/ decisions
makeDecisions <- function(summ, 
                          delta1 = c(0.80,0.90,0.95),
                          delta0 = c(0.90,0.95,0.99),
                          filepath = NULL,
                          savepath = NULL,
                          fname = NULL){
  ## summ: default should be a dataframe
  ##       if character string, then it's a fname to load RDS from
  ## filepath: the path to look for the summary files to process
  ## savepath: path to save the processed file
  ## fname: customize a fname; default "summary_decisions.rds"
  
  if(is.character(summ)){
    if(is.null(savepath)){
      stop('No summary file path provided! You need to give me one.\n')
    }
    summ = readRDS(file.path(filepath,summ))
  }
  
  ## decide on elevated risk signals (H1 true)
  for(d1 in delta1){
    d1_perc = d1 * 100 %>% as.integer()
    summ = summ %>%
      mutate("signal{d1_perc}" := P1 > d1,
             "adjustedSignal{d1_perc}" := adjustedP1 > d1)
  }
  
  ## decide on safety (H0 true)
  for(d0 in delta0){
    d0_perc = d0 * 100 %>% as.integer()
    summ = summ %>%
      mutate("futility{d0_perc}" := P0 > d0,
             "adjustedFutility{d0_perc}" := adjustedP0 > d0)
  }
  
  ## put important columns up front
  summ = summ %>% 
    relocate(database_id, method, analysis_id, exposure_id, 
             outcome_id, period_id, prior_id) 
  
  ## save if...
  if(!is.null(savepath)){
    if(!dir.exists(savepath)) dir.create(savepath)
    if(!is.null(fname)) fname = "summary_decisions.rds"
    
    saveRDS(summ, file = file.path(savepath, fname))
  }else{
    # if not save, then directly return
    return(summ)
  }
  
}

### example below ---------
# ## try it
# big_comb_res = makeDecisions(big_comb_res)


## function to get earliest time to signal/futility for each analysis across all periods

## (1) helper function to get a table of threshold values
getThresholdTable <- function(delta1 = c(0.80,0.90,0.95),
                              delta0 = c(0.90,0.95,0.99),
                              savepath = NULL){
  ## read from saved table if exists
  if(!is.null(savepath) && file.exists(file.path(savepath, 'thresholdTable.rds'))){
    readRDS(file.path(savepath, 'thresholdTable.rds'))
  } 
  
  ## otherwise, construct it and save it to `savepath`
  nd1s = length(delta1)
  nd0s = length(delta0)
  thresholds = data.frame(p1Threshold = rep(delta1, each=nd0s),
                          p0Threshold = rep(delta0, nd1s))
  thresholds$threshold_ID = c(1:nrow(thresholds))
  
  ## save it
  if(!is.null(savepath)){
    saveRDS(thresholds, file = file.path(savepath, 'thresholdTable.rds'))
  }
  
  ## return it
  thresholds
}

# ## try it and save a local version if possible
# thresholds = getThresholdTable(savepath = './localCache/')


## (2) function to determine decision, earliest-time-to-decision, 
##     if earliest decision contradicts final decision
##     AND if the decision is same with truth (v.s. negative/positive control)
## FOR one particular threshold (either adjusted or un-adjusted) of Signal OR Futility
## AND for one particular (DB, method, analysis, exposure_id, outcome_id, prior_id) combination

## (2.a) little helper function to wrangle any string's first letter to capital
##       no matter what
capitalizeFirst <- function(s){
  sapply(s,
         function(w) {
           ll = stringr::str_split(w, '') %>% unlist()
           ll[1] = toupper(ll[1])
           paste0(ll, collapse = '')
         }) %>%
    as.character()
}

## (2.b) function for one particular (DB, method, analysis, exposure_id, outcome_id, prior_id) combination
##       for signal/futility with one threshold, with/without adjustment

## if the threshold is NEVER crossed, then return "-1", 
## meaning not enough evidence for a clear signal in either way
checkOneDecisionRule <- function(df,
                                 threshold, 
                                 decisionType = 'signal',
                                 adjusted = FALSE){
  ## if somehow df has no rows at all
  ## return NULL
  if(nrow(df) == 0){return(NULL)}
  
  # get the relevant variable name
  if(adjusted){
    decisionVar = paste0('adjusted', 
                         stringr::str_to_title(decisionType), 
                         as.integer(threshold * 100))
  }else{
    decisionVar = paste0(decisionType, 
                         as.integer(threshold * 100))
  }
  df = df %>% select(database_id, method,analysis_id,
                     exposure_id,outcome_id,period_id,
                     prior_id,one_of(decisionVar),
                     negativeControl) %>%
    arrange(period_id)
  
  
  ## if all are FALSE, return -1 as earliest-time (meaning threshold never hit!)
  if(!any(df[[decisionVar]])){
    decision = FALSE
    res = data.frame(timeToSignal = -1,
               contradict = FALSE)
  }else{
  ## otherwise, the earliest period where decision is TRUE
    decision = TRUE
    np = nrow(df)
    res = data.frame(timeToSignal = min(df$period_id[df[[decisionVar]]]),
                     contradict = !(df[[decisionVar]][np]))
  }
  
  ## then compare the decision with the truth
  ## for NC: we need signal=FALSE but futility=TRUE
  ## for IPC: we need signal=TRUE but futility=FALSE
  isNC = df$negativeControl[1]
  if(isNC){
    res$correct = ifelse(decisionType=='signal', !decision, decision)
  }else{
    res$correct = ifelse(decisionType=='signal', decision, !decision)
  }
  
  # ## work on column names (adjusted v.s. un-adjusted)
  # if(adjusted){
  #   names(res) = paste0('adjusted',capitalizeFirst(names(res)))
  # }
  
  ## add columns for decisionType, threshold and adjusted
  res$decisionType = decisionType
  res$threshold = threshold
  res$adjusted = adjusted
  
  res
  
}

## (2.c) function to do everything (combinely process all possible decision rule combos) for each
##       (DB, method, analysis, exposure_id, outcome_id, prior_id) combination)

### !!! NEED to decide how to handle the case of "neither" decision boundary crossed !!!
### RIGHT NOW: only focus on whether or not it signals --- is the "to signal" decision is correct, then 
checkAllDecisionRules <- function(df, 
                                  thresholdTable){
  ## if somehow df has no rows at all
  ## return NULL
  if(nrow(df) == 0){return(NULL)}
  
  ## produce all one-rule decision results first
  delta1 = sort(unique(thresholdTable$p1Threshold))
  delta0 = sort(unique(thresholdTable$p0Threshold))
  
  ## get all one-decision-rule results
  res = NULL
  for(adj in c(FALSE, TRUE)){
    resSignal = foreach(d1=delta1, .combine = 'bine_rows') %dopar% {
      checkOneDecisionRule(df, threshold = d1, decisionType ='signal', adjusted=adj)
    }
    
    # early stopping if returns NULL
    if(is.null(resSignal)){return(NULL)}
    
    resSignal = foreach(d0=delta0, .combine = 'bine_rows') %dopar% {
      checkOneDecisionRule(df, threshold = d0, decisionType ='futility', adjusted=adj)
    }
    
    res = bind_rows(res, resSignal, resSignal)
  }
  
  ## go through them and produce long format results
  threshold_IDs = sort(thresholdTable$threshold_ID)
  
  decisionChecks = 
  foreach(adj=c(FALSE, TRUE), .combine = 'bind_rows') %dopar% {
    foreach(id=threshold_IDs, .combine = 'bind_rows') %dopar% {
      d1 = thresholdTable$p1Threshold[id]
      d0 = thresholdTable$p0Threshold[id]
      thispair = res %>% 
        filter((decisionType=='signal' & threshold==d1) | (decisionType=='futility' & threshold==d0)) %>%
        filter(adjusted == adj)
      
      ## work with the timeToSignal column
      if(all(thispair$timeToSignal < 0)){
        # if neither threshold crossed...
        this.combo = 
          data.frame(decision='neither',
                     timeToSignal = -1,
                     correct = thispair %>% 
                       filter(decisionType=='signal') %>% 
                       select(correct) %>% pull()) ## <- need to rethink how to judge a "neither" decision!!
        
      }else if(sum(thispair$timeToSignal > 0) == 1){
        # if only one threshold crossed
        this.combo = thispair %>% filter(timeToSignal > 0) %>% 
          select(decision = decisionType, timeToSignal, correct)
      }else{
        # if both thresholds crossed
        # go with whichever is the early one
        this.combo = thispair %>% 
          filter(timeToSignal == min(thispair$timeToSignal)) %>% 
          select(decision = decisionType, timeToSignal, correct)
      }
      
      ## return with information about adjustment and decision rule ID
      this.combo %>% mutate(threshold_id = id,
                            adjusted = adj)
    }
    
  }
  
  ## attach with other info for this analysis combination
  infos = df[1,] %>% 
    select(database_id, method, analysis_id, exposure_id, outcome_id, prior_id, negativeControl)
  
  bind_cols(infos, decisionChecks)
  
}


## (2) helper function to do stuff for each (DB, method, analysis, exposure_id, outcome_id, prior_id) combination
getTimeToSignal <- function(decisions){
  decisions %>% 
    group_by(database_id, method, analysis_id, exposure_id, outcome_id, prior_id) %>%
    summarize() # TBD!!!
}


