# functions used to post process cluster run results
# (from using likelihood profiles)

library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel()

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
  chunks = stringr::str_split(fname, pattern = '([1-9]*_period)|(_a)') %>% unlist()
  as.numeric(chunks[2])
}

#### examples tried ------------
# ## try it
# IPCs = readRDS('./localCache/allIPCs.rds')
# comb_res = pullResultsOneExpo(database_id = 'IBM_MDCD',
#                               method = 'SCCS',
#                               exposure_id = 211981,
#                               resultsPath = '~/Documents/results/betterResults/',
#                               IPCtable = IPCs, verbose=TRUE)
# 
# all_expos = sort(unique(IPCs$EXPOSURE_ID))
# big_comb_res = pullResults(database_id = 'IBM_MDCD',
#                            method = 'SCCS',
#                            exposure_id = all_expos,
#                            resultsPath = '~/Documents/results/betterResults/',
#                            IPCpath = './localCache/')



## function to make test decisions using the Bayesian postprob thresholds
makeDecisions <- function(summ, 
                          delta1 = c(0.80,0.90,0.95),
                          delta0 = c(0.90,0.95,0.99),
                          filepath = NULL,
                          savepath = NULL){
  ## summ: default should be a dataframe
  ##       if character string, then it's a fname to load RDS from
  ## filepath: the path to look for the summary files to process
  
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
  
  summ
  
}

### example below ---------
# ## try it
# big_comb_res = makeDecisions(big_comb_res)


## function to get earliest time to signal/futility for each analysis across all periods

## (1) helper function to do stuff for each (DB, method, analysis, exposure_id, outcome_id, prior_id) combination
getTimeToSignal <- function(decisions){
  decisions %>% 
    group_by(database_id, method, analysis_id, exposure_id, outcome_id, prior_id) %>%
    summarize() # TBD!!!å
}


