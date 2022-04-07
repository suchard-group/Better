# functions used to post process cluster run results
# (from using likelihood profiles)

library(tidyverse)
library(foreach)
library(doParallel)
#registerDoParallel(cores = 4)

#####
## function to put together all summaries for each (or a set of) exposure
## and for each method

##also need to produce grouping for negative/positive controls

## update: rearrange column order here..
pullResults <- function(database_id, 
                        method, 
                        exposure_id, 
                        resultsPath, 
                        savePath = NULL,
                        IPCpath='./localCache/'){
  
  # get IPC table
  IPCtable = readRDS(file.path(IPCpath, 'allIPCs.rds'))
  
  res = 
    foreach(expo = exposure_id, .combine = 'bind_rows', 
            .multicombine = TRUE, .export = 'pullResultsOneExpo') %dopar% {
    pullResultsOneExpo(database_id, method, expo, resultsPath, IPCtable)
    } %>%
    mutate(database_id = database_id,
           method = method) %>% 
      relocate(database_id, method, analysis_id, exposure_id,
               outcome_id, period_id, prior_id)
  
  # if savePath provided, save...
  if(!is.null(savePath)){
    ## if exposures are NOT all, then include there ids in file name
    ## otherwise, no exposure ids included
    if(length(exposure_id) == 10){
      fname = sprintf('AllSummary-%s-%s.rds',
                      database_id, method)
    }else{
      fname = paste0(sprintf('AllSummary-%s-%s-',
                             database_id,
                             method),
                     paste(exposure_id, collapse = '+'),
                     '.rds')
    }
    saveRDS(res, file.path(savePath, fname))
    
    cat(sprintf('Results for database %s and method %s saved at %s\n', 
                database_id, method, file.path(savePath, fname)))
  }
  
  # return no matter what
  res
  
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
    foreach(fname = allFiles, .combine = 'bind_rows',
            .multicombine = TRUE) %dopar% {
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

# 04/06/2022: GBS specific functions
pullGBSResultsOneExpo <- function(database_id, 
                                  method, 
                                  exposure_id, 
                                  resultsPath){
  fname = sprintf('AllSummary-%s-%s-%s.rds', 
                  database_id, method, exposure_id)
  fpath = file.path(resultsPath, fname)
  if(!file.exists(fpath)){
    mes = sprintf('Summary file at %s does NOT exist! Skipped', fpath)
    cat(mes)
    return(NULL)
  }
  readRDS(fpath)
}

pullGBSResults <- function(database_id, 
                           methods, 
                           exposure_ids, 
                           resultsPath){
  df = NULL
  for(me in methods){
    for(expo in exposure_ids){
      df = rbind(df, 
                 pullGBSResultsOneExpo(database_id = database_id,
                                       method = me,
                                       exposure_id = expo,
                                       resultsPath = resultsPath))
    }
  }
  
  df
}

#### examples tried ------------
# ## try it
# IPCs = readRDS('./localCache/allIPCs.rds')
# resPath = '~/Documents/Research/betterResults/betterResults-CCAE/' # -MDCD
# comb_res = pullResultsOneExpo(database_id = 'CCAE', # 'IBM_MDCD'
#                               method = 'SCCS',
#                               exposure_id = 211981,
#                               resultsPath = resPath,
#                               IPCtable = IPCs, verbose=TRUE)
# 
# all_expos = sort(unique(IPCs$EXPOSURE_ID))
# big_comb_res = pullResults(database_id = 'CCAE', # 'IBM_MDCD'
#                            method = 'SCCS',
#                            exposure_id = all_expos,
#                            resultsPath = resPath,
#                            IPCpath = './localCache/')


#### Feb 17 updated faster functions to make decisions and judge them---------------
## March 2, update the earliest period_id with a decision

## function to get decisions (adjusted & unadjusted) based on d1 and d0 thresholds
getOverallDecisions <- function(df, d1=0.8, d0=0.9){
  
  nr = nrow(df)
  
  signals = apply(df %>% select(P1, adjustedP1), 2,
                  function(v){
                    c(min(df$period_id[v > d1]), 
                      v[nr] > d1)
                  })
  futilities = apply(df %>% select(P0, adjustedP0), 2,
                     function(v){
                       c(min(df$period_id[v > d0]),
                         v[nr] > d1)
                     })
  res = cbind(signals, futilities)
  res[2,] = as.logical(res[2,]) != (res[1,] < Inf)
  interpretDecision(res) %>% as.data.frame()
}

## utils function that works with `getOverallDecisions`
## March 2 update: add period_id for the decision 
##                 (earliest month with decision, Inf if no decision is made)
interpretDecision <- function(block){
  
  if(block[1,1] < Inf){
    decision = list(decision='signal',
                    contradict = as.logical(block[2,1]),
                    timeToDecision = block[1,1])
  }else if(block[1,3] < Inf){
    decision = list(decision='futility',
                    contradict = as.logical(block[2,3]),
                    timeToDecision = block[1,3])
  }else{
    decision = list(decision='neither',contradict =FALSE,
                    timeToDecision = Inf)
  }
  
  if(block[1,2] < Inf){
    decision$adjustedDecision = 'signal'
    decision$adjustedContradict = as.logical(block[2,2])
    decision$adjustedTimeToDecision = block[1,2]
  }else if(block[1,4] < Inf){
    decision$adjustedDecision = 'futility'
    decision$adjustedContradict = as.logical(block[2,4])
    decision$adjustedTimeToDecision = block[1,4]
  }else{
    decision$adjustedDecision = 'neither'
    decision$adjustedContradict = FALSE
    decision$adjustedTimeToDecision = Inf
  }
  decision
}

## function to make decisions using all combos of thresholds
getAllOverallDecisions <- function(df, thresholdTable){
  lapply(1:nrow(thresholdTable), 
         function(i) {
           getOverallDecisions(df, 
                               d1=thresholdTable$p1Threshold[i],
                               d0=thresholdTable$p0Threshold[i])
         }) %>% 
    bind_rows() %>%
    mutate(threshold_id = thresholdTable$threshold_ID)
}



## helper function to get a table of threshold values
## update: make it possible to have differently named threshold values
getThresholdTable <- function(delta1 = c(0.80,0.90,0.95),
                              delta0 = c(0.90,0.95,0.99),
                              savepath = NULL,
                              tablename = 'thresholdTable'){
  tableFname = paste0(tablename, '.rds')
  ## read from saved table if exists
  if(!is.null(savepath) && file.exists(file.path(savepath, tableFname))){
    readRDS(file.path(savepath, tableFname))
  }

  ## otherwise, construct it and save it to `savepath`
  nd1s = length(delta1)
  nd0s = length(delta0)
  thresholds = data.frame(p1Threshold = rep(delta1, each=nd0s),
                          p0Threshold = rep(delta0, nd1s))
  thresholds$threshold_ID = c(1:nrow(thresholds))

  ## save it
  if(!is.null(savepath)){
    saveRDS(thresholds, file = file.path(savepath, tableFname))
  }

  ## return it
  thresholds
}





#### function to judge if decision is right or wrong ---------
#### and compute Type I/II errors and rates of don't know

## 2 judging styles:
## (1) strict: if a signal isn't a "signal", or a non-signal isn't "futility", then WRONG
## (2) lenient: if a signal isn't a "futility", or a non-signal isn't "signal", then RIGHT

## need to stratify by effect size as well

## helper function to produce a table of effect sizes for all outcome_ids
getEffectSizes <- function(IPCpath = './localCache/', cachePath = './localCache/'){
  ## check if already saved
  fname = 'effectSizes.rds'
  if(file.exists(file.path(cachePath, fname))){
    readRDS(file.path(cachePath, fname))
  }else{
    IPCs = readRDS(file.path(IPCpath, 'allIPCs.rds'))
    NCs = unique(IPCs$NEGATIVE_CONTROL_ID)
    res = IPCs %>% 
      select(outcome_id = OUTCOME_ID, effect_size = EFFECT_SIZE)
    res = bind_rows(res, data.frame(outcome_id = NCs, effect_size = 1))
    saveRDS(res, file.path(cachePath, fname))
    res
  }
}

## main function: update with effect size stratify
## UPDATE: add more judge style: 
##         a H0 futility is correct if no signal (but not vice versa for H1 signal)
judgeDecisions <- function(df, judgeStyle = 'strict', 
                           cachePath = './localCache/'){
  # judgeStyle: strict, H0neither, lenient
  
  # load effect size table first
  effects = getEffectSizes(cachePath, cachePath)
  
  # read in effect size for all outcomes
  df = df %>% left_join(effects)
  
  # then judge
  if(judgeStyle == 'strict'){
    df = df %>% 
      mutate(judge = if_else(negativeControl, 
                             decision == 'futility', 
                             decision == 'signal'),
             adjustedJudge = if_else(negativeControl, 
                                     adjustedDecision == 'futility', 
                                     adjustedDecision == 'signal'))
  }else if(judgeStyle == 'H0neither'){
    # more similar to the frequentist logic:
    # if P1 > ..., then signal
    # otherwise, treat it like H0 true decision
    df = df %>% 
      mutate(judge = if_else(negativeControl, 
                             decision != 'signal', 
                             decision == 'signal'),
             adjustedJudge = if_else(negativeControl, 
                                     adjustedDecision != 'signal', 
                                     adjustedDecision == 'signal'))
  }else{
    df = df %>% 
      mutate(judge = if_else(negativeControl, 
                             decision != 'signal', 
                             decision != 'futility'),
             adjustedJudge = if_else(negativeControl, 
                                     adjustedDecision != 'signal', 
                                     adjustedDecision != 'futility'))
  }
  
  # group and summarize: stratify by effect size here
  # update fix: errorRate = 1-accuracy
  df %>% 
    group_by(database_id, method, analysis_id, exposure_id, prior_id, threshold_id,
             effect_size) %>%
    summarize(errorRate = 1-mean(judge), 
              adjustedErrorRate = 1-mean(adjustedJudge),
              neitherRate = mean(decision == 'neither'),
              adjustedNeitherRate = mean(adjustedDecision == 'neither'),
              sampleSize = n()) %>%
    mutate(negativeControl = (effect_size == 1))
}


## Feb 25 2022 update ------
## function to do weighted average given error rate 
weightedAverage <- function(m,n){
  sum(m * n)/sum(n)
}

## function to get the rank of each character string within a whole vector
getRank <- function(v){
  sorted = sort(unique(v))
  sapply(v, function(x) which(sorted==x)) %>% as.vector()
}

##### Outdated slow functions of making decisions below--------------------
### saved for records....

# ## function to make test decisions using the Bayesian postprob thresholds
# ## re-arrange the columns a little, and allow saving the dataframe w/ decisions
# makeDecisions <- function(summ, 
#                           delta1 = c(0.80,0.90,0.95),
#                           delta0 = c(0.90,0.95,0.99),
#                           filepath = NULL,
#                           savepath = NULL,
#                           fname = NULL){
#   ## summ: default should be a dataframe
#   ##       if character string, then it's a fname to load RDS from
#   ## filepath: the path to look for the summary files to process
#   ## savepath: path to save the processed file
#   ## fname: customize a fname; default "summary_decisions.rds"
#   
#   if(is.character(summ)){
#     if(is.null(savepath)){
#       stop('No summary file path provided! You need to give me one.\n')
#     }
#     summ = readRDS(file.path(filepath,summ))
#   }
#   
#   ## decide on elevated risk signals (H1 true)
#   for(d1 in delta1){
#     d1_perc = d1 * 100 %>% as.integer()
#     summ = summ %>%
#       mutate("signal{d1_perc}" := P1 > d1,
#              "adjustedSignal{d1_perc}" := adjustedP1 > d1)
#   }
#   
#   ## decide on safety (H0 true)
#   for(d0 in delta0){
#     d0_perc = d0 * 100 %>% as.integer()
#     summ = summ %>%
#       mutate("futility{d0_perc}" := P0 > d0,
#              "adjustedFutility{d0_perc}" := adjustedP0 > d0)
#   }
#   
#   ## put important columns up front
#   summ = summ %>% 
#     relocate(database_id, method, analysis_id, exposure_id, 
#              outcome_id, period_id, prior_id) 
#   
#   ## save if...
#   if(!is.null(savepath)){
#     if(!dir.exists(savepath)) dir.create(savepath)
#     if(!is.null(fname)) fname = "summary_decisions.rds"
#     
#     saveRDS(summ, file = file.path(savepath, fname))
#   }else{
#     # if not save, then directly return
#     return(summ)
#   }
#   
# }
# 
# ### example below ---------
# # ## try it
# # big_comb_res = makeDecisions(big_comb_res)
# 
# 
# ## function to get earliest time to signal/futility for each analysis across all periods
# 
# ## (1) helper function to get a table of threshold values
# getThresholdTable <- function(delta1 = c(0.80,0.90,0.95),
#                               delta0 = c(0.90,0.95,0.99),
#                               savepath = NULL){
#   ## read from saved table if exists
#   if(!is.null(savepath) && file.exists(file.path(savepath, 'thresholdTable.rds'))){
#     readRDS(file.path(savepath, 'thresholdTable.rds'))
#   } 
#   
#   ## otherwise, construct it and save it to `savepath`
#   nd1s = length(delta1)
#   nd0s = length(delta0)
#   thresholds = data.frame(p1Threshold = rep(delta1, each=nd0s),
#                           p0Threshold = rep(delta0, nd1s))
#   thresholds$threshold_ID = c(1:nrow(thresholds))
#   
#   ## save it
#   if(!is.null(savepath)){
#     saveRDS(thresholds, file = file.path(savepath, 'thresholdTable.rds'))
#   }
#   
#   ## return it
#   thresholds
# }
# 
# # ## try it and save a local version if possible
# # thresholds = getThresholdTable(savepath = './localCache/')
# 
# 
# ## (2) function to determine decision, earliest-time-to-decision, 
# ##     if earliest decision contradicts final decision
# ##     AND if the decision is same with truth (v.s. negative/positive control)
# ## FOR one particular threshold (either adjusted or un-adjusted) of Signal OR Futility
# ## AND for one particular (DB, method, analysis, exposure_id, outcome_id, prior_id) combination
# 
# ## (2.a) little helper function to wrangle any string's first letter to capital
# ##       no matter what
# capitalizeFirst <- function(s){
#   sapply(s,
#          function(w) {
#            ll = stringr::str_split(w, '') %>% unlist()
#            ll[1] = toupper(ll[1])
#            paste0(ll, collapse = '')
#          }) %>%
#     as.character()
# }
# 
# ## (2.b) function for one particular (DB, method, analysis, exposure_id, outcome_id, prior_id) combination
# ##       for signal/futility with one threshold, with/without adjustment
# 
# ## if the threshold is NEVER crossed, then return "-1", 
# ## meaning not enough evidence for a clear signal in either way
# checkOneDecisionRule <- function(df,
#                                  threshold, 
#                                  decisionType = 'signal',
#                                  adjusted = FALSE){
#   ## if somehow df has no rows at all
#   ## return NULL
#   if(nrow(df) == 0){return(NULL)}
#   
#   # get the relevant variable name
#   if(adjusted){
#     decisionVar = paste0('adjusted', 
#                          stringr::str_to_title(decisionType), 
#                          as.integer(threshold * 100))
#   }else{
#     decisionVar = paste0(decisionType, 
#                          as.integer(threshold * 100))
#   }
#   df = df %>% select(database_id, method,analysis_id,
#                      exposure_id,outcome_id,period_id,
#                      prior_id,one_of(decisionVar),
#                      negativeControl) %>%
#     arrange(period_id)
#   
#   
#   ## if all are FALSE, return -1 as earliest-time (meaning threshold never hit!)
#   if(!any(df[[decisionVar]])){
#     decision = FALSE
#     res = data.frame(timeToSignal = -1,
#                contradict = FALSE)
#   }else{
#   ## otherwise, the earliest period where decision is TRUE
#     decision = TRUE
#     np = nrow(df)
#     res = data.frame(timeToSignal = min(df$period_id[df[[decisionVar]]]),
#                      contradict = !(df[[decisionVar]][np]))
#   }
#   
#   ## then compare the decision with the truth
#   ## for NC: we need signal=FALSE but futility=TRUE
#   ## for IPC: we need signal=TRUE but futility=FALSE
#   isNC = df$negativeControl[1]
#   if(isNC){
#     res$correct = ifelse(decisionType=='signal', !decision, decision)
#   }else{
#     res$correct = ifelse(decisionType=='signal', decision, !decision)
#   }
#   
#   # ## work on column names (adjusted v.s. un-adjusted)
#   # if(adjusted){
#   #   names(res) = paste0('adjusted',capitalizeFirst(names(res)))
#   # }
#   
#   ## add columns for decisionType, threshold and adjusted
#   res$decisionType = decisionType
#   res$threshold = threshold
#   res$adjusted = adjusted
#   
#   res
#   
# }
# 
# # ## try it out
# # checkOneDecisionRule(big_comb_res %>% 
# #                        filter(exposure_id==21184, outcome_id ==73302, analysis_id==1, prior_id==1),
# #                      threshold = 0.8)
# 
# ## (2.c) function to do everything (combinely process all possible decision rule combos) for each
# ##       (DB, method, analysis, exposure_id, outcome_id, prior_id) combination)
# 
# ### !!! NEED to decide how to handle the case of "neither" decision boundary crossed !!!
# ### RIGHT NOW: only focus on whether or not it signals --- is the "to signal" decision is correct, then 
# checkAllDecisionRules <- function(df, 
#                                   thresholdTable){
#   ## if somehow df has no rows at all
#   ## return NULL
#   if(nrow(df) == 0){return(NULL)}
#   
#   ## produce all one-rule decision results first
#   delta1 = sort(unique(thresholdTable$p1Threshold))
#   delta0 = sort(unique(thresholdTable$p0Threshold))
#   
#   ## get all one-decision-rule results
#   res = NULL
#   for(adj in c(FALSE, TRUE)){
#     resSignal1 = foreach(d1=delta1, .combine = 'bind_rows') %do% {
#       checkOneDecisionRule(df, threshold = d1, decisionType ='signal', adjusted=adj)
#     }
#     
#     # early stopping if returns NULL
#     if(is.null(resSignal1)){return(NULL)}
#     
#     resSignal0 = foreach(d0=delta0, .combine = 'bind_rows') %do% {
#       checkOneDecisionRule(df, threshold = d0, decisionType ='futility', adjusted=adj)
#     }
#     
#     res = bind_rows(res, resSignal1, resSignal0)
#   }
#   
#   #print(res)
#   
#   ## go through them and produce long format results
#   threshold_IDs = sort(thresholdTable$threshold_ID)
#   
#   decisionChecks = NULL
#   for(adj in c(FALSE, TRUE)){
#   # decisionChecks = 
#   # foreach(adj=c(FALSE, TRUE), .combine = 'bind_rows') %dopar% {
#     adj.res = NULL
#     for(id in threshold_IDs){
#     # foreach(id=threshold_IDs, .combine = 'bind_rows') %dopar% {
#       d1 = thresholdTable$p1Threshold[id]
#       d0 = thresholdTable$p0Threshold[id]
#       thispair = res %>% 
#         filter((decisionType=='signal' & threshold==d1) | (decisionType=='futility' & threshold==d0)) %>%
#         filter(adjusted == adj)
#       
#       #print(thispair)
#       
#       ## work with the timeToSignal column
#       if(all(thispair$timeToSignal < 0)){
#         # if neither threshold crossed...
#         this.combo = 
#           data.frame(decision='neither',
#                      timeToSignal = -1,
#                      correct = thispair %>% 
#                        filter(decisionType=='signal') %>% 
#                        select(correct) %>% pull()) ## <- need to rethink how to judge a "neither" decision!!
#         
#       }else if(sum(thispair$timeToSignal > 0) == 1){
#         # if only one threshold crossed
#         this.combo = thispair %>% filter(timeToSignal > 0) %>% 
#           select(decision = decisionType, timeToSignal, correct)
#       }else{
#         # if both thresholds crossed
#         # go with whichever is the early one
#         this.combo = thispair %>% 
#           filter(timeToSignal == min(thispair$timeToSignal)) %>% 
#           select(decision = decisionType, timeToSignal, correct)
#       }
#       
#       ## return with information about adjustment and decision rule ID
#       # this.combo %>% mutate(threshold_id = id,
#       #                       adjusted = adj)
#       
#       adj.res = bind_rows(adj.res, this.combo %>% mutate(threshold_id = id,
#                                                          adjusted = adj))
#     }
#     
#     #adj.res
#     
#     decisionChecks = rbind(decisionChecks, adj.res)
#     
#   }
#   
#   ## a wide format with adjusted-named columns
#   adjusted_checks = decisionChecks %>% 
#     filter(adjusted==TRUE) %>%
#     select(adjustedDecision = decision, 
#            adjustedTimeToSignal = timeToSignal, 
#            adjustedCorrect = correct)
#   
#   decisionChecks = decisionChecks %>% 
#     filter(adjusted==FALSE) %>%
#     select(-adjusted) %>%
#     bind_cols(adjusted_checks) %>%
#     relocate(threshold_id)
#   
#   ## attach with other info for this analysis combination
#   infos = df[1,] %>% 
#     select(database_id, method, analysis_id, exposure_id, outcome_id, prior_id, negativeControl)
#   
#   bind_cols(infos, decisionChecks)
#   
# }
# 
# # ## try it out
# # allRules =
# #   checkAllDecisionRules(big_comb_res %>%
# #                        filter(exposure_id==21184, outcome_id ==73302, analysis_id==1, prior_id==1),
# #                      thresholds)
# 
# 
# ## (2) helper function to do stuff for each (DB, method, analysis, exposure_id, outcome_id, prior_id) combination
# getDecisionChecks <- function(decisions, savepath=NULL, fname=NULL){
#   # decisions %>% 
#   #   group_by(database_id, method, analysis_id, exposure_id, outcome_id, prior_id) %>%
#   #   summarize() # TBD!!!
#   
#   # split up by grouping
#   groups = split(
#     decisions,
#     list(
#       decisions$database_id,
#       decisions$method,
#       decisions$analysis_id,
#       decisions$exposure_id,
#       decisions$outcome_id,
#       decisions$prior_id
#     ),
#     drop=TRUE
#   )
#   
#   # for each group, do the decision check stuff
#   allChecks = lapply(groups, 
#                      checkAllDecisionRules,
#                      thresholdTable = getThresholdTable(savepath=savepath))
#   
#   # back to data frame format
#   allChecks = bind_rows(allChecks)
#   
#   if(!is.null(savepath)){
#     if(!is.null(fname)){
#       fpath = file.path(savepath, fname)
#     }else{
#       fpath = file.path(savepath, 'allDecisionChecks.rds')
#     }
#     
#     saveRDS(allChecks, fpath)
#     
#   }
#   
#   allChecks
#   
# }
# 
# ## try it out
# # allChecks = getDecisionChecks(big_comb_res %>% sample_n(200),
# #                               savepath = './localCache/')
# 
