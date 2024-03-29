# functions used to post process cluster run results
# (from using likelihood profiles)

library(tidyverse)
library(ggplot2)
library(ggridges)
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
    ## create the savePath if not exist
    if(!dir.exists(savePath)){
      dir.create(savePath)
    }
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


## 09/06/2022 updates:
## functions for pulling posterior samples for an outcome----
## 09/13/2022: add relevant info to the pulled samples
pullPostSamples <- function(database_id, 
                            method,
                            analysis_id,
                            exposure_id,
                            outcome_id = 343, # default to GBS id
                            prior_id = 3, # default SD=4 prior
                            resultsPath = '~/Documents/Research/betterGBSanalysesResults/GBSsamples/',
                            savePath = NULL){
  if((!is.null(savePath)) && (!dir.exists(savePath))){
    dir.create(savePath)
  }
  
  fnamePattern = sprintf(
    '%s_%s_%s_period[0-9]*_analysis%s_samples.rds',
    database_id,
    method,
    exposure_id,
    analysis_id
  )
  
  sample_files = list.files(path = resultsPath, 
                            pattern = fnamePattern)
  
  if(length(sample_files) == 0){
    cat(sprintf('No samples exist for %s, %s analysis %s, exposure %s. Skipped!\n',
                database_id,
                method, analysis_id,
                exposure_id))
    return(numeric())
  }
  
  postSamps = NULL
  adjustedPostSamps = NULL
  raw_period_ids = NULL
  adj_period_ids = NULL
  
  for(f in sample_files){
    samps = readRDS(file.path(resultsPath, f))[[prior_id]]
    
    # check if outcome_id is available
    allOutcomes = rownames(samps$postSamps) %>% as.numeric()
    if(outcome_id %in% allOutcomes){
      cat(sprintf('Found outcome %s in file %s\n', outcome_id, f))
      
      this.postSamps = samps$postSamps[as.character(outcome_id),] %>% as.vector()
      if(length(this.postSamps) > 0){
        postSamps = c(postSamps, 
                      this.postSamps)
        raw_period_ids = c(raw_period_ids,
                           unlist(str_split(f, 'period|_'))[5] %>% as.numeric())
        
      }else{
        cat(sprintf('No raw samples for %s somehow...\n', outcome_id))
      }
      
      this.adjustedPostSamps = samps$adjustedPostSamps[as.character(outcome_id),] %>% as.vector()
      if(length(this.adjustedPostSamps) > 0){
        adjustedPostSamps = c(adjustedPostSamps, this.adjustedPostSamps)
        adj_period_ids = c(adj_period_ids, unlist(str_split(f, 'period|_'))[5] %>% as.numeric())
      }else{
        cat(sprintf('No adjusted samples for %s somehow...\n', outcome_id))
      }
      
    }else{
      cat(sprintf('Did not find outcome %s in file %s!!\n', outcome_id, f))
    }
    
    if(which(sample_files == f) == 1){
      numsamps = ncol(samps$postSamps)
    }
  }
  
  if(length(postSamps) > 0){
    # construct a longggg-format dataframe of posterior samples
    res = data.frame(posteriorSample = c(postSamps, adjustedPostSamps),
                     method = c(rep('unadjusted', length(postSamps)),
                                rep('adjusted', length(adjustedPostSamps))),
                     period_id = c(rep(raw_period_ids, each = numsamps), 
                                   rep(adj_period_ids, each = numsamps))
    )
  }else{
    cat('No samples available!!\n')
    res = NULL
  }
  
  if(!is.null(savePath)){
    fname = sprintf('AllSamples-%s-%s-%s-analysis%s-prior%s.rds',
                    database_id, method, exposure_id, analysis_id, prior_id)
    saveRDS(res, file.path(savePath, fname))
  }
  
  attr(res, 'info') = 
    list(database_id = database_id, 
         method = method,
         exposure_id = exposure_id,
         analysis_id = analysis_id,
         prior_id = prior_id)
  
  return(res)
  
}

# ## try it
# saveSamplePath = '~/Documents/Research/betterGBSanalysesResults/SamplesDataFrame/'
# allSamps = pullPostSamples(database_id = 'CCAE',
#                            method = 'HistoricalComparator',
#                            analysis_id = 2,
#                            exposure_id = 211981,
#                            savePath = saveSamplePath)

# function to plot the posterior distribution by period_id----
# 09/13/2022: try to overlay with prior densities.....
plotGBSPosteriors <- function(allSamps, 
                              adjust = FALSE, 
                              fillColor = 'gray80',
                              markMedian = TRUE,
                              showPlot = TRUE,
                              valueRange = c(-5,10),
                              textSize = 14,
                              logScale = TRUE,
                              showPrior = FALSE,
                              cachePath = './localCache/'){
  # valueRange: the lower and upper limits of estimates (on log scale)
  # logScale: whether or not to keep the log scale on y-axis;
  #           if FALSE, then use original rate ratio scale on y-axis
  # showPrior: whether or not to overlay prior densities on top of the ridge lines...
  #           default FALSE
  
  if(adjust){
    methodText = 'adjusted'
  }else{
    methodText = 'unadjusted'
  }
  
  
  if(logScale){
    xbreaks = seq(from = valueRange[1], to = valueRange[2], length.out = 5)
    xlabels = as.character(round(xbreaks, 1))
    xname = 'Effect size (log rate ratio)'
  }else{
    xbreaks = c(valueRange[1], 0, log(1.5), log(2), log(4), valueRange[2])
    xlabels = as.character(round(exp(xbreaks),1))
    xname = 'Effect size (rate ratio)'
  }
  
  
  dat = allSamps %>% filter(method == methodText) %>% select(-method)
  
  # prior lines
  if(showPrior){
    priors = readRDS(file.path(cachePath, 'priorTable.rds'))
    pr_id = attr(allSamps, 'info')$prior_id
    priorSd = priors %>% filter(prior_id == pr_id) %>% select(Sd) %>% pull()
    periods = unique(dat$period_id)
    priorSamps = rnorm(50000, mean = 0, sd = priorSd)
    priorDat = data.frame(posteriorSample = rep(priorSamps, length(periods)),
                          period_id = rep(periods, each = 10000),
                          showGroup = 'prior')
    dat = bind_rows(dat %>% mutate(showGroup = 'posterior'),
                    priorDat)
  }else{
    dat = dat %>% mutate(showGroup = 'posterior')
  }
  
  if(markMedian){
    medians = dat %>% group_by(period_id) %>%
      summarize(med = median(posteriorSample)) %>%
      ungroup() %>%
      mutate(period_id = as.factor(period_id))
    
    p = ggplot(dat) +
      geom_density_ridges(scale = 0.9, 
                          mapping = aes(y=as.factor(period_id), 
                                        x = posteriorSample,
                                        linetype = showGroup,
                                        fill = showGroup)) +
      geom_vline(xintercept = 0, size = 0.8, color = 'gray60')+ 
      geom_point(data = medians, 
                 mapping = aes(y = period_id, x = med),
                 shape = 4, 
                 size = 1.5,
                 position = position_nudge(y = 0.3)) +
      scale_x_continuous(limits = valueRange,
                         breaks = xbreaks,
                         labels = xlabels) +
      scale_y_discrete(expand = expansion(add = c(0.5, 1)))+
      scale_fill_manual(values = c(fillColor, 'transparent'))+
      labs(y = 'Analysis time (month)', 
           x = xname)+
      coord_flip() +
      theme_bw(base_size = textSize) +
      theme(legend.position = 'none')
    
  }else{
    p = ggplot(dat, 
               aes(y=as.factor(period_id), x = posteriorSample,
                   group = showGroup)) +
      geom_density_ridges(scale = 0.9, fill = fillColor) +
      geom_vline(xintercept = 0, size = 0.8, color = 'gray60')+
      scale_x_continuous(limits = c(-5,10)) +
      scale_y_discrete(expand = expansion(add = c(0.5, 1)))+
      labs(y = 'Analysis time (month)', 
           x = 'Effect size (log relative rate ratio)')+
      coord_flip() +
      theme_bw(base_size = textSize)
  }
  
  if(showPlot){
    print(p)
  }else{
    attr(p,'data') = dat
    return(p)
  }
  
}

# another plot to show P(H1 | data) and P(H0 | data) by period_id
plotPosteriorProbs <- function(allSamps, 
                               adjust = FALSE, 
                               colors = NULL,
                               showPlot = TRUE,
                               xpaddings = c(0,0),
                               textSize = 14,
                               legendPosition = 'bottom'){
  
  if(adjust){
    methodText = 'adjusted'
  }else{
    methodText = 'unadjusted'
  }
  
  dat = allSamps %>% filter(method == methodText) %>%
    group_by(period_id) %>%
    summarize(P1 = mean(posteriorSample > 0),
              P0 = mean(posteriorSample < 0)) %>%
    ungroup()
  
  dat = bind_rows(dat %>% select(period_id, prob = P1) %>% mutate(label = 'P(H1)'),
                  dat %>% select(period_id, prob = P0) %>% mutate(label = 'P(H0)'))
  
  p = ggplot(dat, 
             aes(x=period_id, y = prob, color = label)) +
    geom_line(size = 1)+
    scale_x_continuous(breaks = seq(from = min(dat$period_id), 
                                    to = max(dat$period_id),
                                    by = 1),
                       expand = expansion(add = xpaddings)) +
    scale_y_continuous(limits = c(0,1))+
    labs(y = 'Posterior probability', 
         x = 'Analysis time (month)', 
         color = '')+
    theme_bw(base_size = textSize)+
    theme(legend.position = legendPosition)
  
  if(!is.null(colors)){
    p = p+scale_color_manual(values = colors)
  }
  
  attr(p, 'data') = dat
  
  if(showPlot){
    print(p)
  }else{
    return(p)
  }
  
}

# 09/07/2022
# try to look at the KL divergence between consecutive posteriors
consecutiveKL <- function(allSamps, 
                          adjust,
                          K = 3){
  if(adjust){
    samps = allSamps %>% filter(method == 'adjusted')
  }else{
    samps = allSamps %>% filter(method == 'unadjusted')
  }
  
  periods = sort(unique(samps$period_id))
  if(length(periods) < 2){
    warning('Less than 2 periods of samples available. Cannot proceed!!')
    return()
  }
  
  all_KLs = NULL
  prev.samps = samps %>% filter(period_id == periods[1]) %>%
    select(posteriorSample) %>% pull()
  for(i in 2:length(periods)){
    this.samps = samps %>% filter(period_id == periods[i]) %>%
      select(posteriorSample) %>% pull()
    
    this.KL = KL_div(prev.samps,this.samps,K)
    all_KLs = c(all_KLs, this.KL)
    
    prev.samps = this.samps
  }
  
  data.frame(period = periods[-1], KL = all_KLs)
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
