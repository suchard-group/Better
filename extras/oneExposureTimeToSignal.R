# 05/04/2022: simplied version of `computeEarliestTimeToSignal.R`
# focusing on ONE exposure ONLY!
# can be used for Zoster vaccine results....

library(tidyverse)
library(ggh4x)
library(wesanderson)

source('./extras/helperFunctions.R')
source('./extras/postProcessUtils.R')

## 1. a helper function to produce earliest period that reached a sensitivity level
periodToSensitivity <- function(ttd, dec, sens = 0.8){
  # if all time to decision entries are Inf...
  min_time = min(ttd)
  if(min_time == Inf){return(Inf)}
  
  # start from min_time and go forward and check sensitivity level one by one
  n = length(dec)
  time_points = sort(unique(ttd[ttd < Inf]))
  for(t in time_points){
    sensitivity = sum(dec[ttd<=t] == 'signal')/n
    if(sensitivity >= sens){
      return(t)
      break
    }
  }
  
  # if sensitivity level is not reached, return Inf again
  return(Inf)
}

## 2. main function to compute earliest time to specified sensitivity level (default to 0.8)
##    to the computation and save results first, NO PLOTTING IN HERE
## 05/04/2022: (1) do this for a subset of exposures (mainly, for single exposure)
##             (2) also change "Calibrated" --> "Bias adjusted"; also re-order to put unadjusted up front
selectExposureEarliestTimeToSignal <- function(database_id,
                                               method,
                                               exposure_id,
                                               resPath,
                                               # where to find the summary files (decisions in them)
                                               cachePath,
                                               # where to find helper tables
                                               savePath = NULL,
                                               # path for saving results
                                               returnResults = FALSE,
                                               # if return results as dataframe
                                               saveResults = FALSE,
                                               # if save the summary results to savePath
                                               sensitivity = 0.8,
                                               posControlOnly = TRUE) {
  # if only do this for positive controls
  # first check if savePath is provided when saveResults=TRUE
  if(saveResults & is.null(savePath)){
    stop('Must provide a `savePath` when saveResults=TRUE!!!\n')
  }
  
  # also check if results files are already saved
  # if so, skip and message
  if(!is.null(savePath)){
    fname = sprintf('TimeToSignal-%s-%s-Sensitivity%.2f.rds', database_id, method, sensitivity)
    if(file.exists(file.path(savePath,fname))){
      cat(sprintf('Results already exist at %s, will use available results directly...\n', file.path(savePath,fname)))
      if(returnResults){
        return(readRDS(file.path(savePath,fname)))
      }else{
        return()
      }
    }
  }
  
  # read in file
  fname = sprintf('AllDecisions-%s-%s.rds', database_id, method)
  res = readRDS(file.path(resPath, fname))
  
  # massage HistoricalComparator results a little...
  # RIGHT NOW: disregard results for "filtered" analyses (analysis_id: 13-24)
  if(method == 'HistoricalComparator'){
    res = res %>% filter(analysis_id < 13)
  }
  
  # if only focus on positive control outcomes, then remove negative control outcomes
  if(posControlOnly){
    res = res %>% filter(!negativeControl)
  }
  
  # subset on exposure_ids if not null
  if(!is.null(exposure_id)){
    res = res %>% filter(exposure_id %in% !!exposure_id)
  }
  
  # read in the effect size table and join with decisions
  effects = getEffectSizes(cachePath, cachePath)
  res = res %>% left_join(effects)
  
  # compute time to specified sensitivity level 
  # for unadjusted results and 
  unadjusted = res %>% 
    group_by(database_id, method, analysis_id, 
             threshold_id, prior_id, 
             exposure_id, effect_size) %>%
    summarise(timeToSignal = periodToSensitivity(timeToDecision, 
                                                 decision, 
                                                 sens = sensitivity)) %>%
    mutate(Type = 'Unadjusted')
  
  adjusted = res %>% 
    group_by(database_id, method, analysis_id, 
             threshold_id, prior_id, 
             exposure_id, effect_size) %>%
    summarise(timeToSignal = periodToSensitivity(adjustedTimeToDecision, 
                                                 adjustedDecision, 
                                                 sens = sensitivity)) %>%
    mutate(Type = 'Bias adjusted')
  
  # combine and return/save if necessary
  res = bind_rows(unadjusted, adjusted) %>%
    mutate(Type = factor(Type, levels = c('Unadjusted','Bias adjusted')))
  
  # save results if ...
  if(saveResults){
    if(!dir.exists(savePath)){dir.create(savePath)}
    fname = sprintf('TimeToSignal-%s-%s-Sensitivity%.2f.rds', 
                    database_id, method, sensitivity)
    saveRDS(res, file.path(savePath, fname))
    cat(sprintf('\nTime to signal summary file saved for %s, %s, with sensitivity level %.2f!\nFile can be found at %s\n\n',
                database_id, method, sensitivity, file.path(savePath, fname)))
  }
  
  # return results if
  if(returnResults){
    return(res)
  }
}


## 3. plotting function
##   (1) plot density of time to signal (those finite times)
##   (2) proportion of finite times to signal
## 05/04/2022: focus on ONE exposure only here...
plotEarliestTimeToSignalOneExposure <- function(database_id,
                                                method,
                                                exposure_id,
                                                analysesToExclude = NULL,
                                                # analysis_ids to exclude
                                                plotType = 'timeDensity',
                                                resPath,
                                                # where to find the summary files (decisions in them)
                                                cachePath,
                                                # where to find helper tables
                                                savePath = NULL,
                                                # path for saving plots
                                                saveResults = FALSE,
                                                # if save the summary results to savePath
                                                sensitivity = 0.8,
                                                posControlOnly = TRUE,
                                                baseExposures = TRUE,
                                                pHeight = 8,
                                                pWidth = 12,
                                                usePalette = wes_palette("Darjeeling2")[2:3]) {
  
  # first check if savePath is provided when saveResults=TRUE
  if(saveResults & is.null(savePath)){
    stop('Must provide a `savePath` when saveResults=TRUE!!!\n')
  }
  
  # get the data frame needed for plotting
  tts = selectExposureEarliestTimeToSignal(
    database_id = database_id,
    method = method,
    exposure_id = exposure_id,
    resPath = resPath,
    cachePath = cachePath,
    savePath = savePath,
    saveResults = saveResults,
    returnResults = TRUE,
    sensitivity = sensitivity,
    posControlOnly = posControlOnly
  )
  # massage HistoricalComparator results a little...
  # RIGHT NOW: disregard results for "filtered" analyses (analysis_id: 13-24)
  if(method == 'HistoricalComparator'){
    tts = tts %>% filter(analysis_id < 13)
  }
  
  # remove unwanted analysis_ids
  if(!is.null(analysesToExclude)){
    tts = tts %>% filter(!analysis_id %in% analysesToExclude)
  }
  
  # read thresholds table
  # and a bit of processing
  thresholds = readRDS(file.path(cachePath, 'thresholdTable.rds')) %>%
    select(threshold_id = threshold_ID, 
           d1 = p1Threshold, 
           d0 = p0Threshold) %>%
    #arrange(threshold_id) %>%                    # make sure it's ordered by threshold_id
    mutate(d1Label = sprintf('delta1=%.2f',d1),
           d0Label = sprintf('delta0=%.2f',d0))  # create the character label here too
  
  # also prior table
  priors = readRDS(file.path(cachePath, 'priorTable.rds')) %>%
    arrange(Sd) %>%
    mutate(priorLabel = sprintf('Mean=%s, SD=%s', Mean, Sd)) %>%
    select(-Mean)
  
  # ... and exposures table as well
  exposures = getExposures(NULL, NULL, savepath = cachePath) %>%
    select(exposure_id, exposure_name, base_exposure_name)
  
  # if only stratify by base exposure names...
  if(baseExposures){
    exposures = exposures %>% select(exposure_id, exposure_name = base_exposure_name)
  }
  
  # join with exposure table and thresholds table
  tts = tts %>% 
    left_join(exposures, by='exposure_id') %>%
    left_join(thresholds, by = 'threshold_id') %>%
    select(-d1, -d0) %>%
    left_join(priors, by = 'prior_id')
  
  
  # open up pdf file if saving plots
  if(saveResults){
    plotName = sprintf('%s-%s-timeToSignal-%s.pdf', 
                       database_id, method,
                       plotType)
    pdf(file = file.path(savePath, plotName),
        height = pHeight, width = pWidth)
  }
  
  # plot by plotType
  if(plotType == 'timeDensity'){
    ## for one exposure: combine priors and effect size + delta1 in one big plot
    tts$priorLabel = factor(tts$priorLabel, levels = priors$priorLabel)
    pg = ggplot(tts, aes(y=timeToSignal, 
                         x=as.factor(effect_size), 
                         fill=Type)) +
      geom_violin(width = 1)+
      geom_vline(xintercept = c(1.5,2.5), color='gray40')+
      scale_y_continuous(breaks = seq(from=3, to=12, by=3)) +
      labs(x='Effect Size', 
           y=sprintf('Earliest time to reach %.0f%% sensitivity', sensitivity * 100),
           fill = '') +
      facet_grid(d1Label~priorLabel,
                 labeller = label_wrap_gen(width=15)) +
      theme_bw(base_size = 14)+
      theme(panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            legend.position = 'bottom')
    
    if(!is.null(usePalette)){
      print(
        pg + 
          scale_fill_manual(values = usePalette)
      )
    }
    
  }else{
    # bar plot for rate of finite earliest times
    tts_counts = tts %>% 
      group_by(exposure_name, effect_size, priorLabel, d1Label, Type) %>%
      summarise(finiteRate = mean(timeToSignal < Inf))
    ## fix order of prior labels by re-leveling
    tts_counts$priorLabel = factor(tts_counts$priorLabel, levels = priors$priorLabel)
    # prior.labs <- priors$priorLabel
    # names(prior.labs) <- as.character(priors$Sd)
    pg = 
      ggplot(tts_counts, aes(y=finiteRate, 
                             x=as.factor(effect_size),
                             fill = Type)) +
      geom_bar(stat = 'identity', 
               position = position_dodge()) +
      geom_vline(xintercept = c(1.5,2.5), color='gray40')+
      scale_y_continuous(breaks = c(0.1,0.5,1)) +
      coord_flip()+
      labs(x='Effect Size', 
           y=sprintf('Rate of reaching %.0f%% sensitivity before end', 
                     sensitivity * 100),
           fill='') +
      facet_grid(d1Label ~ priorLabel, 
                 labeller = label_wrap_gen(width=15))+
      theme_bw(base_size = 14)+
      theme(panel.grid.major.x = element_blank(),
            #axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_blank(),
            #panel.border = element_blank(),
            legend.position = 'bottom')
    if(!is.null(usePalette)){
      print(
        pg + 
          scale_fill_manual(values = usePalette)
      )
    }
  }
  
  if(saveResults){dev.off()}
  
}

resultspath = '~/Documents/Research/betterResults/summary/'
cachepath = './localCache/'
savepath = '~/Documents/Research/betterResults/timeToSignalZoster/'

sensitivity_level = 0.5

db = 'CCAE'
mt = 'HistoricalComparator'

eid = c(211981:211983)

if(mt == 'SCCS'){
  analysesExclude = c(9:12, 15)
}else{
  analysesExclude = c(9:12)
}

## (1) time density
plotEarliestTimeToSignalOneExposure(database_id = db, 
                         method = mt,
                         exposure_id = eid,
                         analysesToExclude = analysesExclude,
                         plotType = 'timeDensity',
                         resPath = resultspath,
                         cachePath = cachepath,
                         savePath = savepath,
                         saveResults = TRUE,
                         sensitivity = sensitivity_level,
                         posControlOnly = TRUE,
                         baseExposures = TRUE,
                         pHeight = 6, pWidth = 9,
                         usePalette = wes_palette("Darjeeling2")[c(2,4)])

## (2) rate of finite times
plotEarliestTimeToSignalOneExposure(database_id = db, 
                         method = mt,
                         exposure_id = eid,
                         analysesToExclude = analysesExclude,
                         plotType = 'finiteRates',
                         resPath = resultspath,
                         cachePath = cachepath,
                         savePath = savepath,
                         saveResults = TRUE,
                         sensitivity = sensitivity_level,
                         posControlOnly = TRUE,
                         baseExposures = TRUE,
                         pHeight = 6, pWidth = 6,
                         usePalette = wes_palette("Darjeeling2")[c(2,4)])
