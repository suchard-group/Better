# March 2: compute earliest time to signal (for positive control only)
#          (stratified by effect size)
#          time to claim signal for at least 80% (or some other number) of the positive controls

# March 7: arrange priors by order of SD
#          & change palette

# March 8: try flipping coords for the finite signal time bar plot

library(ggh4x)

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
computeEarliestTimeToSignal <- function(database_id,
                                        method,
                                        resPath,          # where to find the summary files (decisions in them)
                                        cachePath,        # where to find helper tables
                                        savePath = NULL, # path for saving results
                                        returnResults = FALSE, # if return results as dataframe
                                        saveResults = FALSE, # if save the summary results to savePath
                                        sensitivity = 0.8,
                                        posControlOnly = TRUE){ # if only do this for positive controls
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
  res = bind_rows(unadjusted, adjusted)
  
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
plotEarliestTimeToSignal <- function(database_id,
                                     method,
                                     analysesToExclude = NULL, # analysis_ids to exclude
                                     plotType = 'timeDensity',
                                     resPath,          # where to find the summary files (decisions in them)
                                     cachePath,        # where to find helper tables
                                     savePath = NULL, # path for saving plots
                                     saveResults = FALSE, # if save the summary results to savePath
                                     sensitivity = 0.8,
                                     posControlOnly = TRUE,
                                     baseExposures = TRUE,
                                     pHeight = 8, pWidth = 12,
                                     usePalette = wes_palette("Darjeeling2")[2:3]){
  # first check if savePath is provided when saveResults=TRUE
  if(saveResults & is.null(savePath)){
    stop('Must provide a `savePath` when saveResults=TRUE!!!\n')
  }
  
  # get the data frame needed for plotting
  tts = computeEarliestTimeToSignal(database_id = database_id,
                                    method = method,
                                    resPath = resPath,
                                    cachePath = cachePath,
                                    savePath = savePath,
                                    saveResults = saveResults,
                                    returnResults = TRUE,
                                    sensitivity = sensitivity,
                                    posControlOnly = posControlOnly)
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
    ## violin plot for the earliest periods
    ## one plot for each prior
    pls = priors$priorLabel
    for(pl in pls){
      this.tts = tts %>% filter(priorLabel == pl)
      pg = ggplot(this.tts, aes(y=timeToSignal, 
                                x=as.factor(effect_size), 
                                fill=Type)) +
        geom_violin(width = 1)+
        geom_vline(xintercept = c(1.5,2.5), color='gray40')+
        scale_y_continuous(breaks = seq(from=3, to=12, by=3)) +
        labs(x='Effect Size', 
             y=sprintf('Earliest time to reach %.0f%% sensitivity', sensitivity * 100),
             fill = '',
             caption = pl) +
        facet_grid(d1Label~exposure_name,
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
      
    }
  }else if(plotType == 'timeBars'){
    # 08/05/2022: compromise --- if decision time > 12, substitute with 12 for now!!!
    tts_distribution = tts %>% 
      mutate(timeToSignal = if_else(timeToSignal==Inf, 12, timeToSignal)) %>%
      group_by(exposure_name, effect_size, priorLabel, d1Label, Type) %>%
      summarise(medianTime = median(timeToSignal),
                lb = quantile(timeToSignal, 0.25),
                ub = quantile(timeToSignal, 0.75)) %>%
      ungroup()
    
    tts_distribution$priorLabel = factor(tts_distribution$priorLabel, 
                                         levels = priors$priorLabel)
    
    ## switch order of approach for better presentation
    tts_distribution$Type = factor(tts_distribution$Type, 
                                   levels = c('Bias adjusted','Unadjusted'))
    
    # if(addCommentary){
    #   tts_distribution = tts_distribution %>% 
    #     left_join(priors, by = 'priorLabel') %>%
    #     select(-Sd, prior_id) %>%
    #     mutate(priorLabel = factor(commentary, 
    #                                levels = priors$commentary)) %>%
    #     left_join(thresholds, by = 'd1Label') %>%
    #     mutate(d1Label = factor(d1comment,
    #                             levels = unique(thresholds$d1comment)))
    # }
    
    # # 09/13/2022: add error rates info
    # if(!is.null(errorRatePath)){
    #   #print(names(tts_distribution))
    #   tts_distribution = tts_distribution %>% 
    #     left_join(errorRates, by = c('exposure_id', 'Type', 'prior_id', 'threshold_id')) %>%
    #     # mutate(shade = case_when(avgType1Error <= 0.1 ~ 0,
    #     #                          avgType1Error <= 0.25 ~ 0.3,
    #     #                          avgType1Error <= 0.35 ~ 0.7,
    #     #                          TRUE ~ 0.95))
    #     mutate(shade = 1-avgType1Error)
    #   shadeRange = c(0.05, 0.55)
    #   capt = 'Color shade by specificity'
    # }else{
    #   tts_distribution$shade = 0
    #   shadeRange = c(0.1, 1)
    #   capt = ''
    # }
    
    pg = ggplot(tts_distribution, 
                aes(y=as.factor(effect_size), 
                    x=medianTime, 
                    fill=Type)) +
      geom_bar(stat='identity', position = 'dodge')+
      geom_errorbar(aes(xmax = ub, xmin = lb, color = Type),
                    position = position_dodge(width=1), 
                    width = 0.6, size = 0.8)+
      geom_hline(yintercept = c(1.5,2.5), color='gray40')+
      scale_x_continuous(breaks = seq(from=6, to=12, by=6)) +
      #scale_alpha_continuous(range = shadeRange, guide = 'none')+
      labs(y='Effect Size', 
           x=sprintf('Earliest time to reach %.0f%% sensitivity', sensitivity * 100),
           fill = '', color = '') +
      facet_nested(priorLabel + d1Label ~ exposure_name,
                   labeller = label_wrap_gen(width=15),
                   nest_line = element_line(linetype = 1))+
      theme_bw(base_size = 14)+
      theme(panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            #strip.text.y = element_text(angle = 0),
            legend.position = 'bottom') +
      scale_color_manual(values = c('gray20','gray20'), guide='none')
    
    if(!is.null(usePalette)){
      
      pg = pg + 
        scale_fill_manual(values = usePalette,
                          breaks = c('Unadjusted', 'Bias adjusted'))
      
      print(pg)
    }
    
    #attr(pg, 'data') = tts_distribution
    #attr(pg, 'errorRates') = errorRates
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
        facet_nested(priorLabel + d1Label ~ exposure_name,
                     labeller = label_wrap_gen(width=15),
                     nest_line = element_line(linetype = 1)) +
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




#### Try running using the functions --------
resultspath = '~/Documents/Research/betterResults/summary/'
cachepath = './localCache/'
savepath = '~/Documents/Research/betterResults/timeToSignal/'
#plotpath = '~/Documents/Research/betterResults/plots/'

sensitivity_level = 0.5

databases = c('IBM_MDCD','CCAE','OptumEhr','OptumDod','MDCR')
#databases = 'IBM_MDCD'
methods = c('SCCS','HistoricalComparator')

for(db in databases){
  for(mt in methods){
    # ## generate and save summary only...
    # computeEarliestTimeToSignal(database_id = db, 
    #                             method = mt,
    #                             resPath = resultspath,
    #                             cachePath = cachepath,
    #                             savePath = savepath,
    #                             returnResults = FALSE,
    #                             saveResults = TRUE,
    #                             sensitivity = 0.5,
    #                             posControlOnly = TRUE)
    
    ## actually do the plotting..
    if(mt == 'SCCS'){
      analysesExclude = c(9:12, 15)
    }else{
      analysesExclude = c(9:12)
    }
    
    ## (1) time density
    plotEarliestTimeToSignal(database_id = db, 
                             method = mt,
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
                             usePalette = wes_palette("Darjeeling2")[2:3])
    
    ## (2) rate of finite times
    plotEarliestTimeToSignal(database_id = db, 
                             method = mt,
                             analysesToExclude = analysesExclude,
                             plotType = 'finiteRates',
                             resPath = resultspath,
                             cachePath = cachepath,
                             savePath = savepath,
                             saveResults = TRUE,
                             sensitivity = sensitivity_level,
                             posControlOnly = TRUE,
                             baseExposures = TRUE,
                             pHeight = 12, pWidth = 8,
                             usePalette = wes_palette("Darjeeling2")[2:3])
    
  }
}



#### TESTING CODE BELOW; DON'T RUN --------

# db = 'IBM_MDCD'
# mt = "SCCS"
# if(mt == 'SCCS'){
#   analysesExclude = c(9:12, 15)
# }else{
#   analysesExclude = c(9:12)
# }
# 
# 
# plotEarliestTimeToSignal(database_id = db, 
#                          method = mt,
#                          analysesToExclude = analysesExclude,
#                          plotType = 'timeDensity',
#                          resPath = resultspath,
#                          cachePath = cachepath,
#                          savePath = savepath,
#                          saveResults = TRUE,
#                          sensitivity = 0.5,
#                          posControlOnly = TRUE,
#                          baseExposures = TRUE,
#                          pHeight = 6, pWidth = 9)


# db = 'IBM_MDCD'
# mt = "SCCS"
# 
# etts = computeEarliestTimeToSignal(database_id = db, 
#                                    method = 'SCCS',
#                                    resPath = resultspath,
#                                    cachePath = cachepath,
#                                    savePath = savepath,
#                                    returnResults = FALSE,
#                                    saveResults = TRUE,
#                                    sensitivity = 0.5,
#                                    posControlOnly = TRUE)

savepath = '~/Downloads/'

sensitivity_level = 0.5

db = 'CCAE'
mt = 'SCCS'
analysesToExclude = c(9:12, 15)

plotEarliestTimeToSignal(database_id = db,
                         method = mt, 
                         analysesToExclude = analysesToExclude,
                         plotType = 'timeBars',
                         resPath = resultspath,
                         cachePath = cachepath,
                         savePath = savepath,
                         saveResults = TRUE,
                         sensitivity = sensitivity_level,
                         posControlOnly = TRUE,
                         baseExposures = TRUE,
                         pHeight = 12, pWidth = 8,
                         usePalette = wes_palette("Darjeeling2")[c(2,4)])

