# 08/29/2022
# time to signal for MaxSPRT results
# (using the local saves by Martijn)

library(dplyr)
library(wesanderson)
library(ggplot2)

# 1. function to make decisions and obtain earliest time to decision
makeDecisionsFreq <- function(estimates, 
                              method = 'HistoricalComparator',
                              analysis_ids = 1:4,
                              exposure_id = 211983,
                              posOnly = TRUE){
  estimates = estimates %>%
    select(method:outcomeId, llr, criticalValue, calibratedLlr, effectSize) %>%
    filter(method == !!method,
           exposureId == !!exposure_id,
           analysisId %in% analysis_ids, 
           !is.na(llr), !is.na(calibratedLlr))
  
  if(posOnly){
    estimates = estimates %>% filter(effectSize > 1)
  }
  
  decisions = estimates %>%
    group_by(databaseId, method, analysisId, outcomeId, effectSize) %>%
    arrange(periodId) %>%
    mutate(signal = (llr > criticalValue),
           calibratedSignal = (calibratedLlr > criticalValue)) %>%
    summarise(timeToSignal = min(which(signal)),
              calibratedTimeToSignal = min(which(calibratedSignal))) %>%
    ungroup()
  
  return(decisions)
}

# ## test it ----
# estPath = './localCache/EstimateswithImputedPcs_CCAE.rds'
# estimates = readRDS(estPath)
# decisions = makeDecisionsFreq(estimates)

# 2. time to ?? sensitivity

## (a) helper function to get time to at least % sensitivity from a vector of time-to-signal's
timeToSensVec <- function(timeToSignals, sens){
  ts = sort(timeToSignals)
  quantile(ts, sens, names = FALSE)
}

## (b) function to return tts
timeToSensitivityFreq <- function(decisions, sens = 0.5){
  
  tts = decisions %>% group_by(analysisId, effectSize) %>%
    summarise(timeToSens = timeToSensVec(timeToSignal, sens = sens),
              calibratedTimeToSens = timeToSensVec(calibratedTimeToSignal, sens = sens)) %>%
    ungroup()
  
  return(tts)
}

# ## test it ----
# tts = timeToSensitivityFreq(decisions, 0.25)

# 3. plotting function
plotTimeToSignalMaxSPRTOneExposure <- function(localFile,
                                               method,
                                               sensitivity_level, 
                                               colors = NULL,
                                               analysis_labels = FALSE,
                                               maxTime = 12,
                                               cachePath = './localCache',
                                               plotType = 'separate',
                                               showPlot = TRUE,
                                               ...){
  estimates = readRDS(localFile)
  decisions = makeDecisionsFreq(estimates, method = method, ...)
  tts = timeToSensitivityFreq(decisions, sens = sensitivity_level)
  
  # long format
  tts_long = bind_rows(tts %>% 
                         select(analysis_id = analysisId, 
                                effect_size = effectSize, 
                                time_to_sens = timeToSens) %>%
                         mutate(Type = 'uncalibrated'),
                       tts %>% select(analysis_id = analysisId, 
                                      effect_size = effectSize, 
                                      time_to_sens = calibratedTimeToSens) %>%
                         mutate(Type = 'calibrated'))
  
  ## set the Inf times to maxTime
  tts_long = tts_long %>% 
    mutate(time_to_sens = if_else(time_to_sens < Inf, time_to_sens, maxTime))
  
  # analysis names
  if(analysis_labels){
    # cross reference analysis description
    analyses = readRDS(file.path(cachePath, 'analyses.rds')) %>%
      filter(method == !!method) %>% 
      select(analysis_id, description, time_at_risk)
    
    ## manual fix for SCCS 12 TaR
    if(method == 'SCCS'){
      analyses[analyses$analysis_id == 12,]$time_at_risk = '0-1'
    }
    
    ## generate description text string
    analyses = analyses %>%
      mutate(analysis_text = sprintf('%s, TaR %s days', description, time_at_risk))
    
    # join with data
    tts_long = tts_long %>% left_join(analyses, by = 'analysis_id') %>%
      arrange(analysis_id)
    
    # panel labels for analysis names
    analysis_labs =  unique(tts_long$analysis_text)
    names(analysis_labs) = as.character(unique(tts_long$analysis_id))
  }else{
    analysis_labs =  as.character(unique(tts_long$analysis_id))
    names(analysis_labs) = as.character(unique(tts_long$analysis_id))
  }
  
  
  # plotting
  if(plotType == 'separate'){
    # time bars separate for each analysis id
    pg = ggplot(tts_long, 
                aes(y=as.factor(effect_size), 
                    x=time_to_sens, 
                    fill=Type)) +
      geom_bar(stat='identity', position = 'dodge')+
      geom_hline(yintercept = c(1.5,2.5), color='gray40')+
      scale_x_continuous(breaks = seq(from=0, to=12, by=3)) +
      labs(y='Effect Size', 
           x=sprintf('Time to %.0f%% sensitivity', sensitivity_level * 100),
           fill = '') +
      facet_grid(analysis_id~.,
                 labeller = labeller(analysis_id = analysis_labs)) +
      theme_bw(base_size = 14)+
      theme(panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = 'bottom')+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
    
  }else{
    # time distributions over all analysis
    tts_distribution = tts_long %>%
      group_by(effect_size, Type) %>%
      summarise(medianTime = median(time_to_sens),
                lb = quantile(time_to_sens, 0.25),
                ub = quantile(time_to_sens, 0.75)) %>%
      ungroup()
    
    pg = ggplot(tts_distribution, 
                aes(y=as.factor(effect_size), 
                    x=medianTime, 
                    fill=Type)) +
      geom_bar(stat='identity', position = 'dodge')+
      geom_errorbar(aes(xmax = ub, xmin = lb),
                    color = 'black',
                    position = position_dodge(width=1), 
                    width = 0.6, size = 0.8)+
      geom_hline(yintercept = c(1.5,2.5), color='gray40')+
      scale_x_continuous(breaks = seq(from=0, to=12, by=3)) +
      labs(y='Effect Size', 
           x=sprintf('Time to %.0f%% sensitivity', sensitivity_level * 100),
           fill = '') +
      # facet_grid(analysis_id~.,
      #            labeller = labeller(analysis_id = analysis_labs)) +
      theme_bw(base_size = 14)+
      theme(panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            # strip.background = element_blank(),
            # strip.text.y = element_text(angle = 0),
            legend.position = 'bottom')+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
    
  }
  
  if(!is.null(colors)){
    pg = pg + 
        scale_fill_manual(values = colors,
                          breaks = c('uncalibrated', 'calibrated'),
                          labels = c('uncalibrated MaxSPRT',
                                     'calibrated MaxSPRT'))
  }
  
  attr(pg, 'data') = tts_long %>% 
    mutate(method = method)
  
  if(showPlot){
    print(pg)
  }
  
  return(pg)
  
  
}


## plots----

# p1 = plotTimeToSignalMaxSPRTOneExposure(localFile = './localCache/EstimateswithImputedPcs_CCAE.rds',
#                                         method = 'HistoricalComparator',
#                                         exposure_id = 211983,
#                                         analysis_ids = 1:8,
#                                         sensitivity_level = 0.25,
#                                         colors = wes_palette("Darjeeling2")[c(2,4)],
#                                         #plotType = 'separate')
#                                         plotType = 'together')


# p1 = plotTimeToSignalMaxSPRTOneExposure(localFile = './localCache/EstimateswithImputedPcs_CCAE.rds',
#                                         method = 'HistoricalComparator',
#                                         analysis_ids = 1:4,
#                                         sensitivity_level = 0.25,
#                                         colors = wes_palette("Darjeeling2")[c(2,4)])


# 01/13/2023
# 3. function to cross-ref with MaxSPRT error rates
pullErrorRatesForTTS <- function(connection,
                                 tts_plot,
                                 database_id,
                                 method,
                                 exposure_id,
                                 localEstimatesPath = NULL,
                                 calibration = FALSE,
                                 cachePath = './localCache/'){
  
  # tts dataframe
  tts = attr(tts_plot, 'data') %>%
    mutate(exposure_id = exposure_id, database_id = database_id)
  
  # get error rates dataframe
  aids = unique(tts$analysis_id)
  
  if(!is.null(localEstimatesPath)){
    localEstimates = readRDS(localEstimatesPath)
  }else{
    localEstimates = NULL
  }
  
  end_type1s = 
    foreach(aid = aids, .combine = rbind) %do% 
    {
      resLst = frequentistDecisions(connection,
                                    'eumaeus',
                                    database_id = database_id,
                                    method = method,
                                    exposure_id = exposure_id,
                                    analysis_id = aid,
                                    calibration = calibration,
                                    estimates = localEstimates,
                                    correct_shift = FALSE,
                                    bonferroni_baseline = FALSE,
                                    FDR = FALSE,
                                    cachePath = cachePath)
      resLst$errorRate %>%
        filter(stats == 'type 1', period_id == max(period_id)) %>%
        select(database_id:exposure_id, period_id, errorRate)
    }
  
  # join two tables
  tts_withError = tts %>% left_join(end_type1s) %>%
    select(analysis_id:database_id, type1error = errorRate) %>%
    mutate(approach = if_else(calibration, 'calibrated MaxSPRT', 'MaxSPRT'))
    
  return(tts_withError)
}

# # test it
# ConnectionDetails <- DatabaseConnector::createConnectionDetails(
#   dbms = "postgresql",
#   server = paste(keyring::key_get("eumaeusServer"),
#                  keyring::key_get("eumaeusDatabase"),
#                  sep = "/"),
#   user = keyring::key_get("eumaeusUser"),
#   password = keyring::key_get("eumaeusPassword"))
# 
# # set up the DB connection
# connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)
# 
# tts_error = pullErrorRatesForTTS(connection,
#                                  p1,
#                                  database_id = 'CCAE',
#                                  method = 'HistoricalComparator',
#                                  exposure_id = 211983,
#                                  localEstimatesPath = './localCache/EstimateswithImputedPcs_CCAE.rds',
#                                  cachePath = './localCache/')
# 
# DatabaseConnector::disconnect(connection)
