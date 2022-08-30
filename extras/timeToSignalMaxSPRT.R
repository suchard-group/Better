# 08/29/2022
# time to signal for MaxSPRT results
# (using the local saves by Martijn)

library(dplyr)
library(wesanderson)

# 1. function to make decisions and obtain earliest time to decision
makeDecisionsFreq <- function(estimates, 
                              method = 'HistoricalComparator',
                              exposure_id = 211981,
                              posOnly = TRUE){
  estimates = estimates %>%
    select(method:outcomeId, llr, criticalValue, calibratedLlr, effectSize) %>%
    filter(!is.na(llr), !is.na(calibratedLlr))
  
  if(posOnly){
    estimates = estimates %>% filter(effectSize > 1)
  }
  
  decisions = estimates %>%
    group_by(databaseId, method, analysisId, outcomeId, effectSize) %>%
    arrange(periodId) %>%
    mutate(signal = (llr > criticalValue),
           calibratedSignal = (calibratedLlr > criticalValue)) %>%
    summarise(timeToSignal = if_else(sum(signal) > 0, 
                                     which.min(signal == TRUE), Inf),
              calibratedTimeToSignal = if_else(sum(calibratedSignal) > 0, 
                                               which.min(calibratedSignal == TRUE), 
                                               Inf))
  
  return(decisions)
}

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

# 3. plotting function
plotTimeToSignalMaxSPRTOneExposure <- function(localFile,
                                               method,
                                               sensitivity_level, 
                                               colors = NULL,
                                               cachePath = './localCache',
                                               ...){
  estimates = readRDS(localFile)
  decisions = makeDecisionsFreq(estimates, method = method, ...)
  tts = timeToSensitivityFreq(decisions, sens = sensitivity_level) %>%
    select(analysis_id = analysisId, effect_size = effectSize)
  
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
  
  # analysis names
  
  
  # plotting
  pg = ggplot(tts_distribution, 
              aes(y=as.factor(effect_size), 
                  x=medianTime, 
                  fill=Type)) +
    geom_bar(stat='identity', position = 'dodge')+
    geom_errorbar(aes(xmax = ub, xmin = lb, color = Type),
                  position = position_dodge(width=1), 
                  width = 0.6, size = 0.8)+
    geom_hline(yintercept = c(1.5,2.5), color='gray40')+
    scale_x_continuous(breaks = seq(from=0, to=12, by=3)) +
    labs(y='Effect Size', 
         x=sprintf('Earliest time to reach %.0f%% sensitivity', sensitivity * 100),
         fill = '', color = '') +
    facet_grid(d1Label~priorLabel,
               labeller = label_wrap_gen(width=15)) +
    theme_bw(base_size = 14)+
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position = 'bottom') +
    scale_color_manual(values = c('gray20','gray20'), guide='none')
  
  
}

