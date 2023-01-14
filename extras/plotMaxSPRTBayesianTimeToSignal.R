# 01/13/2023
# for presentation purposes

# plot time-to-signal for both MaxSPRT and Bayesian methods


# 01/13/2023: start with CCAE database here!
db = 'CCAE'
method = 'HistoricalComparator'
aids = 1:8

eid = 211983 # two-doses of Zoster vaccine

sens = 0.25 # desired sensitivity level

## local estimates save for MaxSPRT
maxSPRT_filepath = './localCache/EstimateswithImputedPcs_CCAE.rds'


# 1. ------
# load functions needed for MaxSPRT time-to-signal
source('./extras/timeToSignalMaxSPRT.R')

# functions for Bayesian time-to-signal
source('./extras/oneExposureTimeToSignal_clean.R')

# load functions needed for MaxSPRT and Bayesian error rates
source('./extras/simpleCalibration.R')
source('./extras/frequentistDecisionComparisons.R')


theColors = wes_palette("Darjeeling2")[c(2,4)]

# 2. -----
## get time-to-signal dataframe for MaxSPRT
p1 = plotTimeToSignalMaxSPRTOneExposure(localFile = maxSPRT_filepath,
                                        method = method,
                                        exposure_id = eid,
                                        analysis_ids = aids,
                                        sensitivity_level = sens,
                                        colors = theColors,
                                        #plotType = 'separate')
                                        plotType = 'together',
                                        showPlot = FALSE)

## cross-ref the MaxSPRT TTS dataframe with Type 1 error dataframe to get Type 1 error rate
maxSPRT_tts = pullErrorRatesForTTS(connection, p1,
                                   database_id = db,
                                   method = method, 
                                   exposure_id = eid, 
                                   localEstimatesPath = maxSPRT_filepath,
                                   calibration = FALSE)


## get the time-to-signal dataframe for Bayesian method
bayes_tts_vanilla = foreach(aid = aids, .combine = rbind) %do% {
  
  ## get the "empirical alpha"
  this.alpha = maxSPRT_tts %>% filter(analysis_id == aid) %>%
    slice(1) %>% select(type1error) %>% pull()
  
  if(length(this.alpha) == 1 && is.double(this.alpha)){
    ## use it to "calibrate" Bayesian threshold to produce TTS results
    tts_bayes_raw = getTTSbyEmpiricalAlpha(database_id = db,
                                           method = method,
                                           analysis_id = aid,
                                           exposure_id = eid,
                                           prior_ids = c(1:3), # include all priors for easier query later
                                           alpha = this.alpha,
                                           summaryPath = summarypath,
                                           cachePath = cachepath,
                                           useAdjusted = FALSE,
                                           calibrate = TRUE,
                                           sensitivity_level = sens,
                                           adjust_max_time = TRUE)
    tts_bayes_raw
  }else{
    NULL
  }
  
}

bayes_tts_adjusted = foreach(aid = aids, .combine = rbind) %do% {
  
  ## get the "empirical alpha"
  this.alpha = maxSPRT_tts %>% filter(analysis_id == aid) %>%
    slice(1) %>% select(type1error) %>% pull()
  
  if(length(this.alpha) == 1 && is.double(this.alpha)){
    ## use it to "calibrate" Bayesian threshold to produce TTS results
    tts_bayes_raw = getTTSbyEmpiricalAlpha(database_id = db,
                                           method = method,
                                           analysis_id = aid,
                                           exposure_id = eid,
                                           prior_ids = c(1:3), # include all priors for easier query later
                                           alpha = this.alpha,
                                           summaryPath = summarypath,
                                           cachePath = cachepath,
                                           useAdjusted = TRUE,
                                           calibrate = TRUE,
                                           sensitivity_level = sens,
                                           adjust_max_time = TRUE)
    tts_bayes_raw
  }else{
    NULL
  }
  
}

# try an example plot
## ADJUSTED BAYES
prior_SD = 4

tts_long = rbind(
  maxSPRT_tts %>% 
    select(database_id, exposure_id, method, analysis_id, effect_size, time_to_sens, approach),
  bayes_tts_adjusted %>%
    filter(prior_sd == prior_SD) %>%
    select(database_id, exposure_id, method, analysis_id, effect_size, time_to_sens, approach)
)

tts_long$approach = factor(tts_long$approach,
                           levels = c('MaxSPRT', 'Adjusted Bayes'))

(
  pg = ggplot(tts_long, 
              aes(y=as.factor(effect_size), 
                  x=time_to_sens, 
                  fill=approach)) +
    geom_bar(stat='identity', position = 'dodge')+
    geom_hline(yintercept = c(1.5,2.5), color='gray40')+
    scale_x_continuous(breaks = seq(from=0, to=12, by=3)) +
    labs(y='Effect Size', 
         x=sprintf('Time to %.0f%% sensitivity', sens * 100),
         fill = '') +
    facet_grid(analysis_id~.) +
    theme_bw(base_size = 14)+
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position = 'bottom')+
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    scale_fill_manual(values = theColors,
                      breaks = c('Adjusted Bayes','MaxSPRT'),
                      labels = c('Bayesian','MaxSPRT'))
)
print(pg)
