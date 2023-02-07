# 01/13/2023
# for presentation purposes

# plot time-to-signal for both MaxSPRT and Bayesian methods
library(foreach)

# 01/13/2023: start with CCAE database here!
db = 'CCAE'
method = 'HistoricalComparator'
aids = 1:8

eid = 211983 # two-doses of Zoster vaccine

sens = 0.25 # desired sensitivity level

## local estimates save for MaxSPRT
maxSPRT_filepath = './localCache/EstimateswithImputedPcs_CCAE.rds'

## other basic setup
library(wesanderson)
theColors = wes_palette("Darjeeling2")[c(2,4)]

# connection to public db server for Eumaeus results
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))
# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)


# 1. ------
# load functions needed for MaxSPRT time-to-signal
source('./extras/timeToSignalMaxSPRT.R')

# functions for Bayesian time-to-signal
source('./extras/oneExposureTimeToSignal_clean.R')

# load functions needed for MaxSPRT and Bayesian error rates
source('./extras/simpleCalibration.R')
source('./extras/frequentistDecisionComparisons.R')


# 2.(a) ---
# sub-function to make TTS plot given maxSPRT and Bayesian results tables
ttsComparePlotFromRes <- function(maxSPRT_tts, 
                                  bayes_tts,
                                  prior_SD = 4,
                                  sensitivity_level = 0.25,
                                  bayesAdjusted = TRUE,
                                  showPlot = TRUE,
                                  colors = theColors,
                                  cachePath = './localCache/'){
  tts_long = rbind(
    maxSPRT_tts %>% 
      select(database_id, exposure_id, method, analysis_id, effect_size, time_to_sens, approach),
    bayes_tts %>%
      filter(prior_sd == prior_SD) %>%
      select(database_id, exposure_id, method, analysis_id, effect_size, time_to_sens, approach)
  )
  
  bayesApproachLabel = ifelse(bayesAdjusted, 'Adjusted Bayes', 'Unadjusted Bayes')
  
  tts_long$approach = factor(tts_long$approach,
                             levels = c('MaxSPRT', bayesApproachLabel))
  
  ## analysis name pull
  analyses = readRDS(file.path(cachePath, 'analyses.rds'))
  tts_long = tts_long %>% 
    left_join(analyses, by = c('method', 'analysis_id')) %>%
    rename('analysis' = 'description')
  
  pg = ggplot(tts_long, 
              aes(y=as.factor(effect_size), 
                  x=time_to_sens, 
                  fill=approach)) +
    geom_bar(stat='identity', position = 'dodge')+
    geom_hline(yintercept = c(1.5,2.5), color='gray40')+
    scale_x_continuous(breaks = seq(from=0, to=12, by=3),
                       limits = c(0,12)) +
    labs(y='Effect Size', 
         x=sprintf('Time to %.0f%% sensitivity', 
                   sensitivity_level * 100),
         fill = '') +
    #facet_grid(analysis_id~.) +
    facet_grid(analysis~.,
               labeller = label_wrap_gen(width=25)) +
    theme_bw(base_size = 14)+
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position = 'bottom')+
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    scale_fill_manual(values = colors,
                      breaks = c(bayesApproachLabel,'MaxSPRT'),
                      labels = c('Bayesian','MaxSPRT'))
  
  if(showPlot){
    print(pg)
  }
  
  # attach data
  attr(pg, 'data') = tts_long
  
  return(pg)
}

# # ## try this
# p_tts = ttsComparePlotFromRes(maxSPRT_tts,
#                               bayes_tts_adjusted,
#                               prior_SD = 4,
#                               bayesAdjusted = TRUE,
#                               showPlot = TRUE)


# 2.(b) ---
# big function to do everything: generate TTS results AND make the plots AND add caption
# ONLY WORKS FOR `CCAE` FOR NOW!!!
# (but can change local file path to other database results)
plotMaxSPRTBayesianTimeToSignal <- function(database_id,
                                            exposure_id,
                                            method,
                                            analysis_ids,
                                            sensitivity_level,
                                            prior_SD = 4,
                                            alpha = .05, # default alpha
                                            calibrateToAlpha = TRUE, # if calibrate to MaxSPRT alpha level?
                                            bayesAdjusted = TRUE, # use bias adjusted Bayesian results?
                                            calibration = FALSE, # use calibrated MaxSPRT?
                                            summaryPath = summarypath, # summary results for Bayesian 
                                            localEstimatesPath = maxSPRT_filepath, # local estimate for MaxSPRT
                                            cachePath = './localCache/', # cache for caption info etc.
                                            colors = theColors, # fill colors for time bars
                                            showPlot = TRUE,
                                            addCaption = TRUE){
  # (1) get time-to-signal dataframe for MaxSPRT
  p1 = plotTimeToSignalMaxSPRTOneExposure(localFile = localEstimatesPath,
                                          method = method,
                                          exposure_id = exposure_id,
                                          analysis_ids = analysis_ids,
                                          sensitivity_level = sensitivity_level,
                                          colors = colors,
                                          plotType = 'separate',
                                          showPlot = FALSE)
  
  ## cross-ref the MaxSPRT TTS dataframe with Type 1 error dataframe to get Type 1 error rate
  maxSPRT_tts = pullErrorRatesForTTS(connection, 
                                     p1,
                                     database_id = database_id,
                                     method = method, 
                                     exposure_id = exposure_id, 
                                     localEstimatesPath = localEstimatesPath,
                                     calibration = calibration)
  
  # (2) Bayesian TTS
  ## only focus on needed prior to save time...
  priors = readRDS(file.path(cachePath, 'priorTable.rds'))
  this.prior_id = priors %>% filter(Sd == prior_SD) %>% select(prior_id) %>% pull()
  
  bayes_tts = foreach(aid = analysis_ids, .combine = rbind) %do% {
    
    ## get the "empirical alpha"
    this.alpha = maxSPRT_tts %>% filter(analysis_id == aid) %>%
      slice(1) %>% select(type1error) %>% pull()
    
    if(length(this.alpha) == 1 && is.double(this.alpha)){
      ## use it to "calibrate" Bayesian threshold to produce TTS results
      tts_bayes_raw = getTTSbyEmpiricalAlpha(database_id = database_id,
                                             method = method,
                                             analysis_id = aid,
                                             exposure_id = exposure_id,
                                             prior_ids = this.prior_id, # include all priors for easier query later
                                             alpha = this.alpha,
                                             summaryPath = summaryPath,
                                             cachePath = cachePath,
                                             useAdjusted = bayesAdjusted,
                                             calibrate = calibrateToAlpha,
                                             sensitivity_level = sensitivity_level,
                                             adjust_max_time = TRUE)
      tts_bayes_raw
    }else{
      NULL
    }
    
  }
  
  # (3) combine and generate plot
  p_tts = ttsComparePlotFromRes(maxSPRT_tts,
                                bayes_tts_adjusted,
                                prior_SD = prior_SD,
                                sensitivity_level = sensitivity_level,
                                bayesAdjusted = bayesAdjusted,
                                showPlot = FALSE,
                                cachePath = cachePath)
  
  # (4) get caption
  exposure_name = readRDS(file.path(cachePath,'exposures.rds')) %>%
    filter(exposure_id == !!exposure_id) %>%
    select(exposure_name) %>% pull() %>% as.character()
  
  methodText = ifelse(method == 'SCCS',
                      'Self-Controlled Case Series',
                      'Historical Comparator')
  
  bayesDescription = ifelse(bayesAdjusted, 
                            "Bayesian with bias correction",
                            "Bayesian without bias correction")
  bayesPriorText = sprintf('prior SD = %.1f', prior_SD)
  
  capt = sprintf('Exposure: %s\nDatabase): %s\nMethod: %s\n%s, %s',
                 exposure_name, database_id, methodText,
                 bayesDescription, bayesPriorText)
  
  p_tts = p_tts + labs(caption = capt) +
    theme(plot.caption = element_text(hjust=0))
  
  if(showPlot){
    print(p_tts)
  }
  
  return(p_tts)
  
}

# 2. (c) ------
# alternative function with different calculation of TTS for MaxSPRT 
# (suspecting previous MaxSPRT calc. is bugged somehow...)
# 01/31/2023: add in hypothetical experiment of longer/shorter analysis plan
plotMaxSPRTBayesianTimeToSignal2 <- function(database_id,
                                            exposure_id,
                                            method,
                                            analysis_ids,
                                            sensitivity_level,
                                            prior_SD = 4,
                                            alpha = .05, # default alpha
                                            calibrateToAlpha = TRUE, # if calibrate to MaxSPRT alpha level?
                                            bayesAdjusted = TRUE, # use bias adjusted Bayesian results?
                                            calibration = FALSE, # use calibrated MaxSPRT?
                                            plan_extension_factor = 1, # shorter/longer analysis plan for MaxSPRT?
                                            summaryPath = summarypath, # summary results for Bayesian 
                                            localEstimatesPath = maxSPRT_filepath, # local estimate for MaxSPRT
                                            cachePath = './localCache/', # cache for caption info etc.
                                            colors = theColors, # fill colors for time bars
                                            showPlot = TRUE,
                                            addCaption = TRUE,
                                            maxCores = 8){
  # (1) compile all needed maxSPRT TTS results
  maxSPRT_tts = NULL
  localEstimates = readRDS(localEstimatesPath)
  
  for(aid in analysis_ids){
    # get everything from Type 1 and Type 2 error instead!!
    # for each analysis_id
    maxSPRT_res = frequentistDecisions(connection,
                                       'eumaeus',
                                       database_id = database_id,
                                       method = method, 
                                       exposure_id = exposure_id, 
                                       analysis_id = aid,
                                       estimates = localEstimates,
                                       calibration = calibration,
                                       correct_shift = TRUE,
                                       bonferroni_baseline = FALSE,
                                       FDR = FALSE,
                                       plan_extension_factor = plan_extension_factor,
                                       cachePath = cachePath,
                                       maxCores = maxCores)
    
    # pull TTS (with Type1error) from the result list
    this.maxSPRT_tts = getTTSfromMaxSPRTErrorRates(maxSPRT_res,
                                                   sensitivity_level = sensitivity_level,
                                                   adjust_max_time = TRUE)
    maxSPRT_tts = rbind(maxSPRT_tts, this.maxSPRT_tts)
  }
  
  
  # (2) Bayesian TTS
  ## only focus on needed prior to save time...
  priors = readRDS(file.path(cachePath, 'priorTable.rds'))
  this.prior_id = priors %>% filter(Sd == prior_SD) %>% select(prior_id) %>% pull()
  
  bayes_tts = foreach(aid = analysis_ids, .combine = rbind) %do% {
    
    ## get the "empirical alpha"
    this.alpha = maxSPRT_tts %>% filter(analysis_id == aid) %>%
      slice(1) %>% select(type1error) %>% pull()
    
    if(length(this.alpha) == 1 && is.double(this.alpha)){
      ## use it to "calibrate" Bayesian threshold to produce TTS results
      tts_bayes_raw = getTTSbyEmpiricalAlpha(database_id = database_id,
                                             method = method,
                                             analysis_id = aid,
                                             exposure_id = exposure_id,
                                             prior_ids = this.prior_id, # include all priors for easier query later
                                             alpha = this.alpha,
                                             summaryPath = summaryPath,
                                             cachePath = cachePath,
                                             useAdjusted = bayesAdjusted,
                                             calibrate = calibrateToAlpha,
                                             sensitivity_level = sensitivity_level,
                                             adjust_max_time = TRUE)
      tts_bayes_raw
    }else{
      NULL
    }
    
  }
  
  # (3) combine and generate plot
  p_tts = ttsComparePlotFromRes(maxSPRT_tts,
                                bayes_tts,
                                prior_SD = prior_SD,
                                sensitivity_level = sensitivity_level,
                                bayesAdjusted = bayesAdjusted,
                                showPlot = FALSE,
                                cachePath = cachePath,
                                colors = colors)
  
  # (4) get caption
  exposure_name = readRDS(file.path(cachePath,'exposures.rds')) %>%
    filter(exposure_id == !!exposure_id) %>%
    select(exposure_name) %>% pull() %>% as.character()
  
  methodText = ifelse(method == 'SCCS',
                      'Self-Controlled Case Series',
                      'Historical Comparator')
  
  bayesDescription = ifelse(bayesAdjusted, 
                            "Bayesian with bias correction",
                            "Bayesian without bias correction")
  bayesPriorText = sprintf('prior SD = %.1f', prior_SD)
  
  capt = sprintf('Exposure: %s\nDatabase: %s\nMethod: %s\n%s, %s',
                 exposure_name, database_id, methodText,
                 bayesDescription, bayesPriorText)
  
  p_tts = p_tts + labs(caption = capt) +
    theme(plot.caption = element_text(hjust=0))
  
  if(showPlot){
    print(p_tts)
  }
  
  return(p_tts)
  
}


## try this
# p_tts_full = plotMaxSPRTBayesianTimeToSignal(database_id = db,
#                                              method = method,
#                                              exposure_id = 211983,
#                                              analysis_ids = aids,
#                                              sensitivity_level = .25)


# 3.make a lot of plots ------
## for CCAE right now, and a few vaccines...

## 02/02/2023: try calibrated MaxSPRT and see how we are doing????
calibrationFlag = TRUE

## (a) Historical Comparator-----
exposures_HC= c(211983, 211833, 21184, 21185, 21215)
# 21214 doesn't have unadjusted results???
# 211981 has a weird bug somewhere for Bayesian -- not looking at it for now

plotsPath = '~/Documents/Research/betterResults/timeToSignal-Bayesian-MaxSPRT-calibrated-comparison/'

# pdf(file.path(plotsPath,'HistoricalComparator-CCAE.pdf'), 
#     width = 7, height = 7.1)

pdf(file.path(plotsPath,'HistoricalComparator-CCAE-calibrated.pdf'), 
    width = 7, height = 7.1)


for(eid in exposures_HC){
  p_tts_full = plotMaxSPRTBayesianTimeToSignal2(database_id = db,
                                                method = method,
                                                exposure_id = eid,
                                                analysis_ids = c(1:8),
                                                sensitivity_level = .5,
                                                calibration = calibrationFlag)
}


dev.off()

## 01/31/2023
# plots with shorter/longer analysis plans for MaxSPRT...
p_tts_full = plotMaxSPRTBayesianTimeToSignal2(database_id = 'CCAE',
                                              method = 'HistoricalComparator',
                                              exposure_id = 211833,
                                              analysis_ids = c(1:4),
                                              sensitivity_level = .5,
                                              plan_extension_factor = 0.5,
                                              maxCores = 8)

p_tts_full = plotMaxSPRTBayesianTimeToSignal2(database_id = 'CCAE',
                                              method = 'HistoricalComparator',
                                              exposure_id = 211833,
                                              analysis_ids = c(1:4),
                                              sensitivity_level = .5,
                                              plan_extension_factor = 2,
                                              maxCores = 8)


## (b) SCCS ------
exposures_SCCS = c(211983, 211833, 21184, 21185, 21215)
# 211981 has a weird bug somewhere for Bayesian,too
# performance not-so-different if SCCS

# pdf(file.path(plotsPath,'SCCS-CCAE.pdf'), 
#     width = 7, height = 7.3)

pdf(file.path(plotsPath,'SCCS-CCAE-calibrated.pdf'), 
    width = 7, height = 7.3)

for(eid in exposures_SCCS){
  p_tts_full = plotMaxSPRTBayesianTimeToSignal2(database_id = db,
                                                method = 'SCCS',
                                                exposure_id = eid,
                                                analysis_ids = c(1:8, 13:14),
                                                sensitivity_level = 0.5,
                                                showPlot = TRUE,
                                                calibration = calibrationFlag)
}

dev.off()


# # TEST CODE BELOW ------
# # DO NOT RUN!!!
# # 3. -------
# ## get time-to-signal dataframe for MaxSPRT
# p1 = plotTimeToSignalMaxSPRTOneExposure(localFile = maxSPRT_filepath,
#                                         method = method,
#                                         exposure_id = 211983,
#                                         analysis_ids = aids,
#                                         sensitivity_level = 0.5,
#                                         colors = theColors,
#                                         plotType = 'separate',
#                                         #plotType = 'together',
#                                         showPlot = FALSE)
# 
# ## cross-ref the MaxSPRT TTS dataframe with Type 1 error dataframe to get Type 1 error rate
# maxSPRT_tts = pullErrorRatesForTTS(connection, p1,
#                                    database_id = db,
#                                    method = method, 
#                                    exposure_id = eid, 
#                                    localEstimatesPath = maxSPRT_filepath,
#                                    calibration = FALSE)
# 
# 
# ## get the time-to-signal dataframe for Bayesian method
# bayes_tts_vanilla = foreach(aid = aids, .combine = rbind) %do% {
#   
#   ## get the "empirical alpha"
#   this.alpha = maxSPRT_tts %>% filter(analysis_id == aid) %>%
#     slice(1) %>% select(type1error) %>% pull()
#   
#   if(length(this.alpha) == 1 && is.double(this.alpha)){
#     ## use it to "calibrate" Bayesian threshold to produce TTS results
#     tts_bayes_raw = getTTSbyEmpiricalAlpha(database_id = db,
#                                            method = method,
#                                            analysis_id = aid,
#                                            exposure_id = eid,
#                                            prior_ids = c(1:3), # include all priors for easier query later
#                                            alpha = this.alpha,
#                                            summaryPath = summarypath,
#                                            cachePath = cachepath,
#                                            useAdjusted = FALSE,
#                                            calibrate = TRUE,
#                                            sensitivity_level = sens,
#                                            adjust_max_time = TRUE)
#     tts_bayes_raw
#   }else{
#     NULL
#   }
#   
# }
# 
# bayes_tts_adjusted = foreach(aid = aids, .combine = rbind) %do% {
#   
#   ## get the "empirical alpha"
#   this.alpha = maxSPRT_tts %>% filter(analysis_id == aid) %>%
#     slice(1) %>% select(type1error) %>% pull()
#   
#   if(length(this.alpha) == 1 && is.double(this.alpha)){
#     ## use it to "calibrate" Bayesian threshold to produce TTS results
#     tts_bayes_raw = getTTSbyEmpiricalAlpha(database_id = db,
#                                            method = method,
#                                            analysis_id = aid,
#                                            exposure_id = eid,
#                                            prior_ids = c(1:3), # include all priors for easier query later
#                                            alpha = this.alpha,
#                                            summaryPath = summarypath,
#                                            cachePath = cachepath,
#                                            useAdjusted = TRUE,
#                                            calibrate = TRUE,
#                                            sensitivity_level = sens,
#                                            adjust_max_time = TRUE)
#     tts_bayes_raw
#   }else{
#     NULL
#   }
#   
# }
# 
# # try an example plot
# ## ADJUSTED BAYES
# prior_SD = 4
# 
# tts_long = rbind(
#   maxSPRT_tts %>% 
#     select(database_id, exposure_id, method, analysis_id, effect_size, time_to_sens, approach),
#   bayes_tts_adjusted %>%
#     filter(prior_sd == prior_SD) %>%
#     select(database_id, exposure_id, method, analysis_id, effect_size, time_to_sens, approach)
# )
# 
# tts_long$approach = factor(tts_long$approach,
#                            levels = c('MaxSPRT', 'Adjusted Bayes'))
# 
# ##
# 
# pg = ggplot(tts_long, 
#             aes(y=as.factor(effect_size), 
#                 x=time_to_sens, 
#                 fill=approach)) +
#   geom_bar(stat='identity', position = 'dodge')+
#   geom_hline(yintercept = c(1.5,2.5), color='gray40')+
#   scale_x_continuous(breaks = seq(from=0, to=12, by=3)) +
#   labs(y='Effect Size', 
#        x=sprintf('Time to %.0f%% sensitivity', sens * 100),
#        fill = '') +
#   facet_grid(analysis_id~.) +
#   theme_bw(base_size = 14)+
#   theme(panel.grid.major.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.background = element_blank(),
#         strip.text.y = element_text(angle = 0),
#         legend.position = 'bottom')+
#   guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
#   scale_fill_manual(values = theColors,
#                     breaks = c('Adjusted Bayes','MaxSPRT'),
#                     labels = c('Bayesian','MaxSPRT'))
# 
# print(pg)
