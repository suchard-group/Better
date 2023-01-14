# 01/13/2023

# cleaner script for one-exposure time-to-signal

library(stringr)
source('./extras/simpleCalibration.R')

# 1. get the TTS data frame from calibrated Bayesian error rates
getTTSfromCalibratedErrorRates <- function(res,
                                           sensitivity_level,
                                           adjust_max_time = TRUE){
  res = res %>% 
    select(period_id, prior_sd = Sd, effect_size, error_rate = y,
           stats) %>%
    filter(str_starts(stats,'type 2')) %>%
    select(-stats)
  
  # create fake "last entries" of type2error = 0
  # for better manipulation
  max_period_plus = max(res$period_id) + 1
  res = rbind(res,
              tibble(period_id = max_period_plus, 
                     prior_sd = rep(unique(res$prior_sd), 
                                    length(unique(res$effect_size))),
                     effect_size= rep(unique(res$effect_size),
                                      length(unique(res$prior_sd))),
                     error_rate = 0)
              )
  
  type2_ub = 1 - sensitivity_level
  
  tts = res %>%
    group_by(prior_sd, effect_size) %>%
    filter(error_rate <= type2_ub) %>%
    summarize(time_to_sens = min(period_id)) %>%
    ungroup()
  
  if(adjust_max_time){
    tts = tts %>%
      mutate(time_to_sens = if_else(time_to_sens == max_period_plus,
                                    max(res$period_id),
                                    time_to_sens))
  }
  
  return(tts)
}

## try it here
#
# raw_tts = getTTSfromCalibratedErrorRates(res_raw, 
#                                          sensitivity_level = 0.25,
#                                          adjust_max_time = TRUE)

# 2. function to get calibrated Bayesian results, then generate a "full information" tts table
getTTSbyEmpiricalAlpha <- function(database_id,
                                   method,
                                   analysis_id,
                                   exposure_id,
                                   prior_ids,
                                   alpha,
                                   summaryPath,
                                   cachePath,
                                   useAdjusted,
                                   calibrate = TRUE,
                                   sensitivity_level = 0.25,
                                   adjust_max_time = TRUE){
  # (1) get the Bayesian error rates results
  res_bayes = plotTempDelta1ByPriors(database_id = database_id,
                                     method = method, 
                                     analysis_id = analysis_id,
                                     exposure_id = exposure_id,
                                     prior_ids = prior_ids, # include all priors for easier query later
                                     alpha = alpha, #0.25, #
                                     summaryPath = summaryPath,
                                     cachePath = cachePath,
                                     useAdjusted = useAdjusted,
                                     showPlots = FALSE,
                                     stratifyByEffectSize = TRUE,
                                     calibrate = calibrate, 
                                     outcomesInEstimates = NULL)
  
  # (2) get TTS from the error rates
  tts_bayes = getTTSfromCalibratedErrorRates(res_bayes, 
                                             sensitivity_level = sensitivity_level,
                                             adjust_max_time = adjust_max_time)
  
  # (3) pad with useful info
  tts_bayes = tts_bayes %>%
    mutate(database_id = database_id,
           method = method, 
           analysis_id = analysis_id,
           exposure_id = exposure_id,
           approach = if_else(useAdjusted, 'Adjusted Bayes', 'Unadjusted Bayes'))
  
  return(tts_bayes)
  
}


# # test it!
# summarypath = '~/Documents/Research/betterResults/summary'
# cachepath = './localCache/'
# 
# tts_bayes_raw = getTTSbyEmpiricalAlpha(database_id = 'CCAE',
#                                    method = 'HistoricalComparator', 
#                                    analysis_id = 4,
#                                    exposure_id = 211983,
#                                    prior_ids = c(1:3), # include all priors for easier query later
#                                    alpha =0.25, #0.25, #
#                                    summaryPath = summarypath,
#                                    cachePath = cachepath,
#                                    useAdjusted = FALSE,
#                                    calibrate = TRUE,
#                                    sensitivity_level = 0.25,
#                                    adjust_max_time = TRUE)
# 
# tts_bayes_adj = getTTSbyEmpiricalAlpha(database_id = 'CCAE',
#                                       method = 'HistoricalComparator', 
#                                       analysis_id = 4,
#                                       exposure_id = 211983,
#                                       prior_ids = c(1:3), # include all priors for easier query later
#                                       alpha =0.25, #0.25, #
#                                       summaryPath = summarypath,
#                                       cachePath = cachepath,
#                                       useAdjusted = TRUE,
#                                       calibrate = TRUE,
#                                       sensitivity_level = 0.25,
#                                       adjust_max_time = TRUE)
