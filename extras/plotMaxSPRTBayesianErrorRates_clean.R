# 01/16/2023
# cleaner version of
# `extras/plotMaxSPRTBayesianErrorRates.R`

# 03/29/2023
# add a plotting power function as well!
# with each effect size on one panel

library(wesanderson)

source('./extras/simpleCalibration.R')
source('./extras/frequentistDecisionComparisons.R')

# path for saving intermediate results
summarypath = '~/Documents/Research/betterResults/summary'
#samplepath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/"
cachepath = './localCache/'

maxSPRT_filepath = './localCache/EstimateswithImputedPcs_CCAE.rds'
maxSPRT_estimates = readRDS(maxSPRT_filepath)

# analyses results to check
db = 'CCAE'
#db = 'OptumEhr'
eid = 211983 # Zoster two doses
#eid = 21184 #H1N1
me = 'HistoricalComparator'
#me = 'SCCS'
aid = 4
#pid = 3 # 1: sd=10; 2: sd=1.5; 3: sd=4 
#tolerance = 0.004

## color
theColors = c(wes_palette("Darjeeling1")[5],
              wes_palette("Darjeeling2")[2],
              wes_palette("BottleRocket2")[3])


# (1) ----
# main function to plot Type 1 error comparisons across different design variants----
# 03/30/2023: add functionality to remove caption text....
plotMaxSPRTBayesianType1 <- function(database_id,
                                     method,
                                     exposure_id,
                                     analysis_ids,
                                     prior_id = 3, # default to SD = 4
                                     raw_Bayesian_alpha = 0.05,
                                     adj_Bayesian_alpha = 0.05,
                                     calibrateToAlpha = FALSE,
                                     plan_extension_factor = 1,
                                     cachePath = './localCache/',
                                     maxSPRTestimatePath = './localCache/EstimateswithImputedPcs_CCAE.rds',
                                     summaryPath = summarypath,
                                     colors = theColors,
                                     showPlot = TRUE,
                                     showCaption = TRUE,
                                     maxCores = 8){
  
  
  maxSPRT_estimates = readRDS(maxSPRTestimatePath)
  
  maxSPRT_errorRates = NULL
  for(aid in analysis_ids){
    resLst = frequentistDecisions(NULL,
                                  NULL,
                                  database_id = database_id,
                                  method = method,
                                  exposure_id = exposure_id,
                                  analysis_id = aid,
                                  calibration = FALSE,
                                  correct_shift = TRUE,
                                  plan_extension_factor = plan_extension_factor,
                                  cachePath = cachePath,
                                  estimates = maxSPRT_estimates,
                                  maxCores = maxCores)
    
    this.maxsprt_errors = resLst$errorRate %>%
      select(database_id, method, exposure_id, analysis_id,
             period_id, y = errorRate, effect_size, stats) %>%
      mutate(approach = '1: MaxSPRT')
    
    maxSPRT_errorRates = rbind(maxSPRT_errorRates,
                               this.maxsprt_errors)
    
    
  }
  
  # get raw Bayesian error rates 
  Bayes_raw_errorRates = NULL
  for(aid in analysis_ids){
    res_raw = plotTempDelta1ByPriors(database_id = database_id,
                                     method = method, 
                                     analysis_id = aid,
                                     exposure_id = exposure_id,
                                     prior_ids = prior_id, # include all priors for easier query later
                                     alpha = raw_Bayesian_alpha, #0.25, #
                                     summaryPath = summaryPath,
                                     cachePath = cachePath,
                                     useAdjusted = FALSE,
                                     showPlots = FALSE,
                                     stratifyByEffectSize = TRUE,
                                     calibrate = calibrateToAlpha,
                                     outcomesInEstimates = NULL)
    
    this.Bayes_errorRates = res_raw %>% 
      mutate(database_id = database_id, method = method, 
             exposure_id = exposure_id, analysis_id = aid) %>%
      select(-prior_id, -Sd, -priorLabel) %>%
      mutate(approach = '2: Bayesian w/o correction')
    
    Bayes_raw_errorRates = rbind(Bayes_raw_errorRates,
                                 this.Bayes_errorRates)
  }
  
  # get adjusted Bayesian results
  # get raw Bayesian error rates 
  Bayes_adj_errorRates = NULL
  for(aid in analysis_ids){
    res_adj = plotTempDelta1ByPriors(database_id = database_id,
                                     method = method, 
                                     analysis_id = aid,
                                     exposure_id = exposure_id,
                                     prior_ids = prior_id, # include all priors for easier query later
                                     alpha = adj_Bayesian_alpha, #0.25, #
                                     summaryPath = summaryPath,
                                     cachePath = cachePath,
                                     useAdjusted = TRUE,
                                     showPlots = FALSE,
                                     stratifyByEffectSize = TRUE,
                                     calibrate = calibrateToAlpha,
                                     outcomesInEstimates = NULL)
    
    this.Bayes_errorRates = res_adj %>% 
      mutate(database_id = database_id, method = method, 
             exposure_id = exposure_id, analysis_id = aid) %>%
      select(-prior_id, -Sd, -priorLabel) %>%
      mutate(approach = '3: Bayesian w/ correction')
    
    Bayes_adj_errorRates = rbind(Bayes_adj_errorRates,
                                 this.Bayes_errorRates)
  }
  
  # combine and make plot
  errors_combined = rbind(maxSPRT_errorRates,
                          Bayes_raw_errorRates,
                          Bayes_adj_errorRates)
  
  earlySplit = 4
  alphaLevel = 0.2
  
  errors_combined = errors_combined %>%
    mutate(stage = if_else(period_id <= earlySplit, alphaLevel, 1))
  
  Type1errors = errors_combined %>% filter(stats=='type 1')
  
  # plotting setup
  yinters = 0.05
  ybreaks = c(0,0.05, 0.1, 0.25, 0.5, 0.75,1.0)
  period_breaks = seq(from = min(errors_combined$period_id),
                      to = max(errors_combined$period_id),
                      by = 3)
  period_labels = as.integer(period_breaks)
  type1colors = colors
  
  analyses = readRDS(file.path(cachePath,'analyses.rds'))
  Type1errors = Type1errors %>% 
    left_join(analyses, by = c('method', 'analysis_id')) %>%
    rename('analysis' = 'description')
  
  # figure out caption
  if(showCaption){
    exposure_name = readRDS(file.path(cachePath,'exposures.rds')) %>%
      filter(exposure_id == !!exposure_id) %>%
      select(exposure_name) %>% pull() %>% as.character()
    
    methodText = ifelse(method == 'SCCS',
                        'Self-Controlled Case Series',
                        'Historical Comparator')
    
    #bayesPriorText = sprintf('prior SD = %.1f', prior_SD)
    
    capt = sprintf('Exposure: %s\nDatabase: %s\nMethod: %s',
                   exposure_name, database_id, methodText)
  }else{
    capt = ''
  }
  
  
  # actual plotting
  p = ggplot(Type1errors, 
             aes(x=period_id, y=y, color=approach, alpha = stage))+
    geom_line(size = 1.5) +
    geom_point(size=2)+
    geom_hline(yintercept = yinters, 
               color = 'gray60', 
               size = 1, linetype=2)+
    scale_y_continuous(limits = c(0,1),
                       breaks = ybreaks,
                       trans = 'sqrt'
    )+
    facet_grid(.~analysis,
               labeller = label_wrap_gen(width=28)) +
    scale_x_continuous(breaks = period_breaks, labels = period_labels)+
    labs(x='analysis period (months)', y='Type 1 error rate', 
         caption = capt, color='Type 1 error of:')+
    scale_color_manual(values = type1colors) +
    scale_alpha_continuous(range = c(0.2, 1), guide = 'none')+
    guides(color=guide_legend(nrow=1,byrow=TRUE))+
    #facet_grid(.~approach)+
    theme_bw(base_size = 16)+
    theme(legend.position = 'bottom',
          plot.caption = element_text(hjust=0))# change to bottom legend...
  
  if(showPlot){
    print(p)
  }
  
  attr(p, 'data') = Type1errors
  
  return(p)
  
}

# (2) ------
# another function to power; similar but only do one analysis at a time...
# default: calibrate toward the SAME type 1 error level as MaxSPRT!
plotMaxSPRTBayesianPower <- function(database_id,
                                     method,
                                     exposure_id,
                                     analysis_id,
                                     prior_id = 3, # default to SD = 4
                                     calibrateToAlpha = TRUE,
                                     plan_extension_factor = 1,
                                     cachePath = './localCache/',
                                     maxSPRTestimatePath = './localCache/EstimateswithImputedPcs_CCAE.rds',
                                     summaryPath = summarypath,
                                     colors = theColors,
                                     showPlot = TRUE,
                                     showCaption = TRUE,
                                     maxCores = 8){
  
  
  maxSPRT_estimates = readRDS(maxSPRTestimatePath)
  
  # (a) maxSPRT results
  resLst = frequentistDecisions(NULL,
                                NULL,
                                database_id = database_id,
                                method = method,
                                exposure_id = exposure_id,
                                analysis_id = analysis_id,
                                calibration = FALSE,
                                correct_shift = TRUE,
                                plan_extension_factor = plan_extension_factor,
                                cachePath = cachePath,
                                estimates = maxSPRT_estimates,
                                maxCores = maxCores)
  maxSPRT_errorRates = resLst$errorRate %>%
    select(database_id, method, exposure_id, analysis_id,
           period_id, y = errorRate, effect_size, stats) %>%
    mutate(approach = '1: MaxSPRT')
  
  maxSPRT_end_alpha = maxSPRT_errorRates %>% 
    filter(period_id == max(period_id), stats == 'type 1') %>%
    select(y) %>% pull()
  
  # (b) get raw Bayesian error rates 
  #     (calibrate to MaxSPRT's end-of-analysis alpha level!!!)
  Bayes_raw_errorRates = NULL
  res_raw = plotTempDelta1ByPriors(database_id = database_id,
                                   method = method, 
                                   analysis_id = analysis_id,
                                   exposure_id = exposure_id,
                                   prior_ids = prior_id, 
                                   alpha = maxSPRT_end_alpha,
                                   summaryPath = summaryPath,
                                   cachePath = cachePath,
                                   useAdjusted = FALSE,
                                   showPlots = FALSE,
                                   stratifyByEffectSize = TRUE,
                                   calibrate = calibrateToAlpha,
                                   outcomesInEstimates = NULL)
  Bayes_raw_errorRates = res_raw %>% 
    mutate(database_id = database_id, method = method, 
           exposure_id = exposure_id, analysis_id = aid) %>%
    select(-prior_id, -Sd, -priorLabel) %>%
    mutate(approach = '2: Bayesian w/o correction')
  
  # (c) get adjusted Bayesian results
  #     (calibrate to MaxSPRT's end-of-analysis alpha level!!!)
  res_adj = plotTempDelta1ByPriors(database_id = database_id,
                                   method = method, 
                                   analysis_id = aid,
                                   exposure_id = exposure_id,
                                   prior_ids = prior_id, # include all priors for easier query later
                                   alpha = maxSPRT_end_alpha,
                                   summaryPath = summaryPath,
                                   cachePath = cachePath,
                                   useAdjusted = TRUE,
                                   showPlots = FALSE,
                                   stratifyByEffectSize = TRUE,
                                   calibrate = calibrateToAlpha,
                                   outcomesInEstimates = NULL)
  Bayes_adj_errorRates = res_adj %>% 
    mutate(database_id = database_id, method = method, 
           exposure_id = exposure_id, analysis_id = aid) %>%
    select(-prior_id, -Sd, -priorLabel) %>%
    mutate(approach = '3: Bayesian w/ correction')
  
  # combine and make plot
  errors_combined = rbind(maxSPRT_errorRates,
                          Bayes_raw_errorRates,
                          Bayes_adj_errorRates)
  
  earlySplit = 4
  alphaLevel = 0.2
  
  errors_combined = errors_combined %>%
    mutate(stage = if_else(period_id <= earlySplit, alphaLevel, 1))
  
  powers = errors_combined %>% 
    filter(stringr::str_starts(stats, 'type 2')) %>% 
    filter(stats != 'delta1') %>% 
    mutate(power = 1-y) %>%
    mutate(trueRR = sprintf('RR = %.1f', effect_size))
  
  # plotting setup
  # yinters = 0.05
  # ybreaks = c(0,0.05, 0.1, 0.25, 0.5, 0.75,1.0)
  period_breaks = seq(from = min(errors_combined$period_id),
                      to = max(errors_combined$period_id),
                      by = 2)
  period_labels = as.integer(period_breaks)
  powerColors = colors
  
  # figure out caption
  if(showCaption){
    analyses = readRDS(file.path(cachePath,'analyses.rds'))
    analysisDesc = analyses %>% filter(method == !!method, 
                                       analysis_id == !!analysis_id) %>%
      select(description) %>% pull()
    methodText = ifelse(method == 'SCCS',
                        'Self-Controlled Case Series',
                        'Historical Comparator')
    methodText = paste0(methodText, ', ', analysisDesc)
    
    exposure_name = readRDS(file.path(cachePath,'exposures.rds')) %>%
      filter(exposure_id == !!exposure_id) %>%
      select(exposure_name) %>% pull() %>% as.character()
    
    #bayesPriorText = sprintf('prior SD = %.1f', prior_SD)
    
    capt = sprintf('Exposure: %s\nDatabase: %s\nMethod: %s',
                   exposure_name, database_id, methodText)
  }else{
    capt = ''
  }
  
  
  # actual plotting
  p = ggplot(powers, 
             aes(x=period_id, y=power, 
                 color = approach, 
                 alpha = stage))+
    geom_line(size = 1.5) +
    geom_point(size=2)+
    scale_y_continuous(limits = c(0,1))+
    scale_x_continuous(breaks = period_breaks, labels = period_labels)+
    labs(x='analysis period (months)', 
         y='power', 
         caption = capt, 
         color='Statistical power of:')+
    scale_color_manual(values = powerCols) +
    scale_alpha_continuous(range = c(0.2, 1), guide = 'none')+
    facet_grid(.~trueRR)+
    theme_bw(base_size = 15)+
    theme(legend.position = 'bottom',
          plot.caption = element_text(hjust=0)) # change to bottom legend...
  
  if(showPlot){
    print(p)
  }
  
  attr(p, 'data') = powers
  
  return(p)
  
}

## try Type 1 error plots ----
# 03/30/2023 update: plots without captions
pType1 = plotMaxSPRTBayesianType1(database_id = db,
                                  method = me,
                                  exposure_id = eid,
                                  analysis_ids = 1:4,
                                  raw_Bayesian_alpha = 0.05,
                                  adj_Bayesian_alpha = 0.05,
                                  calibrateToAlpha = FALSE,
                                  showCaption = FALSE)

pType1 = plotMaxSPRTBayesianType1(database_id = db,
                                  method = 'SCCS',
                                  exposure_id = eid,
                                  analysis_ids = c(1:4,13),
                                  raw_Bayesian_alpha = 0.03,
                                  adj_Bayesian_alpha = 0.03,
                                  calibrateToAlpha = FALSE,
                                  showCaption = FALSE)

## try power plots----
pPowers = plotMaxSPRTBayesianPower(database_id = db,
                                   method = me,
                                   exposure_id = eid,
                                   analysis_id = 4,
                                   prior_id = 3,
                                   calibrateToAlpha = TRUE,
                                   showPlot = TRUE)

pPowers = plotMaxSPRTBayesianPower(database_id = db,
                                   method = 'SCCS',
                                   exposure_id = eid,
                                   analysis_id = 13,
                                   prior_id = 3,
                                   calibrateToAlpha = TRUE,
                                   showPlot = TRUE)


## 01/31/2023
# try with shorter analysis 
pType1 = plotMaxSPRTBayesianType1(database_id = db,
                                  method = me,
                                  exposure_id = 211833,
                                  analysis_ids = 1:4,
                                  raw_Bayesian_alpha = 0.05,
                                  adj_Bayesian_alpha = 0.05,
                                  calibrateToAlpha = FALSE,
                                  plan_extension_factor = 0.5)


# AND longer analysis
pType1 = plotMaxSPRTBayesianType1(database_id = db,
                                  method = me,
                                  exposure_id = 211833,
                                  analysis_ids = 1:4,
                                  raw_Bayesian_alpha = 0.05,
                                  adj_Bayesian_alpha = 0.05,
                                  calibrateToAlpha = FALSE,
                                  plan_extension_factor = 2)



# (a) produce plots for Historical Comparator----
exposures_HC= c(211983, 211833, 21184, 21185, 21215)

plotPath = '~/Documents/Research/betterResults/plots'

pdf(file.path(plotPath, 'Type1Error-comparison-HC-variants.pdf'),
    width = 10.5, height = 4.8)

for(eid in exposures_HC){
  pType1 = plotMaxSPRTBayesianType1(database_id = db,
                                    method = 'HistoricalComparator',
                                    exposure_id = eid,
                                    analysis_ids = 1:4,
                                    raw_Bayesian_alpha = 0.05,
                                    adj_Bayesian_alpha = 0.05,
                                    calibrateToAlpha = FALSE)
}


dev.off()


# (a) produce plots for Historical Comparator----
exposures_SCCS = c(211983, 211833, 21184, 21185, 21215)

plotPath = '~/Documents/Research/betterResults/plots'

pdf(file.path(plotPath, 'Type1Error-comparison-SCCS-variants.pdf'),
    width = 13.6, height = 4.8)

for(eid in exposures_SCCS){
  pType1 = plotMaxSPRTBayesianType1(database_id = db,
                                    method = 'SCCS',
                                    exposure_id = eid,
                                    analysis_ids = c(1:4,13),
                                    raw_Bayesian_alpha = 0.03,
                                    adj_Bayesian_alpha = 0.03,
                                    calibrateToAlpha = FALSE)
}


dev.off()
                                  
