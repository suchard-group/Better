# 08/29/2022
# pull frequentist GBS estimates
# and plot estimates and CIs

library(dplyr)

# 1. function to pull estimates first------
pullGBSfreqEstimates <- function(database_id,
                                 method,
                                 exposure_id,
                                 analysis_ids,
                                 resultsPath = '~/Documents/Research/better_gbs/',
                                 removeNA = TRUE){
  
  subpath = sprintf('Results_%s', database_id)
  estimatesFile = file.path(resultsPath, subpath, 'estimate_withCalibration.csv')
  
  estimates = readr::read_csv(estimatesFile)
  
  names(estimates) = SqlRender::camelCaseToSnakeCase(names(estimates))
  
  estimates = estimates %>% 
    filter(method == !!method,
           exposure_id == !!exposure_id,
           analysis_id %in% analysis_ids) %>%
    select(database_id,
           method,
           analysis_id,
           exposure_id,
           outcome_id,
           period_id,
           log_rr,
           ci_95_lb,
           ci_95_ub,
           calibrated_log_rr,
           calibrated_ci_95_lb,
           calibrated_ci_95_ub,
           llr,
           calibrated_llr,
           critical_value)
  
  if(removeNA){
    estimates = estimates %>% 
      filter(!is.na(log_rr), !is.na(calibrated_log_rr))
  }
  
  return(estimates)
  
}



# 2. function to plot estimates and CIs ----
# 09/12/2022 update: transform x-axis labels to rate ratio scale
plotGBSfreqEstimatesCIs <- function(method,
                                    colors = NULL,
                                    logScale = TRUE,
                                    showPlot = TRUE,
                                    cachePath = './localCache',
                                    ...){
  
  estimates = pullGBSfreqEstimates(method = method,...) %>%
    filter(period_id == max(period_id))
  
  uncalibrated = estimates %>% 
    select(database_id,
           method,
           analysis_id,
           exposure_id,
           outcome_id,
           period_id,
           estimate = log_rr,
           lb = ci_95_lb,
           ub = ci_95_ub,
           llr = llr,
           critical_value) %>%
    mutate(approach = 'uncalibrated')
  
  calibrated = estimates %>% 
    select(database_id,
           method,
           analysis_id,
           exposure_id,
           outcome_id,
           period_id,
           estimate = calibrated_log_rr,
           lb = calibrated_ci_95_lb,
           ub = calibrated_ci_95_ub,
           llr = calibrated_llr,
           critical_value) %>%
    mutate(approach = 'calibrated')
  
  dat = rbind(calibrated, uncalibrated) %>%
    mutate(adjust = factor(approach, 
                              levels = c('uncalibrated', 'calibrated'))) %>%
    mutate(lb = log(lb),
           ub = log(ub))
  
  
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
  
  # join with dat
  dat = dat %>% left_join(analyses, by = 'analysis_id') %>%
    arrange(analysis_id)
  
  # plotting
  
  y_breaks = as.character(unique(dat$analysis_id))
  y_labels = unique(dat$analysis_text)
  
  xmin = min(dat$estimate, dat$lb) %>% floor()
  xmax = max(dat$estimate, dat$ub) %>% ceiling()
  x_breaks = seq(from = xmin, to = xmax, by = 1)
  if(logScale){
    xtext = 'Log relative rate ratio (95% CI)'
    x_labels = as.character(x_breaks)
  }else{
    # xmin = min(dat$estimate, dat$lb) %>% exp() %>% floor()
    # xmax = max(dat$estimate, dat$ub) %>% exp() %>% ceiling()
    xtext = 'Relative rate ratio (95% CI)'
    x_labels = as.character(round(exp(x_breaks),1))
  }
  
  
  
  p = ggplot(dat, aes(x=estimate, 
                      y = as.factor(analysis_id),
                      color = adjust)) +
    geom_point(shape = 18, size = 5,
               position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = lb, xmax = ub), 
                   height = 0.1,
                   position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    scale_y_discrete(breaks = y_breaks, labels = y_labels)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    labs(x=xtext,
         y = '',
         color = '') +
    theme_bw() + 
    theme(legend.position = 'bottom')
  
  
  if(!is.null(colors)){
    p = p + scale_color_manual(values = colors)
  }
  
  attr(p, 'data') = dat
  
  if(showPlot){
    print(p)
  }
    
  return(p)
}

# try it ----

resultspath = '~/Documents/Research/better_gbs/'
#db = 'CCAE'
db = 'MDCR'
#me = 'HistoricalComparator'
me = 'SCCS'
eid = 211981

freqPlot = plotGBSfreqEstimatesCIs(database_id = db,
                        method = me,
                        exposure_id = eid,
                        analysis_ids = 1:12,
                        colors = wes_palette("Darjeeling2")[2:3],
                        logScale = FALSE)

# not doing OptumEhr yet...----
# haven't run calibration for it somehow
#
# freqPlot2 = plotGBSfreqEstimatesCIs(database_id = 'OptumEHR',
#                                    method = me,
#                                    exposure_id = eid,
#                                    analysis_ids = 1:12,
#                                    colors = wes_palette("Darjeeling2")[2:3])
# 

