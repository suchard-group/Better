# functions and code to plot
# estiamtes and 95% Bayesian CIs
# for GBS-Zoster

source('./extras/postProcessUtils.R')
source('./R/betterHelperFunctions.R')

library(wesanderson)

## specify paths, database, methods, and exposures details
resultspath = '~/Documents/Research/betterGBSanalysesResults/betterGBSResults/'
summarypath = '~/Documents/Research/betterGBSanalysesResults/resultsSummary/'

#db = 'CCAE'
#db = 'OptumEhr'
db = 'IBM_MDCR'
methods = c('SCCS', 'HistoricalComparator')
exposures = c(211981:211983)

## pull summary dataframe
summ = pullGBSResultsMeta(database_id = db,
                          methods = methods,
                          exposure_ids = exposures,
                          resultsPath = resultspath,
                          savePath = summarypath)


# 1. function to plot effect estimates for GBS-Zoster ----
# (also plot 95% Bayesian CIs)
plotEffectEstimates <- function(summ,
                                database_id,
                                method,
                                exposure_id,
                                analysis_ids,
                                period_id,
                                prior_id = 3, # default prior SD = 4
                                adjust = FALSE,
                                showPlot = TRUE,
                                cachePath = './localCache/'){
  dat = summ %>%
    filter(database_id == !!database_id,
           method == !!method,
           exposure_id == !!exposure_id,
           analysis_id %in% analysis_ids,
           period_id == !!period_id,
           prior_id == !!prior_id)
  
  if(nrow(dat) == 0){
    cat('No results available!!\n')
    return()
  }
  
  if(adjust){
    dat = dat %>% 
      select(estimate = adjustedPostMedian, 
             lb = adjustedCI95_lb,
             ub = adjustedCI95_ub,
             analysis_id)
  }else{
    dat = dat %>% 
      select(estimate = postMedian, 
             lb = CI95_lb,
             ub = CI95_ub,
             analysis_id)
  }
  
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
  
  y_breaks = as.character(unique(dat$analysis_id))
  y_labels = unique(dat$analysis_text)
  
  p = ggplot(dat, aes(x=estimate, y = as.factor(analysis_id))) +
    geom_point(shape = 18, size = 5) +
    geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.1) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    scale_y_discrete(breaks = y_breaks, labels = y_labels)+
    labs(x='Log relative rate ratio (95% CI)',
         y = '') +
    theme_bw()
  
  ## add "adjust" label to data
  dat$adjust = ifelse(adjust, 'adjusted', 'unadjusted')
  
  attr(p, 'data') = dat
  
  if(showPlot){
    print(p)
  }
  
  return(p)
}


# ## do the plotting (a bit rough for now...) -----
# p1 = plotEffectEstimates(summ, 
#                          database_id = 'CCAE', 
#                          method = 'SCCS', 
#                          exposure_id = 211981, 
#                          analysis_ids = c(1:8, 13:14), # exclude 0-1 days 
#                          period_id = 12,
#                          prior_id = 3,
#                          adjust = TRUE)
# 
# p2 = plotEffectEstimates(summ, 
#                          database_id = 'CCAE', 
#                          method = 'HistoricalComparator', 
#                          exposure_id = 211981, 
#                          analysis_ids = 1:12, 
#                          period_id = 12,
#                          prior_id = 2,
#                          adjust = FALSE)


# 2. function to compare adjusted and unadjusted results
compareEffectEstimates <- function(colors = NULL,
                                   logScale = TRUE,
                                   ...){
  adjusted = plotEffectEstimates(...,adjust = TRUE,  
                                 showPlot = FALSE)
  unadjusted = plotEffectEstimates(..., adjust = FALSE,  
                                   showPlot = FALSE)
  
  dat.adj = attr(adjusted, 'data')
  dat.unadj = attr(unadjusted, 'data')
  
  dat = rbind(dat.adj, dat.unadj) %>% 
    arrange(analysis_id) %>%
    mutate(adjust = factor(adjust, 
                           levels = c('unadjusted', 'adjusted')))
  
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
  
  print(p)
}


## make comparison plots ----

me = 'SCCS'
aids = c(5,6,8,14)

me = 'HistoricalComparator'
aids = c(5:8)

eid = 211983


compareEffectEstimates(summ = summ, 
                       database_id = db, 
                       method = me, 
                       exposure_id = eid, 
                       analysis_ids = aids, 
                       period_id = 12,
                       prior_id = 2,
                       colors = wes_palette("Darjeeling2")[2:3],
                       logScale = FALSE)
