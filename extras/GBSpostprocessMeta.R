# Aug 2022
# GBS analyses post processing again 
# with meta analysis for negative control outcomes

source('./extras/postProcessUtils.R')
source('./R/betterHelperFunctions.R')
#source('./extras/simpleCalibration.R')


resultspath = '~/Documents/Research/betterGBSanalysesResults/betterGBSResults/'
summarypath = '~/Documents/Research/betterGBSanalysesResults/resultsSummary/'

db = 'CCAE'
methods = c('SCCS', 'HistoricalComparator')
exposures = c(211981:211983)

summ = pullGBSResultsMeta(database_id = db,
                          methods = methods,
                          exposure_ids = exposures,
                          resultsPath = resultspath,
                          savePath = summarypath)

# function to plot effect estimates for GBS-Zoster ----
# (also plot 95% Bayesian )
plotEffectEstimates <- function(summ,
                                database_id,
                                method,
                                exposure_id,
                                analysis_ids,
                                period_id,
                                prior_id = 3, # default prior SD = 4
                                adjust = FALSE,
                                localCache = './localCache/'){
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
  
  p = ggplot(dat, aes(x=estimate, y = as.factor(analysis_id))) +
    geom_point(shape = 18, size = 5) +
    geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.1) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    labs(x='Log relative rate ratio (95% CI)',
         y = 'Analysis design') +
    theme_bw()
  
  attr(p, 'data') = dat
  
  print(p)
  
  return(p)
}

## try it
p1 = plotEffectEstimates(summ, 
                         database_id = 'CCAE', 
                         method = 'SCCS', 
                         exposure_id = 211981, 
                         analysis_ids = 1:15, 
                         period_id = 12,
                         prior_id = 3,
                         adjust = TRUE)

p2 = plotEffectEstimates(summ, 
                         database_id = 'CCAE', 
                         method = 'HistoricalComparator', 
                         exposure_id = 211981, 
                         analysis_ids = 1:12, 
                         period_id = 12,
                         prior_id = 2,
                         adjust = TRUE)


                          