# July 2022
# compare results with different prior choices for negative control analysis
# (using Bayesian meta analysis)
# Aug 2022
# add t-model for meta analysis

source('./extras/postProcessUtils.R')
source('./extras/simpleCalibration.R')

## be hacky about it for now
## pretending we are dealing with ALL exposures
exposures = readRDS('./localCache/exposures.rds')
allExposure_ids = exposures$exposure_id

## a function to compute error rates and show plots 
summarizeAndPlotTempDelta1ByPrior <- function(database_id, method, exposures, rawPath, summaryPath, ...){
  # pull results first
  summ = pullResults(database_id, method, exposures, resultsPath = rawPath, savePath = summaryPath)
  
  # then get error rates results and plots
  res_summ = plotTempDelta1ByPriors(database_id, method, summaryPath = summaryPath, ...)
  
  return(res_summ)
}


## specify database, method, analysis and exposure
db = 'CCAE'
me = 'SCCS' 
#me = 'HistoricalComparator'
aid = 8
eid = 211981
# summarypath = outputPath
cachepath = './localCache/'


## go through the different prior settings

adjust = TRUE
pid = 3 # choose SD = 4.0 prior from each result

dir_suffixes = c('default3', 'shrinkMu3', 'shrinkBoth3')

errors_combined = NULL

for(dir_suffix in dir_suffixes){
  if(dir_suffix == dir_suffixes[1]){
    title = 'muPriorSd = 2, tauPriorSd = 0.5'
  }else if(dir_suffix == dir_suffixes[2]){
    title = 'muPriorSd = 0.5, tauPriorSd = 0.5'
  }else{
    title = 'muPriorSd = 0.2, tauPriorSd = 0.2'
  }
  
  resultPath = sprintf('~/Documents/Research/betterResults/betterResults-%s', dir_suffix)
  outputPath = sprintf('~/Documents/Research/betterResults/summary-%s/', dir_suffix)
  
  res = summarizeAndPlotTempDelta1ByPrior(database_id = db,
                                          method = me,
                                          exposures = allExposure_ids,
                                          rawPath = resultPath,
                                          summaryPath = outputPath,
                                          analysis_id = aid,
                                          exposure_id = eid,
                                          prior_ids = c(1:3), # include all priors for easier query later
                                          alpha = 0.05,
                                          cachePath = cachepath,
                                          useAdjusted = adjust,
                                          showPlots = FALSE,
                                          stratifyByEffectSize = TRUE,
                                          calibrate = FALSE, 
                                          outcomesInEstimates = NULL)
  
  errors_combined = rbind(errors_combined, 
                          res %>% filter(prior_id == pid) %>%
                            mutate(hyperLabel = title))
  
  # pp = attr(res, 'plot')
  # print(pp + labs(title = title))
}

## the prior choices only make some difference in earlier periods
## but the endpoints look very similar!

## perhaps sdPriors = c(0.5, 0.5) has the best performance
## (in hitting 0.05 Type I and having reasonable power)

# combined error plot ----
yinters = 0.05
type2cols = c(wes_palette("Zissou1")[3:4],wes_palette("Royal1")[4])
othercols =  wes_palette("Royal1")[2]
allCols = c(othercols, type2cols)

period_breaks = seq(from = min(errors_combined$period_id),
                    to = max(errors_combined$period_id),
                    by = 2)
period_labels = as.integer(period_breaks)

capt = sprintf('%s analysis %s, on %s',
               me, aid, db)

## (i) frequentist, raw Bayes, adjusted Bayes
p = ggplot(errors_combined, 
           aes(x=period_id, y=y, color=stats))+
  geom_line(size = 1.5) +
  geom_point(size=2)+
  geom_hline(yintercept = yinters, 
             color = 'gray60', 
             size = 1, linetype=2)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(breaks = period_breaks, labels = period_labels)+
  labs(x='analysis period (months)', y='error rates', 
       caption = capt, color='Error type')+
  scale_color_manual(values = allCols) +
  facet_grid(.~hyperLabel)+
  theme_bw(base_size = 13)+
  theme(legend.position = 'bottom') # change to bottom legend...

print(p)



##### don't run-----
# ## set suffix
# dir_suffix = 'default2'
# 
# ##
# 
# 
# 
# ## default prior setting 
# summDefault = pullResults(database_id = db,
#                           method = me,
#                           exposure_id = allExposure_ids,
#                           resultsPath = resultPath,
#                           savePath = outputPath)
# 
# res_default_raw = plotTempDelta1ByPriors(database_id = db,
#                                          method = me, 
#                                          analysis_id = aid,
#                                          exposure_id = eid,
#                                          prior_ids = c(1:3), # include all priors for easier query later
#                                          alpha = 0.05,
#                                          summaryPath = summarypath,
#                                          cachePath = cachepath,
#                                          useAdjusted = FALSE,
#                                          showPlots = TRUE,
#                                          stratifyByEffectSize = TRUE,
#                                          calibrate = FALSE, 
#                                          outcomesInEstimates = NULL)


