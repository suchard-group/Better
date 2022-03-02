# Feb 28, plot the "neither" decision rates
# Hope to see that it declines over time (on average)

## 1. a very simple function that does aggregated neither rates plotting over periods
simpleNeitherRate <- function(database_id,
                              method,
                              resPath,          # where to find the summary files
                              cachePath,        # where to find all the helper data tables
                              whichThresholds = c(1:9), # the threshold_ids to check out
                              baseExposures = TRUE, # if only stratify by base exposure names
                              savePath = NULL, # path for saving plot file
                              fnameSuffix='', # suffix for pdf plot file
                              returnResults = FALSE, # if return results as dataframe
                              saveResults = FALSE, # if save the results to savePath
                              pHeight = 8, pWidth = 15){
  # first check if savePath is provided when saveResults=TRUE
  if(saveResults & is.null(savePath)){
    stop('Must provide a `savePath` when saveResults=TRUE!!!\n')
  }
  
  # read in file
  fname = sprintf('AllSummary-%s-%s.rds', database_id, method)
  res = readRDS(file.path(resPath, fname))
  
  # massage HistoricalComparator results a little...
  # RIGHT NOW: disregard results for "filtered" analyses (analysis_id: 13-24)
  if(method == 'HistoricalComparator'){
    res = res %>% filter(analysis_id < 13)
  }
  
  # read thresholds table
  # and a bit of processing
  thresholds = readRDS(file.path(cachePath, 'thresholdTable.rds')) %>%
    select(threshold_id = threshold_ID, 
           d1 = p1Threshold, 
           d0 = p0Threshold) %>%
    arrange(threshold_id) %>%                    # make sure it's ordered by threshold_id
    mutate(d1Label = sprintf('delta1=%.2f',d1),
           d0Label = sprintf('delta0=%.2f',d0))  # create the character label here too
  
  # ... and exposures table as well
  exposures = getExposures(NULL, NULL, savepath = cachePath) %>%
    select(exposure_id, exposure_name, base_exposure_name)
  
  # if only stratify by base exposure names...
  if(baseExposures){
    exposures = exposures %>% select(exposure_id, exposure_name = base_exposure_name)
  }
  
  
  # join with exposure table
  res = res %>% left_join(exposures, by='exposure_id')
  
  # compute neither rates and plot with each threshold combo
  ## set up pdf path if want to save
  if(!is.null(savePath)){
    plotPath = file.path(savePath, 'plots/')
    ## make sure plotPath exists
    if(!dir.exists(plotPath)){
      dir.create(plotPath)
    }
    plotName = sprintf('%s-%s-neitherRates-%s.pdf', 
                       database_id, method,
                       fnameSuffix)
    pdf(file = file.path(plotPath, plotName),
        height = pHeight, width = pWidth)
  }
  
  ## loop through the thresholds...
  NRs = NULL
  npr = max(res$period_id)
  for(id in whichThresholds){
    d1 = thresholds$d1[id]
    d0 = thresholds$d0[id]
    
    # need to split by calibrated & uncalibrated results
    this.NRs = res %>% 
      mutate(unknown = (P1 < d1 & P0 < d0),
             adjustedUnknown = (adjustedP1 < d1 & adjustedP0 < d0)) %>%
      group_by(period_id, exposure_name) %>% 
      summarise(neitherRate = mean(unknown),
                adjustedNeitherRate = mean(adjustedUnknown))
    
    this.NRs = bind_rows(this.NRs %>% select(-adjustedNeitherRate) %>% 
                           mutate(Type = 'Uncalibrated'),
                         this.NRs %>% select(period_id, exposure_name, 
                                             neitherRate = adjustedNeitherRate) %>%
                           mutate(Type = 'Calibrated')) %>%
      mutate(d1Label = thresholds$d1Label[id], d0Label = thresholds$d0Label[id])
    
    # plot...
    pg = ggplot(this.NRs, 
                aes(x=period_id, y=neitherRate, color=exposure_name)) +
      geom_line() +
      scale_x_continuous(breaks = seq(from=2, to=npr, by=2)) +
      scale_y_continuous(limits = c(0,1)) +
      labs(x='Months in study', y='Rate of indecision', 
           color='', 
           caption = sprintf('delta1=%.2f, delta0=%.2f', d1, d0)) +
      facet_grid(.~Type) +
      theme_bw(base_size = 14) +
      theme(legend.position = 'bottom')
    
    print(pg)
    
    # combine results to big dataframe
    NRs = bind_rows(NRs, this.NRs)
  }
  
  if(!is.null(savePath)){dev.off()}
  
  # save the results if...
  if(saveResults){
    saveResPath = file.path(savePath, 'neitherRatesByPeriod/')
    if(!dir.exists(saveResPath)){dir.create(saveResPath)}
    saveName = sprintf('%s-%s-neitherRates-%s.rds', 
                       database_id, method,fnameSuffix)
    saveRDS(NRs, file.path(saveResPath, saveName))
    cat(sprintf('\nNeither rates summary file saved for %s, %s!\nFile can be found at %s\n\n',
                database_id, method, file.path(saveResPath, saveName)))
  }
  
  # return results if...
  if(returnResults){
    return(NRs)
  }
}


### RUN CODE USING FUNCTION ------
# resultspath = '~/Documents/Research/betterResults/summary/'
# cachepath = './localCache/'
# savepath = '~/Documents/Research/betterResults/'
# 
# # db = 'IBM_MDCD'
# # mt = "SCCS"
# 
# thresholdsToCheck = c(1:9)
# toSaveFile = TRUE
# baseExpo = TRUE
# 
# 
# databases = c('IBM_MDCD','CCAE','MDCR','OptumDod','OptumEhr')
# methods = c('SCCS','HistoricalComparator')
# 
# for(db in databases){
#   for(mt in methods){
#     simpleNeitherRate(database_id = db, 
#                       method = mt,
#                       resPath = resultspath,
#                       cachePath = cachepath,
#                       savePath = savepath,
#                       fnameSuffix = 'byPeriod',
#                       baseExposures = baseExpo,
#                       whichThresholds = thresholdsToCheck,
#                       returnResults = FALSE,
#                       saveResults = toSaveFile,
#                       pHeight = 5, pWidth = 9)
#   }
# }




### TESTING CODE BELOW, DON't RUN -------------

# ## load an example result summary
# summ = readRDS('~/Documents/Research/betterResults/summary/AllSummary-CCAE-SCCS.rds')
# 
# ## get the "absolute unknown" column, with P1 < 0.8 and P0 < 0.9
# summ = summ %>% mutate(unknown = (P1 < 0.8 & P0 < 0.9))
# 
# ## (1) only stratify by exposure id
# neitherRates = summ %>% 
#   group_by(period_id, exposure_id) %>% 
#   summarise(neitherRate = mean(unknown))
# 
# ggplot(neitherRates, 
#        aes(x=period_id, y=neitherRate, color=as.factor(exposure_id))) +
#   geom_line()
# 
# ## (1.b) stratify by base exposure name and see if different
# exposures = readRDS('./localCache/exposures.rds') %>%
#   select(exposure_id, exposure_name = base_exposure_name)
# 
# summ = summ %>% left_join(exposures, by='exposure_id')
# 
# neitherRates = summ %>% 
#   group_by(period_id, exposure_name) %>% 
#   summarise(neitherRate = mean(unknown))
# 
# ggplot(neitherRates, 
#        aes(x=period_id, y=neitherRate, color=exposure_name)) +
#   geom_line()
# ### some differences between using exposure_id's and exposure_name's
# ### but overall the trend is same: declined unknown rate by period_id
# 
# 
# ## (2) stratify by priors also
# neitherRates = summ %>% 
#   group_by(period_id, exposure_id, prior_id) %>% 
#   summarise(neitherRate = mean(unknown))
# 
# ggplot(neitherRates, 
#        aes(x=period_id, y=neitherRate, color=as.factor(exposure_id))) +
#   geom_line() +
#   facet_grid(prior_id ~.) # there isn't much difference between different priors!
#   