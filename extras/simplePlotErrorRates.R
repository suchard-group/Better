# March 4, simple plot for error rates

library(ggh4x)
library(wesanderson)
source('extras/helperFunctions.R')

# main function
simplePlotErrorRates <- function(database_id,
                                 method,
                                 resPath,          # where to find the judged files
                                 cachePath,        # where to find helper tables
                                 judgeStyle = 'H0neither', # c('strict', 'lenient', 'H0neither')
                                 exposure_ids = NULL, # select exposure_ids; if NULL, then no subsetting
                                 analysis_ids = c(2,4),
                                 delta1 = NULL, # select delta1 thresholds to show
                                 delta0 = 0.9, # select delta0 threshold to show
                                 prior_Sd = 10, # select prior SDs to show
                                 savePath = NULL, # path for saving plot file
                                 fnameSuffix='', # suffix for pdf plot file
                                 returnResults = FALSE, # if return results as dataframe
                                 saveResults = FALSE, # if save the results to savePath
                                 pHeight = 8, pWidth = 15,
                                 usePalette = wes_palette("Darjeeling2")[c(1,4,3,2)]){ # if use other color palettes
  
  # open up saved judged file (using the most fine-grained file for now)
  fname = sprintf('%s-%s-%s-errorRates-%s.rds', 
                  database_id, method, judgeStyle, 'boxplots')
  judged = readRDS(file.path(resPath, fname)) %>% ungroup()
  
  # massage HistoricalComparator results a little...
  # RIGHT NOW: disregard results for "filtered" analyses (analysis_id: 13-24)
  if(method == 'HistoricalComparator'){
    judged = judged %>% filter(analysis_id < 13)
  }
  
  
  # open up helper files
  thresholds = readRDS(file.path(cachePath, 'thresholdTable.rds')) %>%
    select(threshold_id = threshold_ID, 
           d1 = p1Threshold, 
           d0 = p0Threshold)
  priors = readRDS(file.path(cachePath, 'priorTable.rds'))
  exposures = getExposures(NULL, NULL, savepath = cachePath) %>%
    select(exposure_id, exposure_name = base_exposure_name)
  analyses = getAnalyses(NULL, NULL, savepath = cachePath) %>% 
    filter(method == !!method) %>%
    select(analysis_id, method_description = description)
  
  # work on filters
  if(is.null(exposure_ids)){exposure_ids = unique(exposures$exposure_id)}
  if(is.null(analysis_ids)){analysis_ids = unique(analyses$method)}
  if(is.null(delta1)){delta1 = unique(thresholds$d1); flag = 'byDelta1'}
  if(is.null(delta0)){delta0 = unique(thresholds$d0); flag = 'byDelta0'}
  if(is.null(prior_Sd)){prior_Sd = unique(priors$Sd); flag = 'byPrior'}
  
  # select the relevant records here
  judged = judged %>% left_join(priors, by='prior_id') %>%
    filter(analysis_id %in% analysis_ids, exposure_id %in% exposure_ids,
           d1 %in% delta1, d0 %in% delta0, Sd %in% prior_Sd) %>%
    left_join(analyses, by ='analysis_id')
  
  # group by and calculate mean error rates
  judged = judged %>% 
    left_join(exposures, by='exposure_id') %>%
    group_by(database_id, exposure_name, 
             method_description, 
             d1Label, d0Label, threshold_id, 
             Sd, errorLabel, negativeControl) %>%
    summarise(avgErrorRate = weightedAverage(errorRate, sampleSize)) %>%
    ungroup()
  
  # add labels for Type 1/2 and calibrated/uncalibrated
  judged = judged %>% 
    mutate(errorType = if_else(negativeControl, 'Type I', 'Type II'),
           Type = if_else(stringr::str_starts(errorLabel, 'Cali'), 'Calibrated', 'Uncalibrated'))
  
  # ! switch Type I error back to the positive scale...
  judged = judged %>% 
    mutate(avgErrorRate = if_else(negativeControl, 
                                  -avgErrorRate, avgErrorRate))
  
  
  # set up saving
  if(!is.null(savePath)){
    if(!dir.exists(savePath)){dir.create(savePath)}
    if(flag == 'byDelta1'){
      flagword = sprintf('delta0%.2f-priorSD%.1f',
                         delta0, prior_Sd)
    }
    if(flag == 'byDelta0'){
      flagword = sprintf('delta1%.2f-priorSD%.1f',
                         delta1, prior_Sd)
    }
    if(flag == 'byPrior'){
      flagword = sprintf('delta1%.2f-delta0%.2f',
                         delta1, delta0)
    }
    plotName = sprintf("%s-%s-%s-simpleErrorRates-%s-%s.pdf", 
                       database_id, method, judgeStyle,
                       flagword, fnameSuffix)
    
    pdf(file.path(savePath, plotName),
        height = pHeight,
        width = pWidth)
  }
  
  
  # 3 possible plots
  if(flag == 'byDelta1'){
    pg = ggplot(judged, aes(x=d1Label, y=avgErrorRate, 
                            fill=errorLabel)) +
      geom_bar(stat = 'identity', width = 0.7,
               position = position_dodge()) +
      labs(x = 'delta1', y='average error rate', 
           fill = '',
           caption = sprintf('delta0=%.2f, prior SD = %s',
                             delta0, prior_Sd))+
      scale_x_discrete(labels = as.character(sort(delta1))) +
      scale_y_continuous(limits = c(0,1)) +
      #coord_flip()+
      theme_bw(base_size = 14)+
      theme(strip.background = element_blank(),
            panel.border = element_blank(),
            legend.position = 'bottom',
            strip.text.y = element_text(angle = 0)) +
      facet_nested(method_description ~ exposure_name + errorType,
                   labeller = label_wrap_gen(width=15),
                   nest_line = element_line(linetype = 1))
  }
  if(flag == 'byDelta0'){
    pg = ggplot(judged, aes(x=d0Label, y=avgErrorRate, 
                            fill=errorLabel)) +
      geom_bar(stat = 'identity', width = 0.7,
               position = position_dodge()) +
      labs(x = 'delta0', y='average error rate', 
           fill = '',
           caption = sprintf('delta1=%.2f, prior SD = %s',
                             delta1, prior_Sd))+
      scale_x_discrete(labels = as.character(sort(delta0))) +
      scale_y_continuous(limits = c(0,1)) +
      #coord_flip()+
      theme_bw(base_size = 14)+
      theme(strip.background = element_blank(),
            panel.border = element_blank(),
            legend.position = 'bottom',
            strip.text.y = element_text(angle = 0)) +
      facet_nested(method_description ~ exposure_name + errorType,
                   labeller = label_wrap_gen(width=15),
                   nest_line = element_line(linetype = 1))
  }
  if(flag == 'byPrior'){
    pg = ggplot(judged, aes(x=as.factor(Sd), y=avgErrorRate, fill=errorLabel)) +
      geom_bar(stat = 'identity', width = 0.7,
               position = position_dodge()) +
      labs(x = 'Prior SD', y='average error rate', 
           fill = '',
           caption = sprintf('Thresholds: delta1=%.2f, delta0=%.2f',
                             delta1, delta0))+
      scale_x_discrete(breaks = as.character(prior_Sd)) +
      scale_y_continuous(limits = c(0,1)) +
      #coord_flip()+
      theme_bw(base_size = 14)+
      theme(strip.background = element_blank(),
            panel.border = element_blank(),
            legend.position = 'bottom',
            strip.text.y = element_text(angle = 0)) +
      facet_nested(method_description ~ exposure_name + errorType,
                   labeller = label_wrap_gen(width=15),
                   nest_line = element_line(linetype = 1))
  }
  
  # add color palette if
  if(!is.null(usePalette)){
    print(
      pg + 
        scale_fill_manual(values = usePalette)
      )
  }
  
  # close pdf connection
  if(!is.null(savePath)){
    dev.off()
  }
  
  # save results if
  if(saveResults & !is.null(savePath)){
    saveName = sprintf("%s-%s-%s-simpleErrorRates-%s-%s.rds", 
                       database_id, method, judgeStyle,
                       flagword, fnameSuffix)
    saveRDS(judged, file.path(savePath, saveName))
  }
  
  # return results if...
  if(returnResults){
    return(judged)
  }
  
}




### TEST CODE; DON'T RUN -----

resultspath = '~/Documents/Research/betterResults/judged/'
cachepath = './localCache/'
savepath = '~/Documents/Research/betterResults/simpleErrorPlots/'

db = 'CCAE'
mt = "SCCS"

# SCCS: 2 & 4 (= TaR 1~28, age + season adjusted SCCS & SCRI with post control window )

# x-axis: prior sds
simplePlotErrorRates(database_id = db, method = mt, 
                     resPath = resultspath, 
                     cachePath = cachepath,
                     judgeStyle = 'H0neither',
                     exposure_ids = NULL,
                     analysis_ids = c(10,12),
                     delta1 = 0.95,
                     delta0 = 0.9,
                     prior_Sd = NULL,
                     returnResults = FALSE,
                     usePalette = wes_palette("Darjeeling2")[c(1,4,3,2)])

# x-axis: delta1
simplePlotErrorRates(database_id = db, method = mt, 
                     resPath = resultspath, 
                     cachePath = cachepath,
                     judgeStyle = 'H0neither',
                     exposure_ids = NULL,
                     analysis_ids = c(6,8),
                     delta1 = NULL,
                     delta0 = 0.9,
                     prior_Sd = 1.5,
                     returnResults = FALSE,
                     usePalette = wes_palette("Darjeeling2")[c(1,4,3,2)])

# x-axis: delta0
simplePlotErrorRates(database_id = db, method = mt, 
                     resPath = resultspath, 
                     cachePath = cachepath,
                     judgeStyle = 'H0neither',
                     exposure_ids = NULL,
                     analysis_ids = c(10,12),
                     delta1 = 0.95,
                     delta0 = NULL,
                     prior_Sd = 1.5,
                     returnResults = FALSE,
                     usePalette = wes_palette("Darjeeling2")[c(1,4,3,2)]) 
## this one: supposed to have no change by delta0 with 'H0neither'


