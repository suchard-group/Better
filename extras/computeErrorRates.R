## Feb 24, 2022
## compute Type I, Type II and unknown rates

source('./extras/postProcessUtils.R')
source('./extras/helperFunctions.R')

## 1. function to compute Type I and II error rates,
#     given optional judging rules (on what counts as a mistake)
#     rates of unknowns given grouping
#     rates of contradictory decisions


## NOTE: for HC, filtered analyses 13-24 (those filterd things)
##       since I didn't do them for most databases...

computeErrorRates <- function(database_id,
                              method,
                              resPath,          # where to find the summary files (decisions in them)
                              cachePath,        # where to find helper tables
                              judgeStyle = 'H0neither', # c('strict', 'lenient', 'H0neither')
                              stratifyByExposures = FALSE,
                              stratifyByAnalyses = FALSE,
                              stratifyByPriors = TRUE, # if TRUE, average over exposures and analyses
                              stratifyByEffectSize = FALSE, # option not implemented yet
                              savePath = NULL, # path for saving plot file
                              fnameSuffix='', # suffix for pdf plot file
                              returnResults = FALSE, # if return results as dataframe
                              saveResults = FALSE, # if save the judged results to savePath
                              pHeight = 8, pWidth = 15){
  # first check if savePath is provided when saveResults=TRUE
  if(saveResults & is.null(savePath)){
    stop('Must provide a `savePath` when saveResults=TRUE!!!\n')
  }
  
  # read in file
  fname = sprintf('AllDecisions-%s-%s.rds', database_id, method)
  res = readRDS(file.path(resPath, fname))
  
  # massage HistoricalComparator results a little...
  # RIGHT NOW: disregard results for "filtered" analyses (analysis_id: 13-24)
  if(method == 'HistoricalComparator'){
    res = res %>% filter(analysis_id < 13)
  }
  
  # read in prior table and thresholds table
  # and a bit of processing
  thresholds = readRDS(file.path(cachePath, 'thresholdTable.rds')) %>%
    select(threshold_id = threshold_ID, 
           d1 = p1Threshold, 
           d0 = p0Threshold)
  priors = readRDS(file.path(cachePath, 'priorTable.rds')) %>%
    mutate(priorLabel = sprintf('Mean=%s, SD=%s', Mean, Sd))
  
  # judge the resulting decisions
  judged = judgeDecisions(res, judgeStyle)

  # join with helper tables
  judged = judged %>% left_join(thresholds) %>%
    left_join(priors) %>%
    select(-Mean, -Sd)
  
  # long format into adjusted and unadjusted categories
  judged.unadj = judged %>% 
    select(-adjustedErrorRate, -adjustedNeitherRate) %>%
    mutate(errorLabel = ifelse(negativeControl, 
                               'Uncalibrated Type I',
                               'Uncalibrated Type II'))
  judged.adj = judged %>% select(-errorRate, -neitherRate) %>%
    mutate(errorLabel = ifelse(negativeControl, 
                               'Calibrated Type I',
                               'Calibrated Type II'))
  names(judged.adj) = names(judged.unadj) # change the names for consistency
  
  judged = bind_rows(judged.unadj, judged.adj)
  
  
  
  # If not stratified by effect size: make Type I error rate "negative" so it shows up on the other side
  if(!stratifyByEffectSize){
    judged = judged %>% 
      mutate(errorRate = ifelse(negativeControl, -errorRate, errorRate))
  }else{
    # otherwise, not sure what to do yet...
    # TBD...
    invisible()
  }
  
  
  # start plotting
  # set up pdf path if want to save
  if(!is.null(savePath)){
    plotPath = file.path(savePath, 'plots/')
    ## make sure plotPath exists
    if(!dir.exists(plotPath)){
      dir.create(plotPath)
    }
    plotName = sprintf('%s-%s-%s-errorRates-%s.pdf', 
                       database_id, method, judgeStyle,
                       fnameSuffix)
    pdf(file = file.path(plotPath, plotName),
        height = pHeight, width = pWidth)
  }
  
  # (1) stratified by thresholds and priors, only Type I and Type II error rate
  if(stratifyByPriors & !stratifyByEffectSize){
    # d0.labs = sprintf('delta0=%.2f', sort(unique(thresholds$d0)))
    # names(d0.labs) = c('0.9', '0.95', '0.99')
    # d1.labs = sprintf('delta0=%.2f', sort(unique(thresholds$d1)))
    # names(d1.labs) = c('0.8', '0.9', '0.95')
    judged = judged %>%
      mutate(d1Label = sprintf('delta1=%.2f', d1),
             d0Label = sprintf('delta0=%.2f', d0))
    ybrs = c(-1, -0.5, 0, 0.5, 1)
    ylbs = c(1, 0.5, 0, 0.5, 1)
    
    if (judgeStyle == 'H0neither') {
      ## a. if "H0neither", no need to stratify by d0
      print(
        ggplot(judged,
               aes(
                 x = d1Label,
                 y = errorRate,
                 fill = errorLabel
               )) +
          geom_hline(yintercept = 0, size = 1) +
          geom_boxplot(outlier.size = 1) +
          #scale_x_discrete(breaks = d1.labs, labels = names(d1.labs)) +
          scale_y_continuous(breaks = ybrs, labels = ylbs) +
          labs(x = '', y = 'Error Rate', fill = '') +
          coord_flip() +
          annotate(
            geom = 'text',
            y = -0.5,
            x = 0.5,
            label = 'Type I'
          ) +
          annotate(
            geom = 'text',
            y = 0.5,
            x = 0.5,
            label = 'Type II'
          ) +
          facet_grid(. ~ priorLabel) +
          theme_bw(base_size = 14) +
          theme(strip.background = element_blank(),
                panel.border = element_blank(),
                legend.position = 'bottom')
      )
    } else{
      ## b. if "strict" or "lenient", also stratify by d0
      ## make up little dummy table for annotation
      prior_labels = sort(unique(judged$priorLabel))
      d0_labels = sort(unique(judged$d0Label))
      annos = data.frame(priorLabel = rep(prior_labels, length(d0_labels)),
                         exposure_name = rep(d0_labels, each = length(prior_labels)),
                         label1 = rep(c('','','Type I'), length(d0_labels)),
                         label2 = rep(c('','','Type II'), length(d0_labels)))
      print(
        ggplot(judged,
               aes(
                 x = d1Label,
                 y = errorRate
               )) +
          geom_hline(yintercept = 0, size = 1) +
          geom_boxplot(outlier.size = 1, aes(fill = errorLabel)) +
          #scale_x_discrete(breaks = names(d1.labs), labels = d1.labs) +
          scale_y_continuous(breaks = ybrs, labels = ylbs) +
          labs(x = '', y = 'Error Rate', fill = '') +
          coord_flip() +
          # annotate(
          #   geom = 'text',
          #   y = -0.5,
          #   x = 0.6,
          #   label = 'Type I'
          # ) +
          # annotate(
          #   geom = 'text',
          #   y = 0.5,
          #   x = 0.6,
          #   label = 'Type II'
          # ) +
          facet_grid(priorLabel ~ d0Label) +
          theme_bw(base_size = 14)+
          theme(strip.background = element_blank(),
                panel.border = element_blank(),
                legend.position = 'bottom')+
          geom_text(data=annos, x = 0.55, y = -0.5, aes(label=label1), size = 3.5) +
          geom_text(data=annos, x = 0.55, y = 0.5, aes(label=label2), size = 3.5)
      )
    }
  }
  
  # (2) stratified by thresholds and priors, error rates by effect size as well
  if(stratifyByPriors & stratifyByEffectSize){
    # skip this one for now
    invisible()
  }
  
  # (3) stratified by exposures 
  #      and analyses (optional -- analyses within same method don't differ that much)
  #      need multiple plots output
  if(stratifyByExposures){
    ybrs = c(-1, -0.5, 0, 0.5, 1)
    ylbs = c(1, 0.5, 0, 0.5, 1)
    bar_shift = 0.3
    
    ## get exposure table first
    exposures = getExposures(NULL, NULL, savepath = cachePath) %>%
      select(exposure_id, exposure_name = base_exposure_name)
    
    ## group by "base exposure name" here...
    judged = judged %>% left_join(exposures)
    
    ## (a) stratify by analyses as well? 
    if(stratifyByAnalyses){
      ## get analyses table first 
      analyses = getAnalyses(NULL, NULL, savepath = cachePath) %>%
        # mutate(METHOD = method) %>%
        # filter(METHOD == method) %>% 
        select(method, analysis_id, method_description = description)
      judged = judged %>% left_join(analyses, by=c('method','analysis_id'))
      
      ## group by both exposure and analyses
      judged = judged %>% 
        group_by(database_id, exposure_name, method_description, 
                 d1, d0, threshold_id,
                 priorLabel, errorLabel) %>%
        summarise(avgErrorRate = weightedAverage(errorRate, sampleSize),
                  avgNeitherRate = weightedAverage(neitherRate, sampleSize))
      
      ## massage the data frame a bit to get ranking of method description string
      judged$method_rank = getRank(judged$method_description)
      judged = judged %>% 
        #mutate(method_rank = getRank(.data$method_description)) %>%
        mutate(shifted_method_rank = if_else(stringr::str_starts(errorLabel, 'Uncalibrated'),
                                             method_rank + 0.0, # uncalibrated bar on below
                                             method_rank + bar_shift)) # calibrated bar on top
      
      method_ranks = sort(unique(judged$method_rank))
      method_labels = sort(unique(judged$method_description))
      
      ## make a lot of plots
      threshold_ids = sort(unique(judged$threshold_id))
      #priorLabels = sort(unique(judged$priorLabel))
      
      ## make up little dummy table for annotation
      prior_labels = sort(unique(judged$priorLabel))
      exposure_names = sort(unique(judged$exposure_name))
      annos = data.frame(priorLabel = rep(prior_labels, length(exposure_names)),
                         exposure_name = rep(exposure_names, each = length(prior_labels)),
                         label1 = rep(c('Type I','',''), length(exposure_names)),
                         label2 = rep(c('Type II','',''), length(exposure_names)))
      
      #for(pl in priorLabels){
        for(tid in threshold_ids){
          judged.sel = judged %>% filter(threshold_id == tid)
          d1.this = judged.sel$d1[1]
          d0.this = judged.sel$d0[1]
          print(
            ggplot(judged.sel, aes(x = method_rank, 
                                   y = avgErrorRate)) +
              # geom_hline(yintercept = 0, size = 1) +
              geom_segment(x = min(method_ranks)-bar_shift - 0.2, 
                           xend = max(method_ranks)+bar_shift + 0.2,
                           y = 0, yend = 0, size = 1) +
              geom_rect(aes(ymin = pmin(0,avgErrorRate),
                            ymax = pmax(0,avgErrorRate),
                            xmin = shifted_method_rank - bar_shift,
                            xmax = shifted_method_rank,
                            fill = errorLabel)) +
              scale_y_continuous(breaks = ybrs, labels = ylbs, limits = c(-1.0, 1.0)) +
              scale_x_continuous(breaks = method_ranks, labels = method_labels,
                                 limits = c(min(method_ranks)-bar_shift,max(method_ranks)+0.5)) +
              labs(x = '', y = 'Error Rate', fill = '',
                   caption = sprintf('Thresholds: delta1=%s, delta0=%s', 
                                     d1.this, d0.this)) +
              coord_flip() +
              facet_grid(priorLabel ~ exposure_name,
                         labeller = label_wrap_gen(width=15)) +
              theme_bw(base_size = 14) +
              theme(strip.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = 'bottom') +
              geom_text(data=annos, x = max(method_ranks)+0.5, y = -0.5, aes(label=label1), size = 3.5) +
              geom_text(data=annos, x = max(method_ranks)+0.5, y = 0.5, aes(label=label2), size = 3.5)
          )
        }
      #}
        
    }else{
      # otherwise, only group by exposures (and average over analyses)
      
      # depends on if judgeStyle is "H0neither"
      if(judgeStyle == 'H0neither'){
        # no split on d0, only group with d1
        judged = judged %>% 
          group_by(database_id, exposure_name, 
                   d1, 
                   priorLabel, errorLabel) %>%
          summarise(avgErrorRate = weightedAverage(errorRate, sampleSize),
                    avgNeitherRate = weightedAverage(neitherRate, sampleSize))
      }else{
        # group by d0 as well
        judged = judged %>% 
          group_by(database_id, exposure_name, 
                   d1, d0, threshold_id,
                   priorLabel, errorLabel) %>%
          summarise(avgErrorRate = weightedAverage(errorRate, sampleSize),
                    avgNeitherRate = weightedAverage(neitherRate, sampleSize))
      }
      
      # more post-processing on d1
      # give character label to d1 
      judged = judged %>% 
        mutate(d1Label = sprintf('delta1=%.2f',d1))
      judged$d1_rank = getRank(judged$d1Label)
      judged = judged %>% 
        mutate(shifted_d1_rank = if_else(stringr::str_starts(errorLabel, 'Uncalibrated'),
                                         d1_rank + 0.0, # uncalibrated bar on below
                                         d1_rank + bar_shift))
      
      d1_ranks = sort(unique(judged$d1_rank))
      d1_labels = sort(unique(judged$d1Label))
      
      # diff. plotting given judgeStyle
      if(judgeStyle == 'H0neither'){
        # y axis: d1; panels: prior
        ## make up little dummy table for annotation
        prior_labels = sort(unique(judged$priorLabel))
        exposure_names = sort(unique(judged$exposure_name))
        annos = data.frame(priorLabel = rep(prior_labels, length(exposure_names)),
                           exposure_name = rep(exposure_names, each = length(prior_labels )),
                           label1 = rep(c('Type I','',''), length(exposure_names)),
                           label2 = rep(c('Type II','',''), length(exposure_names)))
        
        print(
          ggplot(judged, aes(x = d1_rank, 
                             y = avgErrorRate)) +
            # geom_hline(yintercept = 0, size = 1) +
            geom_segment(x = min(d1_ranks)-bar_shift - 0.2, 
                         xend = max(d1_ranks)+bar_shift + 0.2,
                         y = 0, yend = 0, size = 1) +
            geom_rect(aes(ymin = pmin(0,avgErrorRate),
                          ymax = pmax(0,avgErrorRate),
                          xmin = shifted_d1_rank - bar_shift,
                          xmax = shifted_d1_rank,
                          fill = errorLabel)) +
            scale_y_continuous(breaks = ybrs, labels = ylbs, limits = c(-1.0, 1.0)) +
            scale_x_continuous(breaks = d1_ranks, labels = d1_labels, 
                               limits = c(min(d1_ranks)-bar_shift,max(d1_ranks)+0.5)) +
            labs(x = '', y = 'Error Rate', fill = '') +
            coord_flip() +
            # annotate(
            #   geom = 'text',
            #   y = -0.5,
            #   x = 0.5,
            #   label = 'Type I'
            # ) +
            # annotate(
            #   geom = 'text',
            #   y = 0.5,
            #   x = 0.5,
            #   label = 'Type II'
            # ) +
            facet_grid(priorLabel ~ exposure_name,
                       labeller = label_wrap_gen(width=15)) +
            theme_bw(base_size = 14) +
            theme(strip.background = element_blank(),
                  panel.border = element_blank(),
                  legend.position = 'bottom')  +
            geom_text(data=annos, x = max(d1_ranks)+0.5, y = -0.5, aes(label=label1), size = 3.5) +
            geom_text(data=annos, x = max(d1_ranks)+0.5, y = 0.5, aes(label=label2), size = 3.5)
        )
      }else{
        # y axis: d1; panels: d0; diff. plots: prior
        ## give d0 values label as well
        judged = judged %>% mutate(d0Label = sprintf('delta0=%.2f',d0))
        ## make up little dummy table for annotation
        d0_labels = sort(unique(judged$d0Label))
        exposure_names = sort(unique(judged$exposure_name))
        annos = data.frame(d0Label = rep(d0_labels, length(exposure_names)),
                           exposure_name = rep(exposure_names, each = length(d0_labels)),
                           label1 = rep(c('Type I','',''), length(exposure_names)),
                           label2 = rep(c('Type II','',''), length(exposure_names)))
        for(pl in sort(priors$priorLabel)){
          judged.sel = judged %>% filter(priorLabel == pl)
          print(
            ggplot(judged.sel, aes(x = d1_rank, 
                                   y = avgErrorRate)) +
              # geom_hline(yintercept = 0, size = 1) +
              geom_segment(x = min(d1_ranks)-bar_shift - 0.2, 
                           xend = max(d1_ranks)+bar_shift + 0.2,
                           y = 0, yend = 0, size = 1) +
              geom_rect(aes(ymin = pmin(0,avgErrorRate),
                            ymax = pmax(0,avgErrorRate),
                            xmin = shifted_d1_rank - bar_shift,
                            xmax = shifted_d1_rank,
                            fill = errorLabel)) +
              scale_y_continuous(breaks = ybrs, labels = ylbs, limits = c(-1.0, 1.0)) +
              scale_x_continuous(breaks = d1_ranks, labels = d1_labels, 
                                 limits = c(min(d1_ranks)-bar_shift,max(d1_ranks)+0.5)) +
              labs(x = '', y = 'Error Rate', fill = '',
                   caption = sprintf('Prior: %s', pl)) +
              coord_flip() +
              # annotate(
              #   geom = 'text',
              #   y = -0.5,
              #   x = 0.5,
              #   label = 'Type I'
              # ) +
              # annotate(
              #   geom = 'text',
              #   y = 0.5,
              #   x = 0.5,
              #   label = 'Type II'
              # ) +
              facet_grid(d0Label ~ exposure_name,
                         labeller = label_wrap_gen(width=15)) +
              theme_bw(base_size = 14) +
              theme(strip.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = 'bottom') +
              geom_text(data=annos, x = max(d1_ranks)+0.5, y = -0.5, aes(label=label1), size = 3.5) +
              geom_text(data=annos, x = max(d1_ranks)+0.5, y = 0.5, aes(label=label2), size = 3.5)
          )
        }
      }
    }
    
  }
  
  
  # turn off plotting device if save the plots
  if(!is.null(savePath)){dev.off()}
  
  
  # finally, return the produced judged file if...
  if(returnResults){
    return(judged)
  }
  
  # save results if ...
  if(saveResults){
    saveResPath = file.path(savePath, 'judged/')
    if(!dir.exists(saveResPath)){dir.create(saveResPath)}
    saveName = sprintf('%s-%s-%s-errorRates-%s.rds', 
                       database_id, method, judgeStyle,
                       fnameSuffix)
    saveRDS(judged, file.path(saveResPath, saveName))
    cat(sprintf('\nError rates summary file saved for %s, %s, with judge style %s!\nFile can be found at %s\n\n',
                database_id, method, judgeStyle, file.path(saveResPath, saveName)))
  }
}



#### Start making some plots--------

resultspath = '~/Documents/Research/betterResults/summary/'
cachepath = './localCache/'
savepath = '~/Documents/Research/betterResults/'

db = 'CCAE'
# mt = "SCCS" #"HistoricalComparator" 

for(mt in c('SCCS',"HistoricalComparator")){
judgements = c('strict','lenient','H0neither')

toSaveFile = TRUE

for(js in judgements){
  cat('Plotting with judging standard', js, '...\n')
  # (1) boxplots by thresholds and priors
  computeErrorRates(database_id = db,
                    method = mt,
                    resPath = resultspath,
                    cachePath = cachepath,
                    judgeStyle = js,
                    savePath = savepath,
                    fnameSuffix = 'boxplots',
                    stratifyByEffectSize = FALSE,
                    saveResults = toSaveFile,
                    pWidth = 12)
  
  # (2) stratify by both analyses and exposures, multiple plots over diff. thresholds
  computeErrorRates(database_id = db,
                    method = mt,
                    resPath = resultspath,
                    cachePath = cachepath,
                    judgeStyle = js,
                    savePath = savepath,
                    fnameSuffix = 'stratifyByAnalyses',
                    stratifyByExposures = TRUE,
                    stratifyByAnalyses = TRUE,
                    stratifyByPriors = FALSE,
                    saveResults = toSaveFile,
                    pWidth = 15)
  
  # (3) stratify by exposures only, (maybe) mutiple plots over priors
  computeErrorRates(database_id = db,
                    method = mt,
                    resPath = resultspath,
                    cachePath = cachepath,
                    judgeStyle = js,
                    savePath = savepath,
                    fnameSuffix = 'stratifyByExposures',
                    stratifyByExposures = TRUE,
                    stratifyByAnalyses = FALSE,
                    stratifyByPriors = FALSE,
                    saveResults = toSaveFile,
                    pWidth = 12)
}

}




##### TESTING CODE, NOT TO RUN------------

# computeErrorRates(database_id = 'IBM_MDCD',
#                   method = "SCCS",
#                   resPath = '~/Documents/Research/betterResults/summary/',
#                   cachePath = './localCache/',
#                   judgeStyle = 'strict',
#                   savePath = '~/Documents/Research/betterResults/plots/',
#                   fnameSuffix = 'boxplots',
#                   stratifyByEffectSize = FALSE,
#                   pWidth = 12)


# pdf('~/Documents/Research/betterResults/ErrorRates.pdf', 
#     height = 8, width = 15)
# computeErrorRates(database_id = 'IBM_MDCD',
#                             method = "SCCS",
#                             resPath = '~/Documents/Research/betterResults/summary/',
#                             cachePath = './localCache/',
#                             stratifyByExposures = TRUE,
#                             stratifyByAnalyses = TRUE,
#                             stratifyByPriors = FALSE)
# dev.off()


# pdf('~/Documents/Research/betterResults/ErrorRates1.pdf', 
#     height = 8, width = 12)
# computeErrorRates(database_id = 'IBM_MDCD',
#                   method = "SCCS",
#                   judgeStyle = 'H0neither',
#                   resPath = '~/Documents/Research/betterResults/summary/',
#                   cachePath = './localCache/',
#                   stratifyByExposures = TRUE,
#                   stratifyByAnalyses = FALSE,
#                   stratifyByPriors = FALSE)
# dev.off()

  