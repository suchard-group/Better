## Feb 24, 2022
## compute Type I, Type II and unknown rates

source('./extras/postProcessUtils.R')
source('./extras/helperFunctions.R')

## 1. function to compute Type I and II error rates,
#     given optional judging rules (on what counts as a mistake)
#     rates of unknowns given grouping
#     rates of contradictory decisions


## TO DO
## 1. enable returnResults option
## 2. enable saving the plots to a pdf file with prefix like "database-method-judgeStyle-......pdf"
## 3. sort out what to do with only exposure stratification (variation within a method is small)
##    y-var has delta 1 on it, and facet by delta0 (if NOT "H0neither") OR by prior (if "H0neither")
##    print multiple plots for diff. priors (if NOT "H0neither")

computeErrorRates <- function(database_id,
                              method,
                              resPath,
                              cachePath,
                              savePath,
                              fnameSuffix='errorRates',
                              judgeStyle = 'H0neither',
                              stratifyByExposures = FALSE,
                              stratifyByAnalyses = FALSE,
                              stratifyByPriors = TRUE,
                              stratifyByEffectSize = FALSE,
                              returnResults = FALSE){
  # read in file
  fname = sprintf('AllDecisions-%s-%s.rds', database_id, method)
  res = readRDS(file.path(resPath, fname))
  
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
  
  
  
  # (1) stratified by thresholds and priors, only Type I and Type II error rate
  if(stratifyByPriors & !stratifyByEffectSize){
    d0.labs = paste0('d0=', c('0.9', '0.95', '0.99'))
    names(d0.labs) = c('0.9', '0.95', '0.99')
    d1.labs = paste0('d1=', c('0.8', '0.9', '0.95'))
    names(d1.labs) = c('0.8', '0.9', '0.95')
    ybrs = c(-1, -0.5, 0, 0.5, 1)
    ylbs = c(1, 0.5, 0, 0.5, 1)
    
    if (judgeStyle == 'H0neither') {
      ## a. if "H0neither", no need to stratify by d0
      print(
        ggplot(judged,
               aes(
                 x = as.factor(d1),
                 y = errorRate,
                 fill = as.factor(errorLabel)
               )) +
          geom_hline(yintercept = 0, size = 1) +
          geom_boxplot(outlier.size = 1) +
          scale_x_discrete(breaks = d1.labs, labels = names(d1.labs)) +
          scale_y_continuous(breaks = ybrs, labels = ylbs) +
          labs(x = 'delta1', y = 'Error Rate', fill = '') +
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
      print(
        ggplot(judged,
               aes(
                 x = as.factor(d1),
                 y = errorRate,
                 fill = as.factor(errorLabel)
               )) +
          geom_hline(yintercept = 0, size = 1) +
          geom_boxplot(outlier.size = 1) +
          scale_x_discrete(breaks = names(d1.labs), labels = d1.labs) +
          scale_y_continuous(breaks = ybrs, labels = ylbs) +
          labs(x = 'delta1', y = 'Error Rate', fill = '') +
          coord_flip() +
          annotate(
            geom = 'text',
            y = -0.5,
            x = 0.6,
            label = 'Type I'
          ) +
          annotate(
            geom = 'text',
            y = 0.5,
            x = 0.6,
            label = 'Type II'
          ) +
          facet_grid(priorLabel ~ d0,
                     labeller = labeller(d0 = d0.labs)) +
          theme_bw(base_size = 14)+
          theme(strip.background = element_blank(),
                panel.border = element_blank(),
                legend.position = 'bottom')
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
    
    ## get exposure table first
    exposures = getExposures(NULL, NULL, savepath = cachePath) %>%
      select(exposure_id, exposure_name = base_exposure_name)
    
    ## group by "base exposure name" here...
    judged = judged %>% left_join(exposures)
    
    ## (a) stratify by analyses as well? 
    if(stratifyByAnalyses){
      ## get analyses table first 
      analyses = getAnalyses(NULL, NULL, savepath = cachePath) %>%
        filter(method == method) %>% 
        select(method, method_description = description)
      judged = judged %>% left_join(analyses)
      
      ## group by both exposure and analyses
      judged = judged %>% 
        group_by(database_id, exposure_name, method_description, 
                 d1, d0, threshold_id,
                 priorLabel, errorLabel) %>%
        summarise(avgErrorRate = weightedAverage(errorRate, sampleSize),
                  avgNeitherRate = weightedAverage(neitherRate, sampleSize))
      
      ## massage the data frame a bit to get ranking of method description string
      bar_shift = 0.3
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
      #for(pl in priorLabels){
        for(tid in threshold_ids){
          judged.sel = judged %>% filter(threshold_id == tid)
          d1.this = judged.sel$d1[1]
          d0.this = judged.sel$d0[1]
          print(
            ggplot(judged.sel, aes(x = method_rank, 
                                   y = avgErrorRate,
                                   fill = errorLabel)) +
              geom_hline(yintercept = 0, size = 1) +
              # geom_bar(stat='identity', position='dodge2') +
              geom_rect(aes(ymin = pmin(0,avgErrorRate),
                            ymax = pmax(0,avgErrorRate),
                            xmin = shifted_method_rank - bar_shift,
                            xmax = shifted_method_rank,
                            fill = errorLabel)) +
              scale_y_continuous(breaks = ybrs, labels = ylbs, limits = c(-1.0, 1.0)) +
              scale_x_continuous(breaks = method_ranks, labels = method_labels) +
              labs(x = '', y = 'Error Rate', fill = '',
                   caption = sprintf('Thresholds: delta1=%s, delta0=%s', 
                                     d1.this, d0.this)) +
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
              facet_grid(priorLabel ~ exposure_name,
                         labeller = label_wrap_gen(width=15)) +
              theme_bw(base_size = 14) +
              theme(strip.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = 'bottom')
          )
        }
      #}
        
    }
    
  }
  
  
  
}




computeErrorRates(database_id = 'IBM_MDCD',
                  method = "SCCS",
                  resPath = '~/Documents/Research/betterResults/summary/',
                  cachePath = './localCache/',
                  judgeStyle = 'lenient',
                  stratifyByEffectSize = FALSE)


pdf('~/Documents/Research/betterResults/ErrorRates.pdf', 
    height = 8, width = 15)
computeErrorRates(database_id = 'IBM_MDCD',
                            method = "SCCS",
                            resPath = '~/Documents/Research/betterResults/summary/',
                            cachePath = './localCache/',
                            stratifyByExposures = TRUE,
                            stratifyByAnalyses = TRUE,
                            stratifyByPriors = FALSE)
dev.off()


  