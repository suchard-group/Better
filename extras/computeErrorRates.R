## Feb 24, 2022
## compute Type I, Type II and unknown rates

source('./extras/postProcessUtils.R')

## 1. function to compute Type I and II error rates,
#     given optional judging rules (on what counts as a mistake)
#     rates of unknowns given grouping
#     rates of contradictory decisions

computeErrorRates <- function(database_id,
                              method,
                              resPath,
                              cachePath,
                              judgeStyle = 'H0neither',
                              stratifyByExposures = FALSE,
                              stratifyByAnalyses = FALSE,
                              stratifyByPriors = TRUE,
                              stratifyByEffectSize = TRUE,
                              plot = TRUE){
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
  }
  # otherwise, not sure what to do yet...
  
  
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
    return()
  }
  
  # (3) stratified by exposures and analyses (similar to EUMAEUS plot)
  #     need multiple ones 
  if()
  
  
  
}

computeErrorRates(database_id = 'IBM_MDCD',
                  method = "SCCS",
                  resPath = '~/Documents/Research/betterResults/summary/',
                  cachePath = './localCache/',
                  judgeStyle = 'lenient',
                  stratifyByEffectSize = FALSE)


  