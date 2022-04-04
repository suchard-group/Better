# April 2022
# bayesian analyses using likelihood profiles for GBS outcome
# (similar to `bayesianWithLikelihoodProfiles.R`)

library(EvidenceSynthesis)
source('./extras/getLikelihoodProfile.R')
source('./extras/fitNegativeControlDistribution.R')
source('./extras/helperFunctions.R')


# 1. run one analysis (database, exposure, method, analysis, etc.) for GBS
oneGBSAnalysis <- function(connection,
                           schema,
                           database_id,
                           method,
                           exposure_id,
                           analysis_id, 
                           period_id,
                           outcome_id = 343, # GBS id
                           # IPCtable = NULL,
                           # savedEstimates = NULL,
                           priorMean = 0,
                           priorSd = 1, 
                           numsamps = 10000,
                           thin = 10,
                           minNCs = 5,
                           cachePath = NULL,
                           savedLPs = NULL){
  
  # check that pre-saved LPs are provided as a dataframe
  if(is.null(savedLPs)){
    stop('Pre-saved likelihood profiles much be provided!\n')
  }
  
  # select relevant entry of LP
  # if(nrow(savedLPs) == 0){
  #   return(numeric(0))
  # }
  lik = selectLikelihoodProfileEntry(savedLPs, 
                                     database_id,
                                     method,
                                     exposure_id,
                                     period_id,
                                     outcome_id,
                                     analysis_id)
  if(length(lik) == 0){
    return(numeric(0))
  }
  
  # get NCs and null distribution
  if(!is.null(cachePath)){
    savedEstimates = readRDS(file.path(cachePath, 'CompNegControls.rds'))
  }else{
    savedEstimates = NULL
  }
  null = fitNegativeControlDistribution(connection, schema, database_id, 
                                        method, exposure_id, 
                                        analysis_id, period_id, 
                                        savedEstimates = savedEstimates,
                                        outcomeToExclude = NULL,
                                        numsamps = numsamps, thin = thin,
                                        minNCs = minNCs)
  
  # if(length(null) > 0){
  #   biases = null$bias
  # }else{
  #   biases = NA # set sampled biases to NA if null is unavailable
  # }
  
  # posterior sampling
  mcmc = tryCatch(
    expr = {approximateSimplePosterior(
      lik,
      chainLength = numsamps * thin + 1e5,
      burnIn = 1e5, # (use default burn-in = 1e5)
      subSampleFrequency = thin,
      priorMean = priorMean,
      priorSd = priorSd
    )},
    error = function(e){
      ParallelLogger::logInfo('Error occurred while trying to run MCMC! Skipped...\n\n')
      'error'
    }
  )
  if(length(mcmc) == 1 && mcmc == 'error') return(numeric(0))
  
  samps = mcmc$theta1
  
  # de-biasing
  if(length(null) > 0){
    adjSamps =  samps - null$bias
  }else{
    adjSamps = NA
  }
  
  # return result list
  # update: add outcome_id here
  res = list(outcome_id = outcome_id,
             postSamps = samps, adjustedPostSamps = adjSamps,
             postMean = mean(samps), postMAP = getMAP(samps), 
             postMedian = median(samps),
             adjustedPostMean = mean(adjSamps), 
             adjustedPostMAP = getMAP(adjSamps),
             adjustedPostMedian = median(adjSamps))
  
  # message
  ParallelLogger::logInfo(sprintf('Finished Bayesian analysis for database %s, exposure %s, in period %s, using %s analysis %s, with priorMean=%s, priorSd=%s\n',
          database_id, exposure_id, period_id, method, analysis_id, priorMean, priorSd))

  return(res)
}


# try it...
# res = oneGBSAnalysis(NULL, NULL, 
#                      'CCAE', 'HistoricalComparator',
#                      211981,
#                      2, 10, 
#                      cachePath = './localCache/', 
#                      savedLPs = LPs)


# 2. functions to process one single result (for single outcome GBS)
# 2.a. summarize
summarizeOneOutcome <- function(res, getCI = TRUE) {
  
  # get CIs and posterior probabilities of H1 and H0
  for (s in c('postSamps', 'adjustedPostSamps')) {
    prefix = ifelse(s == 'postSamps', '', 'adjusted')
    # 95% credible intervals
    if (getCI) {
      res[[paste0(prefix, 'CI95_lb')]] = quantile(res[[s]], 0.025) %>% as.numeric()
      res[[paste0(prefix, 'CI95_ub')]] = quantile(res[[s]], 0.975) %>% as.numeric()
    }
    # posterior hypothesis probs
    res[[paste0(prefix, 'P1')]] = mean(res[[s]] > 0)
    res[[paste0(prefix, 'P0')]] = mean(res[[s]] <= 0)
    
    # remove the samples
    res[[s]] = NULL
  }
  
  # put together summary into a little 1-row dataframe
  as.data.frame(res)
}

# 2.b. retain posterior samples
getSamplesFromOneAnalysis <- function(resls){
  
  # if empty result, return NULL
  if(length(resls) == 0){
    NULL
  }else{
    list(postSamps = resls$postSamps, 
         adjustedPostSamps = resls$adjPostSamps)
  }
}




# 3. large-scale function to run multiple analyses
multiGBSAnalyses <- function(connection,
                             schema,
                             database_id,
                             method,
                             exposure_id,
                             analysis_ids,
                             period_ids,
                             LPpath,
                             LPfname,
                             cachePath = './localCache/',
                             IPCtable = NULL,
                             priors = list(Mean = rep(0, 3),
                                           Sd = c(10, 1.5, 4)),
                             returnCIs = TRUE,
                             numsamps = 10000,
                             thin = 10,
                             preLearnNull = TRUE,
                             preSaveEstimates = TRUE,
                             savepath = 'localCache/testResults/',
                             sampspath = 'localCache/sampleSaves/',
                             removeTempSummary = TRUE){
  
  
}