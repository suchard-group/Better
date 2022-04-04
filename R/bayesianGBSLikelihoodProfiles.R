# April 2022
# bayesian analyses using likelihood profiles for GBS outcome
# (similar to `bayesianWithLikelihoodProfiles.R`)

library(EvidenceSynthesis)
# source('./extras/getLikelihoodProfile.R')
# source('./extras/fitNegativeControlDistribution.R')
# source('./extras/helperFunctions.R')


# 1. run one analysis (database, exposure, method, analysis, etc.) for GBS
oneGBSAnalysis <- function(database_id,
                           method,
                           exposure_id,
                           analysis_id, 
                           period_id,
                           outcome_id = 343, # GBS id
                           priorMeans = rep(0,3),
                           priorSds = c(10,1.5,4.0), 
                           numsamps = 10000,
                           thin = 10,
                           minNCs = 5,
                           savedEstimates = NULL,
                           savedLPs = NULL){
  
  # check that pre-saved LPs are provided as a dataframe
  if(is.null(savedLPs)){
    stop('Pre-saved likelihood profiles must be provided!\n')
  }
  # check that pre-saved NC estimates are provided as well
  if(is.null(savedEstimates)){
    stop('Pre-saved negative control estimates must be provided!\n')
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
  
  # get null distribution
  null = fitNegativeControlDistribution(NULL, NULL, database_id, 
                                        method, exposure_id, 
                                        analysis_id, period_id, 
                                        savedEstimates = savedEstimates,
                                        outcomeToExclude = NULL,
                                        numsamps = numsamps, thin = thin,
                                        minNCs = minNCs)
  
  # posterior sampling for each different prior
  res = list()
  for(pr in 1:length(priorMeans)){
    priorMean = priorMeans[pr]
    priorSd = priorSds[pr]
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
    if(length(mcmc) == 1 && mcmc == 'error') next
    
    samps = mcmc$theta1
    
    # de-biasing
    if(length(null) > 0){
      adjSamps =  samps - null$bias
    }else{
      adjSamps = NA
    }
    
    # save to resls
    res[[as.character(pr)]] = list(outcome_id = outcome_id,
                                   postSamps = samps, adjustedPostSamps = adjSamps,
                                   postMean = mean(samps), postMAP = getMAP(samps), 
                                   postMedian = median(samps),
                                   adjustedPostMean = mean(adjSamps), 
                                   adjustedPostMAP = getMAP(adjSamps),
                                   adjustedPostMedian = median(adjSamps))
  }
  
  # message
  ParallelLogger::logInfo(sprintf('Finished Bayesian analysis for database %s, exposure %s, outcome %s, in period %s, using %s analysis %s.\n',
          database_id, exposure_id, outcome_id, period_id, method, analysis_id))

  # return result list
  return(res)
}


# # try it...
# LPpath = '~/Documents/Research/better_gbs/Results_CCAE/'
# LPfname = 'likelihood_profile.csv'
# LPs = getGBSLikelihoodProfiles('CCAE', 'HistoricalComparator',
#                                LPpath, LPfname)
# NCestimates = readRDS('./localCache/CompNegControls.rds')
# res = oneGBSAnalysis('CCAE', 'HistoricalComparator',
#                      211981,
#                      2, 10,
#                      savedLPs = LPs,
#                      savedEstimates = NCestimates)


# 2. functions to process one single result (for single outcome GBS)
# 2.a. summarize one outcome for one prior_id
summarizeOnePrior <- function(res, getCI = TRUE) {
  
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

# 2.b summarize over all priors
summarizeAllPriors <- function(resls, getCI = TRUE){
  if(length(resls) == 0){
    return(data.frame())
  }
  
  df = NULL
  prs = names(resls)
  for(pr in prs){
    this.pr = summarizeOnePrior(resls[[pr]], getCI)
    this.pr$prior_id = as.numeric(pr)
    df = rbind(df, this.pr)
  }
  df
}

# 2.b. retain posterior samples
getSamplesFromOneAnalysis <- function(resls){
  
  # if empty result, return NULL
  if(length(resls) == 0){
    NULL
  }else{
    samps = list()
    for(pr in names(resls)){
      samps[[pr]] = list(postSamps = resls[[pr]]$postSamps, 
                         adjustedPostSamps = resls[[pr]]$adjustedPostSamps)
    }
    samps
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
                             priors = list(Mean = rep(0, 3),
                                           Sd = c(10, 1.5, 4)),
                             returnCIs = TRUE,
                             numsamps = 10000,
                             thin = 10,
                             minNCs = 5,
                             savepath = 'localCache/testResults/',
                             sampspath = 'localCache/sampleSaves/',
                             maxCores = 4){
  # create paths if...
  if(!dir.exists(savepath)){
    dir.create(savepath)
  }
  if(!dir.exists(sampspath)){
    dir.create(sampspath)
  }
  
  # generate prior table
  priorTable = getPriorTable(priors = priors, default = TRUE)
  #prior_ids = priorTable$prior_id
  
  # save prior table if it doesn't already exist in the results folder
  if (!file.exists(file.path(savepath, 'priorTable.rds'))) {
    saveRDS(priorTable, file.path(savepath, 'priorTable.rds'))
  }
  
  
  
  # obtain LPs
  LPs = getGBSLikelihoodProfiles(database_id, method, 
                                 LPpath, LPfname,
                                 exposures = exposure_id, 
                                 analyses = analysis_ids,
                                 periods = period_ids)
  
  # get NC estimates
  ncfname = 'CompNegControls.rds'
  if(file.exists(file.path(cachePath,ncfname))){
    cat(file.path(cachePath,ncfname))
    NCestimates = readRDS(file.path(cachePath,ncfname))
  }else{
    sql <-  "SELECT database_id,
    method,
    exposure_id,
    analysis_id,
    period_id,
    estimate.outcome_id AS outcome_id,
    log_rr,
    se_log_rr
  FROM @schema.estimate
  INNER JOIN @schema.negative_control_outcome
    ON estimate.outcome_id = negative_control_outcome.outcome_id
  WHERE database_id = '@database_id'
    AND method = '@method'
    AND exposure_id = @exposure_id"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             exposure_id = exposure_id)
    NCestimates <- DatabaseConnector::querySql(connection, sql)
    names(NCestimates) = tolower(names(NCestimates))
    NCestimates = NCestimates %>% 
      filter(period_id %in% period_ids, analysis_id %in% analysis_ids) %>%
      filter(!is.na(log_rr) & !is.na(se_log_rr))
    names(NCestimates) = toupper(names(NCestimates))
  }
  
  # go through analyses and then go through periods
  ## function to process a bunch of periods for one analysis type
  runPeriods <- function(aid, saveResults = TRUE){
    
    # first check if results are already available
    if(saveResults){
      sampsFname = sprintf("%s_%s_%s_analysis%s_samples.rds",
                           database_id, method, exposure_id, aid)
      summsFname = sprintf("%s_%s_%s_analysis%s_summary.rds",
                           database_id, method, exposure_id, aid)
      if(file.exists(file.path(sampspath, sampsFname)) & file.exists(file.path(savepath, summsFname))){
        summs = readRDS(file.path(savepath, summsFname))
        return(summs)
      }
    }
    
    # if results not saved in local path, need to run analyses
    summs = NULL
    samps = list()
    for(pid in period_ids){
      # run analysis for this period_id
      res.p = oneGBSAnalysis(database_id = database_id,
                             method = method,
                             exposure_id = exposure_id,
                             analysis_id = aid,
                             period_id = pid,
                             priorMeans = priors$Mean,
                             priorSds = priors$Sd,
                             numsamps = numsamps,
                             thin = thin,
                             minNCs = minNCs,
                             savedEstimates = NCestimates,
                             savedLPs = LPs)
      # get summary and add period_id info
      summ.p = summarizeAllPriors(res.p, getCI = returnCIs)
      if(nrow(summ.p) > 0){
        summ.p$period_id = pid
      }
      summs = rbind(summs, summ.p)
      # get posterior samples and save it to the list
      samps.p = getSamplesFromOneAnalysis(res.p)
      if(!is.null(samps.p)){
        samps[[as.character(pid)]] = samps.p
      }
    }
    
    # add analysis_id info
    if(nrow(summs) > 0){
      summs$analysis_id = aid
    }
    
    # save results to local file
    if(saveResults){
      #sampsFname = sprintf("%s_%s_%s_analysis%s_samples.rds")
      saveRDS(samps, file.path(sampspath, sampsFname))
      saveRDS(summs, file.path(savepath, summsFname))
    }
    
    return(summs)
  }
  
  ## run things in parallel...
  cluster = ParallelLogger::makeCluster(min(4, maxCores))
  ParallelLogger::clusterRequire(cluster, 'better')
  allSummary = ParallelLogger::clusterApply(cluster,
                                            x = analysis_ids,
                                            fun = runPeriods,
                                            saveResults = TRUE)
  ParallelLogger::stopCluster(cluster)
  allSummary = bind_rows(allSummary)
  
  return(allSummary)
}

# try it...
# LPpath = '~/Documents/Research/better_gbs/Results_CCAE/'
# LPfname = 'likelihood_profile.csv'
# savepath = '~/Documents/Research/better_gbs/summary'
# sampspath = '~/Documents/Research/better_gbs/samples'
# allRes = multiGBSAnalyses(connection, 'eumaeus',
#                           database_id = 'CCAE', method = 'HistoricalComparator',
#                           exposure_id = 211981, 
#                           analysis_ids = c(1:12),
#                           period_ids = c(1:12),
#                           LPpath = LPpath,
#                           LPfname = LPfname,
#                           savepath = savepath,
#                           sampspath = sampspath,
#                           maxCores = 4)
