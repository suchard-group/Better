# Aug 2022 update
# bayesian inference for GBS with 
# meta analysis for fitting NC distribution
# from likelihood profiles
# instead of using estimates & SEs

library(EvidenceSynthesis)
source('./R/getLikelihoodProfile.R')
source('./R/fitNegativeControlDistribution.R')
source('./R/betterHelperFunctions.R')

library(foreach)
library(doParallel)
registerDoParallel()

# 1. function to carry out one single analysis (exposure-design-database) ------
# using meta analysis for negative control analysis
# Aug 2022 update: add robust option for t-model
oneGBSAnalysisMeta <- function(connection,
                               schema,
                               database_id,
                               method,
                               exposure_id,
                               analysis_id, 
                               period_id,
                               IPCtable = NULL,
                               noSyntheticPos=TRUE,
                               priorMean = 0,
                               priorSd = 1, 
                               nullPriorSds = c(0.5,0.5),
                               robustMetaAnalysis = FALSE,
                               numsamps = 10000,
                               thin = 10,
                               minNCs = 5,
                               outcomeLPPath = '~/Documents/Research/better_gbs/',
                               outcomeLPs = NULL,
                               NCLPs = NULL,
                               preSaveNull = NULL,
                               negControls = NULL){
  
  # access list of negative control outcome ids
  if(is.null(negControls)){
    sql <- "SELECT outcome_id from eumaeus.NEGATIVE_CONTROL_OUTCOME"
    NCs = DatabaseConnector::querySql(connection, sql)$OUTCOME_ID
  }else{
    NCs = negControls
  }
  
  ### LP for the outcome:
  # if the LPs for the special outcomes (GBS) are provided as dataframe, then use them
  # otherwise load from "better_gbs" directory
  if(is.null(outcomeLPs)){
    LP_path = file.path(outcomeLPPath, sprintf("Results_%s",database_id), 'likelihood_profile.csv')
    outcomeLPs = read_csv(LP_path)
  }
  # query relevant LP row
  theLP = outcomeLPs %>% filter(exposure_id == !!exposure_id,
                                method == !!method,
                                analysis_id == !!analysis_id,
                                period_id == !!period_id)
  if(nrow(theLP) == 0){
    cat('No likelihood profile for the special outcome of interest!\n')
    return(numeric(0))
  }
  

  ### LP for NCs
  # if NCLPs is provided as a dataframe, then use pre-pulled LPs, otherwise query
  if(is.null(NCLPs)){
    LPs = getMultiLikelihoodProfiles(connection, schema,
                                     database_id,
                                     exposure_id, 
                                     analysis_id,
                                     period_id,
                                     method = method)
  }else{
    LPs = NCLPs
  }
  
  # if no relevant likelihood profiles are available 
  # do nothing and return numeric(0)
  if(length(LPs) == 0){
    return(numeric(0))
  }
  
  # get likelihood profiles for NCs only
  NC_LPs = LPs %>% filter(OUTCOME_ID %in% NCs)
  
  
  # get a subset IPC table if available
  if(!is.null(IPCtable)){
    IPCtable = IPCtable %>% filter(EXPOSURE_ID == exposure_id)
  }
  
  noAdjust = FALSE
  
  # learn the negative control distribution
  if(is.null(preSaveNull)){
    null = fitNegativeControlDistributionLikelihood(connection, schema, database_id, 
                                                    method, exposure_id, 
                                                    analysis_id, period_id, 
                                                    savedLPs = NC_LPs, 
                                                    outcomeToExclude = NULL,
                                                    priorSds = nullPriorSds, 
                                                    minNCs = minNCs,
                                                    numsamps = numsamps, 
                                                    thin = thin,
                                                    robust = robustMetaAnalysis)
  }else{
    null = preSaveNull
  }
  # if no negative controls are available at all
  # set a flag so that we don't return anything adjusted
  if(length(null)==0){
    noAdjust = TRUE
    biases = NULL
  }else{
    biases = null$bias
  }
  
  # get outcomes to do (should be length-1 for now for GBS only)
  outcomesToDo = unique(theLP$outcome_id)
  res = list()
  
  for(outcome in outcomesToDo){
    
    cat('\n\nAnalysis for outcome', outcome, 'underway...\n')
    
    # get the likelihood profile for this outcome
    this.outcomeLP = theLP %>% filter(outcome_id == outcome)
    
    if(nrow(this.outcomeLP) == 0){
      cat(sprintf('No likelihood profile available for outcome %s! Skipped.\n', outcome))
      next
    }
    
    # post-process this LP row
    points = str_split(this.outcomeLP$point,';') %>% 
      unlist() %>% as.numeric()
    values = str_split(this.outcomeLP$value,';') %>% 
      unlist() %>% as.numeric()
    lik = values
    names(lik) = as.character(points)
    
    
    ## run MCMC for posterior inference
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
        cat('\n')
        print(e)
        print(traceback())
        cat('Error occurred while trying to run MCMC! Skipped...\n\n')
        'error'
      }
    )
    ## if there is error, skip this outcome and move on to next
    if(length(mcmc) == 1 && mcmc == 'error') next
    
    samps = mcmc$theta1
    # with calibration
    if(!is.null(biases) && !noAdjust){
      adjSamps =  samps - biases
    }else{
      adjSamps = NA
    }
      
      # get result to list
      # key: outcome id
      res[[as.character(outcome)]] = 
        #tempRes[[as.character(this.outcome)]] = 
        list(postSamps = samps, adjustedPostSamps = adjSamps,
             postMean = mean(samps), postMAP = getMAP(samps), 
             postMedian = median(samps),
             adjustedPostMean = mean(adjSamps), 
             adjustedPostMAP = getMAP(adjSamps),
             adjustedPostMedian = median(adjSamps))
      #cat('\nResult written to list!\n\n')
    }
  
  # return final result
  cat(sprintf('Finished Bayesian analysis for database %s, exposure %s, in period %s, using %s analysis %s, with priorMean=%s, priorSd=%s\n',
              database_id, exposure_id, period_id, method, analysis_id, priorMean, priorSd))
  res
}

# ## try it
# IPCs = readRDS('localCache/allIPCs.rds')
# outcomepath = '~/Documents/Research/better_gbs'
# resLst = oneGBSAnalysisMeta(connection, 'eumaeus',
#                             'CCAE', 'HistoricalComparator',
#                             exposure_id = 211981,
#                             analysis_id = 2,
#                             period_id = 12,
#                             IPCtable = IPCs,
#                             priorSd = 4,
#                             nullPriorSds = c(0.5,0.5),
#                             outcomeLPPath = outcomepath)


# 2. summary functions to post process analysis results----
## (a) get the summary AND calculate posterior P1 and P0 for ONE outcome only
##     (with option to get 95% posterior credible interval)
##    returns a 1-row dataframe
summarizeOneOutcome <- function(res, getCI = TRUE) {
  
  # get CIs and posterior probabilities of H1 and H0
  for (s in c('postSamps', 'adjustedPostSamps')) {
    prefix = ifelse(s == 'postSamps', '', 'adjusted')
    # 95% credible intervals
    if (getCI) {
      res[[paste0(prefix, 'CI95_lb')]] = 
        quantile(res[[s]], 0.025, na.rm = TRUE) %>% 
        as.numeric()
      res[[paste0(prefix, 'CI95_ub')]] = 
        quantile(res[[s]], 0.975, na.rm = TRUE) %>% 
        as.numeric()
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

## (b) get the summary for a list of outcomes
##       return multiple-row dataframe of summary
summarizeOneAnalysis <- function(resls, getCI=TRUE){
  
  # if no resuls are returned, return NULL
  if(length(resls) == 0){
    NULL
  }else{
    outcomes = names(resls)
    foreach(o = outcomes, .combine = 'bind_rows') %do% {
      res.o = summarizeOneOutcome(resls[[o]], getCI = getCI)
      res.o$outcome_id = as.numeric(o)
      res.o
    }
  }
  
}


## (c) get posterior samples
getSamplesFromOneAnalysis <- function(resls){
  
  # if empty result, return NULL
  if(length(resls) == 0){
    NULL
  }else{
    outcomes = names(resls)
    postSamps = 
      foreach(o=outcomes, .combine = 'rbind') %dopar% {
        res.o = matrix(resls[[o]]$postSamps, nrow=1)
        row.names(res.o) = o
        res.o
      }
    adjPostSamps = 
      foreach(o=outcomes, .combine = 'rbind') %dopar% {
        res.o = matrix(resls[[o]]$adjustedPostSamps, nrow=1)
        row.names(res.o) = o
        res.o
      }
    
    list(postSamps = postSamps, adjustedPostSamps = adjPostSamps)
  }
  
}


# 3.larger scale function --------
# for a (database, method, exposure) combo
# can loop through analysis_ids and/or period_ids
# OR loop through prior choices
# return a dataframe of result summary 
# AND a large list of posterior samples

multiGBSAnalysesMeta <- function(connection,
                                 schema,
                                 database_id,
                                 method,
                                 exposure_id,
                                 analysis_ids,
                                 period_ids,
                                 IPCtable = NULL,
                                 priors = list(Mean = rep(0, 3),
                                               Sd = c(10, 1.5, 4)),
                                 returnCIs = TRUE,
                                 numsamps = 10000,
                                 thin = 10,
                                 nullPriorSds = c(2,0.5),
                                 robustMetaAnalysis = FALSE,
                                 minNCs = 5,
                                 preLearnNull = TRUE,
                                 negControls = NULL,
                                 includeSyntheticPos = FALSE,
                                 LPPath = '~/Documents/Research/better_gbs/',
                                 savepath = 'localCache/testResults/',
                                 sampspath = 'localCache/sampleSaves/',
                                 removeTempSummary = TRUE) {
  
  # create folders if not existent
  if(!dir.exists(savepath)){
    dir.create(savepath)
  }
  if(!dir.exists(sampspath)){
    dir.create(sampspath)
  }
  
  # generate prior table
  priorTable = getPriorTable(priors = priors, default = TRUE)
  prior_ids = priorTable$prior_id
  
  # save prior table if it doesn't already exist in the results folder
  if (!file.exists(file.path(savepath, 'priorTable.rds'))) {
    saveRDS(priorTable, file.path(savepath, 'priorTable.rds'))
  }
  
  # get a list of negative controls to use across all runs
  if (is.null(negControls)) {
    sql <- "SELECT outcome_id from eumaeus.NEGATIVE_CONTROL_OUTCOME"
    NCs = DatabaseConnector::querySql(connection, sql)$OUTCOME_ID
  } else{
    NCs = negControls
  }
  
  # whether or not to include synthetic positive controls
  if (includeSyntheticPos) {
    oids = NULL
  } else{
    oids = NCs
  }
  
  # check if to include imputed positive controls
  # and filter it down to the relevant exposure ids
  if (!is.null(IPCtable)) {
    IPCtable = IPCtable %>% filter(EXPOSURE_ID == exposure_id)
  }
  
  
  # pull likelihood profiles for the outcome of interest
  if(database_id == 'OptumEhr'){
    dbname = 'OptumEHR'
  }else if(str_starts(database_id, pattern = 'IBM')){
    dbname = str_split(database_id, pattern = '_')[[1]][2]
  }else{
    dbname = database_id
  }
  LPfilePath = file.path(LPPath, sprintf('Results_%s', dbname), 'likelihood_profile.csv')
  outcomeLPs = read_csv(LPfilePath)
  
  
  # go through period_ids and combine the summary data tables
  summary_dat_list = list()
  for (p in period_ids) {
    ## check if period summary already exists in the savepath
    if(!is.null(savepath) && dir.exists(savepath)){
      flag = tryCatch(
        expr = {
          checkPeriodSummary(savepath, database_id,
                             method, exposure_id,
                             p, analysis_ids)
        },
        error = function(e){
          cat('Got error when checking period summary file. Treat like results are not available..\n\n')
          FALSE
        }
      )
      
      # if period summary already exists, move on...
      if(flag) {
        cat(sprintf('Analysis already done for database %s, method %s, exposure %s, period %s and analysis %s-%s! Skipped...\n\n',
                    database_id, method, exposure_id, 
                    p, min(analysis_ids), max(analysis_ids)))
        next
      }
    }
    
    ## get likelihood profiles to use
    ## (if no subsetting on outcomes, will include old synthetic positive controls)
    LPs = getMultiLikelihoodProfiles(
      connection,
      schema,
      database_id,
      exposure_id,
      analysis_id = NULL,
      period_id = p,
      outcome_ids = oids,
      method = method
    )
    # if no LPs returned, skip this one
    if (nrow(LPs) == 0) {
      cat('No pre-saved negative control likelihood profiles available! Skipped.\n')
      next
    }
    
    # go through analysis_ids and combine the summary data tables
    analysis_dat_list = list()
    
    #analysis_dat = 
    #  foreach(a = analysis_ids, .combine = 'bind_rows') %dopar% {
    for (a in analysis_ids) {
      
      ## FIRST check if this analysis has already been done
      ## do so by checking if posterior samples are already saved
      if(!is.null(sampspath)){
        fname = sprintf(
          '%s_%s_%s_period%s_analysis%s_samples.rds',
          database_id,
          method,
          exposure_id,
          p,
          a
        )
        if(file.exists(file.path(sampspath,fname))) next
      }
      
      ## get relevant rows for this analysis and period
      # of the special outcome
      # if unavailable, skip!
      outcomeLPs.a = outcomeLPs %>% 
        filter(analysis_id == a, 
               period_id == p)
      if(nrow(outcomeLPs.a) == 0){
        cat(sprintf('No likelihood profiles for the special outcome(s) available for analysis %s! Skipped.\n', a))
        next
      }
      
      ## check NC LPs
      LPs.a = LPs %>% filter(ANALYSIS_ID == a)
      if(nrow(LPs.a) == 0){
        cat(sprintf('No negative control likelihood profiles available for analysis %s! Skipped.\n', a))
        next
      } 
      
      ## get NC LPs 
      NC_LPs.a = LPs.a %>% filter(OUTCOME_ID %in% NCs)
      
      
      ## pre-learn null distribution if...
      if (preLearnNull) {
        learnedNull = fitNegativeControlDistributionLikelihood(
          connection,
          schema,
          database_id,
          method,
          exposure_id,
          analysis_id = a,
          period_id = p,
          savedLPs = NC_LPs.a,
          outcomeToExclude = NULL,
          priorSds = nullPriorSds,
          minNCs = minNCs,
          numsamps = numsamps,
          thin = thin,
          robust = robustMetaAnalysis
        )
      } else{
        ## otherwise, just set it to NULL
        learnedNull = NULL
      }
      
      
      ## change the organization of the list for a little bit
      ## easier for post-processing after the for loop
      ## NOT TO USE for the foreach parallel run !!! #
      big_list = list(samples = list(), summary = list())
      for (pr in prior_ids) {
        ## output message
        cat(
          sprintf(
            '\n\n\n Analysis for period %s, analysis %s and prior %s......\n',
            p,
            a,
            pr
          )
        )
        
        ## get prior mean and sd
        this.prior = priorTable %>% filter(prior_id == pr) %>% select(Mean, Sd)
        pMean = this.prior$Mean
        pSd = this.prior$Sd
        
        ## run analysis to get results
        this.res = oneGBSAnalysisMeta(connection, 
                                      schema,
                                      database_id,
                                      method,
                                      exposure_id,
                                      analysis_id = a,
                                      period_id = p,
                                      IPCtable = IPCtable,
                                      priorMean = pMean,
                                      priorSd = pSd,
                                      nullPriorSds = nullPriorSds,
                                      minNCs = minNCs, 
                                      numsamps = numsamps,
                                      thin = thin,
                                      outcomeLPs = outcomeLPs, 
                                      outcomeLPPath = LPPath,
                                      NCLPs = NC_LPs.a,
                                      preSaveNull = learnedNull,
                                      negControls = NCs,
                                      robustMetaAnalysis = robustMetaAnalysis)
        
        ## if no results returned, skip
        if (length(this.res) == 0) {
          cat(sprintf('No results returned for prior %s!\n', pr))
          next
        }
        
        ## summarize it
        this.summ = summarizeOneAnalysis(this.res, getCI = returnCIs)
        # add column for prior_id
        this.summ$prior_id = pr
        
        ## pull post. samples
        this.samps = getSamplesFromOneAnalysis(this.res)
        
        ### change organization of the big_list for easier post-processing
        ### FOR THE FOR LOOP ONLY ###
        ### NOT TO USE for the foreach parallel run !!! ###
        big_list$summary[[which(prior_ids == pr)]] = this.summ
        big_list$samples[[which(prior_ids == pr)]] = this.samps
        
        # cat(sprintf('\nObject big_list$samples now has space %s\n\n',
        #             object.size(big_list$samples)))
      }
      
      # if NO results are returned for this analysis_id
      # then move on to the next one
      # if (length(big_list) == 0) {
      #   cat(sprintf('No results available for analysis %s!\n', a))
      #   return()
      # }
      
      # ONLY for the for loop! #
      # check if big_list$summary is empty
      if (('summary' %in% names(big_list)) &
          (length(big_list$summary) == 0)) {
        cat(sprintf('No results available for analysis %s!\n', a))
        #return()
        next
      }
      
      # save the posterior samples in the big_list
      # IF a "sampspath" is provided
      # otherwise just skip posterior sample saving
      if (!is.null(sampspath)) {
        # create folder
        if (!dir.exists(file.path(sampspath)))
          dir.create(file.path(sampspath))
        
        # glue together file name for saving
        fname = sprintf(
          '%s_%s_%s_period%s_analysis%s_samples.rds',
          database_id,
          method,
          exposure_id,
          p,
          a
        )
        
        # give the samples a name for saving
        samplesToSave = big_list$samples
        
        # save as rds
        saveRDS(samplesToSave, file = file.path(sampspath, fname))
      }
      
      # save the little summary data table as well
      summaryToSave = bind_rows(big_list$summary) %>%
        mutate(analysis_id = a)
      
      if (!is.null(savepath)) {
        # create folder
        if (!dir.exists(file.path(savepath)))
          dir.create(file.path(savepath))
        
        # glue together file name for saving
        fname = sprintf(
          '%s_%s_%s_period%s_analysis%s_summary.rds',
          database_id,
          method,
          exposure_id,
          p,
          a
        )
        
        # give the summary a name
        # summaryToSave = analysis_dat_list[[as.character(a)]]
        # save as rds
        saveRDS(summaryToSave,
                file = file.path(savepath, fname))
      }
      
      # return the summary dataframe part
      # analysis_dat_list[[as.character(a)]] =
      #   bind_rows(big_list$summary) %>%
      #   mutate(analysis_id = a)
      
      #return(summaryToSave)
      
      analysis_dat_list[[as.character(a)]] = summaryToSave
      
    }
    
    # check if any results are returned
    if (length(analysis_dat_list) == 0) {
      #if (nrow(analysis_dat) == 0) {
      cat(sprintf(
        '\nNo results available for all analyses in period %s!\n\n',
        p
      ))
      next
    }
    
    # save result for this whole period
    period_summary = bind_rows(analysis_dat_list) %>%
      mutate(period_id = p)
    
    if(!is.null(savepath) && nrow(period_summary) > 0){
      # create folder if...
      if (!dir.exists(file.path(savepath)))
        dir.create(file.path(savepath))
      
      # glue together file name for saving
      ## also include analysis_id range
      ## for splitting runs
      
      ## Feb 9 change: include *actual* analysis_id range in the summary table
      
      analysis_range = range(period_summary$analysis_id)
      
      fname = sprintf(
        'period_summary_%s_%s_%s_period%s_analysis%s-%s.rds',
        database_id,
        method,
        exposure_id,
        p,
        analysis_range[1],#min(analysis_ids),
        analysis_range[2]#max(analysis_ids)
      )
      
      # save as rds
      saveRDS(period_summary,
              file = file.path(savepath, fname))
      
      # remove the temporary summary files of each analysis
      if(removeTempSummary){
        fnamePattern = sprintf(
          '%s_%s_%s_period%s_analysis[1-9]*_summary.rds',
          database_id,
          method,
          exposure_id,
          p
        )
        filesToRM = list.files(path = savepath, 
                               pattern = fnamePattern)
        if(length(filesToRM) > 0) {file.remove(file.path(savepath,filesToRM))}
      }
    }
    
    # include it in the multi-period big list
    summary_dat_list[[as.character(p)]] = period_summary
  }
  
  # return result
  # update: return summary datatable only
  
  # update: only return a big summary datatable if there are >1 periods...
  if(length(period_ids) > 1){
    final_summary = bind_rows(summary_dat_list)
    
    # also save it if savepath is provided
    if(!is.null(savepath) && nrow(final_summary) > 0){
      # create folder if...
      if (!dir.exists(file.path(savepath)))
        dir.create(file.path(savepath))
      
      # glue together file name for saving
      ## also include analysis_id range
      ## for splitting runs
      
      ## Feb 9 change: 
      ## include *actual* analysis_id range in the summary table
      ## also do that for period_id ranage
      analysis_range = range(final_summary$analysis_id)
      period_range = range(final_summary$period_id)
      
      fname = sprintf(
        'Summary_%s_%s_%s_period%s-%s_analysis%s-%s.rds',
        database_id,
        method,
        exposure_id,
        period_range[1],
        period_range[2],
        analysis_range[1],#min(analysis_ids),
        analysis_range[2]#max(analysis_ids)
      )
      
      # save as rds
      saveRDS(final_summary,
              file = file.path(savepath, fname))
    }
  }
  
  cat(sprintf('All analysis finished for database %s, exposure %s, using method %s analysis %s-%s for periods %s-%s...\n\n',
              database_id,
              exposure_id,
              method,
              min(analysis_ids),
              max(analysis_ids),
              min(period_ids),
              max(period_ids)))
  return()
  
  
  
  #summary_dat = bind_rows(summary_dat_list)
  # (a summary datatable and a prior table)
  # list(summary_dat = summary_dat,
  #      priors = priorTable)
}



# ## try it
# IPCs = getIPCs(connection, 'eumaeus', 'localCache/')
# multiRes = multiGBSAnalysesMeta(connection,
#                                 'eumaeus',
#                                 'OptumDod',
#                                 'HistoricalComparator',
#                                 exposure_id = 211981,
#                                 analysis_ids = c(6,8),
#                                 period_ids = 12,
#                                 IPCtable = IPCs,
#                                 priors = list(Mean = 0,
#                                               Sd = 4),
#                                 nullPriorSds = c(.5,.5),
#                                 preLearnNull = TRUE,
#                                 LPPath = '~/Documents/Research/better_gbs/',
#                                 robustMetaAnalysis = FALSE)



