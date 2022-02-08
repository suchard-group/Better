# Jan 2022
# function to run Bayesian analysis based on likelihood profiles
# with or without calibration (using negative control analysis)

library(EvidenceSynthesis)
source('./extras/getLikelihoodProfile.R')
source('./extras/fitNegativeControlDistribution.R')
source('./extras/helperFunctions.R')

library(foreach)
library(doParallel)
registerDoParallel()


# small module function for
# one (database, method, exposure, analysis, period) combo (with all available outcomes)
# can use pre-saved raw table of likelihood profiles

## return numeric(0) if no likelihood profiles are available 
## return empty list() if no negative control results are available 

## Feb 2022 update: 
## add imputed positive control outcomes
## AND try to parallelize over the outcomes
oneBayesianAnalysis <- function(connection,
                                schema,
                                database_id,
                                method,
                                exposure_id,
                                analysis_id, 
                                period_id,
                                IPCtable = NULL,
                                noSyntheticPos=TRUE,
                                savedEstimates = NULL,
                                priorMean = 0,
                                priorSd = 1, 
                                numsamps = 10000,
                                thin = 10,
                                preSaveLPs = NULL,
                                preSaveNull = NULL,
                                negControls = NULL){
  
  # access list of negative control outcome ids
  if(is.null(negControls)){
    sql <- "SELECT outcome_id from eumaeus.NEGATIVE_CONTROL_OUTCOME"
    NCs = DatabaseConnector::querySql(connection, sql)$OUTCOME_ID
  }else{
    NCs = negControls
  }
  
  # if preSaveLPs is provided as a dataframe, then use pre-pulled LPs, otherwise query
  if(is.null(preSaveLPs)){
    LPs = getMultiLikelihoodProfiles(connection, schema,
                                     database_id,
                                     exposure_id, 
                                     analysis_id,
                                     period_id,
                                     method = method)
  }else{
    LPs = preSaveLPs
  }
  
  # if no relevant likelihood profiles are available 
  # do nothing and return numeric(0)
  if(length(LPs) == 0){
    return(numeric(0))
  }
  
  # get a subset IPC table if available
  if(!is.null(IPCtable)){
    IPCtable = IPCtable %>% filter(EXPOSURE_ID == exposure_id)
  }
  
  # get the list of all outcomes that have LPs available
  outcomesToDo = unique(LPs$OUTCOME_ID)
  
  # filter outcomes to negative control outcomes only (NO synthetic pos controls) if...
  if(noSyntheticPos){
    outcomesToDo = outcomesToDo[outcomesToDo %in% NCs]
  }
  
  # if there is pre-saved learned null distribution 
  # AND there are positive control outcomes
  # use the pre-learned one
  # AND need to learn one if there are positive control outcomes
  if(any(!outcomesToDo %in% NCs)){
    if(is.null(preSaveNull)){
      null = fitNegativeControlDistribution(connection, schema, database_id, 
                                            method, exposure_id, 
                                            analysis_id, period_id, 
                                            savedEstimates = savedEstimates,
                                            outcomeToExclude = NULL,
                                            numsamps = numsamps, thin = thin)
    }else{
      null = preSaveNull
    }
    # if no negative controls are available at all, return numeric(0)
    if(length(null)==0){
      return(numeric(0))
    }
  }
  
  # go through all the outcomes
  
  res = list()
  for(outcome in outcomesToDo){
  
  ## try using foreach instead here
  ## failed.. will try later on
  # res =
  #   foreach(outcome = outcomesToDo, .combine='c') %dopar% {
    # message
    cat('\n\nAnalysis for outcome', outcome, 'underway...\n')
    
    # get the likelihood profile
    lik = selectLikelihoodProfileEntry(LPs,
                                       database_id,
                                       method,
                                       exposure_id, 
                                       period_id,
                                       outcome,
                                       analysis_id)
    # if there isn't an entry
    if(length(lik)==0) next
    
    # check if it's a negative control outcome
    if(outcome %in% NCs){
      # negative control: LOO null distribution needed
      this.null = fitNegativeControlDistribution(connection, schema, database_id, 
                                                 method, exposure_id, analysis_id, 
                                                 period_id, 
                                                 savedEstimates = savedEstimates,
                                                 outcomeToExclude = outcome,
                                                 numsamps = numsamps, thin = thin)
      
      # if can't fit null (no other negative control results available)
      # move on to next one
      if(length(this.null) == 0) next()
      
      biases = this.null$bias
    }else{
      # positive control: use null directly
      biases = null$bias
    }
    
    # run MCMC
    ## need to do imputed positive controls as well if this one is NC
    if(!is.null(IPCtable) && outcome %in% NCs){
      nudgeValues = c(1, 1.5, 2, 4)
      outcomeVector = c(outcome, getThisIPCs(exposure_id, outcome, IPCtable))
    }else{
      nudgeValues = c(1)
      outcomeVector = outcome
    }
    names(outcomeVector) = as.character(nudgeValues)
    
    #tempRes = list()
    original_lik = lik
    for(nv in nudgeValues){
      this.outcome = outcomeVector[[as.character(nv)]]
      if(nv != 1){
        cat('\nAlso analyze imputed positive control outcome', 
            this.outcome, '...\n')
      }
      ## nudge the likelihood profiles given the effect size
      lik = nudgeLikelihood(original_lik, nv)
      mcmc = approximateSimplePosterior(lik, 
                                        chainLength = numsamps * thin + 1e5,
                                        burnIn = 1e5, # (use default burn-in = 1e5)
                                        subSampleFrequency = thin,
                                        priorMean = priorMean, priorSd = priorSd)
      
      samps = mcmc$theta1
      # with calibration
      adjSamps =  samps - biases
      
      # get result to list
      # key: outcome id
      res[[as.character(this.outcome)]] = 
      #tempRes[[as.character(this.outcome)]] = 
        list(postSamps = samps, adjustedPostSamps = adjSamps,
             postMean = mean(samps), postMAP = getMAP(samps), 
             postMedian = median(samps),
             adjustedPostMean = mean(adjSamps), 
             adjustedPostMAP = getMAP(adjSamps),
             adjustedPostMedian = median(adjSamps))
      #cat('\nResult written to list!\n\n')
    }
   #tempRes 
  }

  # return final result
  cat(sprintf('Finished Bayesian analysis for database %s, exposure %s, in period %s, using %s analysis %s\n',
          database_id, exposure_id, period_id, method, analysis_id))
  res
}


# try it
# IPCs = getIPCs(connection, 'eumaeus', './localCache/')
# ## example of returning a lot of results
# bayesRes = oneBayesianAnalysis(connection,
#                                'eumaeus',
#                                database_id = 'IBM_MDCD',
#                                method = 'SCCS',
#                                exposure_id = 21184,
#                                analysis_id = 1,
#                                period_id = 9,
#                                IPCtable = IPCs)
# ## another example of returning less or nothing (list())
# bayesRes = oneBayesianAnalysis(connection,
#                                'eumaeus',
#                                database_id = 'IBM_MDCD',
#                                method = 'SCCS',
#                                exposure_id = 21184,
#                                analysis_id = 15,
#                                period_id = 1,
#                                IPCtable = IPCs)

# ## try to use pre-saved estimates for fitting null distribution
# estimates = getNegControlEstimates(connection, 'eumaeus', 'IBM_MDCD', 'SCCS',
#                                    21184, period_id = 9)
# ## estimates %>% filter(ANALYSIS_ID==1) %>% select(OUTCOME_ID) %>% pull() %>% unique()
# bayesRes = oneBayesianAnalysis(connection,
#                                'eumaeus',
#                                database_id = 'IBM_MDCD',
#                                method = 'SCCS',
#                                exposure_id = 21184,
#                                analysis_id = 1,
#                                period_id = 9,
#                                IPCtable = IPCs,
#                                savedEstimates = estimates,
#                                negControls = c(196044, 196347)) # only use two NCs and see if that works
# bayesRes = oneBayesianAnalysis(connection,
#                                'eumaeus',
#                                database_id = 'IBM_MDCD',
#                                method = 'SCCS',
#                                exposure_id = 21184,
#                                analysis_id = 15,
#                                period_id = 3,
#                                IPCtable = IPCs,
#                                savedEstimates = estimates)


## helper functions to extract results from the result list
## (1.a) get the summary AND calculate posterior P1 and P0 for ONE outcome only
##     (with option to get 95% posterior credible interval)
##    returns a 1-row dataframe
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

# # try it
# summarizeOneOutcome(bayesRes$`10109`, getCI=TRUE)

## (1.b) get the summary for a list of outcomes
##       return multiple-row dataframe of summary
summarizeOneAnalysis <- function(resls, getCI=TRUE){
  
  # if no resuls are returned, return NULL
  if(length(resls) == 0){
    NULL
  }else{
    outcomes = names(resls)
    foreach(o = outcomes, .combine = 'bind_rows') %dopar% {
      res.o = summarizeOneOutcome(resls[[o]], getCI = getCI)
      res.o$outcome_id = as.numeric(o)
      res.o
    }
  }
  
}

# # try it
# BRdat = summarizeOneAnalysis(bayesRes)


## (2.a) get the post. samples only
## return a list of two matrices 
## ("postSamps" & "adjustedPostSamps")
## (rownames are outcome ids)
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

# # try it
# allSamps = getSamplesFromOneAnalysis(bayesRes)




#####
# larger scale function
# for a (database, method, exposure) combo
# can loop through analysis_ids and/or period_ids
# OR loop through prior choices
# return a dataframe of result summary 
# AND a large list of posterior samples

multiBayesianAnalyses <- function(connection,
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
                                  preLearnNull = TRUE,
                                  preSaveEstimates = TRUE,
                                  negControls = NULL,
                                  includeSyntheticPos = FALSE,
                                  savepath = 'localCache/testResults/',
                                  sampspath = 'localCache/sampleSaves/',
                                  removeTempSummary = TRUE) {
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
  
  # go through period_ids and combine the summary data tables
  summary_dat_list = list()
  for (p in period_ids) {
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
      #cat('No pre-saved likelihood profiles available! Skipped.\n')
      next
    }
    
    # pre-pull negative control Estimates
    if(preSaveEstimates){
      savedEstimates = getNegControlEstimates(connection, schema, 
                                              database_id, method, 
                                              exposure_id,
                                              period_id = p)
      
      ## check if estimates exist at all
      if(nrow(savedEstimates) == 0){
        cat('No negative controls estimates available! Skipped.\n')
        next
      }
    }else{
      savedEstimates = NULL
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
          '%s_%s_%s_period%s_analysis%s_samples.RData',
          database_id,
          method,
          exposure_id,
          p,
          a
        )
        if(file.exists(file.path(sampspath,fname))) next
      }
      
      ## check LPs
      LPs.a = LPs %>% filter(ANALYSIS_ID == a)
      if(nrow(LPs.a) == 0){
        return()
      }
      
      ## pre-learn null distribution if...
      if (preLearnNull) {
        learnedNull = fitNegativeControlDistribution(
          connection,
          schema,
          database_id,
          method,
          exposure_id,
          analysis_id = a,
          period_id = p,
          savedEstimates = savedEstimates,
          outcomeToExclude = NULL,
          numsamps = numsamps,
          thin = thin
        )
      } else{
        ## otherwise, just set it to NULL
        learnedNull = NULL
      }
      
      # go through each prior setting
      # (default: combine as a list)
      # big_list = list()
      
      # NOT doing the parallel thing for now...
      # big_list =
      #   foreach(pr = prior_ids) %dopar% {
      #
      #   }
      
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
        this.res = oneBayesianAnalysis(
          connection,
          schema,
          database_id,
          method,
          exposure_id,
          analysis_id = a,
          period_id = p,
          IPCtable = IPCtable,
          savedEstimates = savedEstimates,
          priorMean = pMean,
          priorSd = pSd,
          numsamps = numsamps,
          thin = thin,
          preSaveLPs = LPs,
          preSaveNull = learnedNull,
          negControls = NCs
        )
        
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
        # this.samps = list(this.samps)
        # attr(this.samps, 'names') = as.character(pr)
        
        ## return a big list with summary and samples
        # big_list[[which(prior_ids == pr)]] =
        #   list(summary = this.summ, samples = this.samps)
        
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
        return()
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
          '%s_%s_%s_period%s_analysis%s_samples.RData',
          database_id,
          method,
          exposure_id,
          p,
          a
        )
        
        # give the samples a name for saving
        samplesToSave = big_list$samples
        
        # save as RData
        save(samplesToSave, file = file.path(sampspath, fname))
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
          '%s_%s_%s_period%s_analysis%s_summary.RData',
          database_id,
          method,
          exposure_id,
          p,
          a
        )
        
        # give the summary a name
        # summaryToSave = analysis_dat_list[[as.character(a)]]
        # save as RData
        save(summaryToSave,
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
    #if (length(analysis_dat_list) == 0) {
    if (nrow(analysis_dat) == 0) {
      cat(sprintf(
        '\nNo results available for all analyses in period %s!\n\n',
        p
      ))
      next
    }
    
    summary_dat_list[[as.character(p)]] =
      bind_rows(analysis_dat_list) %>%
      #analysis_dat %>%
      mutate(period_id = p)
  }
  
  # return result
  # update: return summary datatable only
  
  final_summary = bind_rows(summary_dat_list)
  
  # also save it if savepath is provided
  if(!is.null(savepath)){
    # create folder if...
    if (!dir.exists(file.path(savepath)))
      dir.create(file.path(savepath))
    
    # glue together file name for saving
    fname = sprintf(
      'period_summary_%s_%s_%s_period%s.RData',
      database_id,
      method,
      exposure_id,
      p
    )
    
    # save as RData
    save(final_summary,
         file = file.path(savepath, fname))
    
    # remove the temporary summary files of each analysis
    if(removeTempSummary){
      fnamePattern = sprintf(
        '%s_%s_%s_period%s_analysis[1-9]*_summary.RData',
        database_id,
        method,
        exposure_id,
        p
      )
      filesToRM = list.files(path = savepath, 
                             pattern = fnamePattern)
      if(length(filesToRM) > 0) {file.remove(filesToRM)}
    }
  }
  
  #summary_dat = bind_rows(summary_dat_list)
  # (a summary datatable and a prior table)
  # list(summary_dat = summary_dat,
  #      priors = priorTable)
}

# # try it
IPCs = getIPCs(connection, 'eumaeus', './localCache/')
# multiRes = multiBayesianAnalyses(connection,
#                                  'eumaeus',
#                                  database_id = 'IBM_MDCD',
#                                  method = 'SCCS',
#                                  exposure_id = 21184,
#                                  analysis_ids = c(15),
#                                  period_ids = c(5),
#                                  includeSyntheticPos = FALSE)
# 
# multiRes2 = multiBayesianAnalyses(connection,
#                                  'eumaeus',
#                                  database_id = 'IBM_MDCD',
#                                  method = 'SCCS',
#                                  exposure_id = 21184,
#                                  analysis_ids = c(15),
#                                  period_ids = c(5),
#                                  includeSyntheticPos = TRUE)

multiRes3 = multiBayesianAnalyses(connection,
                                 'eumaeus',
                                 database_id = 'IBM_MDCD',
                                 method = 'SCCS',
                                 exposure_id = 21184,
                                 analysis_ids = c(15),
                                 period_ids = c(5),
                                 includeSyntheticPos = FALSE,
                                 IPCtable = IPCs,
                                 preLearnNull = FALSE,
                                 savepath = './localCache/testResults')
# try a very small parallel run example
multiRes4 = multiBayesianAnalyses(connection,
                                  'eumaeus',
                                  database_id = 'IBM_MDCD',
                                  method = 'SCCS',
                                  exposure_id = 21184,
                                  analysis_ids = c(12,14),
                                  period_ids = c(5),
                                  includeSyntheticPos = FALSE,
                                  IPCtable = IPCs,
                                  priors = list(Mean = 0, Sd = 1.5),
                                  preLearnNull = FALSE,
                                  negControls = c(443421, 196347),
                                  savepath = './localCache/testResults')
