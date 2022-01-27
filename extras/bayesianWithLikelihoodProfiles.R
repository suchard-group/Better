# Jan 2022
# function to run Bayesian analysis based on likelihood profiles
# with or without calibration (using negative control analysis)

library(EvidenceSynthesis)
source('./extras/getLikelihoodProfile.R')
source('./extras/fitNegativeControlDistribution.R')

library(foreach)
library(doParallel)
registerDoParallel()


## a little helper function to get maximum density estimate from samples
getMAP <- function(x){
  dens = density(x)
  dens$x[which.max(dens$y)]
}

## another little helper function to generate a table of prior choices
getPriorTable <- function(priors = list(Mean = c(0,0,0),
                                        Sd = c(10, 1.5, 4)),
                          default=TRUE){
  # default: whether or not to use default choices in priors
  #          (FALSE: use historical rates for historical comparator - TBD)
  if(default){
    res = as.data.frame(priors)
    res$prior_id = seq_along(res$Mean)
  }else{
    # the other option not yet implemented
    stop('Non-default prior not yet implemented!\n')
  }
  res
}

# small module function for
# one (database, method, exposure, analysis, period) combo (with all available outcomes)
# can use pre-saved raw table of likelihood profiles

## return numeric(0) if no likelihood profiles are available 
## return empty list() if no negative control results are available 
oneBayesianAnalysis <- function(connection,
                                schema,
                                database_id,
                                method,
                                exposure_id,
                                analysis_id, 
                                period_id,
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
  
  # get the list of all outcomes that have LPs available
  outcomesToDo = unique(LPs$OUTCOME_ID)
  
  # if there is pre-saved learned null distribution 
  # AND there are positive control outcomes
  # use the pre-learned one
  # AND need to learn one if there are positive control outcomes
  if(any(!outcomesToDo %in% NCs)){
    if(is.null(preSaveNull)){
      null = fitNegativeControlDistribution(connection, schema, database_id, 
                                            method, exposure_id, analysis_id, 
                                            period_id, outcomeToExclude = NULL,
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
    # message
    cat('Analysis for outcome', outcome, 'underway...\n')
    
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
                                                 period_id, outcomeToExclude = outcome,
                                                 numsamps = numsamps, thin = thin)
      
      # if can't fit null (no other negative control results available)
      # move on to next one
      if(length(this.null) == 0) next
      
      biases = this.null$bias
    }else{
      # positive control: use null directly
      biases = null$bias
    }
    
    # run MCMC
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
    res[[as.character(outcome)]] = 
      list(postSamps = samps, adjustedPostSamps = adjSamps,
           postMean = mean(samps), postMAP = getMAP(samps), 
           postMedian = median(samps),
           adjustedPostMean = mean(adjSamps), 
           adjustedPostMAP = getMAP(adjSamps),
           adjustedPostMedian = median(adjSamps))
  }

  # return final result
  cat(sprintf('Finished Bayesian analysis for database %s, exposure %s, in period %s, using %s analysis %s\n',
          database_id, exposure_id, period_id, method, analysis_id))
  res
}


# try it
# ## example of returning a lot of results
# bayesRes = oneBayesianAnalysis(connection,
#                                'eumaeus',
#                                database_id = 'IBM_MDCD',
#                                method = 'SCCS',
#                                exposure_id = 21184,
#                                analysis_id = 1,
#                                period_id = 9)
# ## another example of returning less or nothing (list())
# bayesRes = oneBayesianAnalysis(connection,
#                                'eumaeus',
#                                database_id = 'IBM_MDCD',
#                                method = 'SCCS',
#                                exposure_id = 21184,
#                                analysis_id = 15,
#                                period_id = 1)


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
                                  priors = list(Mean = rep(0, 3),
                                                Sd = c(10, 1.5, 4)), 
                                  returnCIs = TRUE,
                                  numsamps = 10000,
                                  thin = 10,
                                  preLearnNull = FALSE,
                                  negControls = NULL,
                                  savepath = 'localCache/testResults/'){
  
  # generate prior table
  priorTable = getPriorTable(priors = priors, default=TRUE)
  prior_ids = priorTable$prior_id
  
  # get a list of negative controls to use across all runs
  if(is.null(negControls)){
    sql <- "SELECT outcome_id from eumaeus.NEGATIVE_CONTROL_OUTCOME"
    NCs = DatabaseConnector::querySql(connection, sql)$OUTCOME_ID
  }else{
    NCs = negControls
  }
  
  # go through period_ids and combine the summary data tables
  summary_dat_list = list()
  for(p in period_ids){
    #foreach(p = period_ids, .combine = 'bind_rows') %dopar% {
      
      # go through analysis_ids and combine the summary data tables
      #foreach(a = analysis_ids, .combine = 'bind_rows') %dopar% {
      analysis_dat_list = list()
      for(a in analysis_ids){
        
        # go through each prior setting
        # (default: combine as a list)
        big_list = list()
        for(pr in prior_ids){
          ## output message
          cat(sprintf(
            '\n\n\n Analysis for period %s, analysis %s and prior %s......\n',
            p,
            a,
            pr
          ))
          
          ## get prior mean and sd
          this.prior = priorTable %>% filter(prior_id == pr) %>% select(Mean, Sd)
          pMean = this.prior$Mean
          pSd = this.prior$Sd
          
          ## get likelihood profiles to use
          LPs = getMultiLikelihoodProfiles(
            connection,
            schema,
            database_id,
            exposure_id,
            analysis_id = a,
            period_id = p,
            method = method
          )
          
          # if no LPs returned, skip this one
          if(nrow(LPs) == 0){
            cat('No pre-saved likelihood profiles available! Skipped.\n')
            next
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
              outcomeToExclude = NULL,
              numsamps = numsamps,
              thin = thin
            )
          } else{
            ## otherwise, just set it to NULL
            learnedNull = NULL
          }
          
          ## run analysis to get results
          this.res = oneBayesianAnalysis(
            connection,
            schema,
            database_id,
            method,
            exposure_id,
            analysis_id = a,
            period_id = p,
            priorMean = pMean,
            priorSd = pSd,
            numsamps = numsamps,
            thin = thin,
            preSaveLPs = LPs,
            preSaveNull = learnedNull,
            negControls = NCs
          )
          
          ## if no results returned, skip
          if(length(this.res) == 0){
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
          big_list[[which(prior_ids == pr)]] = 
            list(summary = this.summ, samples = this.samps)
        }
        
        # NOT doing the parallel thing for now...
        # big_list =
        #   foreach(pr = prior_ids) %dopar% {
        #     
        #   }
        
        # if NO results are returned for this analysis_id
        # then move on to the next one
        if(length(big_list) == 0){
          cat(sprintf('No results available for analysis %s!\n', a))
          next
        }
        
        # save the posterior samples in the big_list
        # IF a savepath is provided
        # otherwise just skip posterior sample saving 
        if(!is.null(savepath)){
          # create folder 
          if(!dir.exists(file.path(savepath))) dir.create(file.path(savepath))
          
          # glue together file name for saving
          fname = sprintf('%s_%s_%s_period%s_analysis%s_samples.RData', 
                          database_id, method, exposure_id, p, a)
          # save as RData
          save(big_list$samples, file = file.path(savepath,fname))
        }
        
        # return the summary dataframe part
        analysis_dat_list[[as.character(a)]] = 
          big_list$summary %>% mutate(analysis_id = a)
      }
      
      # check if any results are returned
      if(length(analysis_dat_list) == 0){
        cat(sprintf('\nNo results available for all analyses in period %s!\n\n', p))
        next
      }
      
      summary_dat_list[[as.character(p)]] = 
        bind_rows(analysis_dat_list) %>% 
        mutate(period_id = p)
    
  }
  
  summary_dat = bind_rows(summary_dat_list)
  
  # return result
  # (a summary datatable and a prior table)
  list(summary_dat = summary_dat,
       priors = priorTable)
}

# try it
multiRes = multiBayesianAnalyses(connection,
                                 'eumaeus',
                                 database_id = 'IBM_MDCD',
                                 method = 'SCCS',
                                 exposure_id = 21184,
                                 analysis_ids = c(14:15),
                                 period_ids = c(2))
