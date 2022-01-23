# Jan 2022
# function to run Bayesian analysis based on likelihood profiles
# with or without calibration (using negative control analysis)

library(EvidenceSynthesis)
source('./extras/getLikelihoodProfile.R')
source('./extras/fitNegativeControlDistribution.R')


## a little helper function to get maximum density estimate from samples
getMAP <- function(x){
  dens = density(x)
  dens$x[which.max(dens$y)]
}

# small module function for
# one (database, method, exposure, analysis, period) combo (with all available outcomes)
# can use pre-saved raw table of likelihood profiles

## need to go through all outcomes!!!!
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
  sprintf('Finished Bayesian analysis for database %s, exposure %s, in period %s, using %s analysis %s',
          database_id, exposure_id, period_id, method, analysis_id)
  res
}


# try it
bayesRes = oneBayesianAnalysis(connection,
                               'eumaeus',
                               database_id = 'IBM_MDCD',
                               method = 'SCCS',
                               exposure_id = 21184,
                               analysis_id = 1, 
                               period_id = 9)
