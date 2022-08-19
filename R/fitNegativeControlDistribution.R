# Jan 2022
# functions for computing negative control distribution
# for a (database, method, exposure, analysis, period) combo
# with optional function to exclude an outcome 
# (if doing LOO for negative controls)

# returns: a list of posterior samples for the normal mean and sd,
# as well as equal number of predictive samples of the systematic error
# (return numeric(0) if no negative controls results exist for required analysis)

#source('./R/getLikelihoodProfile.R')

# Feb 2022 update: make use of pre-pulled estimates if available 
fitNegativeControlDistribution <- function(connection,
                                           schema,
                                           database_id,
                                           method,
                                           exposure_id,
                                           analysis_id, 
                                           period_id,
                                           savedEstimates = NULL,
                                           outcomeToExclude=NULL,
                                           numsamps = 10000,
                                           thin = 10,
                                           plot = FALSE){
  # outcomeToExclude: one or more (negative) outcomes to NOT include
  # numsamps: total num of posterior samples to acquire (default = 10)
  # thin: thinning iters (default = 10)
  
  # get relevant data first
  if(is.null(savedEstimates)){
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
    AND exposure_id = @exposure_id
    AND analysis_id = @analysis_id
    AND period_id = @period_id;"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             exposure_id = exposure_id,
                             period_id = period_id,
                             analysis_id = analysis_id)
    estimates <- DatabaseConnector::querySql(connection, sql)
    cat('Negative controls estimates pulled...\n')
  }else{
    estimates = savedEstimates %>%
      filter(ANALYSIS_ID == analysis_id, 
             PERIOD_ID == period_id)
  }
  
  
  # check if anything returned
  if(nrow(estimates) == 0){
    res = numeric(0)
  }else{
    if(!is.null(outcomeToExclude)){
      estimates = estimates %>% filter(!OUTCOME_ID %in% outcomeToExclude)
    }
    if(nrow(estimates) == 0){
      # if after exclusion, nothing, return empty
      res = numeric(0)
    }else{
      # otherwise, fit a systematic error (normal) distribution with MCMC
      # get mean and precision parameters samples back
      names(estimates) = tolower(names(estimates))
      
      ## run MCMC to get post samples for mean and precision
      null <- EmpiricalCalibration::fitMcmcNull(logRr = estimates$log_rr, 
                                                seLogRr = estimates$se_log_rr,
                                                iter = numsamps * thin)
      cat('Systematic error distribution fitted...\n')
      mcmc <- attr(null, "mcmc")
      sel <- seq(from = thin, to = numsamps * thin, by = thin)
      means = mcmc$chain[sel, 1]
      precs = mcmc$chain[sel, 2]
      sds = 1/sqrt(precs)
      
      ## posterior predictive samples of the bias term
      biases = rnorm(numsamps, mean = means, sd = sds)
      
      ## return result list
      res = list(mean = means, sd = sds, bias = biases)
    }
    
  }
  
  # plotting
  if(plot & length(res) > 0){
    dat = data.frame(x=res$bias)
    print(
      ggplot(dat, aes(x=x)) +
        geom_density() +
        labs(x='log effect size', y=NULL, 
             title='Estimated systematic error bias distribution') +
        theme_bw(base_size = 14)
    )
  }
  
  res
}

# # try it
# NCs = fitNegativeControlDistribution(connection, 'eumaeus',
#                                      database_id = "IBM_MDCD",
#                                      method = "SCCS",
#                                      exposure_id = 21184,
#                                      period_id = 9,
#                                      analysis_id = 1,
#                                      outcomeToExclude = 43020446, # (Sedative withdrawal)
#                                      numsamps = 10000,
#                                      plot = TRUE)


# July 2022 update: 
# use likelihood profiles to fit null distribution
# Aug 2022 update:
# add "robust" option to allow a t-model
fitNegativeControlDistributionLikelihood <- function(connection,
                                                     schema,
                                                     database_id,
                                                     method,
                                                     exposure_id,
                                                     analysis_id, 
                                                     period_id,
                                                     savedLPs = NULL,
                                                     outcomeToExclude=NULL,
                                                     priorSds = c(2,0.5),
                                                     numsamps = 10000,
                                                     thin = 10,
                                                     minNCs = 5,
                                                     robust = FALSE,
                                                     plot = FALSE){
  # outcomeToExclude: one or more (negative) outcomes to NOT include
  # numsamps: total num of posterior samples to acquire (default = 10)
  # thin: thinning iters (default = 10)
  
  # get relevant data first
  if(is.null(savedLPs)){
    sql <- "SELECT database_id, 
    method, 
    exposure_id, 
    analysis_id, 
    period_id,
    likelihood_profile.outcome_id AS outcome_id,
    point,
    value
    FROM @schema.likelihood_profile
    INNER JOIN @schema.NEGATIVE_CONTROL_OUTCOME
    ON likelihood_profile.outcome_id = NEGATIVE_CONTROL_OUTCOME.outcome_id
    WHERE database_id = '@database_id'
          AND method = '@method'
          AND exposure_id = @exposure_id
          AND analysis_id = @analysis_id
          AND period_id = @period_id"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             exposure_id = exposure_id,
                             period_id = period_id,
                             analysis_id = analysis_id)
    LPs <- DatabaseConnector::querySql(connection, sql)
    cat('Negative controls estimates pulled...\n')
  }else{
    LPs = savedLPs
  }
  
  # exclude outcomes
  if(!is.null(outcomeToExclude)){
    LPs = LPs %>% filter(!OUTCOME_ID %in% outcomeToExclude)
  }
  
  # check if there are at least minNCs rows 
  if(nrow(LPs) < minNCs){
    message(sprintf('Less than %s NC likelihood profiles available! Not fitting null distribution.\n',
                    minNCs))
    return(numeric(0))
  }
  
  # post-process to a list of dataframes
  LPlist = postProcessLPs(LPs)
  
  # fit null distribution
  null = EvidenceSynthesis::computeBayesianMetaAnalysis(data = LPlist,
                                                        chainLength = numsamps + numsamps * thin,
                                                        burnIn = numsamps,
                                                        subSampleFrequency = thin,
                                                        priorSd = priorSds,
                                                        robust = robust)
  traces = attr(null, "traces")
  means = traces[,1]
  sds = traces[,2]
  
  biases = rnorm(numsamps, means, sds)
  
  res = list(mean = means, sd = sds, bias = biases)
  
  # plotting
  if(plot & length(res) > 0){
    dat = data.frame(x=res$bias)
    print(
      ggplot(dat, aes(x=x)) +
        geom_density() +
        labs(x='log effect size', y=NULL, 
             title='Estimated systematic error bias distribution') +
        theme_bw(base_size = 14)
    )
  }
  
  res
}

# # try it
# NCs = fitNegativeControlDistributionLikelihood(connection, 'eumaeus',
#                                      database_id = "CCAE",
#                                      method = "HistoricalComparator",
#                                      exposure_id = 211981,
#                                      period_id = 12,
#                                      analysis_id = 2,
#                                      outcomeToExclude = 43020446, # (Sedative withdrawal)
#                                      numsamps = 10000,
#                                      priorSds = c(0.2,0.2),
#                                      plot = TRUE)


# Aug 2022 function to return the fitted null only
fitNegativeControlNullOnly <- function(connection,
                                       schema,
                                       database_id,
                                       method,
                                       exposure_id,
                                       analysis_id, 
                                       period_id,
                                       savedEstimates = NULL,
                                       outcomeToExclude=NULL,
                                       numsamps = 10000,
                                       thin = 10,
                                       minNCs = 5){
  # outcomeToExclude: one or more (negative) outcomes to NOT include
  # numsamps: total num of posterior samples to acquire (default = 10)
  # thin: thinning iters (default = 10)
  
  # get relevant data first
  if(is.null(savedEstimates)){
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
    AND exposure_id = @exposure_id
    AND analysis_id = @analysis_id
    AND period_id = @period_id;"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             exposure_id = exposure_id,
                             period_id = period_id,
                             analysis_id = analysis_id)
    estimates <- DatabaseConnector::querySql(connection, sql)
    cat('Negative controls estimates pulled...\n')
  }else{
    estimates = savedEstimates %>%
      filter(EXPOSURE_ID == exposure_id,
             METHOD == method,
             ANALYSIS_ID == analysis_id, 
             PERIOD_ID == period_id)
  }
  
  
  # check if at least `minNCs` results are available
  if(nrow(estimates) < minNCs){
    cat(sprintf('Available negative control estimates < %s! Skipped.\n', minNCs))
    return(numeric(0))
  }
  if(!is.null(outcomeToExclude)){
    estimates = estimates %>% filter(!OUTCOME_ID %in% outcomeToExclude)
    if(nrow(estimates) < minNCs){
      cat(sprintf('Available negative control estimates < %s! Skipped.\n', minNCs))
      return(numeric(0))
    }
  }
  
  # fit a systematic error (normal) distribution with MCMC
  names(estimates) = tolower(names(estimates))
  
  ## run MCMC to get post samples for mean and precision
  null <- EmpiricalCalibration::fitMcmcNull(logRr = estimates$log_rr, 
                                            seLogRr = estimates$se_log_rr,
                                            iter = numsamps * thin)
  cat('Systematic error distribution fitted...\n')
  
  return(null)
}

# # try it
# null = fitNegativeControlNullOnly(connection, 'eumaeus',
#                                      database_id = "IBM_MDCD",
#                                      method = "SCCS",
#                                      exposure_id = 21184,
#                                      period_id = 9,
#                                      analysis_id = 1,
#                                      outcomeToExclude = 43020446, # (Sedative withdrawal)
#                                      numsamps = 10000,
#                                   savedEstimates = CompNegControls %>% filter(DATABASE_ID == 'IBM_MDCD'))

