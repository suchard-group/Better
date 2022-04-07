# Jan 2022
# functions for computing negative control distribution
# for a (database, method, exposure, analysis, period) combo
# with optional function to exclude an outcome 
# (if doing LOO for negative controls)

# returns: a list of posterior samples for the normal mean and sd,
# as well as equal number of predictive samples of the systematic error
# (return numeric(0) if no negative controls results exist for required analysis)

# Feb 2022 update: make use of pre-pulled estimates if available 
# April 2022: (1) use pre-saved NC estimates and filter on db, method, exposure, etc.
#             (2) force a `minNCs` -- if num. of NCs < minNCs, then don't fit null at all!
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
                                           plot = FALSE,
                                           minNCs = 0){
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
      filter(DATABASE_ID == database_id,
             METHOD == method,
             EXPOSURE_ID == exposure_id,
             ANALYSIS_ID == analysis_id, 
             PERIOD_ID == period_id)
  }
  
  
  # check if anything returned
  # check if minNCs number of estimates is available
  if(nrow(estimates) <= minNCs){
    ParallelLogger::logInfo(sprintf('Minimum %s of negative control estimates not available!\n',
                                    minNCs))
    res = numeric(0)
  }else{
    if(!is.null(outcomeToExclude)){
      estimates = estimates %>% filter(!OUTCOME_ID %in% outcomeToExclude)
    }
    if(nrow(estimates) <= minNCs){
      # if after exclusion, nothing, return empty
      ParallelLogger::logInfo(sprintf('Minimum %s of negative control estimates not available!\n',
                                      minNCs))
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
