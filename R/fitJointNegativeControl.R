# 04/29/2022: fit negative control (systematic error) distribution
# using all likelihood profiles of NCs within one analysis
# using a hierarchical normal model 
# (instead of just a normal using NC estimates and SDs)

fitJointNegativeControl <- function(connection,
                                    schema,
                                    database_id,
                                    method,
                                    exposure_id,
                                    analysis_id, 
                                    period_id,
                                    savedLikehoods = NULL,
                                    outcomeToExclude=NULL,
                                    numsamps = 10000,
                                    thin = 10,
                                    minNCs = 0,
                                    randomSeed = 42,
                                    plot = FALSE){
  # outcomeToExclude: one or more (negative) outcomes to NOT include
  # numsamps: total num of posterior samples to acquire (default = 10)
  # thin: thinning iters (default = 10)
  # minNCs: minimum number of NCs before fitting the null distribution
  
  if(is.null(savedLikehoods)){

    sql <- 
  "SELECT database_id, 
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
    AND period_id = @period_id;"
    
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             exposure_id = exposure_id,
                             period_id = period_id,
                             analysis_id = analysis_id)
    likelihoods <- DatabaseConnector::querySql(connection, sql)
    cat('Negative controls lieklihood profiles pulled...\n')
  }else{
    likelihoods = savedLikehoods %>%
      filter(DATABASE_ID == database_id,
             METHOD == method,
             EXPOSURE_ID == exposure_id,
             ANALYSIS_ID == analysis_id, 
             PERIOD_ID == period_id)
  }
  
  
  # check if anything returned
  # check if minNCs number of likelihoods is available
  if(nrow(likelihoods) <= minNCs){
    ParallelLogger::logInfo(sprintf('Minimum %s of negative control likelihood profiles not available!\n',
                                    minNCs))
    res = numeric(0)
  }else{
    if(!is.null(outcomeToExclude)){
      likelihoods = likelihoods %>% filter(!OUTCOME_ID %in% outcomeToExclude)
    }
    if(nrow(likelihoods) <= minNCs){
      # if after exclusion, nothing, return empty
      ParallelLogger::logInfo(sprintf('Minimum %s of negative control likelihoods not available!\n',
                                      minNCs))
      res = numeric(0)
    }else{
      # otherwise, fit a systematic error (hierarchical normal) distribution 
      # with MCMC
      # on those likelihoods
      # get mean and precision parameters samples back
      names(likelihoods) = tolower(names(likelihoods))
      
      # transform it to a list of one-row dataframes
      likelihoods = split(likelihoods, seq(nrow(likelihoods)), drop = TRUE)
      
      cat(length(likelihoods), '\n', names(likelihoods[[1]]), '\n')
      
      ## run MCMC to get post samples for mean and precision
      null <- EvidenceSynthesis::computeBayesianMetaAnalysis(data = likelihoods,
                                                             chainLength = numsamps * thin,
                                                             burnIn = 1e+4,
                                                             subSampleFrequency = thin,
                                                             seed = randomSeed)
      
      cat('Systematic error distribution fitted...\n')
      traces <- attr(null, "traces")
      means = traces[, 1]
      sds = traces[, 2] #???
      #sds = 1/sqrt(precs)
      
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


# try it
null = fitJointNegativeControl(connection = connection,
                               schema = 'eumaeus',
                               database_id = "IBM_MDCD",
                               method = "SCCS",
                               exposure_id = 21184,
                               period_id = 9,
                               analysis_id = 1,
                               outcomeToExclude = NULL, # (Sedative withdrawal)
                               numsamps = 10000,
                               thin = 100,
                               minNCs = 5,
                               plot = TRUE)
