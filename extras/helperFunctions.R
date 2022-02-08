# helper function for Bayesian analyses using likelihood profiles

## function to get all imputed positive control outcomes
getIPCs <- function(connection, schema, 
                    cacheFolder, fname='allIPCs.rds'){
  # if not exist, then query it
  if(file.exists(file.path(cacheFolder,fname))){
    #cat('IPCs already saved. Skipping.\n')
    IPCs = readRDS(file.path(cacheFolder,fname))
  }else{
    sql <- "SELECT outcome_id,exposure_id, negative_control_id, effect_size 
            FROM @schema.IMPUTED_POSITIVE_CONTROL_OUTCOME"
    sql <- SqlRender::render(sql, 
                             schema = schema)
    IPCs <- DatabaseConnector::querySql(connection, sql)
    
    # save it
    saveRDS(IPCs, file = file.path(cacheFolder, fname))
  }
  # return the table
  IPCs
}

## try it
# IPCs = getIPCs(connection, 'eumaeus', './localCache/')

# IPClist = split(IPCs, 
#                 list(IPCs$EXPOSURE_ID, IPCs$NEGATIVE_CONTROL_ID),
#                 drop = TRUE)


## function to return list of IPCs for any (exposure, negControl) combo
## as a named vector (name=effectSize)
getThisIPCs <- function(exposure_id, negControl_id, IPCtable){
  thisIPCs = IPCtable %>% 
    filter(EXPOSURE_ID == exposure_id, NEGATIVE_CONTROL_ID == negControl_id) %>%
    arrange(EFFECT_SIZE)
  
  # if somehow imputed pos controls don't exist for this one... make up outcome_ids
  if(nrow(thisIPCs) == 0){
    res = negControl_id * 100 + c(15,20,40)
    names(res) = as.character(c(1.5,2.0,4.0))
  }else{
    res = thisIPCs$OUTCOME_ID
    names(res) = as.character(thisIPCs$EFFECT_SIZE)
  }
  
  # return named vector
  res
  
}

## try it
# getThisIPCs(exposure_id = 21184, negControl_id = 443421, IPCtable = IPCs)


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


## helper function to "nudge" likelihood profile given the imputed effect size
nudgeLikelihood <- function(lik, effectSize){
  if(effectSize == 1){
    lik
  }else{
    points = names(lik) %>% as.numeric()
    names(lik) = as.character(points + log(effectSize))
    lik
  }
}


## helper function to get the relevant negative control estimates
## to use for fitting negative control (systematic error) distributions
getNegControlEstimates <- function(connection,
                                   schema,
                                   database_id,
                                   method,
                                   exposure_id,
                                   analysis_id = NULL,
                                   period_id = NULL) {
  # query relevant data
  
  ## if analysis_id not specified, get all analyses for this period_id
  if (is.null(analysis_id)) {
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
    AND period_id = @period_id;"
    sql <- SqlRender::render(
      sql,
      schema = schema,
      database_id = database_id,
      method = method,
      exposure_id = exposure_id,
      period_id = period_id
    )
  } else{
    ## vice versa for period_id
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
    AND analysis_id = @analysis_id;"
    sql <- SqlRender::render(
      sql,
      schema = schema,
      database_id = database_id,
      method = method,
      exposure_id = exposure_id,
      analysis_id = analysis_id
    )
  }
  
  estimates <- DatabaseConnector::querySql(connection, sql)
  estimates
}

estimates = getNegControlEstimates(connection, 'eumaeus', 'IBM_MDCD', 'SCCS', 
                                   21184, period_id = 9)
