# helper function for Bayesian analyses using likelihood profiles

## function to get all imputed positive control outcomes
getIPCs <- function(connection, schema, 
                    cacheFolder, fname='allIPCs.rds'){
  # if not exist, then query it
  if(file.exists(file.path(cacheFolder,fname))){
    #cat('IPCs already saved. Skipping.\n')
    return()
  }else{
    sql <- "SELECT outcome_id,exposure_id, negative_control_id, effect_size 
            FROM @schema.IMPUTED_POSITIVE_CONTROL_OUTCOME"
    sql <- SqlRender::render(sql, 
                             schema = schema)
    IPCs <- DatabaseConnector::querySql(connection, sql)
  }
  
  # save it
  saveRDS(IPCs, file = file.path(cacheFolder, fname))
}

## try it
# getIPCs(connection, 'eumaeus', './localCache/')

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
  points = names(lik) %>% as.numeric()
  names(lik) = as.character(points + log(effectSize))
  lik
}
