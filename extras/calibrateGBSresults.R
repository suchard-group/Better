# Aug 2022:
# script to perform empirical calibration for GBS analyses
# (to get frequentist results)

library(ParallelLogger)

source('R/fitNegativeControlDistribution.R')

calibrateGBSresults <- function(connection,
                                schema,
                                database_id,
                                numsamps = 10000,
                                thin = 10,
                                minNCs = 5,
                                resultsPath = '~/Documents/Research/better_gbs/',
                                localEstimates = NULL,
                                maxCores = 4,
                                saveResults = TRUE){
  # load GBS estimates
  EumaeusDataname = c("IBM_MDCR", "CCAE", "OptumDod", "IBM_MDCD", "OptumEhr")
  names(EumaeusDataname) = c('MDCR', 'CCAE', 'OptumDod', 'MDCD', 'OptumEHR')
  
  fpath = file.path(resultsPath, sprintf('Results_%s', database_id), 'estimate.csv')
  estimates = readr::read_csv(fpath)
  
  # load/query NC estimates
  db_id = EumaeusDataname[database_id]
  if(!is.null(localEstimates)){
    savedEstimates = localEstimates %>%
      filter(DATABASE_ID == db_id)
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
    AND (method = 'SCCS' OR method = 'HistoricalComparator')';"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = db_id)
    savedEstimates = DatabaseConnector::querySql(connection, sql)
    savedEstimates = savedEstimates %>% 
      select(DATABASE_ID, METHOD, ANALYSIS_ID,
      EXPOSURE_ID, OUTCOME_ID, PERIOD_ID, LOG_RR, SE_LOG_RR) %>%
      filter(!is.na(LOG_RR), !is.na(SE_LOG_RR))
    cat('Relevant NC estimates pulled!\n')
  }
  
  # subset and calibrate
  cluster <- ParallelLogger::makeCluster(min(10, maxCores))
  ParallelLogger::clusterRequire(cluster, "dplyr")
  subsets <- split(estimates, paste(estimates$databaseId, 
                                    estimates$method, 
                                    estimates$analysisId, 
                                    estimates$periodId, 
                                    estimates$exposureId))
  message("Computing calibrated one-sided p-values and LLRs")
  estimates <- ParallelLogger::clusterApply(cluster, subsets, calibrateGBS, 
                                            connection, schema, 
                                            numsamps, thin, minNCs,
                                            savedEstimates)
  estimates <- bind_rows(estimates)
  ParallelLogger::stopCluster(cluster)
  
  if(saveResults){
    savePath = file.path(resultsPath, sprintf('Results_%s', database_id), 'estimate_withCalibration.csv')
    write.csv(estimates, savePath)
    cat(sprintf('Empirical calibration for results of database %s is done!\n Results saved at: %s\n',
                database_id, savePath))
  }
  
  return(estimates)
  
}


calibrateGBS <- function(subset, 
                         connection, 
                         schema, 
                         numsamps = 10000, 
                         thin = 10, 
                         minNCs = 5,
                         savedEsimates = NULL){
  
  # deal with NA cases first
  if(all(is.na(subset$logRr))){
    subset$calibratedLogRr = NA
    subset$calibratedRr = NA
  }
  if(all(is.na(subset$seLogRr))){
    subset$calibratedSeLogRr <- NA
    subset$calibratedCi95Lb <- NA
    subset$calibratedCi95Ub <- NA
    subset$calibratedLlr = NA
  }
  if(all(is.na(subset$p))){
    subset$calibratedP = NA
  }
  if(all(is.na(subset$oneSidedP))){
    subset$calibratedOneSidedP = NA
  }
  
  # the non-NA cases
  if(any(!is.na(subset$logRr)) || any(!is.na(subset$seLogRr)) || 
     any(!is.na(subset$p)) || any(!is.na(subset$oneSidedP))){
    
    # fit null first
    null = fitNegativeControlNullOnly(connection, 
                                      schema,
                                      database_id = subset$databaseId[1],
                                      method = subset$method[1],
                                      exposure_id = subset$exposureId[1],
                                      analysis_id = subset$analysisId[1], 
                                      period_id = subset$periodId[1],
                                      savedEstimates = savedEsimates,
                                      numsamps = numsamps,
                                      thin = thin,
                                      minNCs = minNCs)
    
    if(length(null) == 0){
      # if no null dist. is fitted... calibrated results all NA
      subset$calibratedLogRr = NA
      subset$calibratedRr = NA
      subset$calibratedSeLogRr <- NA
      subset$calibratedCi95Lb <- NA
      subset$calibratedCi95Ub <- NA
      subset$calibratedLlr = NA
      subset$calibratedP = NA
      subset$calibratedOneSidedP = NA
    }else{
      calibratedP <- EmpiricalCalibration::calibrateP(null, subset$logRr, subset$seLogRr, 
                                                      twoSided = TRUE)
      subset$calibratedP <- calibratedP$p
      
      calibratedP <- EmpiricalCalibration::calibrateP(null, subset$logRr, subset$seLogRr, 
                                                      twoSided = FALSE, upper = TRUE)
      subset$calibratedOneSidedP <- calibratedP$p
      
      model <- EmpiricalCalibration::convertNullToErrorModel(null)
      calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(logRr = subset$logRr, 
                                                                        seLogRr = subset$seLogRr,
                                                                        model = model)
      subset$calibratedRr <- exp(calibratedCi$logRr)
      subset$calibratedLogRr <- calibratedCi$logRr
      subset$calibratedSeLogRr <- calibratedCi$seLogRr
      subset$calibratedCi95Lb <- exp(calibratedCi$logLb95Rr)
      subset$calibratedCi95Ub <- exp(calibratedCi$logUb95Rr)
      
      # LLR...
      subset$calibratedLlr = NA
      validIdx <- which(!is.na(subset$seLogRr)) 
      if (any(validIdx)) {
        null <- c(null[1], 1/sqrt(null[2]))
        names(null) <- c("mean", "sd")
        class(null) <- "null"
        calibratedLlr <- EmpiricalCalibration::calibrateLlr(null, subset[validIdx, ])
        calibratedLlr[is.infinite(calibratedLlr)] <- 9999
        subset$calibratedLlr[validIdx] <- calibratedLlr
      }
      
    }
    
  }
  
  return(subset)
}




## RUN CODE BELOW ------
CompNegControls = readRDS('./localCache/CompNegControls.rds')
GBSresultsPath = '~/Documents/Research/better_gbs/'

database = 'MDCD'
maxCores = 4

calibratedGBS = calibrateGBSresults(connection, 
                                    'eumaeus',
                                    database_id = database,
                                    localEstimates = CompNegControls,
                                    resultsPath = GBSresultsPath,
                                    maxCores = maxCores)
# write.csv(calibratedGBS, 
#           file.path(GBSresultsPath, 'Results_MDCR', 'estimate_withCalibration.csv'))
  
  
  
  
  
  