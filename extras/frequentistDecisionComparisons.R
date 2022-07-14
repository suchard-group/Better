# 06/16/2022
# add frequentist comparisons

library(tidyverse)
library(ggpubr)
library(xtable)
library(wesanderson)

# set up EUMAEUS results query connection-----
## connection details
# ConnectionDetails <- DatabaseConnector::createConnectionDetails(
#   dbms = "postgresql",
#   server = paste(keyring::key_get("eumaeusServer"),
#                  keyring::key_get("eumaeusDatabase"),
#                  sep = "/"),
#   user = keyring::key_get("eumaeusUser"),
#   password = keyring::key_get("eumaeusPassword"))
# 
# ## set up the DB connection
# connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)


# 1 function to extract results from EUMAEUS data server and make decisions using MaxSPRT-------
# do this for one analysis and all outcomes
frequentistDecisions <- function(connection,
                                 schema,
                                 database_id,
                                 method,
                                 exposure_id,
                                 analysis_id,
                                 estimates = NULL,
                                 calibration = FALSE,
                                 correct_shift = FALSE,
                                 cachePath = './localCache/'){
  # pull estimates
  if(is.null(estimates)){
    # sql <- "SELECT estimate.*
    # FROM @schema.ESTIMATE estimate
    # WHERE database_id = '@database_id'
    #       AND method = '@method'
    #       AND analysis_id = @analysis_id
    #       AND exposure_id = @exposure_id"
    # sql <- SqlRender::render(sql, 
    #                          schema = schema,
    #                          database_id = database_id,
    #                          method = method,
    #                          analysis_id = analysis_id,
    #                          exposure_id = exposure_id)
    # sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
    # estimatesNC <- DatabaseConnector::querySql(connection, sql)
    
    sql <- "SELECT estimate.*
    FROM @schema.ESTIMATE_IMPUTED_PCS estimate
    WHERE database_id = '@database_id'
          AND method = '@method'
          AND analysis_id = @analysis_id
          AND exposure_id = @exposure_id"
    sql <- SqlRender::render(sql, 
                             schema = schema,
                             database_id = database_id,
                             method = method,
                             analysis_id = analysis_id,
                             exposure_id = exposure_id)
    sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
    estimatesPC <- DatabaseConnector::querySql(connection, sql)
  }
  
  # names(estimatesNC) = tolower(names(estimatesNC))
  # estimatesNC = estimatesNC %>%
  #   filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
  #   select(database_id, method, analysis_id, 
  #          exposure_id, outcome_id, period_id, 
  #          p, log_rr, se_log_rr, llr, critical_value,
  #          calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
  #          calibrated_llr)
  
  names(estimatesPC) = tolower(names(estimatesPC))
  estimatesPC = estimatesPC %>%
    filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr)
  
  ## combine NCs and PCs 
  #estimates = rbind(estimatesNC, estimatesPC)
  estimates = estimatesPC
  
  # get effect sizes
  IPCs = readRDS(file.path(cachePath, 'allIPCs.rds'))
  names(IPCs) = tolower(names(IPCs))
  
  estimates = estimates %>%
    left_join(IPCs) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr, effect_size) %>%
    mutate(effect_size = if_else(is.na(effect_size), 1, effect_size),
           negativeControl = (effect_size == 1))
  
  # 07/12/2022: correct for the amount of shift for IPCs
  if(correct_shift){
    estimates = correctShift(estimates, log_CI = FALSE, cachePath = cachePath)$shifted
  }
  
  
  # make decisions
  if(calibration){
    decisions = estimates %>%
      group_by(outcome_id) %>%
      arrange(period_id) %>%
      mutate(yes = (calibrated_llr > critical_value)) %>%
      mutate(reject = cumsum(yes)) %>%
      mutate(reject = (reject > 0)) %>%
      select(-yes) %>%
      ungroup()
  }else{
    decisions = estimates %>%
      group_by(outcome_id) %>%
      arrange(period_id) %>%
      mutate(yes = (llr > critical_value)) %>%
      mutate(reject = cumsum(yes)) %>%
      mutate(reject = (reject > 0)) %>%
      select(-yes) %>%
      ungroup()
  }
  
  
  # summarize Type I and Type II errors
  errorRate = decisions %>%
    group_by(database_id, method, analysis_id, 
             exposure_id, negativeControl,
             effect_size, period_id) %>%
    summarize(rejectRate = mean(reject, na.rm =TRUE)) %>%
    mutate(errorRate = if_else(negativeControl, rejectRate, 1-rejectRate),
           stats = if_else(negativeControl, 'type 1',
                           sprintf('type 2 (effect=%.1f)', effect_size))) %>%
    ungroup()
    
             
  # return
  return(list(estimates = estimates, calibrate = calibration,
              decisions = decisions, errorRate = errorRate))
  
}


# try it-----
# db = 'CCAE'
# me = 'HistoricalComparator'
# eid = 211981
# aid = 2
# 
# 
# ConnectionDetails <- DatabaseConnector::createConnectionDetails(
#   dbms = "postgresql",
#   server = paste(keyring::key_get("eumaeusServer"),
#                  keyring::key_get("eumaeusDatabase"),
#                  sep = "/"),
#   user = keyring::key_get("eumaeusUser"),
#   password = keyring::key_get("eumaeusPassword"))
# 
# # set up the DB connection
# connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)
# 
# resLst = frequentistDecisions(connection, 
#                               'eumaeus',
#                               database_id = db,
#                               method = me,
#                               exposure_id = eid,
#                               analysis_id = aid)
#   
# DatabaseConnector::disconnect(connection)

## 07/12/2022: correctiong shift version (with or without calibration)
# resLst = frequentistDecisions(connection,
#                               'eumaeus',
#                               database_id = db,
#                               method = me,
#                               exposure_id = eid,
#                               analysis_id = aid,
#                               calibration = TRUE, #FALSE
#                               correct_shift = TRUE)


# 07/12/2022
# correct the IPC shift in frequentist estimates -----
## note: there may or may not be CIs columns in the df
## BUT will be log_rr and calibrated_rr and NO NA's in the estimates columns
correctShift <- function(estimates, 
                         log_CI = FALSE,
                         cachePath = './localCache/'){
  
  # get IPC table
  IPCs = readRDS(file.path(cachePath, 'allIPCs.rds'))
  names(IPCs) = tolower(names(IPCs))
  
  # take log scale on the CIs if not on log scale
  if(log_CI & ('ci_95_lb' %in% names(estimates))){
    estimates = estimates %>%
      mutate(ci_95_lb = log(ci_95_lb),
             ci_95_ub = log(ci_95_ub),
             calibrated_ci_95_lb = log(calibrated_ci_95_lb),
             calibrated_ci_95_ub = log(calibrated_ci_95_ub))
  }
  
  # get NC ID for each outcome_id
  estimates = estimates %>% 
    left_join(IPCs %>% select(-effect_size), 
              by = c('outcome_id', 'exposure_id')) %>%
    mutate(negative_control_id = if_else(is.na(negative_control_id), 
                                         outcome_id,
                                         negative_control_id))
  
  # group by NC id to check shift
  shifts = estimates %>% 
    group_by(negative_control_id, period_id) %>%
    arrange(effect_size) %>%
    mutate(correct_shift = log(effect_size),
           actual_shift = log_rr - log_rr[1],
           actual_shift_calibrated = calibrated_log_rr - calibrated_log_rr[1]) %>%
    mutate(correction = correct_shift - actual_shift,
           correction_calibrated = correct_shift - actual_shift_calibrated) %>%
    ungroup()
  
  # make the corrections in terms of shift
  # LLR computation using Normal approximation, following EUMAEUS code (by Martijn)
  # here: NO shifting in terms of p-values!!!! only work on LLR's for now (don't really know how to deal with the p-value yet)
  shifts = shifts %>%
    mutate(log_rr = log_rr + correction,
           calibrated_log_rr = calibrated_log_rr + correction_calibrated) %>%
    mutate(llr = dnorm(log_rr, log_rr, se_log_rr, log = TRUE) - dnorm(0, log_rr, se_log_rr, log = TRUE),
           calirated_llr = dnorm(calibrated_log_rr, calibrated_log_rr, calibrated_se_log_rr, log = TRUE) - 
             dnorm(0, calibrated_log_rr, calibrated_se_log_rr, log = TRUE))
  
  # work on the CIs as well if there are CIs in the df
  if('ci_95_lb' %in% names(estimates)){
    shifts = shifts %>%
      mutate(ci_95_lb = ci_95_lb + correction,
             ci_95_ub = ci_95_ub + correction,
             calibrated_ci_95_lb = calibrated_ci_95_lb + correction_calibrated,
             calibrated_ci_95_ub = calibrated_ci_95_ub + correction_calibrated)
  }
  
  # de-select intermediate columns 
  shifts = shifts %>% select(-(correct_shift:correction_calibrated))
  
  # return results
  list(original_estimates = estimates, shifted = shifts)
}