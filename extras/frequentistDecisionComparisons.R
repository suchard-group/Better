# 06/16/2022
# add frequentist comparisons
# 07/27/2022
# add a naive Bonferroni correction sequential test version
# and also include a list of outcome_ids with available estimates

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


# (i) helper func to do robust division
# a/b = 0 if b=0 else a/b
robustDivide <- function(a, b){
  if(b==0){
    0
  }else{
    a/b
  }
}

# (ii) helper func to compute FDR given rejectCounts vector 
#      for effect sizes 1, 1.5, 2, 4
NAdivide <- function(a, b){
  if(is.na(a) || is.na(b) || b==0){
    NA
  }else{
    a/b
  }
}
computeFDR1 <- function(counts, effects){
  if(length(effects) == 4 || any(effects == 1)){
    res = sapply((counts[1] + counts[-1]), 
                 function(v) NAdivide(counts[1], v))
    #res = NAdivide(counts[1], (counts[1] + counts[-1]))
    res = c(NAdivide(counts[1], sum(counts)),res)
  }else{
    res = rep(NA, length(counts))
  }
  res
}

# 1 function to extract results from EUMAEUS data server and make decisions using MaxSPRT-------
# do this for one analysis and all outcomes
# 09/21/2022: add false discovery rate computation
frequentistDecisions <- function(connection,
                                 schema,
                                 database_id,
                                 method,
                                 exposure_id,
                                 analysis_id,
                                 estimates = NULL,
                                 calibration = FALSE,
                                 correct_shift = FALSE,
                                 bonferroni_baseline = TRUE,
                                 bonferroni_adjust_factor = 1,
                                 FDR = TRUE,
                                 alpha = 0.05,
                                 cachePath = './localCache/'){
  # pull estimates
  if(is.null(estimates)){
    
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
    
    names(estimatesPC) = tolower(names(estimatesPC))
    estimatesPC = estimatesPC %>%
      filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
      select(database_id, method, analysis_id, 
             exposure_id, outcome_id, period_id, 
             p, one_sided_p, log_rr, se_log_rr, llr, critical_value,
             calibrated_p, calibrated_one_sided_p, 
             calibrated_log_rr, calibrated_se_log_rr,
             calibrated_llr)
    
    estimates = estimatesPC
  }else{
    # use local data "estimates"
    names(estimates) = SqlRender::camelCaseToSnakeCase(names(estimates))
    
    estimates = estimates %>% 
      filter(database_id == !!database_id,
             method == !!method,
             analysis_id == !!analysis_id,
             exposure_id == !!exposure_id) %>%
      filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
      select(database_id, method, analysis_id, 
             exposure_id, outcome_id, period_id, 
             p, one_sided_p, log_rr, se_log_rr, llr, critical_value,
             calibrated_p, calibrated_one_sided_p, 
             calibrated_log_rr, calibrated_se_log_rr,
             calibrated_llr)
  }
  
  # get effect sizes
  IPCs = readRDS(file.path(cachePath, 'allIPCs.rds'))
  names(IPCs) = tolower(names(IPCs))
  
  estimates = estimates %>%
    left_join(IPCs) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, one_sided_p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_one_sided_p, 
           calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr, effect_size) %>%
    mutate(effect_size = if_else(is.na(effect_size), 1, effect_size),
           negativeControl = (effect_size == 1))
  
  # 07/12/2022: correct for the amount of shift for IPCs
  if(correct_shift){
    estimates = correctShift(estimates, log_CI = FALSE, cachePath = cachePath)$shifted
  }
  
  # make decisions
  if(bonferroni_baseline){
    bon_alpha = alpha/(max(estimates$period_id) * bonferroni_adjust_factor)
  }
  if(calibration){
    decisions = estimates %>%
      group_by(outcome_id) %>%
      arrange(period_id) %>%
      mutate(yes = (calibrated_llr > critical_value)) %>%
      mutate(reject = cumsum(yes)) %>%
      mutate(reject = (reject > 0)) %>%
      select(-yes) %>%
      ungroup()
    # also with the bonferroni correction...
    if(bonferroni_baseline){
      decisions_bonferroni = estimates %>%
        group_by(outcome_id) %>%
        arrange(period_id) %>%
        mutate(yes = (calibrated_one_sided_p < bon_alpha)) %>%
        mutate(reject = cumsum(yes)) %>%
        mutate(reject = (reject > 0)) %>%
        select(-yes) %>%
        ungroup()
    }
  }else{
    decisions = estimates %>%
      group_by(outcome_id) %>%
      arrange(period_id) %>%
      mutate(yes = (llr > critical_value)) %>%
      mutate(reject = cumsum(yes)) %>%
      mutate(reject = (reject > 0)) %>%
      select(-yes) %>%
      ungroup()
    if(bonferroni_baseline){
      decisions_bonferroni = estimates %>%
        group_by(outcome_id) %>%
        arrange(period_id) %>%
        mutate(yes = (one_sided_p < bon_alpha)) %>%
        mutate(reject = cumsum(yes)) %>%
        mutate(reject = (reject > 0)) %>%
        select(-yes) %>%
        ungroup()
    }
  }
  
  
  # summarize Type I and Type II errors
  ## 09/05/2022 change/debug:
  ## use the TOTAL number of controls as denominator
  ## NOT just the available outcomes up to time t
  ## otherwise, results don't make sense!!
  
  outcomeTotal = estimates %>% 
    filter(negativeControl == TRUE) %>% 
    pull(outcome_id) %>% 
    unique() %>% length()
  ## assuming that for each original negative control, its availability is SAME as the
  ## derived imputed positive control for each effect size
  
  errorRate = decisions %>%
    group_by(database_id, method, analysis_id, 
             exposure_id, negativeControl,
             effect_size, period_id) %>%
    #summarize(rejectRate = mean(reject, na.rm =TRUE)) %>%
    summarize(rejectCounts = sum(reject, na.rm =TRUE)) %>%
    mutate(errorRate = if_else(negativeControl, 
                               robustDivide(rejectCounts, outcomeTotal), 
                               1-robustDivide(rejectCounts, outcomeTotal)),
           stats = if_else(negativeControl, 'type 1',
                           sprintf('type 2 (effect=%.1f)', effect_size))) %>%
    ungroup()
  ## also with Bonferroni correction baseline
  if(bonferroni_baseline){
    errorRate_bonferroni = decisions_bonferroni %>%
      group_by(database_id, method, analysis_id, 
               exposure_id, negativeControl,
               effect_size, period_id) %>%
      summarize(rejectCounts = sum(reject, na.rm =TRUE)) %>%
      mutate(errorRate = if_else(negativeControl, 
                                 robustDivide(rejectCounts, outcomeTotal), 
                                 1-robustDivide(rejectCounts, outcomeTotal)),
             stats = if_else(negativeControl, 'type 1',
                             sprintf('type 2 (effect=%.1f)', effect_size))) %>%
      ungroup()
  }else{
    errorRate_bonferroni = NULL
  }
  
  # 09/21/2022: add FDR computation
  if(FDR){
    FDRs = decisions %>%
      group_by(database_id, method, analysis_id, 
               exposure_id, negativeControl,
               effect_size, period_id) %>%
      summarize(rejectCounts = sum(reject, na.rm =TRUE)) %>% 
      ungroup() %>%
      group_by(database_id, method, analysis_id, 
               exposure_id, period_id) %>%
      arrange(effect_size) %>%
      mutate(errorRate = computeFDR1(rejectCounts, effect_size),
             stats = if_else(negativeControl, 'FDR (effect=1.0)',
                             sprintf('FDR (effect=%.1f)', effect_size))) %>%
      ungroup()
  }else{
    FDRs = NULL
  }
  

  # return
  return(list(estimates = estimates, calibrate = calibration,
              decisions = decisions, errorRate = errorRate,
              errorRate_bonferroni = errorRate_bonferroni,
              FDRs = FDRs))
  
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