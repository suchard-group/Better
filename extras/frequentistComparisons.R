# 06/16/2022
# add frequentist comparisons

library(tidyverse)
library(ggpubr)
library(xtable)

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
                                 cachePath = './localCache/'){
  # pull estimates
  if(is.null(estimates)){
    sql <- "SELECT estimate.*
    FROM @schema.ESTIMATE estimate
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
    estimatesNC <- DatabaseConnector::querySql(connection, sql)
    
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
  
  names(estimatesNC) = tolower(names(estimatesNC))
  estimatesNC = estimatesNC %>%
    filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr)
  
  names(estimatesPC) = tolower(names(estimatesPC))
  estimatesPC = estimatesPC %>%
    filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr)
  
  ## combine NCs and PCs 
  estimates = rbind(estimatesNC, estimatesPC)
  
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
    summarize(rejectRate = mean(reject)) %>%
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


# 2. compute estimation error (MSE?)-----
frequentistMSE <- function(connection,
                           schema,
                           database_id,
                           method,
                           exposure_id,
                           analysis_id,
                           estimates = NULL,
                           calibration = FALSE,
                           cachePath = './localCache/'){
  # pull estimates
  if(is.null(estimates)){
    sql <- "SELECT estimate.*
    FROM @schema.ESTIMATE estimate
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
    estimatesNC <- DatabaseConnector::querySql(connection, sql)
    
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
  
  names(estimatesNC) = tolower(names(estimatesNC))
  estimatesNC = estimatesNC %>%
    filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr)
  
  names(estimatesPC) = tolower(names(estimatesPC))
  estimatesPC = estimatesPC %>%
    filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr)
  
  ## combine NCs and PCs 
  estimates = rbind(estimatesNC, estimatesPC)
  
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
  
  # calculate estimation error and MSE across effect sizes
  # do this for both calibrated and uncalibrated results
  MSEs = estimates %>%
    mutate(truth = log(effect_size)) %>%
    group_by(outcome_id) %>%
    filter(period_id == max(period_id)) %>%
    mutate(error = log_rr - truth,
           calibrated_error = calibrated_log_rr - truth) %>%
    ungroup() %>%
    group_by(database_id, method, analysis_id, 
             exposure_id, effect_size) %>%
    summarise(mse = mean(error^2), calibrated_mse = mean(calibrated_error^2)) %>%
    ungroup()
  
  return(MSEs)

}

## try it
db = 'CCAE'
me = 'HistoricalComparator'
eid = 211981
aid = 8

(
freqMSEs = frequentistMSE(connection,
                              'eumaeus',
                              database_id = db,
                              method = me,
                              exposure_id = eid,
                              analysis_id = aid)
)


# 2.b MSEs for Bayesian version as well
BayesMSE <- function(summaryPath, 
                     database,
                     method,
                     exposure_id,
                     analysis_id,
                     prior = 3,
                     cachePath = './localCache/'){
  
  fname = sprintf('%s_%s_%s_period12_analysis%s_summary.rds',
                  database, method, exposure_id, analysis_id)
  
  fpath = file.path(summaryPath, fname)
  
  summ = readRDS(fpath)
  
  estimates = summ %>% 
    filter(prior_id == prior) %>%
    select(postMedian, adjustedPostMedian, outcome_id, prior_id)
  
  # get effect sizes
  IPCs = readRDS(file.path(cachePath, 'allIPCs.rds'))
  names(IPCs) = tolower(names(IPCs))
  
  estimates = estimates %>%
    left_join(IPCs) %>%
    select(outcome_id, prior_id, 
           postMedian, adjustedPostMedian, 
           effect_size) %>%
    mutate(effect_size = if_else(is.na(effect_size), 1, effect_size),
           negativeControl = (effect_size == 1)) %>%
    mutate(truth = log(effect_size))
  
  MSEs = estimates %>%
    group_by(outcome_id) %>%
    mutate(error = postMedian - truth,
           adjusted_error = adjustedPostMedian - truth) %>%
    ungroup() %>%
    group_by(effect_size) %>%
    summarise(mse = mean(error^2), adjusted_mse = mean(adjusted_error^2)) %>%
    ungroup()
  
  return(list(estimates = estimates, MSEs=MSEs))
}

# # try it
summarypath = '~/Documents/Research/betterResults/betterResults-CCAE/'
#fname = 'CCAE_HistoricalComparator_211981_period12_analysis2_summary.rds'

# db = 'CCAE'
# me = 'HistoricalComparator'
# eid = 211981
# aid = 6

bMSE = BayesMSE(summaryPath = summarypath,
                database = db,
                method = me,
                exposure_id = eid,
                analysis_id = aid,
                prior =3)
bMSE$MSEs

## combine frequentist and Bayesian results and output a table
combined_mses = cbind(freqMSEs %>% select(effect_size, mse, calibrated_mse),
                      bMSE$MSEs %>% select(mse, adjusted_mse))
print(xtable(combined_mses, digits = 3), include.rownames = FALSE)

## some unexpected findings
# 1. sometimes with Bayesian adjustments, MSE increases (compared to doing nothing)
# 2. Bayesian method NOT necessarily has less MSE compared to the frequentist method

# ## try to plot the error terms
# 
# bestimates = bMSE$estimates %>%
#   mutate(error = postMedian - truth,
#          adjusted_error = adjustedPostMedian - truth) %>%
#   group_by(effect_size) %>%
#   mutate(index = seq_along(error)) %>%
#   ungroup() %>%
#   mutate(effect_size_label = factor(effect_size,
#                                     levels = as.character(sort(unique(effect_size)))))
# 
# type2cols = c(wes_palette("Zissou1")[3:4],wes_palette("Royal1")[4])
# othercols =  wes_palette("Royal1")[2]
# allCols = c(othercols, type2cols)
# 
# p1 = 
# ggplot(bestimates, aes(y=error,x=index,
#                        color=effect_size_label)) +
#   geom_hline(yintercept = 0, size = 1, color='gray50')+
#   geom_point() +
#   scale_x_continuous(breaks = NULL)+
#   scale_y_continuous(limits = c(-2.5,2.5))+
#   scale_color_manual(values = allCols)+
#   labs(y='error term', x='', color='effect size',
#        caption='unadjusted estimates')+
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'none')
# 
# p2 = 
#   ggplot(bestimates, aes(y=adjusted_error,x=index,
#                          color=effect_size_label)) +
#   geom_hline(yintercept = 0, size = 1, color='gray50')+
#   geom_point() +
#   scale_x_continuous(breaks = NULL)+
#   scale_y_continuous(limits = c(-2.5,2.5))+
#   scale_color_manual(values = allCols)+
#   labs(y='error term', x='', color='effect size',
#        caption='adjusted estimates')+
#   theme_bw(base_size = 14)
# 
# ggarrange(p1, p2, ncol=2, labels = c('raw', 'adj'),
#           widths = c(3,3.8))

  