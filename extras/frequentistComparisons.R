# 06/16/2022
# add frequentist comparisons

library(tidyverse)
library(ggpubr)
library(xtable)

# set up EUMAEUS results query connection-----
## connection details
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

## set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)


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
# 06/28/2022: add 95% CIs
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
    filter(period_id == max(period_id)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr,
           calibrated_ci_95_lb, calibrated_ci_95_ub,
           ci_95_lb, ci_95_ub)
  
  names(estimatesPC) = tolower(names(estimatesPC))
  estimatesPC = estimatesPC %>%
    filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
    filter(period_id == max(period_id)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr,
           calibrated_ci_95_lb, calibrated_ci_95_ub,
           ci_95_lb, ci_95_ub)
  
  ## combine NCs and PCs 
  ## try: not include queries from "estimates"...
  #estimates = rbind(estimatesNC, estimatesPC)
  estimates = estimatesPC
  
  # get effect sizes
  IPCs = readRDS(file.path(cachePath, 'allIPCs.rds'))
  names(IPCs) = tolower(names(IPCs))
  
  all_outcome_ids = unique(c(IPCs$outcome_id, IPCs$negative_control_id))
  
  estimates = estimates %>%
    filter(outcome_id %in% all_outcome_ids) %>%
    distinct() %>%
    left_join(IPCs) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr, 
           ci_95_lb, ci_95_ub,
           calibrated_ci_95_lb, calibrated_ci_95_ub,
           effect_size) %>%
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
  
  return(list(MSEs = MSEs, estimates = estimates))

}

## try it
db = 'CCAE'
me = 'HistoricalComparator'
#me = 'SCCS'
eid = 211981
aid = 8
#aid = 6

#(
freqMSEs = frequentistMSE(connection,
                              'eumaeus',
                              database_id = db,
                              method = me,
                              exposure_id = eid,
                              analysis_id = aid)
#)


# 2.b MSEs for Bayesian version as well
# 06/28/2022: add 95% credible intervals
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
    select(postMedian, adjustedPostMedian, 
           CI95_lb, CI95_ub,
           adjustedCI95_lb, adjustedCI95_ub,
           outcome_id, prior_id)
  
  # get effect sizes
  IPCs = readRDS(file.path(cachePath, 'allIPCs.rds'))
  names(IPCs) = tolower(names(IPCs))
  
  estimates = estimates %>%
    left_join(IPCs) %>%
    select(outcome_id, prior_id, 
           postMedian, adjustedPostMedian, 
           CI95_lb, CI95_ub,
           adjustedCI95_lb, adjustedCI95_ub, 
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
# bMSE$MSEs

## combine frequentist and Bayesian results and output a table
# combined_mses = cbind(freqMSEs %>% select(effect_size, mse, calibrated_mse),
#                       bMSE$MSEs %>% select(mse, adjusted_mse))
# print(xtable(combined_mses, digits = 3), include.rownames = FALSE)

## some unexpected findings----
# 1. sometimes with Bayesian adjustments, MSE increases (compared to doing nothing)
# 2. Bayesian method NOT necessarily has less MSE compared to the frequentist method

# ## try to plot the error terms----
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


# 06/28/2022: 
# try to plot Bayesian estimates (post.median) and 95% CIs
# and compare with frequentist estimates (MLEs) and 95% CIs

festimates = freqMSEs$estimates %>%
  mutate(ci_95_lb = log(ci_95_lb),
         ci_95_ub = log(ci_95_ub),
         calibrated_ci_95_lb = log(calibrated_ci_95_lb),
         calibrated_ci_95_ub = log(calibrated_ci_95_ub)) %>%
  group_by(effect_size) %>%
  mutate(index = seq_along(log_rr)) %>%
  ungroup() %>%
  mutate(effect_size_label = factor(effect_size,
                                    levels = as.character(sort(unique(effect_size)))))

# only look at outcomes with MLEs in frequentist approach...
bestimates = bMSE$estimates %>%
  filter(outcome_id %in% unique(festimates$outcome_id)) %>%
  mutate(error = postMedian - truth,
         adjusted_error = adjustedPostMedian - truth) %>%
  group_by(effect_size) %>%
  mutate(index = seq_along(error)) %>%
  ungroup() %>%
  mutate(effect_size_label = factor(effect_size,
                                    levels = as.character(sort(unique(effect_size)))))

type2cols = c(wes_palette("Zissou1")[3:4],wes_palette("Royal1")[4])
othercols =  wes_palette("Royal1")[2]
allCols = c(othercols, type2cols)

## option to show only some effect sizes
allEffects = c(1,1.5,2,4)
effect_set = c(1,1.5,2,4)
#effect_set = c(1,4)

festimates = festimates %>% filter(effect_size %in% effect_set)
bestimates = bestimates %>% filter(effect_size %in% effect_set)
allCols = allCols[allEffects %in% effect_set]

#ylims = c(-8, 5)
ylims = c(-5, 5)


## (1) unadjusted estimates
(
p1 =
ggplot(bestimates, aes(y=postMedian,x=index,
                       color=effect_size_label)) +
  geom_hline(yintercept = log(effect_set), 
             size = 1, color='gray50')+
  geom_point() +
  geom_errorbar(aes(ymax = CI95_ub, ymin = CI95_lb))+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(limits = ylims)+
  scale_color_manual(values = allCols)+
  labs(y='estimate (log scale)', x='', color='effect size',
       caption='unadjusted inference results (Bayesian)')+
  theme_bw(base_size = 14) +
  theme(legend.position = 'none')
)

## (1.b) uncalibrated estimates (frequentists)
(
p1b =
  ggplot(festimates, aes(y=log_rr, x=index,
                         color=effect_size_label)) +
  geom_hline(yintercept = log(effect_set), 
             size = 1, color='gray50')+
  geom_point() +
  geom_errorbar(aes(ymax = ci_95_ub, ymin = ci_95_lb))+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(limits = ylims)+
  scale_color_manual(values = allCols)+
  labs(y='estimate (log scale)', x='', color='effect size',
       caption='uncalibrated estimates (frequentist)')+
  theme_bw(base_size = 14) +
  theme(legend.position = 'none')
)

## (2) adjusted estimates
(
  p2 =
    ggplot(bestimates, aes(y=adjustedPostMedian,x=index,
                           color=effect_size_label)) +
    geom_hline(yintercept = log(effect_set), 
               size = 1, color='gray50')+
    geom_point() +
    geom_errorbar(aes(ymax = adjustedCI95_ub, 
                      ymin = adjustedCI95_lb))+
    scale_x_continuous(breaks = NULL)+
    scale_y_continuous(limits = ylims)+
    scale_color_manual(values = allCols)+
    labs(y='estimate (log scale)', x='', color='effect size',
         caption='adjusted inference results (Bayesian)')+
    theme_bw(base_size = 14) +
    theme(legend.position = 'none')
)

## (2.b) calibrated estimates (frequentist)
(
  p2b =
    ggplot(festimates, aes(y=calibrated_log_rr, x=index,
                           color=effect_size_label)) +
    geom_hline(yintercept = log(effect_set), 
               size = 1, color='gray50')+
    geom_point() +
    geom_errorbar(aes(ymax = calibrated_ci_95_ub, 
                      ymin = calibrated_ci_95_lb))+
    scale_x_continuous(breaks = NULL)+
    scale_y_continuous(limits = ylims)+
    scale_color_manual(values = allCols)+
    labs(y='estimate (log scale)', x='', color='effect size',
         caption='calibrated estimates (frequentist)')+
    theme_bw(base_size = 14) +
    theme(legend.position = 'none')
)


## save plots in a document
ppath = '~/Documents/Research/betterResults/estimationResults/'
pdf(file.path(ppath, 'estimationChecks-SCCS1.pdf'),
    width = 15, height = 8)

ggarrange(ggarrange(p1b, p2b, ncol=2),
          ggarrange(p1, p2, ncol=2),
          nrow = 2)

dev.off()

## look at confidence interval/credible interval widths before/after calibration/adjustment
freq_CI_widths = festimates %>%
  mutate(CI_width = ci_95_ub - ci_95_lb,
         calibrated_CI_width = calibrated_ci_95_ub - calibrated_ci_95_lb) %>%
  group_by(effect_size) %>%
  summarise(width_avg = mean(CI_width, na.rm =TRUE),
            width_90_lb = quantile(CI_width, 0.05, na.rm =TRUE),
            width_90_ub = quantile(CI_width, 0.95, na.rm =TRUE),
            calibrated_width_avg = mean(calibrated_CI_width, na.rm =TRUE),
            calibrated_width_90_lb = quantile(calibrated_CI_width, 0.05, na.rm =TRUE),
            calibrated_width_90_ub = quantile(calibrated_CI_width, 0.95, na.rm =TRUE))

Bayes_CI_widths = bestimates %>%
  mutate(CI_width = CI95_ub - CI95_lb,
         adjusted_CI_width = adjustedCI95_ub - adjustedCI95_lb) %>%
  group_by(effect_size) %>%
  summarise(width_avg = mean(CI_width, na.rm =TRUE),
            width_90_lb = quantile(CI_width, 0.05, na.rm =TRUE),
            width_90_ub = quantile(CI_width, 0.95, na.rm =TRUE),
            adjusted_width_avg = mean(adjusted_CI_width, na.rm =TRUE),
            adjusted_width_90_lb = quantile(adjusted_CI_width, 0.05, na.rm =TRUE),
            adjusted_width_90_ub = quantile(adjusted_CI_width, 0.95, na.rm =TRUE))

print(
  xtable(cbind(freq_CI_widths, Bayes_CI_widths[,-1]), digits = 3),
  include.rownames = FALSE
)


## 06/30/2022: look at distance between imputed PCs and NCs
##             is the shift somewhat off? 

(
effect_dist = festimates %>% group_by(index) %>%
  arrange(effect_size) %>%
  summarize(log_rr_dist = diff(log_rr),
            cali_log_rr_dist = diff(calibrated_log_rr),
            effect_diff = effect_size[-1]) %>%
  ungroup()
)

## is the shift all the same?
ggplot(effect_dist, aes(x=index, y=log_rr_dist, color = as.factor(effect_diff))) +
  geom_point()
  