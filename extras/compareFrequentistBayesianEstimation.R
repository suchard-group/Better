# 07/12/2022: move the estimation checking code here

library(tidyverse)
library(ggpubr)
library(xtable)
library(wesanderson)

# set up EUMAEUS results query connection-----
# ## connection details
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

# 2. compute estimation error (MSE?)-----
# 06/28/2022: add 95% CIs
# 07/12/2022: add option to correct for the shifts for IPCs
# 03/03/2023: add CI coverage rates
# 03/03/2023: add an option to include pre-loaded results to save loading time
frequentistMSE <- function(connection,
                           schema,
                           database_id,
                           method,
                           exposure_id,
                           analysis_id,
                           calibration = FALSE,
                           correct_shift = FALSE,
                           localEstimatesPath = NULL,
                           localEstimates = NULL,
                           cachePath = './localCache/'){
  # pull estimates
  if(is.null(localEstimates)){
    if(is.null(localEstimatesPath)){
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
    }else{
      estimatesPC = readRDS(localEstimatesPath)
      names(estimatesPC) = SqlRender::camelCaseToSnakeCase(names(estimatesPC))
      
      estimatesPC = estimatesPC %>%
        filter(database_id == !!database_id,
               method == !!method,
               analysis_id == !!analysis_id,
               exposure_id == !!exposure_id)
    }
  }else{
    estimatesPC = localEstimates %>%
      filter(database_id == !!database_id,
             method == !!method,
             analysis_id == !!analysis_id,
             exposure_id == !!exposure_id)
  }
  
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
  rm(estimatesPC)
  
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
  
  # check if estimates are empty
  if(nrow(estimates) < 1){
    cat(sprintf('No viable MaxSPRT estimates for %s, %s, exposure %s, analysis %s! Skipped.\n',
                database_id, method, exposure_id, analysis_id))
    return(list(estimates = NULL, MSEs = NULL))
  }
  
  # do shift corrections here...
  if(correct_shift){
    estimates = correctShift(estimates, log_CI = TRUE, cachePath = cachePath)$shifted
  }
  
  # calculate estimation error and MSE across effect sizes
  # do this for both calibrated and uncalibrated results
  MSEs = estimates %>%
    mutate(truth = log(effect_size)) %>%
    group_by(outcome_id) %>%
    filter(period_id == max(period_id)) %>%
    mutate(error = log_rr - truth,
           calibrated_error = calibrated_log_rr - truth,
           cover = (ci_95_lb <= truth & ci_95_ub >= truth),
           calibrated_cover = (calibrated_ci_95_lb <= truth & calibrated_ci_95_ub >= truth)) %>%
    ungroup() %>%
    group_by(database_id, method, analysis_id, 
             exposure_id, effect_size) %>%
    summarise(mse = mean(error^2), 
              calibrated_mse = mean(calibrated_error^2),
              coverage = mean(cover, na.rm = TRUE),
              calibrated_coverage = mean(calibrated_cover, na.rm = TRUE)) %>%
    ungroup()
  
  return(list(MSEs = MSEs, estimates = estimates))
  
}

## try it
db = 'CCAE'
me = 'HistoricalComparator'
#me = 'SCCS'
eid = 211833
aid = 1
#aid = 6

#(
# freqMSEs = frequentistMSE(connection,
#                           'eumaeus',
#                           database_id = db,
#                           method = me,
#                           exposure_id = eid,
#                           analysis_id = aid)
#)

# # 07/12/2022: try a version with shift correction
# freqMSEs2 = frequentistMSE(connection,
#                            'eumaeus',
#                            database_id = db,
#                            method = me,
#                            exposure_id = eid,
#                            analysis_id = aid,
#                            correct_shift = TRUE)

# 02/07/2023: try using local saved estimates (Postgre server too slow!)

estimatesPath = './localCache/EstimateswithImputedPcs_CCAE.rds'

freqMSE2 = frequentistMSE(NULL,
                          'eumaeus',
                          database_id = db,
                          method = me,
                          exposure_id = eid,
                          analysis_id = aid,
                          correct_shift = TRUE,
                          localEstimatesPath = estimatesPath)
freqMSE2$MSEs

# 03/03/2023: try with loaded local estimates...
allIpcEstimates = readRDS('./localCache/allIpcEstimates.rds')

freqMSE3 = frequentistMSE(NULL,
                          'eumaeus',
                          database_id = db,
                          method = me,
                          exposure_id = eid,
                          analysis_id = aid,
                          correct_shift = TRUE,
                          localEstimatesPath = NULL,
                          localEstimates = allIpcEstimates)
freqMSE3$MSEs

# 2.b MSEs for Bayesian version as well
# 06/28/2022: add 95% credible intervals
# 03/03/2023: subset on outcome_ids (to focus on MaxSPRT estimable outcomes only)
# 03/03/2023: add CI coverage summary
BayesMSE <- function(summaryPath, 
                     database,
                     method,
                     exposure_id,
                     analysis_id,
                     outcomeSubset = NULL,
                     prior = 3,
                     cachePath = './localCache/'){
  
  fname = sprintf('%s_%s_%s_period12_analysis%s_summary.rds',
                  database, method, exposure_id, analysis_id)
  
  if(file.exists(file.path(summaryPath, fname))){
    # try with period12 first
    summ = readRDS(file.path(summaryPath, fname))
  }else{
    # if doesn't exist, try period9
    fname = sprintf('%s_%s_%s_period9_analysis%s_summary.rds',
                    database, method, exposure_id, analysis_id)
    if(file.exists(file.path(summaryPath, fname))){
      summ = readRDS(file.path(summaryPath, fname))
    }else{
      # if still doesn't exist, output a message and move on
      cat(sprintf('Results summary doesn\'t exist for %s, %s, exposure %s and analysis %s! Skipped.\n',
                  database, method, exposure_id, analysis_id))
      return(list(estimates = NULL, MSEs = NULL))
    }
  }
  
  #summ = readRDS(fpath)
  
  estimates = summ %>% 
    filter(prior_id == prior) %>%
    select(postMedian, adjustedPostMedian, 
           CI95_lb, CI95_ub,
           adjustedCI95_lb, adjustedCI95_ub,
           outcome_id, prior_id)
  
  # subset on outcome_id if...
  if(!is.null(outcomeSubset)){
    estimates = estimates %>%
      filter(outcome_id %in% outcomeSubset)
  }
  
  # check if estimates are empty...
  if(nrow(estimates) < 1){
    cat(sprintf('No viable Bayesian estimates for %s, %s, exposure %s, analysis %s! Skipped.\n',
                database, method, exposure_id, analysis_id))
    return(list(estimates = NULL, MSEs=NULL))
  }
  
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
           adjusted_error = adjustedPostMedian - truth,
           cover = (CI95_lb <= truth & CI95_ub >= truth),
           adjusted_cover = (adjustedCI95_lb <= truth & adjustedCI95_ub >= truth)) %>%
    ungroup() %>%
    group_by(effect_size) %>%
    summarise(mse = mean(error^2), 
              adjusted_mse = mean(adjusted_error^2),
              coverage = mean(cover, na.rm = TRUE),
              adjusted_covearge = mean(adjusted_cover, na.rm = TRUE)) %>%
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

## subset on outcomes estimable by MaxSPRT only...
bMSE = BayesMSE(summaryPath = summarypath,
                database = db,
                method = me,
                exposure_id = eid,
                analysis_id = aid,
                #outcomeSubset = unique(freqMSE2$estimates$outcome_id),
                prior =2)
bMSE$MSEs


# 03/03/2023 UPDATE -------
# pool all MSEs (and coverage results for MaxSPRT and Bayesian results)
# and save!

#summarypath = '~/Documents/Research/betterResults/betterResults-CCAE/'
#estimatesPath = './localCache/EstimateswithImputedPcs_CCAE.rds'

allIpcEstimates = readRDS('./localCache/allIpcEstimates.rds')

exposure_ids = readRDS('./localCache/exposures.rds')$exposure_id
methods = c(rep('HistoricalComparator', 12), rep('SCCS', 15))
analysis_ids = c(1:12, 1:15)
nAnalyses = length(analysis_ids)

databases = c('CCAE', 'MDCD', 'MDCR', 'OptumEhr', 'OptumDod')

allMSEs = NULL

for(db in databases){
  if(db %in% c('MDCD', 'MDCR')){
    maxSPRTdb = paste0('IBM_', db)
  }else{
    maxSPRTdb = db
  }
  summarypath = sprintf('~/Documents/Research/betterResults/betterResults-%s/',
                        db)
  cat(sprintf('\n\nWorking on database: %s...\n\n', db))
  
  for(eid in exposure_ids){
    for(i in 1:nAnalyses){
      me = methods[i]; aid = analysis_ids[i]
      
      # MaxSPRT
      freqMSE = frequentistMSE(NULL,
                               'eumaeus',
                               database_id = maxSPRTdb,
                               method = me,
                               exposure_id = eid,
                               analysis_id = aid,
                               correct_shift = TRUE,
                               localEstimates = allIpcEstimates)
      #localEstimatesPath = estimatesPath)
      
      # Bayesian
      bMSE = BayesMSE(summaryPath = summarypath,
                      database = maxSPRTdb,
                      method = me,
                      exposure_id = eid,
                      analysis_id = aid,
                      outcomeSubset = unique(freqMSE$estimates$outcome_id),
                      prior = 2)
      
      # combine column-wise
      if(!is.null(freqMSE$MSEs) && !is.null(bMSE$MSEs)){
        this.chunk = cbind(freqMSE$MSEs %>% rename(maxSPRT_mse = mse,
                                                   maxSPRT_coverage = coverage),
                           bMSE$MSEs %>% select(-effect_size) %>%
                             rename(Bayesian_mse = mse,
                                    Bayesian_coverage = coverage))
      }else{
        cat(sprintf('Results not available for database %s, %s, exposure %s, analysis %s...\n',
                    db, me, eid, aid))
        this.chunk = NULL
      }
      
      allMSEs = rbind(allMSEs, this.chunk)
    }
  }
}

# save to local 
#saveRDS(allCCAE_MSEs, './localCache/allCCAE_MSEs.rds')
saveRDS(allMSEs, './localCache/allMSEs.rds')

# 03/23/2023: try to save a possibly (hopefully) de-bugged version
saveRDS(allMSEs, './localCache/allMSEs-2.rds')


# 04/18/2023 update -----
# add MSEs for CUIMC and HistoricalComparator
allIpcEstimatesCuimc = readRDS('./localCache/allIpcEstimates-CUIMC.rds')

exposure_ids = readRDS('./localCache/exposures.rds')$exposure_id
methods = c(rep('HistoricalComparator', 12))
analysis_ids = c(1:12)
nAnalyses = length(analysis_ids)

databases = c('CUIMC')

allMSEs = NULL

for(db in databases){
  if(db %in% c('MDCD', 'MDCR')){
    maxSPRTdb = paste0('IBM_', db)
  }else{
    maxSPRTdb = db
  }
  summarypath = sprintf('~/Documents/Research/betterResults/betterResults-%s/',
                        db)
  cat(sprintf('\n\nWorking on database: %s...\n\n', db))
  
  for(eid in exposure_ids){
    for(i in 1:nAnalyses){
      me = methods[i]; aid = analysis_ids[i]
      
      # MaxSPRT
      freqMSE = frequentistMSE(NULL,
                               'eumaeus',
                               database_id = maxSPRTdb,
                               method = me,
                               exposure_id = eid,
                               analysis_id = aid,
                               correct_shift = TRUE,
                               localEstimates = allIpcEstimatesCuimc)
      #localEstimatesPath = estimatesPath)
      
      # Bayesian
      bMSE = BayesMSE(summaryPath = summarypath,
                      database = maxSPRTdb,
                      method = me,
                      exposure_id = eid,
                      analysis_id = aid,
                      outcomeSubset = unique(freqMSE$estimates$outcome_id),
                      prior = 2)
      
      # combine column-wise
      if(!is.null(freqMSE$MSEs) && !is.null(bMSE$MSEs)){
        this.chunk = cbind(freqMSE$MSEs %>% rename(maxSPRT_mse = mse,
                                                   maxSPRT_coverage = coverage),
                           bMSE$MSEs %>% select(-effect_size) %>%
                             rename(Bayesian_mse = mse,
                                    Bayesian_coverage = coverage))
      }else{
        cat(sprintf('Results not available for database %s, %s, exposure %s, analysis %s...\n',
                    db, me, eid, aid))
        this.chunk = NULL
      }
      
      allMSEs = rbind(allMSEs, this.chunk)
    }
  }
}

# save to local 
#saveRDS(allCCAE_MSEs, './localCache/allCCAE_MSEs.rds')
saveRDS(allMSEs, './localCache/allMSEs-cuimc.rds')



# 03/23/2023: plot MSEs and coverage between maxSPRT and Bayesian----
## filter out some very bad/weird rows

# 04/19/2023: update MSE and coverage plots with CUIMC results
allMSEs = bind_rows(readRDS('./localCache/allMSEs-2.rds'),
                    readRDS('./localCache/allMSEs-cuimc.rds'))

MSEs = allMSEs %>% 
  filter(maxSPRT_mse < 100, adjusted_mse < 100)

effect_sizes = sort(unique(MSEs$effect_size))
effect_labs = paste0('effect size = ', effect_sizes)
names(effect_labs) = effect_sizes

# (1) compare MSEs 
ggplot(MSEs, aes(x=maxSPRT_mse, y = adjusted_mse)) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = 2, color = 'gray50', size = 0.8) +
  geom_point(size = 1.2) +
  scale_y_continuous(limits = c(min(MSEs$maxSPRT_mse), max(MSEs$maxSPRT_mse)))+
  labs(x = 'MSE by MaxSPRT', y = 'MSE by Bayesian correction') +
  facet_grid(.~effect_size, labeller = labeller(effect_size = effect_labs)) +
  theme_bw(base_size = 15)

# (2) compare coverage rates
ggplot(MSEs, aes(x=maxSPRT_coverage, 
                 y= adjusted_covearge)) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = 2, color = 'gray50', size = 0.8) +
  geom_hline(yintercept = 0.95, 
             linetype = 2, color = 'gray50', size = 0.8)+
  geom_vline(xintercept = 0.95, 
             linetype = 2, color = 'gray50', size = 0.8)+
  geom_point(size = 1.2) +
  scale_y_continuous(limits = c(min(MSEs$maxSPRT_coverage), max(MSEs$maxSPRT_coverage)))+
  labs(x = 'Coverage rate by MaxSPRT', 
       y = 'Coverage rate by \nBayesian correction') +
  facet_grid(.~effect_size, labeller = labeller(effect_size = effect_labs)) +
  theme_bw(base_size = 15)

# (3) produce summary statistics on results in (1) and (2)----
MSEs_long = tibble(mse = c(MSEs$maxSPRT_mse, MSEs$adjusted_mse),
                   method = rep(c('MaxSPRT', 'Bayesian'), 
                                each = nrow(MSEs)))

## try plotting histograms: still no good!
ggplot(MSEs_long, aes(x=mse,fill = method)) +
  geom_histogram(position = position_dodge()) +
  scale_x_continuous('Mean squared errors (MSEs) for log RR estimation') +
  theme_bw(base_size = 15)


## MSEs:

## utility function to produce summary statistics
get_summ_stats <- function(v){
  res = c(mean(v), quantile(v, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))
  return(res)
}

MSEs_summary_table = as.data.frame(rbind(get_summ_stats(MSEs$adjusted_mse),
                           get_summ_stats(MSEs$maxSPRT_mse)))

names(MSEs_summary_table) = c('Average', 
                              '10%', '25%', 'median', '75%', '90%')
MSEs_summary_table = cbind(c('Bayesian','MaxSPRT'),
                           MSEs_summary_table)
names(MSEs_summary_table)[1] = ''

print(xtable(MSEs_summary_table, digits = 3), include.rownames=FALSE)

## Coverage rates:
coverage_long = tibble(coverage = c(MSEs$maxSPRT_coverage, 
                               MSEs$adjusted_covearge),
                       effect_size = rep(MSEs$effect_size,2),
                       method = rep(c('MaxSPRT', 'Bayesian'), 
                                    each = nrow(MSEs)))

# 04/19/2023: update formatting of the coverage rates table
coverage_summary = coverage_long %>% 
  group_by(effect_size, method) %>%
  summarize(average = mean(coverage),
            Q1 = quantile(coverage, .25),
            median = median(coverage),
            Q3 = quantile(coverage, .75)) %>%
  ungroup() %>% 
  mutate(method = if_else(method == 'Bayesian', 'BBC', 'MLE')) %>%
  select(method, average, Q1, median, Q3)

print(xtable(coverage_summary, digits = 3), include.rownames=FALSE)

# formatted_coverage_summary = NULL
# for(es in unique(coverage_summary$effect_size)){
#   
# }



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



# 02/07/2023 new plot -----
# Make forrest-style (?) plots to compare MaxSPRT & Bayesian estimates
# for select negative control outcomes...
# to showcase de-biasing power

theColors = c(wes_palette("GrandBudapest1")[c(1,3)])

NCs = readRDS('./localCache/allNegativeControls.rds')
names(NCs) = tolower(names(NCs))

#select_outcomes = sample(NCs$outcome_id)
#select_outcomes = sample(NCs$outcome_id, 60, replace = FALSE)
select_outcomes = c(23731, 73302, 74719, 133424, 137977, 
                    195500, 4145627, 4284982)

# extract uncalibrated estimate for RR, for NCs only
maxsprt_estimates = freqMSE2$estimates %>%
  filter(effect_size == 1, 
         outcome_id %in% select_outcomes) %>%
  mutate(estimate = exp(log_rr),
         ci_95_lb = exp(ci_95_lb),
         ci_95_ub = exp(ci_95_ub)) %>%
  select(outcome_id, estimate, ci_95_lb, ci_95_ub) %>%
  mutate(approach = 'A') %>%
  filter(estimate > 1, ci_95_ub <= 5)
  
bayesian_estimates = bMSE$estimates %>%
  filter(outcome_id %in% unique(maxsprt_estimates$outcome_id)) %>%
  mutate(estimate = exp(adjustedPostMedian), 
         ci_95_lb = exp(adjustedCI95_lb),
         ci_95_ub = exp(adjustedCI95_ub)) %>%
  select(outcome_id, estimate, ci_95_lb, ci_95_ub) %>%
  mutate(approach = 'B')

combined_estimates = bind_rows(maxsprt_estimates, bayesian_estimates) %>%
  left_join(NCs, by = 'outcome_id')

## forrest-style plots for estimates comparison

dodgeWidth = 0.6

ggplot(combined_estimates, aes(x=outcome_name, 
                               y = estimate, 
                               color=approach)) +
  geom_hline(yintercept = 1, color = 'gray70',
             size = 2) +
  geom_point(size = 2, position = position_dodge(width = dodgeWidth))+
  geom_errorbar(aes(ymin = ci_95_lb, ymax = ci_95_ub), 
                size = 1, width = 0.5,
                position = position_dodge(width = dodgeWidth)) +
  scale_y_continuous(limits = c(0, 6), breaks = c(0,1,2,4,6))+
  scale_color_manual(values = theColors,
                     labels = c('MaxSPRT', 'Bayesian')) +
  labs(x = '', y = 'RR estimate for negative control outcomes',
       color = 'Estimates by:') +
  coord_flip()+
  theme_bw(base_size = 16)


# 06/28/2022: -----
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


## 06/30/2022: look at distance between imputed PCs and NCs--------
##             is the shift somewhat off? 

# (
# effect_dist = festimates %>% group_by(index) %>%
#   arrange(effect_size) %>%
#   summarize(log_rr_dist = log_rr[-1] - log_rr[1],
#             cali_log_rr_dist = calibrated_log_rr[-1] - calibrated_log_rr[1],
#             effect_diff = effect_size[-1]) %>%
#   ungroup()
# )
# 
# ## is the shift all the same?
# ggplot(effect_dist, aes(x=index, y=log_rr_dist, color = as.factor(effect_diff))) +
#   geom_point()
# 
# ## what's happening with the imputation shifting??
# exp(effect_dist$log_rr_dist)
# # [1]  2.25  4.00 16.00  2.25  4.00 16.00 ....
# seems that the shifting was done twice!!! that's why it was off!


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


#######
# # don't run
# # check synthetic PC results
# frequentistMSE_synPC <- function(connection,
#                            schema,
#                            database_id,
#                            method,
#                            exposure_id,
#                            analysis_id,
#                            estimates = NULL,
#                            calibration = FALSE,
#                            cachePath = './localCache/'){
#   # pull estimates
#   if(is.null(estimates)){
#     sql <- "SELECT estimate.*
#     FROM @schema.ESTIMATE estimate
#     WHERE database_id = '@database_id'
#           AND method = '@method'
#           AND analysis_id = @analysis_id
#           AND exposure_id = @exposure_id"
#     sql <- SqlRender::render(sql, 
#                              schema = schema,
#                              database_id = database_id,
#                              method = method,
#                              analysis_id = analysis_id,
#                              exposure_id = exposure_id)
#     sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
#     estimates<- DatabaseConnector::querySql(connection, sql)
#   }
#   
#   names(estimates) = tolower(names(estimates))
#   estimates = estimates %>%
#     filter(!is.na(log_rr) & !is.na(se_log_rr)) %>%
#     filter(period_id == max(period_id)) %>%
#     select(database_id, method, analysis_id, 
#            exposure_id, outcome_id, period_id, 
#            p, log_rr, se_log_rr, llr, critical_value,
#            calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
#            calibrated_llr,
#            calibrated_ci_95_lb, calibrated_ci_95_ub,
#            ci_95_lb, ci_95_ub)
#   
#   
#   # get effect sizes
#   sPCs = readRDS(file.path(cachePath, 'allsPCs.rds'))
#   names(sPCs) = tolower(names(sPCs))
#   sPCs = sPCs %>% select(-exposure_id)
# 
#   estimates = estimates %>%
#     #filter(outcome_id %in% all_outcome_ids) %>%
#     #distinct() %>%
#     left_join(sPCs, by='outcome_id') %>%
#     select(database_id, method, analysis_id, 
#            exposure_id, outcome_id, period_id,
#            p, log_rr, se_log_rr, llr, critical_value,
#            calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
#            calibrated_llr, 
#            ci_95_lb, ci_95_ub,
#            calibrated_ci_95_lb, calibrated_ci_95_ub,
#            effect_size, negative_control_id) %>%
#     mutate(effect_size = if_else(is.na(effect_size), 1, effect_size),
#            negativeControl = (effect_size == 1),
#            negative_control_id = if_else(is.na(negative_control_id), 
#                                          outcome_id, negative_control_id))
#   
#   # calculate estimation error and MSE across effect sizes
#   # do this for both calibrated and uncalibrated results
#   MSEs = estimates %>%
#     mutate(truth = log(effect_size)) %>%
#     group_by(outcome_id) %>%
#     filter(period_id == max(period_id)) %>%
#     mutate(error = log_rr - truth,
#            calibrated_error = calibrated_log_rr - truth) %>%
#     ungroup() %>%
#     group_by(database_id, method, analysis_id, 
#              exposure_id, effect_size) %>%
#     summarise(mse = mean(error^2), calibrated_mse = mean(calibrated_error^2)) %>%
#     ungroup()
#   
#   return(list(MSEs = MSEs, estimates = estimates))
#   
# }
# 
# ## try it
# db = 'CCAE'
# me = 'HistoricalComparator'
# #me = 'SCCS'
# eid = 211981
# aid = 5
# #aid = 6
# 
# #(
# freqMSEs_synPC = frequentistMSE_synPC(connection,
#                           'eumaeus',
#                           database_id = db,
#                           method = me,
#                           exposure_id = eid,
#                           analysis_id = aid)
# #)
# 
# # look at estimates and CIs with synthetic PCs instead
# festimates_synPC = freqMSEs_synPC$estimates %>%
#   mutate(ci_95_lb = log(ci_95_lb),
#          ci_95_ub = log(ci_95_ub),
#          calibrated_ci_95_lb = log(calibrated_ci_95_lb),
#          calibrated_ci_95_ub = log(calibrated_ci_95_ub)) %>%
#   mutate(effect_size_label = factor(effect_size,
#                                     levels = as.character(sort(unique(effect_size)))))
# 
# ## make plots of estimates and CIs
# type2cols = c(wes_palette("Zissou1")[3:4],wes_palette("Royal1")[4])
# othercols =  wes_palette("Royal1")[2]
# allCols = c(othercols, type2cols)
# 
# allEffects = c(1,1.5,2,4)
# effect_set = c(1,1.5,2,4)
# 
# festimates_synPC = festimates_synPC %>% filter(effect_size %in% effect_set)
# #bestimates = bestimates %>% filter(effect_size %in% effect_set)
# allCols = allCols[allEffects %in% effect_set]
# 
# ylims = c(-8,5)
# 
# (
#   p1b_synPC =
#     ggplot(festimates_synPC, 
#            aes(y=log_rr, x=as.factor(negative_control_id),
#                color=effect_size_label)) +
#     geom_hline(yintercept = log(effect_set), 
#                size = 1, color='gray50')+
#     geom_point() +
#     geom_errorbar(aes(ymax = ci_95_ub, ymin = ci_95_lb))+
#     scale_x_discrete(breaks = NULL)+
#     scale_y_continuous(limits = ylims)+
#     scale_color_manual(values = allCols)+
#     labs(y='estimate (log scale)', x='', color='effect size',
#          caption='uncalibrated estimates (frequentist)')+
#     theme_bw(base_size = 14) +
#     theme(legend.position = 'none')
# )
# 
# (
#   p2b_synPC =
#     ggplot(festimates_synPC, 
#            aes(y=calibrated_log_rr, 
#                x=as.factor(negative_control_id),
#                color=effect_size_label)) +
#     geom_hline(yintercept = log(effect_set), 
#                size = 1, color='gray50')+
#     geom_point() +
#     geom_errorbar(aes(ymax = ci_95_ub, ymin = ci_95_lb))+
#     scale_x_discrete(breaks = NULL)+
#     scale_y_continuous(limits = ylims)+
#     scale_color_manual(values = allCols)+
#     labs(y='estimate (log scale)', x='', color='effect size',
#          caption='calibrated estimates (frequentist)')+
#     theme_bw(base_size = 14) +
#     theme(legend.position = 'none')
# )
# 
# # the synthetic PC results don't look wayyy off
# # BUT they are only available for cases with a lot of data
# # for cases of insufficient data for the NC, there are no estimates produced