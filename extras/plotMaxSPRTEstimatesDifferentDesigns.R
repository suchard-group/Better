# 02/09/2023
# plot to compare RR estimates by different designs
# HC v.s. SCCS

library(tidyverse)
library(ggpubr)
library(xtable)
library(wesanderson)

# copied from "extras/compareFrequentistBayesianEstimation.R"
# (1) extract MaxSPRT estimates
frequentistMSE <- function(connection,
                           schema,
                           database_id,
                           method,
                           exposure_id,
                           analysis_id,
                           calibration = FALSE,
                           correct_shift = FALSE,
                           localEstimatesPath = NULL,
                           cachePath = './localCache/'){
  # pull estimates
  if(is.null(localEstimatesPath)){
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
  
  # names(estimatesNC) = tolower(names(estimatesNC))
  # estimatesNC = estimatesNC %>%
  #   filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
  #   filter(period_id == max(period_id)) %>%
  #   select(database_id, method, analysis_id, 
  #          exposure_id, outcome_id, period_id, 
  #          p, log_rr, se_log_rr, llr, critical_value,
  #          calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
  #          calibrated_llr,
  #          calibrated_ci_95_lb, calibrated_ci_95_ub,
  #          ci_95_lb, ci_95_ub)
  
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
           calibrated_error = calibrated_log_rr - truth) %>%
    ungroup() %>%
    group_by(database_id, method, analysis_id, 
             exposure_id, effect_size) %>%
    summarise(mse = mean(error^2), calibrated_mse = mean(calibrated_error^2)) %>%
    ungroup()
  
  return(list(MSEs = MSEs, estimates = estimates))
  
}


# (a) pull estimates
estimatesPath = './localCache/EstimateswithImputedPcs_CCAE.rds'

eid = 211833 # HPV vaccine (both doses)

HCestimates = frequentistMSE(NULL,
                             'eumaeus',
                             database_id = 'CCAE',
                             method = 'HistoricalComparator',
                             exposure_id = eid,
                             analysis_id = 2,
                             correct_shift = TRUE,
                             localEstimatesPath = estimatesPath)

SCCSestimates = frequentistMSE(NULL,
                               'eumaeus',
                               database_id = 'CCAE',
                               method = 'SCCS',
                               exposure_id = eid,
                               analysis_id = 2,
                               correct_shift = TRUE,
                               localEstimatesPath = estimatesPath)

# (b) clean up and combine a bit

NCs = readRDS('./localCache/allNegativeControls.rds')
names(NCs) = tolower(names(NCs))

select_outcomes = NCs$outcome_id
#select_outcomes = sample(NCs$outcome_id, 60, replace = FALSE)
# select_outcomes = c(23731, 73302, 74719, 133424, 137977, 
#                     195500, 4145627, 4284982)

# extract uncalibrated estimate for RR, for NCs only
maxsprt_SCCS = SCCSestimates$estimates %>%
  filter(effect_size == 1, 
         outcome_id %in% select_outcomes) %>%
  mutate(estimate = exp(log_rr),
         ci_95_lb = exp(ci_95_lb),
         ci_95_ub = exp(ci_95_ub)) %>%
  select(outcome_id, estimate, ci_95_lb, ci_95_ub) %>%
  mutate(approach = 'Self-controlled') %>%
  filter(estimate > 1, ci_95_ub <= 5)

maxsprt_HC = HCestimates$estimates %>%
  filter(effect_size == 1, 
         outcome_id %in% unique(maxsprt_SCCS$outcome_id)) %>%
  mutate(estimate = exp(log_rr),
         ci_95_lb = exp(ci_95_lb),
         ci_95_ub = exp(ci_95_ub)) %>%
  select(outcome_id, estimate, ci_95_lb, ci_95_ub) %>%
  mutate(approach = 'Historical rates')

combined_estimates = bind_rows(maxsprt_HC, maxsprt_SCCS) %>%
  left_join(NCs, by = 'outcome_id')


# (c) make plots
theColors = wes_palette("Darjeeling2")[2:3]

dodgeWidth = -0.6

ggplot(combined_estimates, aes(#x=as.factor(outcome_id), 
                               x = outcome_name,
                               y = estimate, 
                               color=approach)) +
  geom_hline(yintercept = 1, color = 'gray70',
             size = 2) +
  geom_point(size = 2, position = position_dodge(width = dodgeWidth))+
  geom_errorbar(aes(ymin = ci_95_lb, ymax = ci_95_ub), 
                size = 1, width = 0.5,
                position = position_dodge(width = dodgeWidth)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0,1,2,4,6,8))+
  scale_color_manual(values = theColors) +
  labs(x = '', y = 'RR estimate for negative control outcomes',
       color = 'Estimates (95% CI) by:') +
  coord_flip()+
  theme_bw(base_size = 16)

