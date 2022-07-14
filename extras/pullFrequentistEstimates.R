# pull frequentist estimates from EUMAEUS results server

library(tidyverse)

frequentistEst <- function(connection,
                           schema,
                           database_id,
                           method,
                           exposure_id,
                           analysis_id){
  # pull estimates
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
  estimates <- DatabaseConnector::querySql(connection, sql)
  
  # get imputed PCs
  sql <- "SELECT outcome_id, exposure_id, negative_control_id, effect_size
        FROM @schema.IMPUTED_POSITIVE_CONTROL_OUTCOME"
  sql <- SqlRender::render(sql, 
                           schema = schema)
  IPCs <- DatabaseConnector::querySql(connection, sql)
  names(IPCs) = tolower(names(IPCs))
  
  # processing
  names(estimates) = tolower(names(estimates))
  estimates = estimates %>%
    filter(!is.na(log_rr) & !is.na(se_log_rr) & !is.na(llr) & !is.na(critical_value)) %>%
    filter(period_id == max(period_id)) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr,
           calibrated_ci_95_lb, calibrated_ci_95_ub,
           ci_95_lb, ci_95_ub)

  # get effect size and NC ids
  estimates = estimates %>%
    distinct() %>%
    left_join(IPCs) %>%
    select(database_id, method, analysis_id, 
           exposure_id, outcome_id, period_id, 
           p, log_rr, se_log_rr, llr, critical_value,
           calibrated_p, calibrated_log_rr, calibrated_se_log_rr,
           calibrated_llr, 
           ci_95_lb, ci_95_ub,
           calibrated_ci_95_lb, calibrated_ci_95_ub,
           effect_size, negative_control_id) %>%
    mutate(effect_size = if_else(is.na(effect_size), 1, effect_size),
           negativeControl = (effect_size == 1),
           negative_control_id = if_else(is.na(negative_control_id), 
                                         outcome_id, 
                                         negative_control_id))

  # look at the log_rr and calibrated_log_rr difference between PCs and NCs
  # (to check PC imputation)
  shifts = estimates %>% 
    group_by(negative_control_id, period_id) %>%
    arrange(effect_size) %>%
    mutate(correct_shift = log(effect_size),
           actual_shift = log_rr - log_rr[1],
           actual_shift_calibrated = calibrated_log_rr - calibrated_log_rr[1]) %>%
    mutate(actual_shift_exp = exp(actual_shift),
           actual_shift_calibrated_exp = exp(actual_shift_calibrated)) %>%
    ungroup()
  
  
  return(list(estimates = estimates, imputed_shifts = shifts))
  
}



# check an example---

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

## pull estimates for an analysis

db = 'CCAE'
me = 'HistoricalComparator'
#me = 'SCCS'
eid = 211981
aid = 6


estimates = frequentistEst(connection, 'eumaeus',
                           database_id = db, 
                           method = me,
                           exposure_id = eid,
                           analysis_id = aid)

## check the imputation 
# look at the min, mean, and max amount of estimation shift between NC and PC with different effect sizes
peek <- function(v) c(mean(v), min(v), max(v))

estimates$imputed_shifts %>% 
  group_by(effect_size) %>%
  summarize(log_shift = peek(actual_shift),
            shift = peek(actual_shift_exp),
            log_shift_calibrated = peek(actual_shift_calibrated),
            shift_calibrated = peek(actual_shift_calibrated_exp))
# can see that the imputation is "doubled"
# e.g., for effect = 1.5, estimate = NC_estimate * 1.5^2


DatabaseConnector::disconnect(connection)