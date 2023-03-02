# March 1 2023
# generate database "characterization" table as data overview
# "Table 1" in results section

library(dplyr)
library(xtable)

# set up the DB connection----
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

# 1a. pull negative control estimates with info on----
# exposed subjects, exposed days, exposed outcomes
# comparing different designs, exposures, databases
sql = "SELECT estimate.*
      FROM eumaeus.ESTIMATE_IMPUTED_PCS estimate
      INNER JOIN eumaeus.NEGATIVE_CONTROL_OUTCOME
      ON estimate.outcome_id = NEGATIVE_CONTROL_OUTCOME.outcome_id
      WHERE (method = 'SCCS' OR method = 'HistoricalComparator')"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
allNegControls <- DatabaseConnector::querySql(connection, sql)
names(allNegControls) = tolower(names(allNegControls))

allNegControls = allNegControls %>% 
  select(database_id, method, analysis_id, 
         exposure_id, outcome_id, period_id, 
         exposure_subjects, counterfactual_subjects,
         exposure_days, counterfactual_days,
         exposure_outcomes, counterfactual_outcomes,
         rr) # add RR estimate here as well

saveRDS(allNegControls, './localCache/negativeControlEstimatesCharacteristics.rds')


# 1b pull database basic characteristics info----
# for possible later use
sql = 'SELECT * from eumaeus.DATABASE_CHARACTERIZATION'
databaseInfo = DatabaseConnector::querySql(connection, sql)
glimpse(databaseInfo)

names(databaseInfo) = tolower(names(databaseInfo))

saveRDS(databaseInfo, './localCache/database_characterization.rds')


DatabaseConnector::disconnect(connection)


# 2. summarize characteristics for Historical Comparator first
characterizeDataHC <- function(estimates, 
                               analysis_id = 1,
                               exposureSubset = NULL,
                               hasPositiveExposure = TRUE,
                               makeWideTable = TRUE,
                               cachePath = './localCache/'){
  
  if(is.null(exposureSubset)){
    exposureSet = unique(estimates$exposure_id)
  }else{
    exposureSet = exposureSubset
  }
  
  if(hasPositiveExposure){
    exposureLB = 0
  }else{
    exposureLB = -1
  }
  
  dataInfo = estimates %>% 
    filter(analysis_id == !!analysis_id,
           exposure_id %in% exposureSet,
           method == 'HistoricalComparator') %>%
    group_by(database_id, exposure_id) %>%
    filter(period_id == max(period_id)) %>%
    filter(exposure_subjects > exposureLB) %>%
    # summarise(average_exposure_subjects = mean(exposure_subjects, na.rm = TRUE),
    #           #se_exposure_subjects = sd(exposure_subjects, na.rm = TRUE),
    #           min_exposure_subjects = min(exposure_subjects, na.rm = TRUE),
    #           max_exposure_subjects = max(exposure_subjects, na.rm = TRUE),
    #           average_exposure_days = mean(exposure_days, na.rm = TRUE),
    #           #se_exposure_days = sd(exposure_days, na.rm = TRUE),
    #           min_exposure_days = min(exposure_days, na.rm = TRUE),
    #           max_exposure_days = max(exposure_days, na.rm = TRUE)) %>%
    summarize(total_exposure_subjects = max(exposure_subjects, na.rm = TRUE),
              total_exposure_days = max(exposure_days, na.rm = TRUE),
              median_exposure_outcomes = median(exposure_outcomes, na.rm = TRUE),
              quartile1_exposure_outcomes = quantile(exposure_outcomes, 
                                                     probs =.25, na.rm = TRUE),
              quartile3_exposure_outcomes = quantile(exposure_outcomes, 
                                                     probs =.75, na.rm = TRUE)) %>%
    ungroup()
  
  # reformat
  databases = unique(dataInfo$database_id)
  
  dataInfo = dataInfo %>% 
    #arrange(exposure_id) %>%
    mutate(exposure_outcomes = 
             sprintf("%.1f [%.1f, %.1f]", 
                     median_exposure_outcomes,
                     quartile1_exposure_outcomes,
                     quartile3_exposure_outcomes)) %>%
    mutate(subject_count = format(total_exposure_subjects, big.mark = ','),
           exposure_days = format(total_exposure_days, big.mark = ',')) %>%
    mutate(database = case_when(
      database_id == 'IBM_MDCD' ~ 'MDCD',
      database_id == 'IBM_MDCR' ~ 'MDCR',
      database_id == 'OptumEhr' ~ 'Optum EHR',
      database_id == 'OptumDod' ~ 'Optum',
      TRUE ~ 'CCAE')) %>%
    select(exposure_id, database,
           subject_count,
           exposure_days,
           exposure_outcomes)
  
  # check if for any exposure/database there is an entry missing
  # and fill them in...
  template = tibble(exposure_id = rep(unique(dataInfo$exposure_id), each = length(databases)),
                    database = rep(unique(dataInfo$database), length(unique(dataInfo$exposure_id))),
                    subject_count  = format(0, big.mark = ','),
                    exposure_days = format(0, big.mark = ','),
                    exposure_outcomes = 'N.A.')
  leftOverRows = template %>% anti_join(dataInfo, by = c('exposure_id', 'database'))
  dataInfo = dataInfo %>% bind_rows(leftOverRows)
                                   
  # get exposure names instead of id
  exposures = readRDS(file.path(cachePath, 'exposures.rds'))
  
  #colNames = c(names(dataInfo), 'exposure')
  dataInfo = dataInfo %>% inner_join(exposures, by = 'exposure_id') %>%
    rename(exposure = exposure_name) %>%
    arrange(exposure_id) %>%
    mutate(exposure = if_else(database == 'CCAE', exposure, '')) %>%
    select(exposure, database,
           subject_count,
           exposure_days,
           exposure_outcomes)
  
  # cut the long table into two halves...
  if(makeWideTable && nrow(dataInfo) %% 2 != 0){
    cat('Number of rows not even! Will not cut the table into two columns...\n')
  }else if(makeWideTable){
    halfIndex = nrow(dataInfo)/2
    dataInfo = bind_cols(dataInfo[1:halfIndex,],
                         dataInfo[(halfIndex+1):nrow(dataInfo),])
  }
  
  # make exposure name strings shorter...
  # move the "(xxx)" to the second row to save horizontal space
  exposureNameRows = seq(from = 1, to = nrow(dataInfo), by = nData)
  for(i in exposureNameRows){
    cuts = stringr::str_split(dataInfo[i,1], '\\(') %>% unlist()
    dataInfo[i,1] = cuts[1]
    if(length(cuts) > 1){
      dataInfo[i+1,1] = paste0("(",cuts[2])
    }
  }
  
  # add \midrule in between exposure blocks
  nData = length(databases)
  addRuleIndex = seq(from = nData + 1, to = nrow(dataInfo), by = nData)
  if(length(addRuleIndex) > 0){
    for(i in addRuleIndex){
      dataInfo[i,1] = paste0("\\midrule ", dataInfo[i,1])
    }
  }
  
  
  return(dataInfo)
}


## use it to generate formatted table
## CONTENT only; 
## still need latex headers etc.

estimates = readRDS('./localCache/negativeControlEstimatesCharacteristics.rds')
dataInfo = characterizeDataHC(estimates,
                              analysis_id = 1,
                              hasPositiveExposure = TRUE,
                              makeWideTable = FALSE)

print(xtable(dataInfo, format = "latex"),
      include.rownames = FALSE,
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE,
      sanitize.text.function = identity)
