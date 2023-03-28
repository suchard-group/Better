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
## 03/24/2023 update: ----
# change to person-years AND add incidence rates column (divided by person-years)
characterizeDataHC <- function(estimates, 
                               analysis_id = 1,
                               exposureSubset = NULL,
                               hasPositiveExposure = TRUE,
                               makeWideTable = TRUE,
                               incidence_rate_factor = 1e5,
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
              total_person_years = max(exposure_days, na.rm = TRUE)/365.25,
              median_exposure_outcomes = median(exposure_outcomes, na.rm = TRUE),
              quartile1_exposure_outcomes = quantile(exposure_outcomes, 
                                                     probs =.25, na.rm = TRUE),
              quartile3_exposure_outcomes = quantile(exposure_outcomes, 
                                                     probs =.75, na.rm = TRUE)) %>%
    mutate(median_incidence_rates = median_exposure_outcomes/total_person_years * incidence_rate_factor,
           quartile1_incidence_rates = quartile1_exposure_outcomes/total_person_years * incidence_rate_factor,
           quartile3_incidence_rates = quartile3_exposure_outcomes/total_person_years * incidence_rate_factor) %>%
    ungroup()
  
  # reformat
  databases = unique(dataInfo$database_id)
  
  dataInfo = dataInfo %>% 
    #arrange(exposure_id) %>%
    mutate(exposure_outcomes = 
             sprintf("%.1f [%.1f, %.1f]", 
                     median_exposure_outcomes,
                     quartile1_exposure_outcomes,
                     quartile3_exposure_outcomes),
           incidence_rates = 
             sprintf("%.2f [%.2f, %.2f]", 
                     median_incidence_rates,
                     quartile1_incidence_rates,
                     quartile3_incidence_rates)) %>%
    mutate(subject_count = format(total_exposure_subjects, big.mark = ','),
           person_years = format(round(total_person_years,2), 
                                 big.mark = ',', nsmall = 2)) %>%
    mutate(database = case_when(
      database_id == 'IBM_MDCD' ~ 'MDCD',
      database_id == 'IBM_MDCR' ~ 'MDCR',
      database_id == 'OptumEhr' ~ 'Optum EHR',
      database_id == 'OptumDod' ~ 'Optum',
      TRUE ~ 'CCAE')) %>%
    select(exposure_id, 
           database,
           subject_count,
           person_years,
           exposure_outcomes,
           incidence_rates)
  
  # check if for any exposure/database there is an entry missing
  # and fill them in...
  template = tibble(exposure_id = rep(unique(dataInfo$exposure_id), each = length(databases)),
                    database = rep(unique(dataInfo$database), length(unique(dataInfo$exposure_id))),
                    subject_count  = format(0, big.mark = ','),
                    person_years = format(0, big.mark = ','),
                    exposure_outcomes = 'N.A.',
                    incidence_rates = 'N.A.')
  leftOverRows = template %>% anti_join(dataInfo, by = c('exposure_id', 'database'))
  dataInfo = dataInfo %>% bind_rows(leftOverRows)
                                   
  # get exposure names instead of id
  exposures = readRDS(file.path(cachePath, 'exposures.rds'))
  
  # construct a table with exposure names alone
  dataExposures = dataInfo %>% inner_join(exposures, by = 'exposure_id') %>%
    arrange(exposure_id) %>%
    select(exposure_id, exposure_name)
  
  #extract columns WITHOUT the exposure name column...
  dataInfo = dataInfo %>% 
    arrange(exposure_id) %>%
    #mutate(person_years = format(person_years, big.mark = ',', nsmall = 2)) %>%
    select(database,
           subject_count,
           person_years,
           exposure_outcomes,
           incidence_rates)
  
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
  # exposureNameRows = seq(from = 1, to = nrow(dataInfo), by = nData)
  # for(i in exposureNameRows){
  #   cuts = stringr::str_split(dataInfo[i,1], '\\(') %>% unlist()
  #   dataInfo[i,1] = cuts[1]
  #   if(length(cuts) > 1){
  #     dataInfo[i+1,1] = paste0("(",cuts[2])
  #   }
  # }
  
  # add \midrule in between exposure blocks
  # also: add exposure names at top of each block; not as a separate column
  
  ## deal with the first row first
  # i = 1
  # exposureString = sprintf("\\midrule \\multicolumn{5}{l}{\\textbf{%s}} \\\\ [0.25em] ",
  #                          dataExposures$exposure_name[i])
  # dataInfo[i,1] = paste0(exposureString, dataInfo[i,1])
  
  ## then also all the other blocks...
  nData = length(databases)
  addRuleIndex = seq(from = 1, to = nrow(dataInfo), by = nData)
  if(length(addRuleIndex) > 0){
    for(i in addRuleIndex){
      exposureString = sprintf("\\midrule \\multicolumn{5}{l}{\\textbf{%s}} \\\\ [0.25em] ",
                               dataExposures$exposure_name[i])
      dataInfo[i,1] = paste0(exposureString, dataInfo[i,1])
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
                              makeWideTable = FALSE,
                              incidence_rate_factor = 1e4)

print(xtable(dataInfo, format = "latex", digits = 2),
      include.rownames = FALSE,
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE,
      sanitize.text.function = identity)
