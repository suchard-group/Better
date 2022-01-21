# Jan 2022
# some extra stuff I run to save certain results locally

library(tidyverse)

# connection details
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)


## 1. get a table of all available negative control analysis results
sql = "SELECT estimate.*
      FROM eumaeus.ESTIMATE estimate
      INNER JOIN eumaeus.NEGATIVE_CONTROL_OUTCOME
      ON estimate.outcome_id = NEGATIVE_CONTROL_OUTCOME.outcome_id
      WHERE (method = 'SCCS' OR method = 'HistoricalComparator')"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
allNegControls <- DatabaseConnector::querySql(connection, sql)

names(allNegControls) = tolower(names(allNegControls))
### save it
### only need database, method, exposure, outcome, analysis, period and 
### log_rr & se_log_rr
allNegControls = allNegControls %>% 
  select(database_id, method, analysis_id, exposure_id, outcome_id, period_id, 
         log_rr, se_log_rr)
#saveRDS(allNegControls, './localCache/NegControlsSCCSHistComparator.rds')

### save a summary too
### count how many results are available for each existing analysis 
### (that there is ANY neg control result at all)
# NC_summary = allNegControls %>% 
#   group_by(database_id, method, exposure_id, analysis_id, period_id) %>%
#   count()
# saveRDS(NC_summary, './localCache/NegControlsSummarySCCSHistComparator.rds')

### try out one particular analysis
exNC = allNegControls %>% filter(method == 'SCCS', period_id == 5, 
                                 analysis_id ==1, database_id == 'IBM_MDCD',
                                 exposure_id == 21184)

DatabaseConnector::disconnect(connection)
