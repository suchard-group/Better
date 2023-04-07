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

## 0. get a table of all databases
sql = 'SELECT * from eumaeus.DATABASE'
databases = DatabaseConnector::querySql(connection, sql)
databases$DATABASE_ID

## save a local copy of database info
names(databases) = tolower(names(databases))
saveRDS(databases, './localCache/database.rds')

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
### April 14: also pull ci_95_lb and ci_95_ub
allNegControls = allNegControls %>% 
  select(database_id, method, analysis_id, 
         exposure_id, outcome_id, period_id, 
         log_rr, se_log_rr,
         ci_95_lb, ci_95_ub)
#saveRDS(allNegControls, './localCache/NegControlsSCCSHistComparator.rds')

allNegControls %>% filter(!is.na(log_rr) & !is.na(se_log_rr)) %>% count()
# 738460 with estimates

## save complete obs to local
allNegControls = allNegControls %>% 
  filter(!is.na(log_rr) & !is.na(se_log_rr))
names(allNegControls) = toupper(names(allNegControls))

# saveRDS(allNegControls, './localCache/CompNegControls.rds')

### save a summary too
### count how many results are available for each existing analysis 
### (that there is ANY neg control result at all)
### filter out the NA results first
NC_summary = allNegControls %>%
  filter(!is.na(log_rr) & !is.na(se_log_rr)) %>%
  group_by(database_id, method, exposure_id, analysis_id, period_id) %>%
  count()
saveRDS(NC_summary, './localCache/NegControlsSummarySCCSHistComparator.rds')

### try out one particular analysis
exNC = allNegControls %>% filter(method == 'SCCS', period_id == 5, 
                                 analysis_id ==1, database_id == 'IBM_MDCD',
                                 exposure_id == 21184)

# 1b: pull imputed positive control outcomes...
sql = "SELECT estimateipc.*
      FROM eumaeus.ESTIMATE_IMPUTED_PCS estimateipc
      WHERE (method = 'SCCS' OR method = 'HistoricalComparator')"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
allIpcEstimates <- DatabaseConnector::querySql(connection, sql)

names(allIpcEstimates) = tolower(names(allIpcEstimates))

saveRDS(allIpcEstimates, './localCache/allIpcEstimates.rds')

# 2. a summary table of all analyses
sql = "SELECT * from eumaeus.ANALYSIS
       WHERE method = 'HistoricalComparator'
       OR method = 'SCCS'"
analyses = DatabaseConnector::querySql(connection, sql)


# 3. some kind of summary for likelihood profiles

## Negative control results
sql <- "SELECT database_id, method, 
        exposure_id, analysis_id, period_id,
        likelihood_profile.outcome_id AS outcome_id
        FROM eumaeus.likelihood_profile
        INNER JOIN eumaeus.NEGATIVE_CONTROL_OUTCOME
        ON likelihood_profile.outcome_id = NEGATIVE_CONTROL_OUTCOME.outcome_id
        WHERE (method = 'SCCS' OR method = 'HistoricalComparator')"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
NegControlsLik <- DatabaseConnector::querySql(connection, sql)
### 837860 rows with LPs saved
names(NegControlsLik) = tolower(names(NegControlsLik))
### arrange it a little bit and give it an ID column
NegControlsLik = NegControlsLik %>%
  arrange(database_id, method, exposure_id, outcome_id, 
          analysis_id, period_id)
IDs = seq(from = 1, to = nrow(NegControlsLik), by = 1)
NegControlsLik$ID = paste0('NC', IDs)
### save it
saveRDS(NegControlsLik, './localCache/LikProfilesNegControlsSCCSHistComparator.rds')

## Imputed positive controls
sql <- "SELECT database_id, method, 
          likelihood_profile.exposure_id AS exposure_id,
          analysis_id, period_id,
          likelihood_profile.outcome_id AS outcome_id
        FROM eumaeus.likelihood_profile
        INNER JOIN eumaeus.IMPUTED_POSITIVE_CONTROL_OUTCOME
          ON likelihood_profile.outcome_id = IMPUTED_POSITIVE_CONTROL_OUTCOME.outcome_id
          AND likelihood_profile.exposure_id = IMPUTED_POSITIVE_CONTROL_OUTCOME.exposure_id
        WHERE (method = 'SCCS' OR method = 'HistoricalComparator')"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
PosControlsLik <- DatabaseConnector::querySql(connection, sql)

### ??? no likelihood profiles saved for imputed positive controls ???

sql <- "SELECT database_id, method, 
          exposure_id,
          analysis_id, period_id,
          outcome_id
        FROM eumaeus.likelihood_profile
        WHERE (method = 'SCCS' OR method = 'HistoricalComparator')
          AND outcome_id = 219728539"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
PosControlsLik <- DatabaseConnector::querySql(connection, sql)

### try to get ALL outcome results from likelihood profile and see if any imputed PCs exist in there
sql <- "SELECT database_id, method, exposure_id, outcome_id,
        MAX(period_id) max_period
        FROM eumaeus.likelihood_profile
        WHERE (method = 'SCCS' OR method = 'HistoricalComparator')
        GROUP BY database_id, method, exposure_id, outcome_id"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
max_periods_LP <- DatabaseConnector::querySql(connection, sql)

### check if any imputed pos controls have likelihood profiles at all
PosControls = DatabaseConnector::querySql(connection, 'SELECT outcome_id from eumaeus.IMPUTED_POSITIVE_CONTROL_OUTCOME')
PC_ids = unique(PosControls$OUTCOME_ID)
LP_ids = unique(max_periods_LP$OUTCOME_ID)
sum(sapply(PC_ids, function(x) x %in% LP_ids))
### !!!! NONE of the imputed positive control outcomes have likelihood profiles!!!!


# 4. experiment with likelihood profile pull
## single period, multiple analyses --> ~10sec on Mac, okay
sql <- "SELECT * 
          FROM eumaeus.LIKELIHOOD_PROFILE
          WHERE database_id = 'IBM_MDCD'
          AND method = 'SCCS'
          AND exposure_id = 21184
          AND period_id = 9"

## single analysis, multiple periods --> slightly less time
sql <- "SELECT * 
          FROM eumaeus.LIKELIHOOD_PROFILE
          WHERE database_id = 'IBM_MDCD'
          AND method = 'SCCS'
          AND exposure_id = 21184
          AND analysis_id = 1"

## multiple analyses + multiple periods?
sql <- "SELECT * 
          FROM eumaeus.LIKELIHOOD_PROFILE
          WHERE database_id = 'IBM_MDCD'
          AND method = 'SCCS'
          AND exposure_id = 21184"

## time it
t1 = Sys.time()
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
LPs <- DatabaseConnector::querySql(connection, sql)
Sys.time() - t1

# 5. look at imputed positive control outcomes
sql <- "SELECT outcome_id, exposure_id, negative_control_id, effect_size
        FROM eumaeus.IMPUTED_POSITIVE_CONTROL_OUTCOME"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
IPCs <- DatabaseConnector::querySql(connection, sql)

# 6. pull synthetic positive controls 
sql <- "SELECT outcome_id, exposure_id, negative_control_id, effect_size
        FROM eumaeus.POSITIVE_CONTROL_OUTCOME"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
sPCs <- DatabaseConnector::querySql(connection, sql)

DatabaseConnector::disconnect(connection)
