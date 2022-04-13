# April 13 2022
# compare a few negative control test cases
# for GBS run

library(dplyr)

# load test results
respath = '~/Documents/Research/better_gbs/betterTest_mdcd/'
fname = 'hcSummary.csv'

summ = read_csv(file.path(respath, fname))


# load all negative control estimates from EUMAEUS
allNCs = readRDS('./localCache/CompNegControls.rds')
names(allNCs) = tolower(names(allNCs))

# compare
db = 'IBM_MDCD'
method = 'HistoricalComparator'

outcomes = unique(summ$outcomeId); outcomes = outcomes[outcomes != 343]
exposures = 211831 #unique(summ$exposureId)[1]
thePeriod = c(1) #max(summ$seqId)
analysis_ids = 1 #unique(summ$analysisId)[1]

# query rows from EUMAEUS estimates 
(NCs = allNCs %>% 
  filter(database_id == db, method == !!method,
         outcome_id %in% outcomes, 
         exposure_id %in% exposures,
         period_id %in% thePeriod,
         analysis_id %in% analysis_ids) %>%
    arrange(period_id,outcome_id))

# query rows from test results
(selsumm = summ %>% 
  filter(outcomeId %in% outcomes,
         analysisId %in% analysis_ids, 
         exposureId %in% exposures, 
         seqId %in% thePeriod) %>%
    select(outcomeId, logRr, seqId, targetOutcomes) %>%
    arrange(seqId, outcomeId))
# period=1, analysis 1 of HC, outcome 197494 had 0 target outcomes;
# original EUMAEUS analysis wouldn't produce any estimates at all


## also check if likelihood profiles are provided.....-------
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

# pull LP for the 0-count one
LP = getLikelihoodProfile(connection, 'eumaeus',
                          db, exposures, 197494, 1, 12, 
                          method, plot = TRUE)
# we do have an LP!! monotone and decreasing, just like every other guy with zero counts

DatabaseConnector::disconnect(connection)
