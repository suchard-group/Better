# 04/21/2023: prepare example profile likelihood data for `EvidenceSynthesis` vignette

#### 1. extract some example LPs from eumaeus results schema-----
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

# extract some likelihood profiles as example
sql <- "SELECT * FROM eumaeus.likelihood_profile
        WHERE method = 'HistoricalComparator'
              AND database_id = 'CCAE'
              AND analysis_id = 1
              AND exposure_id = 211983;"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
exampleLPs <- DatabaseConnector::querySql(connection, sql)
names(exampleLPs) = tolower(names(exampleLPs))

## keep LPs negative controls...
NCids = readRDS("./localCache/allIPCs.rds")$NEGATIVE_CONTROL_ID %>% unique()
ncLPs = exampleLPs %>% filter(outcome_id %in% NCids)

## BUT! also save one example for an old synthetic positive control
pcLPs = exampleLPs %>% filter(!outcome_id %in% NCids)
saveRDS(pcLPs, './localCache/syntheticPClikelihoodProfilesExample.rds')

## pick one positive outcome that I like...
pcLPs = pcLPs %>% filter(outcome_id == 10205)

DatabaseConnector::disconnect(connection)


#### 2. make the LPs into LP list
source('extras/getLikelihoodProfile.R')

postProcessPeriodLPs <- function(LPdata){
  periods = sort(unique(LPdata$period_id))
  
  LPlist = list()
  
  # for(p in periods){
  #   LPlist[[as.character(p)]] = list()
  # }
  
  for(p in periods){
    this.list = postProcessLPs(LPdata %>% filter(period_id == p),
                               name_by_outcome = FALSE)
    
    LPlist[[as.character(p)]] = this.list
  }
  
  return(LPlist)
}

ooiLikelihoods = postProcessPeriodLPs(pcLPs)
ncLikelihoods = postProcessPeriodLPs(ncLPs)

## save to local RData file
save(ooiLikelihoods, file = './localCache/ooiLikelihoods.RData')
save(ncLikelihoods, file = './localCache/ncLikelihoods.RData')
