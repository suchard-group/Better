# 04/04/2022
# Try running multiple analyses for GBS

library(tidyverse)
library(better)
library(stringr)

# connection details
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

# test code ------
# code to run
LPpath = '~/Documents/Research/better_gbs/Results_CCAE/'
LPfname = 'likelihood_profile.csv'
savepath = '~/Documents/Research/better_gbs/summary'
sampspath = '~/Documents/Research/better_gbs/samples'
allRes = better::multiGBSAnalyses(connection, 'eumaeus',
                          database_id = 'CCAE', 
                          method = 'HistoricalComparator',
                          exposure_id = 211981,
                          analysis_ids = c(1:12),
                          period_ids = c(1:12),
                          LPpath = LPpath,
                          LPfname = LPfname,
                          savepath = savepath,
                          sampspath = sampspath,
                          maxCores = 4)
