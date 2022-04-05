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

# run for multiple methods & exposures----------
database = 'CCAE'
expos = c(211981, 211982, 211983) # Zoster 1st, 2nd and either dose
methods = c('HistoricalComparator', 'SCCS')

LPpath = sprintf('~/Documents/Research/better_gbs/Results_%s/', database)
LPfname = 'likelihood_profile.csv'
savepath = '~/Documents/Research/better_gbs/summary'
sampspath = '~/Documents/Research/better_gbs/samples'

logFile = sprintf('log-%s.txt', database)
errorFile = sprintf('errorReport-%s.txt', database)
ParallelLogger::addDefaultFileLogger(file.path(savepath, logFile))
ParallelLogger::addDefaultErrorReportLogger(file.path(savepath, errorFile))

for(method in methods){
  for(exposure in expos){
    ParallelLogger::logInfo(sprintf("\n\nRunning analyses on %s for %s using %s...\n\n",
                database, exposure, method))
    allRes = better::multiGBSAnalyses(connection, 'eumaeus',
                                      database_id = database, 
                                      method = method,
                                      exposure_id = exposure,
                                      analysis_ids = c(1:12),
                                      period_ids = c(1:12),
                                      LPpath = LPpath,
                                      LPfname = LPfname,
                                      savepath = savepath,
                                      sampspath = sampspath,
                                      maxCores = 4)
    summName = sprintf('AllSummary-%s-%s-%s.rds', 
                       database, method, exposure)
    saveRDS(allRes, file.path(savepath, summName))
    
    ParallelLogger::logInfo(sprintf("\nAnalyses results saved on %s for %s using %s...\n\n\n",
                                    database, exposure, method))
  }
}


# close connection -----
DatabaseConnector::disconnect(connection)
