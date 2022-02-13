# Jan 2022
# to run on Hoffman2

# (Feb 9: re-run aborted jobs, set removeTempSummary = FALSE)

## setup
source('./extras/bayesianWithLikelihoodProfiles.R')

## on cluster setup------
cache_dir = '/u/scratch/f/fanbu42/samples'
output_dir = '/u/home/f/fanbu42/betterResults'

## do one database for now ----------------
database_id = 'IBM_MDCD'

## try CCAE, Feb 12
database_id = 'CCAE'

## process array_id to figure out what to run in this task ---------
array_id = Sys.getenv('SGE_TASK_ID') %>% as.numeric()
# SCCS: analysis_id 1-15; period_id 1-12
# HC: analysis_id 1-24;  period_id 1-12

### Feb 9 update: split it further up to smaller chunks
### with 300 runs

NR = 300
runs = c(1:NR)
expo_groups = (runs-1) %/% 60 + 1
analysis_groups = rep(1:60, NR%/%60)

exid = expo_groups[array_id]
aid = analysis_groups[array_id]

if(aid <= 24){
  ## first 24: SCCS 
  method = 'SCCS'
  period_id = (aid + 1) %/% 2
  if(aid %% 2 == 1){
    # odd number: analysis 1-7
    analysis_ids = c(1:7)
  }else{
    analysis_ids = c(8:15)
  }
}else{
  ## last 36: HC
  method = "HistoricalComparator"
  period_id = (aid-24 + 2) %/% 3
  #cat('aid =', aid, '; period_id =', period_id, '\n')
  if(aid %% 3 == 1){
    analysis_ids = c(1:8)
  }else if(aid %% 3 == 2){
    analysis_ids = c(9:16)
  }else{
    analysis_ids = c(17:24)
  }
}



cat(sprintf('Task %s: run analysis for %s, %s, period %s, and analysis %s to %s\n',
            array_id, database_id, 
            method, period_id, 
            min(analysis_ids), max(analysis_ids)))


## connection details (without using keyring)--------
ConnectionDetails = DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste("shinydb.cqnqzwtn5s1q.us-east-1.rds.amazonaws.com",
                 "shinydb",
                 sep = "/"),
  user = "eumaeus_readonly",
  password = "9BA3DFEE23C176",
  pathToDriver = '~/')

connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

## check if big IPC table is available in cache and get it if not saved-----
IPCs = getIPCs(connection, 'eumaeus', output_dir)

## get exposure IDs and negative controls outcomes IDs
exposures = sort(unique(IPCs$EXPOSURE_ID))
NCs = sort(unique(IPCs$NEGATIVE_CONTROL_ID))

## Feb 9 2022: only subset on exposures
exposures = exposures[c(exid*2-1,exid*2)]
cat('Running analyses for exposures', paste(exposures, collapse = ', '), '\n\n')

## run the analyses corresponding to this array id---------
## go through each exposure one by one
for(expo in exposures){
  cat('\n\nRunning analysis for exposure_id', expo, '...\n\n')
  multiBayesianAnalyses(connection,
                        'eumaeus',
                        database_id = database_id,
                        method = method,
                        exposure_id = expo,
                        analysis_ids = analysis_ids,
                        period_ids = period_id,
                        includeSyntheticPos = FALSE,
                        IPCtable = IPCs,
                        preLearnNull = FALSE,
                        negControls = NCs,
                        savepath = output_dir,
                        sampspath = cache_dir,
                        removeTempSummary = FALSE)
}


## close connection--------
DatabaseConnector::disconnect(connection)
