# pre-run on Hoffman2
# to save local results

## setup
source('./extras/bayesianWithLikelihoodProfiles.R')

## on cluster setup------
cache_dir = '/u/scratch/f/fanbu42/LPs' # different from cache_dir in serverRun.R
output_dir = '/u/home/f/fanbu42/betterResults'
resources_dir = '/u/home/f/fanbu42/betterSaves'

## do one database for now ----------------
database_id = 'IBM_MDCD'

## get IPCs ----------
IPCs = readRDS(file.path(resources_dir, 'allIPCs.rds'))

exposures = sort(unique(IPCs$EXPOSURE_ID))
NCs = sort(unique(IPCs$NEGATIVE_CONTROL_ID))

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


## iteratively save LPs to cache
runs = c(1:300)
expo_groups = (runs-1) %/% 60 + 1
analysis_groups = runs %% 60 + 1
for(i in runs){
  aid = analysis_groups[i]
  exid = expo_groups[i]
  
  ## process the aid to get method, period_id and analysis_ids
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
    period_id = (aid + 2) %/% 3
    if(aid %% 3 == 1){
      analysis_ids = c(1:8)
    }else if(aid %% 3 == 2){
      analysis_ids = c(9:16)
    }else{
      analysis_ids = c(17:24)
    }
  }
  
  for(ex in exposures[c(exid*2-1,exid*2)]){
    fname = sprintf('LPs_%s_%s_%s_period%s.rds',
                    database_id, 
                    method,
                    ex,
                    period_id)
    fpath = file.path(cache_dir, fname)
    # if file already exists, skip
    if(file.exists(fpath)) next
    
    # otherwise, query from server and save
    LPs = getMultiLikelihoodProfiles(connection, 'eumaeus', 
                                     database_id, exposure_id = ex, 
                                     analysis_id = NULL, 
                                     period_id = period_id,
                                     outcome_ids = NCs,
                                     method = method)
    
    saveRDS(LPs, file = fpath)
    cat('LPs saved at', fpath, '.\n')
  }
  
}
