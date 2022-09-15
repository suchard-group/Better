# Aug 2022
# run GBS analyses on Hoffman2

## setup
source('./R/bayesianGBSLikelihoodProfilesMeta.R')

## null prior for NC meta analysis
nullPriorSds = c(0.5,0.5)

## on cluster setup------
cache_dir = '/u/scratch/f/fanbu42/GBSsamples'
output_dir = '/u/home/f/fanbu42/betterGBSResults'
GBSdata_dir = '/u/home/f/fanbu42/GBSdata'


## process array_id to figure out what to run in this task ---------
array_id = Sys.getenv('SGE_TASK_ID') %>% as.numeric()
# SCCS: analysis_id 1-15; period_id 1-12
# HC: analysis_id 1-24;  period_id 1-12
# Feb 14 update: only do analysis 1-12 for HC (13-24 has filtered outcomes)


# Aug 2022: 200 total runs
# 5 databases; 10 expos
# all analyses split in 4 big groups
NR = 200
db_labels = rep(c('CCAE','IBM_MDCR', 'IBM_MDCD','OptumEhr', 'OptumDod'),
                each = NR/5)
expo_labels = rep(rep(c(1:10), each = 4), 5)

if(array_id %% 4 == 1){
  method = 'SCCS'
  period_ids = c(1:6)
  analysis_ids = c(1:15)
}else if(array_id %% 4 == 2){
  method = 'SCCS'
  period_ids = c(7:12)
  analysis_ids = c(1:15)
}else if(array_id %% 4 == 3){
  method = "HistoricalComparator"
  period_ids = c(1:6)
  analysis_ids = c(1:12)
}else{
  method = "HistoricalComparator"
  period_ids = c(7:12)
  analysis_ids = c(1:12)
}

database_id = db_labels[array_id]
expoInd = expo_labels[array_id]

cat(sprintf('Task %s: run analysis for %s, %s, period %s to %s, and analysis %s to %s\n',
            array_id, database_id, 
            method, 
            min(period_ids), max(period_ids),
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
exposures = exposures[expoInd]
cat('Running analyses for exposures', paste(exposures, collapse = ', '), '\n\n')


## run the analyses corresponding to this array id---------
## go through each exposure one by one
for(expo in exposures){
  cat('\n\nRunning GBS analysis for exposure_id', expo, '...\n\n')
  multiGBSAnalysesMeta(connection,
                       'eumaeus',
                       database_id = database_id,
                       method = method,
                       exposure_id = expo,
                       analysis_ids = analysis_ids,
                       period_ids = period_ids,
                       includeSyntheticPos = FALSE,
                       IPCtable = IPCs,
                       nullPriorSds = nullPriorSds,
                       preLearnNull = TRUE,
                       negControls = NCs,
                       LPPath = GBSdata_dir,
                       savepath = output_dir,
                       sampspath = cache_dir,
                       robustMetaAnalysis = FALSE)
}


## close connection--------
DatabaseConnector::disconnect(connection)
