# July 2022
# re-run on Hoffman2
# to experiment on prior settings for negative control meta analysis 
# (on-server debug run again...)

# Aug 2022
# run with "robust" meta analysis (t-distribution for NC meta analysis)
robustFlag = TRUE # TRUE for t-model, FALSE for normal model

## setup
#source('./extras/bayesianWithLikelihoodProfiles.R')
source('./extras/bayesianWithLikelihoodProfilesMeta.R')

## process array_id to figure out what to run in this task ---------
array_id = Sys.getenv('SGE_TASK_ID') %>% as.numeric()

if(array_id %% 3 == 0){
  nullPriorSds = c(2, 0.5)
  dir_suffix = 'default3'
}else if(array_id %% 3 == 1){
  nullPriorSds = c(0.5, 0.5)
  dir_suffix = 'shrinkMu3'
}else{
  nullPriorSds = c(0.2, 0.2)
  dir_suffix = 'shrinkBoth3'
}

## directory setup------
cache_dir = sprintf('/u/scratch/f/fanbu42/samples-%s', dir_suffix)
output_dir = sprintf('/u/home/f/fanbu42/betterResults-%s', dir_suffix)

if(!dir.exists(cache_dir)){
  dir.create(cache_dir)
}
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

## database ----------------
database_id = 'CCAE'

## fix exposure to Zoster vaccine
exposures = 211981

## Aug 1 2022 debug:
## method and analysis
## only run: 
# SCCS & HC
# analysis_ids: 2,4,6,8
q = array_id %% 8
if(q %/% 4 == 1){
  me = 'SCCS'
}else{
  me = 'HistoricalComparator'
}

if(q <= 3){
  aid = (q+1) * 2
}else{
  aid = (q-3) * 2
}

## run periods 1-12
period_ids = c(1:12)


cat(sprintf('Task %s: run analysis for %s, %s, and analysis %s...\n',
            array_id, database_id, 
            me, aid))


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

cat('EUMAEUS server connection established!\n\n')

## check if big IPC table is available in cache and get it if not saved-----
IPCs = getIPCs(connection, 'eumaeus', output_dir)

## get negative controls outcomes IDs
#exposures = sort(unique(IPCs$EXPOSURE_ID))
NCs = sort(unique(IPCs$NEGATIVE_CONTROL_ID))

#exposures = exposures[c(exid*2-1,exid*2)]
cat('Running analyses for exposures', paste(exposures, collapse = ', '), '\n\n')

## run the analyses corresponding to this array id---------
## go through each exposure one by one
for(expo in exposures){
  cat('\n\nRunning analysis for exposure_id', expo, '...\n\n')
  multiBayesianAnalysesMeta(connection,
                            'eumaeus',
                            database_id = database_id,
                            method = me,
                            exposure_id = expo,
                            analysis_ids = aid,
                            period_ids = period_ids,
                            includeSyntheticPos = FALSE,
                            IPCtable = IPCs,
                            preLearnNull = FALSE,
                            negControls = NCs,
                            savepath = output_dir,
                            sampspath = cache_dir,
                            removeTempSummary = FALSE,
                            nullPriorSds = nullPriorSds,
                            minNCs = 5,
                            robustMetaAnalysis = robustFlag)
}


## close connection--------
DatabaseConnector::disconnect(connection)
