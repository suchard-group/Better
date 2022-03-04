# Feb 17, post process main functions

#library(stringr)
source('./extras/postProcessUtils.R')


## main function for post processing
## (1) pull results for each database and method (optional choice on exposure_ids) 
##     --> save the big summary file
## (2) also make decisions based on posterior probs --> no judgement of right/wrong yet!

## NOTE: must already have IPC table saved locally before running this!!!

## NOTE: `thresholdTableName` has to differ from default if delta1 or delta0 does differ!!

postProcess <- function(database_id, 
                        method,
                        exposure_id = NULL,
                        delta1 = c(0.80,0.90,0.95),
                        delta0 = c(0.90,0.95,0.99),
                        thresholdTableName = 'thresholdTable', 
                        resultsPath, 
                        savePath = NULL,
                        IPCpath='./localCache/',
                        maxCores = 10){
  # 0. create savePath if not exist
  if(!is.null(savePath) && !dir.exists(savePath)){
    dir.create(savePath)
  }
  
  # 1. get IPC table and set exposures 
  IPCs = readRDS(file.path(IPCpath, 'allIPCs.rds'))
  if(is.null(exposure_id)){exposure_id = sort(unique(IPCs$EXPOSURE_ID))}
  
  # 2. get threshold table
  thresholds = getThresholdTable(delta1 = delta1, delta0 = delta0, 
                                 savepath = savePath, tablename = thresholdTableName)
  
  # 3. pull needed results
  registerDoParallel(cores = 4)
  summ = pullResults(database_id = database_id, 
                     method = method, 
                     exposure_id = exposure_id,
                     resultsPath = resultsPath,
                     savePath = savePath,
                     IPCpath = IPCpath)
  stopImplicitCluster()
  
  # 4. make decisions
  groups = split(summ,
                 list(
                   summ$database_id,
                   summ$method,
                   summ$analysis_id,
                   summ$exposure_id,
                   summ$outcome_id,
                   summ$prior_id
                 ), drop=TRUE)
  
  cat('\nSummary table subsetted...\n')
  
  ngr = length(groups)
  
  cl <- makeCluster(max(4, maxCores))
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = ngr, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cat('\nStart decision making...\n')
  t1 = Sys.time()
  res <- foreach(i = 1:ngr, .combine = 'bind_rows',
                 .options.snow = opts,
                 .multicombine = TRUE,
                 .errorhandling = 'remove') %dopar%
    {
      source('./extras/postProcessUtils.R')
      #this.name = names(groups)[i]
      getAllOverallDecisions(groups[[i]], thresholdTable = thresholds) %>%
        mutate(analysis_id = groups[[i]]$analysis_id[1],
               exposure_id = groups[[i]]$exposure_id[1],
               outcome_id = groups[[i]]$outcome_id[1], 
               prior_id = groups[[i]]$prior_id[1],
               negativeControl = groups[[i]]$negativeControl[1])
    }
  cat('\nProcessing finished!\n')
  print(Sys.time()-t1)
  cat('\n')
  close(pb)
  stopCluster(cl) 
  
  # 5. reformat the result table and save it
  res = res %>% 
    mutate(database_id = database_id, method =method) %>% 
    relocate(database_id, method, analysis_id, 
             exposure_id, outcome_id, prior_id, threshold_id)
  
  if(!is.null(savePath)){
    if(length(exposure_id) == 10){
      fname = sprintf('AllDecisions-%s-%s.rds',
                      database_id, method)
    }else{
      fname = paste0(sprintf('AllDecisions-%s-%s-',
                             database_id,
                             method),
                     paste(exposure_id, collapse = '+'),
                     '.rds')
    }
    saveRDS(res, file.path(savePath, fname))
    cat(sprintf('Decision results for database %s method %s saved at %s\n\n',
                database_id, method, file.path(savePath, fname)))
  }
}



#### RUN command to post process -------------
library(doSNOW)

## run setup
savepath = '~/Documents/Research/betterResults/summary'

database_id = 'MDCR'

if(stringr::str_starts(database_id, 'IBM')){
  resultspath = sprintf('~/Documents/Research/betterResults/betterResults-%s/',
                        unlist(stringr::str_split(database_id, '_'))[2])
}else{
  resultspath = sprintf('~/Documents/Research/betterResults/betterResults-%s/',
                        database_id)
}


for(method in c('SCCS','HistoricalComparator')){

postProcess(database_id = database_id, 
            method = method,
            resultsPath = resultspath,
            savePath = savepath,
            maxCores = 8)

}


