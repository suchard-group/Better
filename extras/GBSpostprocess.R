# 04/06/2022
# post process GBS analyses

source('./extras/postProcessUtils.R')
source('./R/helperFunctions.R')
source('./extras/simpleCalibration.R')


resultspath = '~/Documents/Research/better_gbs/summary/'
db = 'CCAE'
methods = c('SCCS', 'HistoricalComparator')
expos = c(211981:211983)

allRes = pullGBSResults(database_id = db,
                        methods = methods,
                        exposure_ids = expos,
                        resultsPath = resultspath)

# produce some examples first -----
# exRes1 = allRes %>% 
#   filter(method == 'HistoricalComparator',
#          analysis_id == 2,
#          exposure_id == 211981,
#          prior_id == 1)
# decision1.1 = getOverallDecisions(exRes1, d1 = 0.954, d0 = 0.9, withEstimates = TRUE)
# decision1.2 = getOverallDecisions(exRes1, d1 = 0.999, d0 = 0.9, withEstimates = TRUE)
# 
# exRes2 = allRes %>% 
#   filter(method == 'SCCS',
#          analysis_id == 2,
#          exposure_id == 211981,
#          prior_id == 1)
# decision2.1 = getOverallDecisions(exRes1, d1 = 0.977, d0 = 0.9, withEstimates = TRUE)
# decision2.2 = getOverallDecisions(exRes1, d1 = 0.994, d0 = 0.9, withEstimates = TRUE)
## none of the results are significant!!


# function to mass-produce decisions
getAllDecisionsWithCalibration <- function(res,
                                           savePath = NULL,
                                           calibrate = TRUE,
                                           estimatesPath = NULL,
                                           alpha = 0.05,
                                           adjust = TRUE,
                                           defaults = list(d1=0.95, d0=0.9),
                                           estimateType = 'Mean',
                                           maxCores = 4){
  if(calibrate & is.null(estimatesPath)){
    stop('estaimtesPath must be provided if calibrate!!')
  }
  
  d1 = defaults$d1; d0 = defaults$d0
  
  # sub-function to process 
  decideChunk <- function(chunk){
    db = chunk$database_id[1]
    me = chunk$method[1]
    aid = chunk$analysis_id[1]
    expo = chunk$exposure_id[1]
    prid = chunk$prior_id[1]
    
    # get calibrated delta1
    if(calibrate){
      calibrateRes = calibrateByDelta1(db, me, aid, expo, prid, estimatesPath, 
                                       alpha = alpha, useAdjusted = adjust, evalType2 = FALSE)
      if(is.null(calibrateRes)){
        cat('No calibrated threshold produced!\n')
        return(NULL)
      }
      d1 = calibrateRes$calibratedDelta1
    }
    
    # decide
    chunk.decision = getOverallDecisions(chunk, d1, d0, 
                                         withEstimates = TRUE, 
                                         estimateType = estimateType)
    
    sel.decision = chunk.decision %>% select(!starts_with('adjusted')) 
    if(adjust){
      decision = chunk.decision %>% select(starts_with('adjusted'))
      names(decision) = names(sel.decision)
    }else{
      decision= sel.decision
    }
    
    decision = decision %>% 
      mutate(database_id = db,
             method = me, 
             analysis_id = aid,
             exposure_id = expo, 
             prior_id = prid,
             delta1 = d1)
    
    return(decision)
  }
  
  # run in parallel
  chunks = split(res,
                 list(
                   res$database_id,
                   res$method,
                   res$analysis_id,
                   res$exposure_id,
                   res$prior_id
                 ), drop=TRUE)
  
  cluster = ParallelLogger::makeCluster(min(maxCores, 4))
  resls = ParallelLogger::clusterApply(cluster, chunks, decideChunk,
                               stopOnError = FALSE, progressBar = TRUE)
  ParallelLogger::stopCluster(cluster)
  AllDec = bind_rows(resls) %>% 
    arrange(exposure_id, method, analysis_id) %>%
    relocate(exposure_id, method, analysis_id, prior_id, decision, timeToDecision, 
             estimate, finalEstimate, everything())
  
  # save results to file
  if(!is.null(savePath)){
    if(!dir.exists(savePath)){
      dir.create(savePath)
    }
    saveName = sprintf('%s-cablibrate%s-adjust%s-%s.csv',
                       database_id, calibrate, adjust, estimateType)
    #saveRDS(AllDec, file.path(savePath, saveName))
    write_csv(AllDec, file.path(savePath, saveName))
  }
  
  return(AllDec)
}


## try batch processing
estimatespath = '~/Documents/Research/betterResults/summary/'
savepath = '~/Documents/Research/better_gbs/test/'

AllDecAdj = getAllDecisionsWithCalibration(allRes %>% filter(!is.na(adjustedPostMean)),
                                           savePath = savepath,
                                           calibrate = TRUE,
                                           estimatesPath = estimatespath,
                                           estimateType = 'Median',
                                           maxCores = 1)


