# 04/06/2022
# post process GBS analyses

source('./extras/postProcessUtils.R')
source('./R/helperFunctions.R')


resultspath = '~/Documents/Research/better_gbs/summary/'
db = 'CCAE'
methods = c('SCCS', 'HistoricalComparator')
expos = c(211981:211983)

allRes = pullGBSResults(database_id = db,
                        methods = methods,
                        exposure_ids = expos,
                        resultsPath = resultspath)

# produce some examples first
exRes1 = allRes %>% 
  filter(method == 'HistoricalComparator',
         analysis_id == 2,
         exposure_id == 211981,
         prior_id == 1)
decision1.1 = getOverallDecisions(exRes1, d1 = 0.954, d0 = 0.9, withEstimates = TRUE)
decision1.2 = getOverallDecisions(exRes1, d1 = 0.999, d0 = 0.9, withEstimates = TRUE)

exRes2 = allRes %>% 
  filter(method == 'SCCS',
         analysis_id == 2,
         exposure_id == 211981,
         prior_id == 1)
decision2.1 = getOverallDecisions(exRes1, d1 = 0.977, d0 = 0.9, withEstimates = TRUE)
decision2.2 = getOverallDecisions(exRes1, d1 = 0.994, d0 = 0.9, withEstimates = TRUE)
