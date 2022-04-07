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