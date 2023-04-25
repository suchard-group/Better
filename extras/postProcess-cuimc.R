# 04/18/2023: post process results for CUIMC database

source('./extras/postProcessUtils.R')

resultspath = '~/Documents/Research/betterResults/betterResults-CUIMC/'
savepath = '~/Documents/Research/betterResults/summary/'
db = 'CUIMC'
me = 'HistoricalComparator'

exposure_ids = readRDS('./localCache/exposures.rds')$exposure_id

summ_cuimc = 
pullResults(db, me, exposure_ids, 
            resultsPath = resultspath,
            savePath = savepath)

# No SCCS results!!!
# me = 'SCCS'
# pullResults(db, me, exposure_ids, 
#             resultsPath = '~/Documents/Research/betterResults-CUIMC/',
#             savePath = savepath)

