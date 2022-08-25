# Aug 25 2022: plot posterior distribution for GBS analyses
# by analysis periods

source('./extras/postProcessUtils.R')


samplepath = '~/Documents/Research/betterGBSanalysesResults/GBSsamples/'
saveSamplePath = '~/Documents/Research/betterGBSanalysesResults/SamplesDataFrame/'

db = 'CCAE'
aid = 2
me = 'HistoricalComparator'
eid = 211981

allSamps = pullPostSamples(database_id = db,
                           method = me,
                           analysis_id = aid,
                           exposure_id = eid,
                           resultsPath = samplepath,
                           savePath = saveSamplePath)


plotGBSPosteriors(allSamps, adjust = FALSE)
plotGBSPosteriors(allSamps, adjust = TRUE)
