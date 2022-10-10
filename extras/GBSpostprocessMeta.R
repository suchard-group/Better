# Aug 2022
# GBS analyses post processing again 
# with meta analysis for negative control outcomes

source('./extras/postProcessUtils.R')
source('./R/betterHelperFunctions.R')
#source('./extras/simpleCalibration.R')

#library(wesanderson)


resultspath = '~/Documents/Research/betterGBSanalysesResults/betterGBSResults/'
summarypath = '~/Documents/Research/betterGBSanalysesResults/resultsSummary/'

## try one database
db = 'CCAE'
methods = c('SCCS', 'HistoricalComparator')
exposures = c(211981:211983)

summ = pullGBSResultsMeta(database_id = db,
                          methods = methods,
                          exposure_ids = exposures,
                          resultsPath = resultspath,
                          savePath = summarypath)


# try post processing all databases
for(db in c('IBM_MDCD', 'IBM_MDCR', 'OptumDod', 'OptumEhr')){
  summ = pullGBSResultsMeta(database_id = db,
                            methods = methods,
                            exposure_ids = exposures,
                            resultsPath = resultspath,
                            savePath = summarypath)
}

## 09/20/2022:
## try post processing the re-run GBS results on MDCR
## 10/05/2022:
## try with GBS-flu (but not fluzone, only any flu vacc)
## 10/10/2022: GBS-flu on CCAE
resultspath = './localCache/testResults/'
summarypath = './localCache/GBSsummary/'

#db = 'IBM_MDCR'
db = 'CCAE'
methods = 'SCCS'
#exposures = c(211981:211983)
exposures = c(21215)

summ = pullGBSResultsMeta(database_id = db,
                          methods = methods,
                          exposure_ids = exposures,
                          resultsPath = resultspath,
                          savePath = summarypath)
