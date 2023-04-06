# 04/05/2023
# compile results into csv files to upload to OHDSI server

library(dplyr)

outputFolder = '~/Documents/Research/betterOutput/export/'
# if(!dir.exists(outputFolder)){
#   dir.create(outputFolder)
# }

# (1) things that were in Eumaeus that I need to upload too ------
analyses = readRDS('./localCache/analyses.rds')
readr::write_csv(analyses, file.path(outputFolder, 'analysis.csv'))

exposures = readRDS('./localCache/exposures.rds')
readr::write_csv(exposures, file.path(outputFolder, 'exposures.csv'))

priors = readRDS('./localCache/priorTable.rds')
readr::write_csv(priors, file.path(outputFolder, 'priors.csv'))

all_ipcs = readRDS('./localCache/allIPCs.rds')
names(all_ipcs) = tolower(names(all_ipcs))
readr::write_csv(all_ipcs, file.path(outputFolder, 'all_ipcs.csv'))

# (2) concatenate all results summary files into one single long table and save----
summaryPath = '~/Documents/Research/betterResults/summary/'
databases = c('CCAE','MDCR','IBM_MDCD','OptumDod', 'OptumEhr')
methods = c('HistoricalComparator','SCCS')

allSummary = NULL

for(db in databases){
  for(me in methods){
    fileName = sprintf("AllSummary-%s-%s.rds", db, me)
    this.summary = readRDS(file.path(summaryPath, fileName))
    allSummary = bind_rows(allSummary, this.summary)
  }
}

## change `MDCR` database name to `IBM_MDCR` to be consistent with Eumaeus
allSummary = allSummary %>% 
  mutate(database_id = if_else(database_id == 'MDCR', 'IBM_MDCR', database_id))

readr::write_csv(allSummary, file.path(outputFolder, 'all_summary.csv'))


# (3) save useful intermediate results for faster plotting and summarization on ShinyApp
all_mses = readRDS('./localCache/allMSEs-2.rds')
readr::write_csv(all_mses, file.path(outputFolder, 'all_mses.csv'))

