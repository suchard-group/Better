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
readr::write_csv(exposures, file.path(outputFolder, 'exposure.csv'))

priors = readRDS('./localCache/priorTable.rds')
readr::write_csv(priors, file.path(outputFolder, 'priors.csv'))

all_ipcs = readRDS('./localCache/allIPCs.rds')
names(all_ipcs) = tolower(names(all_ipcs))
readr::write_csv(all_ipcs, file.path(outputFolder, 'all_ipcs.csv'))

database = readRDS('./localCache/database.rds')
names(database) = tolower(names(database))
readr::write_csv(database, file.path(outputFolder, 'database.csv'))

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

## 04/06/2023 update: massage column names a little for better snake case transform
allSummary = allSummary %>% 
  rename(postMap = postMAP,
         adjustedPostMap = adjustedPostMAP,
         ci95_lb = CI95_lb,
         ci95_ub = CI95_ub,
         p1 = P1,
         p0 = P0,
         adjustedCi95_lb = adjustedCI95_lb,
         adjustedCi95_ub = adjustedCI95_ub)

names(allSummary) = SqlRender::camelCaseToSnakeCase(names(allSummary))

readr::write_csv(allSummary, file.path(outputFolder, 'summary.csv'))


# (3) save useful intermediate results for faster plotting and summarization on ShinyApp
all_mses = readRDS('./localCache/allMSEs-2.rds')
readr::write_csv(all_mses, file.path(outputFolder, 'mses.csv'))

all_type1s = readRDS('./localCache/all_type1s_95threshold.rds')
readr::write_csv(all_type1s, file.path(outputFolder, 'type1s.csv'))

all_powers = readRDS('./localCache/all_powers_calibrated.rds')
readr::write_csv(all_powers, file.path(outputFolder, 'powers.csv'))

all_tts = readRDS('./localCache/all_tts_sens50.rds')
#all_tts = readRDS('./localCache/all_tts_sens25.rds')
readr::write_csv(all_tts, file.path(outputFolder, 'time_to_signal.csv'))
