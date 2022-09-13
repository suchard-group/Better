# 09/13/2022
# look at some subject/case counts for SCCS results
# and see if they add up at all

## from EUMAEUS results
sql = "SELECT estimate.*
      FROM eumaeus.ESTIMATE 
      WHERE (method = 'SCCS' AND outcome_id = 443421)"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)

selectNCs <- DatabaseConnector::querySql(connection, sql)

names(selectNCs) = tolower(names(selectNCs))

## check out a few examples here

small_SCCS_example = 
selectNCs %>% filter(database_id == 'IBM_MDCR', exposure_id == 211981,
                     period_id == 12) %>%
  select(database_id:period_id, exposure_subjects:counterfactual_outcomes, rr)

write_csv(small_SCCS_example, '~/Downloads/tiny_EUMAEUS_SCCS_example.csv')


## the GBS analyses
GBS_example = read_csv('~/Documents/Research/better_gbs/Results_MDCR/estimate.csv')
(GBS_example = GBS_example %>% filter(databaseId == 'MDCR', method == 'SCCS', 
                                     exposureId == 211981, periodId == 12) %>%
  select(databaseId:periodId, exposureSubjects:counterfactualOutcomes, rr))

write_csv(GBS_example, '~/Downloads/GBS_example.csv')
