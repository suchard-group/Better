# 10/10/2022
# check GBS analysis estimates
# hunting for "significant" results

# path and file name
est_path = '~/Documents/Research/better_gbs_flu84/Results_MDCR/'
fname = 'estimate.csv'

fpath = file.path(est_path, fname)

# load estimates
res = read_csv(fpath)

# filter GBS and non-NA results
res = res %>% filter(outcomeId == 343) %>%
  filter(!is.na(rr), !is.na(ci95Lb), !is.na(ci95Ub)) %>%
  select(method, analysisId, periodId, exposureId, outcomeId, rr, ci95Lb, ci95Ub)

# check final results at max.period
res_final = res %>% group_by(method, analysisId, exposureId) %>%
  filter(periodId == max(periodId)) %>%
  ungroup()

# look at "significant" results with CI > 1
res_final %>% filter(ci95Lb >= 1)
