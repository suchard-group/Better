# 10/10/2022
# check GBS analysis estimates
# hunting for "significant" results

# path and file name
est_path = '~/Documents/Research/better_gbs_flu84/Results_MDCR/'
est_path = '~/Documents/Research/better_gbs/Results_MDCR/'
est_path = '~/Documents/Research/better_gbs/Results_MDCR-old/'

est_path = '~/Documents/Research/better_gbs/Results_CCAE/'
est_path = '~/Documents/Research/better_gbs/Results_CCAE-old/'

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

# (1) CCAE, 43-183 days for post control
#     all flu vaccine, SCCS
# method analysisId periodId exposureId outcomeId    rr ci95Lb ci95Ub
# <chr>       <dbl>    <dbl>      <dbl>     <dbl> <dbl>  <dbl>  <dbl>
# 1 SCCS            1        9      21215       343  2.03   1.02   3.74
# 2 SCCS            2        9      21215       343  2.25   1.12   4.20
# 3 SCCS            3        9      21215       343  9.73   1.87 178.  
# 4 SCCS            7        9      21215       343  7.70   1.52 140.  

# (2) CCAE, 43-365 (73?) days for post control
# the "significant" results are all the same as in (1)...
