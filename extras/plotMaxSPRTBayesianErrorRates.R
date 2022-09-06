# 07/13/2022
# a separate, cleaner script for plotting error rates
# comparing MaxSPRT and Bayesian testing results
# 07/27/2022
# add a baseline with naive bonferroni correction
# 08/03/2022
# very quickly check Type 1 error detail (with 'sqrt' transformation)
## pretty much the same thing as those under original scale
# and see what happens with Bonferroni correction if adjust on 2 years
## difference is minimal!!
# 08/11/2022
# check a different database-exposure-analysis pair

# 09/05/2022
# try a different scale for plotting 


# 1. plot Type 1 and 2 error rates by time ------

source('./extras/simpleCalibration.R')
source('./extras/frequentistDecisionComparisons.R')

# path for saving intermediate results
summarypath = '~/Documents/Research/betterResults/summary'
#samplepath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/"
cachepath = './localCache/'

# analyses results to check
db = 'CCAE'
#db = 'OptumEhr'
eid = 211981 # Zoster 1st dose
#eid = 21184 #H1N1
me = 'HistoricalComparator'
#me = 'SCCS'
aid = 4
#pid = 3 # 1: sd=10; 2: sd=1.5; 3: sd=4 
#tolerance = 0.004


# get MaxSPRT results----
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

# here: use the raw, uncalibrated results
resLst = frequentistDecisions(connection,
                              'eumaeus',
                              database_id = db,
                              method = me,
                              exposure_id = eid,
                              analysis_id = aid,
                              calibration = FALSE,
                              correct_shift = TRUE,
                              #bonferroni_adjust_factor = 2,
                              cachePath = cachepath)
maxsprt_errors = resLst$errorRate %>%
  select(period_id, y = errorRate, effect_size, stats) %>%
  mutate(approach = '1: MaxSPRT')

## add bonferroni correction result (without MaxSPRT)
bonferroni_errors = resLst$errorRate_bonferroni %>%
  select(period_id, y = errorRate, effect_size, stats) %>%
  mutate(approach = '0: Bonferroni')

# then also pull the calibrated results
resLst_cali = frequentistDecisions(connection,
                                   'eumaeus',
                                   database_id = db,
                                   method = me,
                                   exposure_id = eid,
                                   analysis_id = aid,
                                   calibration = TRUE,
                                   correct_shift = TRUE,
                                   #bonferroni_adjust_factor = 2,
                                   cachePath = cachepath)
maxsprt_errors_cali = resLst_cali$errorRate %>%
  select(period_id, y = errorRate, effect_size, stats) %>%
  mutate(approach = '1: MaxSPRT w/ calibration')

## add bonferroni correction result (without MaxSPRT)
bonferroni_errors_cali = resLst_cali$errorRate_bonferroni %>%
  select(period_id, y = errorRate, effect_size, stats) %>%
  mutate(approach = '0: Bonferroni w/ calibration')

DatabaseConnector::disconnect(connection)


# get Bayesian temporal results----
# (a) the raw, unadjusted Bayesian method
# 07/28/2022: restrict to MaxSPRT-estimable outcomes
res_raw = plotTempDelta1ByPriors(database_id = db,
                                 method = me, 
                                 analysis_id = aid,
                                 exposure_id = eid,
                                 prior_ids = c(1:3), # include all priors for easier query later
                                 alpha = 0.05,
                                 summaryPath = summarypath,
                                 cachePath = cachepath,
                                 useAdjusted = FALSE,
                                 showPlots = TRUE,
                                 stratifyByEffectSize = TRUE,
                                 calibrate = FALSE, 
                                 outcomesInEstimates = resLst$estimates)

# pick a prior as example
#pid = 2 # use SD = 1.5 results for this Bayesian example
pid = 3 # use SD = 4 result
#pid = 1 # SD = 10 results
res_Bayes_raw = res_raw %>% 
  filter(prior_id == pid) %>% 
  select(period_id, y, effect_size, stats) %>%
  mutate(approach = '2: Bayesian, unadjusted')


# (b) bias adjusted Bayesian method
res_adj = plotTempDelta1ByPriors(database_id = db,
                                 method = me, 
                                 analysis_id = aid,
                                 exposure_id = eid,
                                 prior_ids = c(1:3),
                                 summaryPath = summarypath,
                                 cachePath = cachepath,
                                 alpha = 0.04, 
                                 useAdjusted = TRUE,
                                 showPlots = TRUE,
                                 stratifyByEffectSize = TRUE,
                                 calibrate = FALSE, 
                                 outcomesInEstimates = resLst$estimates)

pid = 3 # use SD = 4 results for this Bayesian example
#pid = 1 # use SD = 10
#pid = 2 # use SD = 1.5
res_Bayes = res_adj %>% 
  filter(prior_id == pid) %>% 
  select(period_id, y, effect_size, stats) %>%
  mutate(approach = '3: Bayesian, bias adjusted')



# combine and plot-----
## add bonferroni correction plot
errors_combined = rbind(res_Bayes, 
                        res_Bayes_raw,
                        maxsprt_errors,
                        bonferroni_errors)

yinters = 0.05
type2cols = c(wes_palette("Zissou1")[3:4],wes_palette("Royal1")[4])
othercols =  wes_palette("Royal1")[2]
allCols = c(othercols, type2cols)

period_breaks = seq(from = min(errors_combined$period_id),
                    to = max(errors_combined$period_id),
                    by = 2)
period_labels = as.integer(period_breaks)

capt = '' # no caption for now...

## (i) frequentist, raw Bayes, adjusted Bayes
p = ggplot(errors_combined, 
           aes(x=period_id, y=y, color=stats))+
  geom_line(size = 1.5) +
  geom_point(size=2)+
  geom_hline(yintercept = yinters, 
             color = 'gray60', 
             size = 1, linetype=2)+
  scale_y_continuous(limits = c(0,1),
                     trans = 'sqrt'
                     )+
  scale_x_continuous(breaks = period_breaks, labels = period_labels)+
  labs(x='analysis period (months)', y='error rates', 
       caption = capt, color='Error type')+
  scale_color_manual(values = allCols) +
  facet_grid(.~approach)+
  theme_bw(base_size = 13)+
  theme(legend.position = 'bottom') # change to bottom legend...

print(p)


# (ii) another: calibrated MaxSPRT + adjusted Bayes
errors_combined = rbind(res_Bayes, 
                        maxsprt_errors_cali,
                        bonferroni_errors_cali)

p = ggplot(errors_combined, 
           aes(x=period_id, y=y, color=stats))+
  geom_line(size = 1.5) +
  geom_point(size=2)+
  geom_hline(yintercept = yinters, 
             color = 'gray60', 
             size = 1, linetype=2)+
  scale_y_continuous(limits = c(0,1),
                     trans = 'sqrt'
                     )+
  scale_x_continuous(breaks = period_breaks, labels = period_labels)+
  labs(x='analysis period (months)', y='error rates', 
       caption = capt, color='Error type')+
  scale_color_manual(values = allCols) +
  facet_grid(.~approach)+
  theme_bw(base_size = 13)+
  theme(legend.position = 'bottom') # change to bottom legend...

print(p)


## usually don't run; check that estimates are correct ---


## 09/05/2022-----
## add plots with Type 1 and Type 2 separate-----
## do it without Bonferroni yet...

errors_combined = rbind(res_Bayes, 
                        res_Bayes_raw,
                        maxsprt_errors,
                        maxsprt_errors_cali)

errors_combined$approach[errors_combined$approach == "1: MaxSPRT"] = "1a: MaxSPRT"
errors_combined$approach[errors_combined$approach == "1: MaxSPRT w/ calibration" ] = 
  "1b: MaxSPRT w/ calibration"

## split up periods into early and late
earlySplit = 4
alphaLevel = 0.3

errors_combined = errors_combined %>%
  mutate(stage = if_else(period_id <= earlySplit, alphaLevel, 1))

Type1errors = errors_combined %>% filter(stats=='type 1')
Type2errors = errors_combined %>% filter(stats!='type 1')

yinters = 0.05
period_breaks = seq(from = min(errors_combined$period_id),
                    to = max(errors_combined$period_id),
                    by = 2)
period_labels = as.integer(period_breaks)
capt = '' # no caption for now...

## (1) type 1 error plots, with transformation to highlight near 0.05 regions---

type1colors = wes_palette("GrandBudapest1")[c(1,2,4,3)]

ybreaks = c(0,0.05, 0.1, 0.25, 0.5, 0.75,1.0)

p = ggplot(Type1errors, 
           aes(x=period_id, y=y, color=approach, alpha = stage))+
  geom_line(size = 1.5) +
  geom_point(size=2)+
  geom_hline(yintercept = yinters, 
             color = 'gray60', 
             size = 1, linetype=2)+
  scale_y_continuous(limits = c(0,1),
                     breaks = ybreaks,
                     trans = 'sqrt'
  )+
  scale_x_continuous(breaks = period_breaks, labels = period_labels)+
  labs(x='analysis period (months)', y='error rates', 
       caption = capt, color='Type 1 error of:')+
  scale_color_manual(values = type1colors) +
  scale_alpha_continuous(range = c(0.2, 1), guide = 'none')+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  #facet_grid(.~approach)+
  theme_bw(base_size = 13)+
  theme(legend.position = 'bottom')# change to bottom legend...
  

print(p)


## (2) Type 2 error, with original scale so they don't look weird

type2cols = c(wes_palette("Zissou1")[3:4],wes_palette("Royal1")[4])

p = ggplot(Type2errors, 
           aes(x=period_id, y=y, color=stats, alpha = stage))+
  geom_line(size = 1.5) +
  geom_point(size=2)+
  # geom_hline(yintercept = yinters, 
  #            color = 'gray60', 
  #            size = 1, linetype=2)+
  scale_y_continuous(limits = c(0,1)#,
                     #trans = 'sqrt'
  )+
  scale_x_continuous(breaks = period_breaks, labels = period_labels)+
  labs(x='analysis period (months)', y='error rates', 
       caption = capt, color='Error type')+
  scale_color_manual(values = type2cols) +
  scale_alpha_continuous(range = c(0.2, 1), guide = 'none')+
  facet_grid(.~approach)+
  theme_bw(base_size = 13)+
  theme(legend.position = 'bottom') # change to bottom legend...

print(p)

