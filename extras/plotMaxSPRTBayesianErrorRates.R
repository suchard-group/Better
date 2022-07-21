# 07/13/2022
# a separate, cleaner script for plotting error rates
# comparing MaxSPRT and Bayesian testing results


# 1. plot Type 1 and 2 error rates by time ------

source('./extras/simpleCalibration.R')
source('./extras/frequentistDecisionComparisons.R')

# path for saving intermediate results
summarypath = '~/Documents/Research/betterResults/summary'
#samplepath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/"
cachepath = './localCache/'

# analyses results to check
db = 'CCAE'
eid = 211981
#me = 'HistoricalComparator'
me = 'SCCS'
aid = 4
#pid = 3 # 1: sd=10; 2: sd=1.5; 3: sd=4 
#tolerance = 0.004




# get Bayesian temporal results first
# (a) the raw, unadjusted Bayesian method
res_raw = plotTempDelta1ByPriors(database_id = db,
                                 method = me, 
                                 analysis_id = aid,
                                 exposure_id = eid,
                                 prior_ids = c(1:3), # include all priors for easier query later
                                 alpha = 0.03,
                                 summaryPath = summarypath,
                                 cachePath = cachepath,
                                 useAdjusted = FALSE,
                                 showPlots = TRUE,
                                 stratifyByEffectSize = TRUE,
                                 calibrate = FALSE)

# pick a prior as example
#pid = 1 # use SD = 1.5 results for this Bayesian example
pid = 3 # use SD = 4 result
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
                                 calibrate = FALSE)

pid = 3 # use SD = 4 results for this Bayesian example
res_Bayes = res_adj %>% 
  filter(prior_id == pid) %>% 
  select(period_id, y, effect_size, stats) %>%
  mutate(approach = '3: Bayesian, bias adjusted')


# get MaxSPRT results
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
                              cachePath = cachepath)
maxsprt_errors = resLst$errorRate %>%
  select(period_id, y = errorRate, effect_size, stats) %>%
  mutate(approach = '1: MaxSPRT')

# then also pull the calibrated results
resLst_cali = frequentistDecisions(connection,
                              'eumaeus',
                              database_id = db,
                              method = me,
                              exposure_id = eid,
                              analysis_id = aid,
                              calibration = TRUE,
                              correct_shift = TRUE,
                              cachePath = cachepath)
maxsprt_errors_cali = resLst_cali$errorRate %>%
  select(period_id, y = errorRate, effect_size, stats) %>%
  mutate(approach = '1: MaxSPRT w/ calibration')

DatabaseConnector::disconnect(connection)


# combine and plot
errors_combined = rbind(res_Bayes, 
                        res_Bayes_raw,
                        maxsprt_errors)

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
  scale_y_continuous(limits = c(0,1))+
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
                        maxsprt_errors_cali)

p = ggplot(errors_combined, 
           aes(x=period_id, y=y, color=stats))+
  geom_line(size = 1.5) +
  geom_point(size=2)+
  geom_hline(yintercept = yinters, 
             color = 'gray60', 
             size = 1, linetype=2)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(breaks = period_breaks, labels = period_labels)+
  labs(x='analysis period (months)', y='error rates', 
       caption = capt, color='Error type')+
  scale_color_manual(values = allCols) +
  facet_grid(.~approach)+
  theme_bw(base_size = 13)+
  theme(legend.position = 'bottom') # change to bottom legend...

print(p)


## usually don't run; check that estimates are correct ---
## oh no! Martijn has removed the imputed PCs results...

