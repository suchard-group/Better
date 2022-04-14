# 04/12/2022

library(dplyr)
library(ggplot2)
library(wesanderson)

source('extras/getLikelihoodProfile.R')

# check Bayesian analyses (with LPs) results with very extreme NC estimates
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

db = 'OptumEhr'
method = 'SCCS'
eid = 211981
aid = 2
pid = 3 # mean = 0, sd=4

outcomes = c(23731,196347, 196625, 433716, 440367)


LPs = getMultiLikelihoodProfiles(connection, 'eumaeus', database_id = db, 
                                 exposure_id = eid, analysis_id = aid, 
                                 period_id = 12, outcome_ids = outcomes, 
                                 method = method)

for(o in outcomes){
  selectLikelihoodProfileEntry(LPs, db, method, eid, 12, o, aid, plot=TRUE)
}

# They are all monotone with large density values on very negative numbers 
# (note: point only ranges from -2.3xxx to 2.3xxx)
# SO without very strong priors, estimate is going to be VERY negative...


# close connection
DatabaseConnector::disconnect(connection)


# make contrast plots comparing Bayesian estimates and EUMAEUS estimates----
source('extras/simpleCalibration.R')

summarypath = '~/Documents/Research/betterResults/summary'
samplepath = "/Volumes/WD-Drive/betterResults-likelihoodProfiles/"
cachepath = './localCache/'

db = 'CCAE'
eid = 211981
aid = 2
pid = 1 # 1: sd=10; 2: sd=1.5; 3: sd=4 
tolerance = 0.004

## Bayesian estimates
HCbiases1 = getBiases(database_id = db, method = 'HistoricalComparator',
                      exposure_id = eid, analysis_id = aid, prior_id = pid,
                      resPath = summarypath, estimateType = 'MAP')

SCCSbiases1 = getBiases(database_id = db, method = 'SCCS',
                        exposure_id = eid, analysis_id = aid, prior_id = pid,
                        resPath = summarypath, estimateType = 'MAP')

dat1 = bind_rows(as.data.frame(HCbiases1), as.data.frame(SCCSbiases1)) %>% 
  select(outcome_id, estimates, method, mean, sd, num) %>%
  mutate(label = sprintf('Num. of NCs: %s\nMean: %.3f\nStd: %.3f',
                         num, mean, sd),
         method = if_else(method == 'HistoricalComparator', 
                          'Historical Comparator', 'SCCS'))


## EUMAEUS estimates
HCbiases2 = getBiases(database_id = db, method = 'HistoricalComparator',
                      exposure_id = eid, analysis_id = aid, prior_id = pid,
                      resPath = cachepath, source='EUMAEUS')

SCCSbiases2 = getBiases(database_id = db, method = 'SCCS',
                        exposure_id = eid, analysis_id = aid, prior_id = pid,
                        resPath = cachepath, source='EUMAEUS')

dat2 = bind_rows(as.data.frame(HCbiases2), as.data.frame(SCCSbiases2)) %>% 
  select(outcome_id, estimates, method, mean, sd, num) %>%
  mutate(label = sprintf('Num. of NCs: %s\nMean: %.3f\nStd: %.3f',
                         num, mean, sd),
         method = if_else(method == 'HistoricalComparator', 
                          'Historical Comparator', 'SCCS'))

# outcomes not included in EUMAEUS
extreme_outcomes = unique(dat1$outcome_id)[!unique(dat1$outcome_id) %in% unique(dat2$outcome_id)]


dat1 = dat1 %>% 
  mutate(commentary = if_else(method == 'SCCS', 
                              sprintf('%s more outcomes\nwith zero counts',length(extreme_outcomes)),
                              ''))

xrange = range(c(dat1$estimates, dat2$estimates))

ggplot(dat1, aes(x=estimates)) +
  geom_density(aes(fill = method), alpha = 0.8) +
  facet_grid(method~.)+
  geom_text(x=0.8, y = 0.6, aes(label = label), 
            fontface = 'plain', hjust = 0) +
  geom_text(x=-3.3, y=0.6, aes(label=commentary),
            fontface = 'plain', hjust = 0) +
  scale_fill_manual(values = wes_palette("Chevalier1")[2:3]) +
  scale_x_continuous(limits = xrange) +
  labs(x='Bayesian estimates for negative control outcomes',
       fill = '', y = '',
       caption = 'Bayesian inference using likelihood profiles') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom')
ggplot(dat2, aes(x=estimates)) +
  geom_density(aes(fill = method), alpha = 0.8) +
  facet_grid(method~.)+
  geom_text(x=0.8, y = 0.6, aes(label = label), 
            fontface = 'plain', hjust = 0) +
  scale_fill_manual(values = wes_palette("Chevalier1")[2:3]) +
  scale_x_continuous(limits = xrange) +
  labs(x='Estimates for negative control outcomes',
       fill = '', y = '', 
       caption = 'EUMAEUS results') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom')

## also pull up the extreme outcome estimates
dat1 %>% filter(outcome_id %in% extreme_outcomes) %>% 
  select(outcome_id, estimates)


