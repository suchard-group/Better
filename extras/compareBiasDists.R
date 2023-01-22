# 10/10/2022
# demo plot to show systematic error (negative control distribution)
# across databases/exposure/etc.

# the null fitting function
source('./extras/fitNegativeControlDistribution.R')

# EUMAEUS connection
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)


## nulls across databases ----

#db = 'IBM_MDCD'
#me = 'SCCS'
me = 'HistoricalComparator'
aid = 6
pid = 9
eid = 211983

nullDat = NULL
NCDat = NULL 

databases = c('CCAE', 'IBM_MDCD', 'IBM_MDCR', 'OptumEhr', 'OptumDod')
alterNames = c('CCAE', 'MDCD', 'MDCR', 'OptumEHR', 'OptumDoD')

for(dbId in databases){
  db = dbId
  # fit null dist.
  fittedNull = fitNegativeControlDistributionLikelihood(connection,
                                                        'eumaeus',
                                                        database_id = db,
                                                        method = me,
                                                        analysis_id = aid,
                                                        period_id = pid,
                                                        exposure_id = eid,
                                                        priorSds = c(0.5,0.5))
  # pull MLEs for the NCs
  NCs = fitNegativeControlDistribution(connection,
                                       'eumaeus',
                                       database_id = db,
                                       method = me,
                                       analysis_id = aid,
                                       period_id = pid,
                                       exposure_id = eid,
                                       returnEstimatesOnly = TRUE) %>%
    select(log_rr) %>% pull()
  
  # plot systematic error distribution
  this.db = alterNames[which(databases==dbId)]
  nullDat = rbind(nullDat,
                  data.frame(x=fittedNull$bias,
                             db = this.db))
  NCDat = rbind(NCDat,
                data.frame(x=NCs, y = 0.03,
                           db = this.db))
}


p = ggplot(nullDat, aes(x=x)) +
  geom_density(fill = 'gray80') +
  geom_point(data = NCDat, mapping = aes(x=x,y=y), shape = 4) +
  labs(x='',y='') +
  scale_y_continuous(breaks = NULL) +
  theme_bw(base_size = 14) +
  facet_grid(db~., switch = "y") + 
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0))
p

## nulls across exposures ----

db = 'IBM_MDCD'
me = 'SCCS'
#me = 'HistoricalComparator'
aid = 2
pid = 9
eid = 211983

allExpos = readRDS('./localCache/exposures.rds')

nullDat = NULL
NCDat = NULL 

# databases = c('CCAE', 'IBM_MDCD', 'IBM_MDCR', 'OptumEhr', 'OptumDod')
# alterNames = c('CCAE', 'MDCD', 'MDCR', 'OptumEHR', 'OptumDoD')

for(ex in allExpos$exposure_id){
  eid = ex
  # fit null dist.
  fittedNull = fitNegativeControlDistributionLikelihood(connection,
                                                        'eumaeus',
                                                        database_id = db,
                                                        method = me,
                                                        analysis_id = aid,
                                                        period_id = pid,
                                                        exposure_id = eid,
                                                        priorSds = c(0.5,0.5))
  # pull MLEs for the NCs
  NCs = fitNegativeControlDistribution(connection,
                                       'eumaeus',
                                       database_id = db,
                                       method = me,
                                       analysis_id = aid,
                                       period_id = pid,
                                       exposure_id = eid,
                                       returnEstimatesOnly = TRUE) %>%
    select(log_rr) %>% pull()
  
  # plot systematic error distribution
  #this.db = alterNames[which(databases==dbId)]
  this.ex = allExpos %>% filter(exposure_id == eid) %>% 
    select(base_exposure_name) %>% pull()
  nullDat = rbind(nullDat,
                  data.frame(x=fittedNull$bias,
                             exposure = this.ex))
  NCDat = rbind(NCDat,
                data.frame(x=NCs, y = 0.03,
                           exposure = this.ex))
}

(
p = ggplot(nullDat, aes(x=x)) +
  geom_density(fill = 'gray80') +
  geom_point(data = NCDat, mapping = aes(x=x,y=y), shape = 4) +
  labs(x='',y='') +
  scale_y_continuous(breaks = NULL) +
  theme_bw(base_size = 14) +
  facet_grid(exposure~., switch = "y") + 
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0))
)
p



# -----

# (p1 = ggplot(dat, aes(x=x)) +
#   geom_density(fill = 'gray80') +
#   geom_point(data = datNC, mapping = aes(x=x,y=y), shape = 4) +
#   labs(x='',y='') +
#   scale_y_continuous(breaks = NULL) +
#   theme_bw(base_size = 14)
# )
