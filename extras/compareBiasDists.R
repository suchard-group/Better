# 10/10/2022
# demo plot to show systematic error (negative control distribution)
# across databases/exposure/etc.

# 01/23/2023
# minor update on bias distribution plot
# with same x-axis scale 

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

#db = 'IBM_MDCD'
db = 'CCAE'
#me = 'SCCS'
me = 'HistoricalComparator'
aid = 6
pid = 9
eid = 211983

allExpos = readRDS('./localCache/exposures.rds')
exposures_select = c(211983, 211833, 21184, 21185, 21215)

nullDat = NULL
NCDat = NULL 

# databases = c('CCAE', 'IBM_MDCD', 'IBM_MDCR', 'OptumEhr', 'OptumDod')
# alterNames = c('CCAE', 'MDCD', 'MDCR', 'OptumEHR', 'OptumDoD')

for(ex in exposures_select){
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


## some graphing setups
xlims = c(-3, 3)
RRbreaks = c(0.1, 0.5, 1, 2, 5, 20)
xbreaks = log(RRbreaks)
xlabels = RRbreaks %>% as.character()

(
p = ggplot(nullDat, aes(x=x)) +
  geom_density(fill = 'gray80') +
  geom_point(data = NCDat, 
             mapping = aes(x=x,y=y), shape = 4) +
  scale_x_continuous(limits = xlims, 
                     breaks = xbreaks, 
                     labels = xlabels) +
  labs(x='Rate ratio estimates for negative control outcomes',
       y='') +
  scale_y_continuous(breaks = NULL) +
  theme_bw(base_size = 14) +
  facet_grid(exposure~., switch = "y") + 
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0))
)
p


# 02/13/2023: HC and SCCS overlay-----
## (1) HC
#db = 'IBM_MDCD'
db = 'CCAE'
#me = 'SCCS'
me = 'HistoricalComparator'
aid = 2
pid = 9
eid = 211983

allExpos = readRDS('./localCache/exposures.rds')
exposures_select = c(211983, 211833, 21184, 21185, 21215)

nullDat = NULL
NCDat = NULL 

# databases = c('CCAE', 'IBM_MDCD', 'IBM_MDCR', 'OptumEhr', 'OptumDod')
# alterNames = c('CCAE', 'MDCD', 'MDCR', 'OptumEHR', 'OptumDoD')

for(ex in exposures_select){
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

## (2) SCCS ----
db = 'IBM_MDCD'
db = 'CCAE'
me = 'SCCS'
#me = 'HistoricalComparator'
aid = 2
pid = 9
eid = 211983

allExpos = readRDS('./localCache/exposures.rds')
exposures_select = c(211983, 211833, 21184, 21185, 21215)

nullDat2 = NULL
NCDat2 = NULL 

# databases = c('CCAE', 'IBM_MDCD', 'IBM_MDCR', 'OptumEhr', 'OptumDod')
# alterNames = c('CCAE', 'MDCD', 'MDCR', 'OptumEHR', 'OptumDoD')

for(ex in exposures_select){
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
  nullDat2 = rbind(nullDat2,
                  data.frame(x=fittedNull$bias,
                             exposure = this.ex))
  NCDat2 = rbind(NCDat2,
                data.frame(x=NCs, y = 0.03,
                           exposure = this.ex))
}

## overlay the densities of bias distributions together----
nullDatAll = rbind(nullDat %>% mutate(method = 'HistoricalComparator'), 
                   nullDat2 %>% mutate(method = 'SCCS'))

xlims = c(-2, 2)
RRbreaks = c(0.2, 0.5, 1, 2, 5)
xbreaks = log(RRbreaks)
xlabels = RRbreaks %>% as.character()

(
  p3 = ggplot(nullDatAll, aes(x=x, fill = method)) +
    geom_density(alpha = 0.4) +
    # geom_point(data = NCDat, 
    #            mapping = aes(x=x,y=y), shape = 4) +
    scale_x_continuous(limits = xlims, 
                       breaks = xbreaks, 
                       labels = xlabels) +
    labs(x='Rate ratio estimates for negative control outcomes',
         y='',
         fill = 'Design:') +
    scale_y_continuous(breaks = NULL) +
    scale_fill_manual(values = wes_palette("Darjeeling2")[2:3],
                      labels = c('Historical Comparator', 
                                 'Self-controlled'))+
    theme_bw(base_size = 14) +
    facet_grid(exposure~., switch = "y") + 
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = 'right')
)
p3


# 03/31/2023: plot NC MLEs and estimated null distributions over time  ----
# for one exposure-database combo

db = 'CCAE'
#me = 'SCCS'
me = 'HistoricalComparator'
aid = 5
pid = 9
eid = 211983

allPeriods = 1:12

#allExpos = readRDS('./localCache/exposures.rds')
#exposures_select = c(211983, 211833, 21184, 21185, 21215)

nullDat = NULL
NCDat = NULL 

# databases = c('CCAE', 'IBM_MDCD', 'IBM_MDCR', 'OptumEhr', 'OptumDod')
# alterNames = c('CCAE', 'MDCD', 'MDCR', 'OptumEHR', 'OptumDoD')

for(pid in allPeriods){
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
  nullDat = rbind(nullDat,
                  data.frame(x=fittedNull$bias,
                             period = pid))
  NCDat = rbind(NCDat,
                data.frame(x=NCs, y = 0.03, period = pid))
}


## some graphing setups
xlims = c(-3, 3)
RRbreaks = c(0.1, 0.5, 1, 2, 5, 20)
xbreaks = log(RRbreaks)
xlabels = RRbreaks %>% as.character()

(
  p = ggplot(nullDat, aes(x=x)) +
    geom_vline(xintercept = 0, linetype = 2, 
               linewidth = 1, color = 'gray40')+
    geom_density(fill = 'gray80') +
    geom_point(data = NCDat, 
               mapping = aes(x=x,y=y), shape = 4) +
    scale_x_continuous(limits = xlims, 
                       breaks = xbreaks, 
                       labels = xlabels) +
    labs(x='Rate ratio estimates for negative control outcomes',
         y='',
         caption = 'Analysis time (month)') +
    scale_y_continuous(breaks = NULL) +
    coord_flip() +
    theme_bw(base_size = 14) +
    facet_grid(.~period, switch = "both") + 
    theme(strip.placement = "inside",
          strip.text.y.left = element_text(angle = 0),
          panel.border = element_blank(),
          strip.background = element_blank(),
          plot.caption = element_text(size=14, hjust=0.5))
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
