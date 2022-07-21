# try out the "meta analysis" fitting of negative control distributions
# start with likelihood profiles of negative controls

source('extras/fitNegativeControlDistribution.R')
source('extras/getLikelihoodProfile.R')

library(dplyr)
library(stringr)

# one example
db = 'CCAE'
eid = 211981
#me = 'HistoricalComparator'
me = 'SCCS'
aid = 2
period = 12

IPCs = readRDS('localCache/allIPCs.rds')
NCs = unique(IPCs$NEGATIVE_CONTROL_ID)

## database connection
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

# 1. get the regular fit of negative control distr.----
null1 = fitNegativeControlDistribution(connection = connection,
                                       schema = 'eumaeus',
                                       database_id = db,
                                       method = me, 
                                       exposure_id = eid,
                                       analysis_id = 2,
                                       period_id = period)
### look at distribution for the mean bias and sd
hist(null1$mean)
hist(null1$sd)

# 2 use Bayesian meta analysis to fit NC distr. from likelihood profiles-----
LPs = getMultiLikelihoodProfiles(connection = connection,
                                 schema = 'eumaeus',
                                 database_id = db,
                                 method = me, 
                                 exposure_id = eid,
                                 analysis_id = 2,
                                 period_id = period,
                                 outcome_ids = NCs)
names(LPs) = tolower(names(LPs))
LPs = LPs %>% filter(analysis_id == aid, period_id == period, 
                     outcome_id %in% NCs)
# nrow(LPs)
LPs$index = c(1:nrow(LPs))
LPlist = split(LPs, LPs$index)
for(i in 1:length(LPlist)){
  points = str_split(LPlist[[i]]$point,';') %>% 
    unlist() %>% as.numeric()
  values = str_split(LPlist[[i]]$value,';') %>% 
    unlist() %>% as.numeric()
  LPlist[[i]] = data.frame(point = points, value = values)
}

## (1) with default prior for tau
null2a = EvidenceSynthesis::computeBayesianMetaAnalysis(data = LPlist,
                                                       priorSd = c(2,0.5))
## check summary
null2a 
null2traces = attr(null2a, "traces")
## plot global mean and sd posteriors
hist(null2traces[,1])
hist(null2traces[,2])


## (2) try using more conservative prior for tau
null2b = EvidenceSynthesis::computeBayesianMetaAnalysis(data = LPlist,
                                                        priorSd = c(2,0.1))
## check summary
null2b 
null2btraces = attr(null2b, "traces")
## plot global mean and sd posteriors
hist(null2btraces[,1])
hist(null2btraces[,2])

# there isn't much difference w.r.t. posterior of global bias mu even if we change prior for tau...
