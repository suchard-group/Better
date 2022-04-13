# 04/12/2022

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

db = 'CCAE'
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
