# Fan: 12/15/2021
# try pulling likelihood profiles
# from EUMAEUS results database

library(tidyverse)

# connection details
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("eumaeusServer"),
                 keyring::key_get("eumaeusDatabase"),
                 sep = "/"),
  user = keyring::key_get("eumaeusUser"),
  password = keyring::key_get("eumaeusPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

# try getting something from EUMAEUS database
## analysis table
sql <- "SELECT * from eumaeus.ANALYSIS"
analysisTable <- DatabaseConnector::querySql(connection = connection,
                                              sql = sql)
SCCSanalysis <- analysisTable %>% filter(METHOD == 'SCCS')
# method: "SCCS"; analysis_id ranges from 1 to 15

## database table
sql <- "SELECT * from eumaeus.DATABASE"
databaseTable <- DatabaseConnector::querySql(connection = connection,
                                             sql = sql)
databaseTable$DATABASE_ID
# [1] "CCAE"     "IBM_MDCD" "IBM_MDCR" "OptumEhr" "OptumDod"

## exposure table
sql <- "SELECT * from eumaeus.EXPOSURE"
exposureTable <- DatabaseConnector::querySql(connection = connection,
                                             sql = sql)
exposureTable$EXPOSURE_ID
# [1]  21184  21185  21214  21215 211981 211982 211983 211831 211832
# [10] 211833

## negative control outcome table
sql <- "SELECT * from eumaeus.NEGATIVE_CONTROL_OUTCOME"
negOutTable <- DatabaseConnector::querySql(connection = connection,
                                             sql = sql)
negOutTable$OUTCOME_ID
# [1]   438945   434455   316211   201612   438730   441258   432513
# 93 in total

## positive control outcome table
sql <- "SELECT * from eumaeus.POSITIVE_CONTROL_OUTCOME"
posOutTable <- DatabaseConnector::querySql(connection = connection,
                                           sql = sql)
length(posOutTable$OUTCOME_ID)
# somehow there are 921 in total
# positive control created for each exposure-(neg)outcome pair, with 3 RRs for each?
# BUT only corresponds to 80 negative controls
# Also only 6 exposure ids (not stratified by dose #)

## imputed positive control outcome table
sql <- "SELECT * from eumaeus.IMPUTED_POSITIVE_CONTROL_OUTCOME"
impPosOutTable <- DatabaseConnector::querySql(connection = connection,
                                           sql = sql)
length(impPosOutTable$OUTCOME_ID)
# 2790 in total???
# this, though, corresponds to all 93 negative controls
# with all 10 exposure ids (stratified by vaxx dose)

## time period table
sql <- "SELECT * from eumaeus.TIME_PERIOD"
periodTable <- DatabaseConnector::querySql(connection = connection,
                                           sql = sql)
periodTable$PERIOD_ID
# mostly ranges from 1 to 9

## look at ESTIMATE v.s. ESTIMATE_IMPUTED_PCS
sql <- "SELECT outcome_id from eumaeus.ESTIMATE"
estimateTable <- DatabaseConnector::querySql(connection = connection,
                                          sql = sql)
unique(estimateTable$OUTCOME_ID) # 1014 total outcomes exist in the ESTIMATE table

sql <- "SELECT outcome_id from eumaeus.ESTIMATE_IMPUTED_PCS"
estimateIPCsTable <- DatabaseConnector::querySql(connection = connection,
                                             sql = sql)
length(unique(estimateIPCsTable$OUTCOME_ID))
# 2883 in total (= 2790 imputed pos + 93 neg controls)
# ALTHOUGH, not all of them exist in the estimate table for each database-method combo

### try to see if ESTIMATE_IMPUTED_PCS table corresponds to the likelihood profiles
sql <- "SELECT * from eumaeus.ESTIMATE_IMPUTED_PCS
        WHERE database_id = 'IBM_MDCD'
        AND method = 'SCCS'
        AND analysis_id = 1
        AND exposure_id = 21184
        AND outcome_id = 438945"
exEstimateIPCsTable <- DatabaseConnector::querySql(connection = connection,
                                                   sql = sql)
# 7 rows in total, which matches the likelihood profile query results

### contrast with ESTIMATE
sql <- "SELECT * from eumaeus.ESTIMATE
        WHERE database_id = 'IBM_MDCD'
        AND method = 'SCCS'
        AND analysis_id = 1
        AND exposure_id = 21184
        AND outcome_id = 438945"
exEstimateTable <- DatabaseConnector::querySql(connection = connection,
                                                   sql = sql)
# hmm also 7 rows in total...


# sql <- "SELECT TOP 10 * FROM eumaeus.LIKELIHOOD_PROFILE"
# sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
# likTable <- DatabaseConnector::querySql(connection = connection,
#                                         sql = sql)


sql <- "SELECT * FROM eumaeus.LIKELIHOOD_PROFILE
        WHERE database_id = 'IBM_MDCD'
        AND method = 'SCCS'
        AND analysis_id = 1
        AND exposure_id = 21184
        AND outcome_id = 438945"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
likTable <- DatabaseConnector::querySql(connection = connection,sql = sql)


sql <- "SELECT outcome_id FROM eumaeus.LIKELIHOOD_PROFILE"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
likTableOutcomes <- DatabaseConnector::querySql(connection = connection,sql = sql)


# likelihood profile pull function
# params:
## database_id
## method
## analysis_id (nested within method)
## exposure_id
## outcome_id
## period_id

# need to output a double (names being the grid points) instead
# ALSO need to check for non-existent likelihood profiles
pullLikelihood <- function(database_id, analysis_id, 
                           exposure_id, outcome_id,
                           period_id,
                           method = "'SCCS'", # default to SCCS for our purposes
                           plot=FALSE){
  
  sql <- "SELECT point, value FROM eumaeus.LIKELIHOOD_PROFILE
          WHERE database_id = @database_id
          AND method = @method
          AND analysis_id = @analysis_id
          AND exposure_id = @exposure_id
          AND outcome_id = @outcome_id
          AND period_id = @period_id"
  
  ## using the updated "render" & "translate" functions instead
  ## they return chatacter strings directly; no need to extract from a list
  sql <- SqlRender::render(sql, 
                           database_id = database_id,
                           method = method,
                           analysis_id = analysis_id,
                           exposure_id = exposure_id,
                           outcome_id = outcome_id,
                           period_id = period_id)
  sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
  lik <- DatabaseConnector::querySql(connection, sql)
  
  # check to see if result is empty
  if(nrow(lik) == 0){
    res = numeric(0)
  }else{
    # split string and convert to numbers
    points = str_split(lik$POINT,';') %>% 
      unlist() %>% as.numeric()
    values = str_split(lik$VALUE,';') %>% 
      unlist() %>% as.numeric()
    
    # # create a data frame as result
    # res = data.frame(point = points, value = values)
    
    # UPDATED (for BEAST usage)
    # create a double vector as result
    res = values
    names(res) = as.character(points)
    
    # plot if...
    if(plot){
      # construct dataframe for plotting
      res.dat = data.frame(point = points, value = values)
      g = ggplot(res.dat, aes(x=point, y=value)) +
        geom_line(size = 0.8) +
        labs(x = 'parameter value', y = 'log-likelihood') +
        theme_bw(base_size=14)
      print(g)
    }
  }
  
  # return the likelihood profile as a double vector
  res
}


## try out the function
lik = pullLikelihood(database_id = "'IBM_MDCD'", analysis_id = 1,
                     exposure_id = 21184, outcome_id = 438945,
                     period_id = 3, plot=TRUE)
### (needs a bit of hacking to really work...)



#### Okay there is a small issue though......
### these two lines do NOT work
sql <- "SELECT point, value FROM eumaeus.LIKELIHOOD_PROFILE
        WHERE database_id = @database_id
        AND method = @method
        AND analysis_id = @analysis_id
        AND exposure_id = @exposure_id
        AND outcome_id = @outcome_id
        AND period_id = @period_id"
sql <- SqlRender::render(sql, 
                         database_id = "IBM_MDCD",
                         method = 'SCCS',
                         analysis_id = 1,
                         exposure_id = 21184,
                         outcome_id = 438945,
                         period_id = 3)
## Hmmmmm need to wrap the chars with additional pair of " " to make it work?!?!
sql <- SqlRender::render(sql, 
                         database_id = "\"IBM_MDCD\"",
                         method = "'SCCS'",
                         analysis_id = 1,
                         exposure_id = 21184,
                         outcome_id = 438945,
                         period_id = 3)


### BUT this one works!
sql <- "SELECT point, value FROM eumaeus.LIKELIHOOD_PROFILE
        WHERE database_id = 'IBM_MDCD'
        AND method = 'SCCS'
        AND analysis_id = 1
        AND exposure_id = 21184
        AND outcome_id = 438945
        AND period_id = 3"

sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
lik <- DatabaseConnector::querySql(connection, sql)



# finally need to disconnect
DatabaseConnector::disconnect(connection)
