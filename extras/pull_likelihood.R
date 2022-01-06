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

## time period table
sql <- "SELECT * from eumaeus.TIME_PERIOD"
periodTable <- DatabaseConnector::querySql(connection = connection,
                                           sql = sql)
periodTable$PERIOD_ID
# mostly ranges from 1 to 9


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


# likelihood profile pull function
# TODO: edit this to filter on 
## database_id
## method
## analysis_id (nested within method)
## exposure_id
## outcome_id
## period_id
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
  
  # split string and convert to numbers
  points = str_split(lik$POINT,';') %>% 
    unlist() %>% as.numeric()
  values = str_split(lik$VALUE,';') %>% 
    unlist() %>% as.numeric()
  
  # create a data frame as result
  res = data.frame(point = points, value = values)
  
  # plot if...
  if(plot){
    g = ggplot(res, aes(x=point, y=value)) +
      geom_line(size = 0.8) +
      labs(x = 'parameter value', y = 'log-likelihood') +
      theme_bw(base_size=14)
    print(g)
  }
  
  # return the data frame
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
                         database_id = 'IBM_MDCD',
                         method = 'SCCS',
                         analysis_id = 1,
                         exposure_id = 21184,
                         outcome_id = 438945,
                         period_id = 3)
## Hmmmmm need to wrap the chars with additional pair of " " to make it work?!?!
sql <- SqlRender::render(sql, 
                         database_id = "'IBM_MDCD'",
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

