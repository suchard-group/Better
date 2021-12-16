# Fan: 12/15/2021
# try pulling likelihood profiles
# from EUMAEUS results database

library(tidyverse)

# EUMAEUS results credentials
keyring::key_set_with_value("eumaeusUser", password = "eumaeus_readonly")
keyring::key_set_with_value("eumaeusPassword", password = "9BA3DFEE23C176")
keyring::key_set_with_value("eumaeusServer", password = "shinydb.cqnqzwtn5s1q.us-east-1.rds.amazonaws.com")
keyring::key_set_with_value("eumaeusDatabase", password = "shinydb")

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
sql <- "SELECT * from eumaeus.ANALYSIS"
analysisTable <- DatabaseConnector::querySql(connection = connection,
                                              sql = sql)


sql <- "SELECT TOP 10 * FROM eumaeus.ANALYSIS"
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
likTable <- DatabaseConnector::querySql(connection = connection,
                                        sql = sql)


# likelihood profile pull function
# TODO: edit this to filter on 
## database_id
## method
## analysis_id (nested within method)
## exposure_id
## outcome_id
## period_id
pullLikelihood <- function(target_id, comparator_id,
                           outcome_id, analysis_id, 
                           plot=FALSE){
  
  sql <- "SELECT point, value FROM legendt2dm_class_results.likelihood_profile
          WHERE target_id = @target_id
          AND comparator_id = @comparator_id
          AND outcome_id = @outcome_id
          AND analysis_id = @analysis_id"
  
  ## using the updated "render" & "translate" functions instead
  ## they return chatacter strings directly; no need to extract from a list
  sql <- SqlRender::render(sql, 
                           target_id = target_id,
                           comparator_id = comparator_id,
                           outcome_id = outcome_id,
                           analysis_id = analysis_id)
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




