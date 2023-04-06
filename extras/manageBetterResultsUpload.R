# Upload BETTER results to the OHDSI PostgreSQL public server

# set up connection details
# using my credentials with write permissions.... 
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("betterServer"),
                 keyring::key_get("betterDatabase"),
                 sep = "/"),
  user = keyring::key_get("betterUser"),
  password = keyring::key_get("betterPassword"))
                                            
schema <- "better_results"

# create results data models in the better_results schema
createResultsDataModel(connectionDetails, schema)
