source('extras/betterResultsDataModel.R')

# Upload BETTER results to the OHDSI PostgreSQL public server
exportFolder = '~/Documents/Research/betterOutput/export/'

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
# ONLY DO THIS ONCE!!!
# createResultsDataModel(connectionDetails, schema)

# grant user read-only access to the results schema
# (for pulling results and R ShinyApp)
grantPermissionOnServer(connectionDetails, 
                        'better_results',
                        user = 'legend')
grantPermissionOnServer(connectionDetails, 
                        'better_results',
                        user = 'legendt2dm_readonly')
grantPermissionOnServer(connectionDetails, 
                        'better_results',
                        user = "eumaeus_readonly")

# package up loose local csv's into a zip file for uploading
zipName <- normalizePath(file.path(exportFolder, 'results_better.zip'), mustWork = FALSE)
files <- list.files(exportFolder, pattern = ".*\\.csv$")
oldWd <- getwd()
setwd(exportFolder)
#on.exit(setwd(oldWd))
DatabaseConnector::createZipFile(zipFile = zipName, files = files)
setwd(oldWd)

# No longer needed... (set envvar for using PostgreSQL for bulk load funtionality)
#Sys.setenv(POSTGRES_PATH = "/Library/PostgreSQL/11/bin")

# Try Using legendT2dm functionality instead!
# Not really working :(
# LegendT2dm::uploadResultsToDatabase(
#   connectionDetails = connectionDetails,
#   schema = schema,
#   purgeSiteDataBeforeUploading = FALSE,
#   zipFileName = c(
#     zipName
#   ),
#   specifications = getResultsDataModelSpecifications()
# )

# try uploading to ohdsi shinydb data server
uploadResultsToDatabase(connectionDetails = connectionDetails,
                        schema = schema,
                        zipFileName = zipName,
                        purgeSiteDataBeforeUploading = FALSE)


## manual checking ----
connection = DatabaseConnector::connect(connectionDetails)
sql <- "SELECT COUNT(*) FROM better_results.summary"
up_summary = DatabaseConnector::querySql(connection, sql) # summary wrong....

sql <- "SELECT COUNT(*) FROM better_results.powers"
up_powers = DatabaseConnector::querySql(connection, sql)
DatabaseConnector::disconnect(connection)


## checking for read access -----
ConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("betterServer"),
                 keyring::key_get("betterDatabase"),
                 sep = "/"),
  user = keyring::key_get("betterUser"),
  password = keyring::key_get("betterPassword"))

# set up the DB connection
connection = DatabaseConnector::connect(connectionDetails = ConnectionDetails)

sql = 'SELECT * from better_results.DATABASE'
databases = DatabaseConnector::querySql(connection, sql)
databases$DATABASE_ID

