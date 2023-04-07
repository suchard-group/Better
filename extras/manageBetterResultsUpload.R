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
createResultsDataModel(connectionDetails, schema)

# package up loose local csv's into a zip file for uploading
zipName <- normalizePath(file.path(exportFolder, 'results_better.zip'), mustWork = FALSE)
files <- list.files(exportFolder, pattern = ".*\\.csv$")
oldWd <- getwd()
setwd(exportFolder)
#on.exit(setwd(oldWd))
DatabaseConnector::createZipFile(zipFile = zipName, files = files)
setwd(oldWd)

# set 
#Sys.setenv(POSTGRES_PATH = "/Library/PostgreSQL/11/bin")

# Try Using legendT2dm functionality instead!
LegendT2dm::uploadResultsToDatabase(
  connectionDetails = connectionDetails,
  schema = schema,
  purgeSiteDataBeforeUploading = FALSE,
  zipFileName = c(
    zipName
  ),
  specifications = getResultsDataModelSpecifications()
)

# try uploading to ohdsi shinydb data server
uploadResultsToDatabase(connectionDetails = connectionDetails,
                        schema = schema,
                        zipFileName = zipName,
                        purgeSiteDataBeforeUploading = FALSE)

## manual checking ----
connection = DatabaseConnector::connect(connectionDetails)
sql <- "SELECT COUNT(*) FROM better_results.summary"
up_summary = DatabaseConnector::querySql(connection, sql) # summary wrong....

