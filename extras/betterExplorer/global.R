# NOTE: probably need to replace "keyring::key_get()" with "Sys.getenv()" for deployment!!!

source("dataPulls.R")

connectionPool <- pool::dbPool(drv = DatabaseConnector::DatabaseConnectorDriver(),
                               dbms = "postgresql",
                               server = paste(keyring::key_get("eumaeusServer"),
                                              keyring::key_get("eumaeusDatabase"),
                                              sep = "/"),
                               user = keyring::key_get("eumaeusUser"),
                               password = keyring::key_get("eumaeusPassword"))

onStop(function() {
  if (DBI::dbIsValid(connectionPool)) {
    writeLines("Closing connection pool")
    pool::poolClose(connectionPool)
  }
})

schema <- "eumaeus"
#Sys.getenv("eumaeusdbSchema")

analysis <- loadEntireTable(connectionPool, schema, "analysis")
database <- loadEntireTable(connectionPool, schema, "database")
exposure <- loadEntireTable(connectionPool, schema, "exposure")
negativeControlOutcome <- loadEntireTable(connectionPool, schema, "negative_control_outcome")
positiveControlOutcome <- loadEntireTable(connectionPool, schema, "positive_control_outcome")
timePeriod <- loadEntireTable(connectionPool, schema, "time_period")
databaseCharacterization <- loadEntireTable(connectionPool, schema, "database_characterization")
vaccinations <- getVaccinations(connectionPool, schema)

                                              

# subset <- getEstimates(connection = connectionPool,
#                        schema = schema,
#                        databaseId = "IBM_MDCR",
#                        exposureId = 21184,
#                        timeAtRisk = "1-28")

trueRrs <- c("Any", 1, unique(positiveControlOutcome$effectSize))
#timeAtRisks <- unique(analysis$timeAtRisk)

# calibrationInfoHtml <- readChar("calibration.html", file.info("calibration.html")$size)
# vaccineInfoHtml <- readChar("vaccine.html", file.info("vaccine.html")$size)
# databaseInfoHtml <- readChar("databases.html", file.info("databases.html")$size)
# timeAtRiskInfoHtml <- readChar("timeAtRisk.html", file.info("timeAtRisk.html")$size)
# trueRrInfoHtml <- readChar("trueRr.html", file.info("trueRr.html")$size)
# methodsInfoHtml <- readChar("methods.html", file.info("methods.html")$size)
# periodInfoHtml <- readChar("period.html", file.info("period.html")$size)
# minimumOutcomesInfoHtml <- readChar("minimumOutcomes.html", file.info("minimumOutcomes.html")$size)
# metricInfoHtml <- readChar("metrics.html", file.info("metrics.html")$size)



# BETTER connections..... 
connectionPoolBetter <- pool::dbPool(drv = DatabaseConnector::DatabaseConnectorDriver(),
                                     dbms = "postgresql",
                               server = paste(keyring::key_get("betterServer"),
                                              keyring::key_get("betterDatabase"),
                                              sep = "/"),
                               user = keyring::key_get("betterUser"),
                               password = keyring::key_get("betterPassword"))


onStop(function() {
  if (DBI::dbIsValid(connectionPoolBetter)) {
    writeLines("Closing connection pool")
    pool::poolClose(connectionPoolBetter)
  }
})

schema <- "better_results"

mses = loadEntireTable(connectionPoolBetter, schema, "mses")
# all_ipcs = loadEntireTable(connectionPoolBetter, schema, "all_ipcs") why is this table not uploaded?
priors = loadEntireTable(connectionPoolBetter, schema, "priors")

