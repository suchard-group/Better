# specify results schema data models
# for BETTER

exportFolder = '~/Documents/Research/betterOutput/export/'

# Create initial data model CSV based on CSV files ---------------------------------------------

determineSqlDataType <- function(value) {
  if (is(value, "Date")) {
    return("DATE")
  } else if (is(value, "POSIXct") || is(value, "POSIXt")) {
    return("TIMESTAMP")
  } else if (is(value, "numeric")) {
    if (all(is.na(value) || value == as.integer(value))) {
      if (max(value, na.rm = TRUE) < .Machine$integer.max && min(value, na.rm = TRUE) > -.Machine$integer.max) {
        return("INTEGER")
      } else {
        return("BIGINT")
      }
    } else {
      return("NUMERIC")
    }
  } else if (is(value, "integer")) {
    return("INTEGER")
  } else if (is(value, "character")) {
    fieldCharLength <- max(nchar(value))
    if (is.na(fieldCharLength)) {
      fieldCharLength = 9999
    }
    if (fieldCharLength <= 1) {
      return("VARCHAR(1)")
    } else if (fieldCharLength <= 255) {
      return("VARCHAR(255)")
    } else {
      return("TEXT")
    }
  } else if (is(value, "logical")) {
    return("VARCHAR(1)")
  } else {
    stop("Unkown type:", class(value))
  }
}

analyzeTable <- function(file) {
  tableName <- gsub(".csv", "", file)
  writeLines(sprintf("Analysing table '%s'", tableName))
  data <- readr::read_csv(file.path(exportFolder, file), col_types = readr::cols(), guess_max = 1e5, n_max = 1e5)
  types <- sapply(data, determineSqlDataType)
  return(dplyr::tibble(tableName = tableName,
                       fieldName = colnames(data),
                       type = types))
}

# Code to create database specs (csv)
files <- list.files(exportFolder, "*.csv")
specs <- purrr::map_dfr(files, analyzeTable)

## more customization on the specs
specs$isRequired <- "Yes"
## set up primary key
specs$primaryKey <- "No"
specs$primaryKey[grepl("_id$", specs$fieldName)] <- "Yes"
specs$primaryKey[grepl("(method)|(approach)|(effect_size)", specs$fieldName)] <- "Yes"
specs$primaryKey[specs$tableName == 'time_to_signal' & specs$fieldName == "sensitivity"] = "Yes"
## set `isRequired` for non-primaryKey fields
specs$isRequired[specs$primaryKey != 'Yes'] <- "No"
## make sure `effect_size` is type NUMERIC not INTEGER
specs$type[specs$fieldName == 'effect_size'] <- "NUMERIC"
# # 04/06/2023: add `optional` column so can use LegendT2dm upload functionality
# specs$optional <- "No"
# specs$optional <- "Yes"
# specs$optional[specs$primaryKey == "Yes"] <- "No" # primary key is NOT optional
specs$emptyIsNa <- "Yes"
# specs$fieldName <- gsub("\\.", "_", gsub("-", "min", specs$fieldName))
readr::write_csv(specs, "inst/settings/betterResultsModelSpecs.csv")
# 
# # use LegendT2dm functionality to create sql file for results data model
LegendT2dm::createDataModelSqlFile(specifications = readr::read_csv("inst/settings/betterResultsModelSpecs.csv"),
                                   fileName = "inst/sql/postgresql/CreateResultsTablesBetter.sql")
# 
# # OLD Eumaeus code below --------
# # Generate DDL from specs ----------------------------------------------------------------
# library(dplyr)
# specifications <- readr::read_csv("inst/settings/betterResultsModelSpecs.csv")
# 
# tableNames <- specifications$tableName %>% unique()
# script <- c()
# script <- c(script, "-- Drop old tables if exist")
# script <- c(script, "")
# for (tableName in tableNames) {
#   script <- c(script, paste0("DROP TABLE IF EXISTS ", tableName, ";"))
# }
# script <- c(script, "")
# script <- c(script, "")
# script <- c(script, "-- Create tables")
# for (tableName in tableNames) {
#   script <- c(script, "")
#   script <- c(script, paste("--Table", tableName))
#   script <- c(script, "")
#   table <- specifications %>%
#     dplyr::filter(.data$tableName == !!tableName)
#   
#   script <- c(script, paste0("CREATE TABLE ", tableName, " ("))
#   fieldSql <- c()
#   for (fieldName in table$fieldName) {
#     field <- table %>%
#       filter(.data$fieldName == !!fieldName)
#     
#     if (field$primaryKey == "Yes") {
#       required <- " PRIMARY KEY"
#     }
#     
#     if (field$isRequired == "Yes") {
#       required <- " NOT NULL"
#     } else {
#       required = ""
#     }
#     fieldSql <- c(fieldSql, paste0("\t\t\t",
#                                    fieldName,
#                                    " ",
#                                    toupper(field$type),
#                                    required))
#   }
#   primaryKeys <- table %>%
#     filter(.data$primaryKey == "Yes") %>%
#     select(.data$fieldName) %>%
#     pull()
#   fieldSql <- c(fieldSql, paste0("\t\t\tPRIMARY KEY(", paste(primaryKeys, collapse = ", "), ")"))
#   script <- c(script, paste(fieldSql, collapse = ",\n"))
#   script <- c(script, ");")
# }
# if (!exists("inst/sql/postgresql")) {
#   dir.create("inst/sql/postgresql")
# }
# SqlRender::writeSql(paste(script, collapse = "\n"), "inst/sql/postgresql/CreateResultsTablesBetter.sql")


# 04/12/2023
# function to grant user access to database schema
grantPermissionOnServer <- function(connectionDetails,
                                    schema,
                                    user = "legend") {
  connection <- DatabaseConnector::connect(connectionDetails)
  sql <- paste0("grant usage on schema ", schema, " to ", user, ";")
  DatabaseConnector::executeSql(connection, sql)
  sql <- paste0("grant select on all tables in schema ", schema, " to ", user, ";")
  DatabaseConnector::executeSql(connection, sql)
  sql <- paste0("alter default privileges in schema ", schema, " grant select on tables to ", user, ";")
  DatabaseConnector::executeSql(connection, sql)
  DatabaseConnector::disconnect(connection)
}
