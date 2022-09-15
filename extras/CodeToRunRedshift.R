# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of Better
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library(better)

options(andromedaTempFolder = "E:/andromedaTemp")
options(sqlRenderTempEmulationSchema = NULL)

maxCores <- 4

# For bulk uploading synthetic outcomes:
# Sys.setenv("AWS_OBJECT_KEY" = "bulk")
# Sys.setenv("AWS_ACCESS_KEY_ID" = keyring::key_get("bulkUploadS3Key"))
# Sys.setenv("AWS_SECRET_ACCESS_KEY" = keyring::key_get("bulkUploadS3Secret"))
# Sys.setenv("AWS_BUCKET_NAME" = keyring::key_get("bulkUploadS3Bucket"))
# Sys.setenv("AWS_DEFAULT_REGION" = "us-east-1")
# Sys.setenv("AWS_SSE_TYPE" = "AES256")
# Sys.setenv("DATABASE_CONNECTOR_BULK_UPLOAD" = TRUE)

# specify where the Drivers are
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER='D:/Drivers')

## 1a. Optum DoD------------
cdmDatabaseSchema <- "cdm_optum_extended_dod_v1825"
serverSuffix <-"optum_extended_dod"
cohortDatabaseSchema <- "scratch_fbu2"
databaseId <- "OptumDod"
databaseName <- "Optum Clinformatics Extended Data Mart - Date of Death (DOD)"
databaseDescription <- "Optum Clinformatics Extended DataMart is an adjudicated US administrative health claims database for members of private health insurance, who are fully insured in commercial plans or in administrative services only (ASOs), Legacy Medicare Choice Lives (prior to January 2006), and Medicare Advantage (Medicare Advantage Prescription Drug coverage starting January 2006).  The population is primarily representative of commercial claims patients (0-65 years old) with some Medicare (65+ years old) however ages are capped at 90 years.  It includes data captured from administrative claims processed from inpatient and outpatient medical services and prescriptions as dispensed, as well as results for outpatient lab tests processed by large national lab vendors who participate in data exchange with Optum.  This dataset also provides date of death (month and year only) for members with both medical and pharmacy coverage from the Social Security Death Master File (however after 2011 reporting frequency changed due to changes in reporting requirements) and location information for patients is at the US state level."
tablePrefix <- "legend_monotherapy_OptumDoD"
outputFolder <- "E:/betterGBS_OptumDod2" # DONE # changed the save directory

## 1b. Optum EHR ---------------
# cdmDatabaseSchema <- "cdm_optum_ehr_v1821"
# serverSuffix <- "optum_ehr"
# cohortDatabaseSchema <- "scratch_fbu2"
# databaseId <- "OptumEHR"
# databaseName <- "Optum© de-identified Electronic Health Record Dataset"
# databaseDescription <- "Optum© de-identified Electronic Health Record Dataset represents Humedica’s Electronic Health Record data a medical records database. The medical record data includes clinical information, inclusive of prescriptions as prescribed and administered, lab results, vital signs, body measurements, diagnoses, procedures, and information derived from clinical Notes using Natural Language Processing (NLP)."
# tablePrefix <- "legend_monotherapy_ehr"
# outputFolder <- "E:/betterGBS_OptumEhr2" # DONE

## April 11: try to replicate some EUMAEUS results and see if fixed results are different...
## 2. IBM MDCD ------------------
# #cdmDatabaseSchema <- "cdm_truven_mdcd_v1476" # <- EUMAEUS data version
# cdmDatabaseSchema <- "cdm_truven_mdcd_v1714" 
# serverSuffix <- "truven_mdcd"
# cohortDatabaseSchema <- "scratch_fbu2"
# databaseId<- "MDCD"
# databaseName <- "IBM Health MarketScan® Multi-State Medicaid Database"
# databaseDescription <- "IBM MarketScan® Multi-State Medicaid Database (MDCD) adjudicated US health insurance claims for Medicaid enrollees from multiple states and includes hospital discharge diagnoses, outpatient diagnoses and procedures, and outpatient pharmacy claims as well as ethnicity and Medicare eligibility. Members maintain their same identifier even if they leave the system for a brief period however the dataset lacks lab data."
# tablePrefix <- "legend_monotherapy_mdcd"
# outputFolder <- "E:/betterGBS_mdcd2" # DONE
# #outputFolder <- "E:/betterTest_mdcd" # DONE

## 3. IBM CCAE ------------
# cdmDatabaseSchema <- "cdm_truven_ccae_v1709" #"cdm_idm_ccae_seta"
# serverSuffix <- "truven_ccae" # "ibm"
# cohortDatabaseSchema <- "scratch_fbu2"
# databaseId<- "CCAE"
# databaseName <- "IBM Health MarketScan® Commercial Claims and Encounters"
# databaseDescription <- "IBM MarketScan® Commercial Claims and Encounters (CCAE) adjudicated US health insurance claims for Medicaid enrollees from multiple states and includes hospital discharge diagnoses, outpatient diagnoses and procedures, and outpatient pharmacy claims as well as ethnicity and Medicare eligibility. Members maintain their same identifier even if they leave the system for a brief period however the dataset lacks lab data."
# tablePrefix <- "legend_monotherapy_ccae"
# outputFolder <- "E:/betterGBS_ccae2" # DONE

## 4. IBM MDCR --------------
# cdmDatabaseSchema <- "cdm_truven_mdcr_v1838"
# serverSuffix <- "truven_mdcr"
# cohortDatabaseSchema <- "scratch_fbu2"
# databaseId<- "MDCR"
# databaseName <- "IBM Health MarketScan Medicare Supplemental and Coordination of Benefits Database"
# databaseDescription <- "IBM Health MarketScan® Medicare Supplemental and Coordination of Benefits Database (MDCR) represents health services of retirees in the United States with primary or Medicare supplemental coverage through privately insured fee-for-service, point-of-service, or capitated health plans. These data include adjudicated health insurance claims (e.g. inpatient, outpatient, and outpatient pharmacy). Additionally, it captures laboratory tests for a subset of the covered lives."
# tablePrefix <- "better_gbs_mdcr"
# outputFolder <- "E:/betterGBS_mdcr2" # DONE

## fill out connection details ------------
conn <- DatabaseConnector::createConnectionDetails(
  dbms = "redshift",
  server = paste0(keyring::key_get("epi_server"), "/", !!serverSuffix),
  port = 5439,
  user = keyring::key_get("redshiftUser"),
  password = keyring::key_get("redshiftPassword"),
  extraSettings = "ssl=true&sslfactory=com.amazon.redshift.ssl.NonValidatingFactory",
  pathToDriver = 'D:/Drivers')


#oracleTempSchema <- NULL
cohortTable = 'cohort_fbu2'

# Aug 2022: re-run to compute critical values for GBS

execute(connectionDetails = conn,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cohortDatabaseSchema = cohortDatabaseSchema,
        cohortTable = cohortTable,
        databaseId = databaseId,
        databaseName = databaseName,
        databaseDescription = databaseDescription,
        outputFolder = outputFolder,
        maxCores = maxCores,
        exposureIds = getExposuresOfInterest()$exposureId,
        verifyDependencies = FALSE,
        createCohorts = FALSE,
        createAllControls = FALSE,
        # synthesizePositiveControls = F,
        # runCohortMethod = F,
        runSccs = FALSE,
        # runCaseControl = F,
        runHistoricalComparator = FALSE,
        generateDiagnostics = FALSE,
        computeCriticalValues = TRUE,
        createDbCharacterization =  FALSE,
        exportResults = TRUE)



# uploadResults(outputFolder = outputFolder,
#               privateKeyFileName = "c:/home/keyfiles/study-data-site-covid19.dat",
#               userName = "study-data-site-covid19")

# JnJ specific code to store database version:
# source("extras/GetDatabaseVersion.R")
# version <- getDatabaseVersion(connectionDetails, cdmDatabaseSchema)
# readr::write_csv(version, file.path(outputFolder, "version.csv"))
