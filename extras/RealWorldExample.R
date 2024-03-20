# Implement a real-world example comparing DOACs to warfarin for the outcome of MACE using major bleeds as mediator

# Settings for CCAE
folder <- "d:/mediatorAnalysis/CCAE"
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = keyring::key_get("redShiftConnectionStringOhdaCcae"),
                                                                user = keyring::key_get("redShiftUserName"),
                                                                password = keyring::key_get("redShiftPassword"))
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "cohort_mediation_ccae"
cdmDatabaseSchema <- "cdm_truven_ccae_v2756" 


# DFetch cohort definitions from WebAPI ------------------------------------------------------------
library(dplyr)
ROhdsiWebApi::authorizeWebApi(baseUrl = Sys.getenv("baseUrl"),
                              authMethod = "windows")
cohorts <- tibble(
  cohortId = c(16329,
               16330,
               10870,
               11024,
               11051),
  cohortName =c("DOACs", 
                "Warfarin",
                "Major bleeding",
                "MACE",
                "AMI")
)
cohortDefinitionSet <- ROhdsiWebApi::exportCohortDefinitionSet(baseUrl = Sys.getenv("baseUrl"),
                                                               cohortIds = cohorts$cohortId)
cohortDefinitionSet <- cohortDefinitionSet %>%
  select(-"cohortName") %>%
  inner_join(cohorts, by = join_by("cohortId"))
saveRDS(cohortDefinitionSet, "extras/CohortDefinitionSet.rds")

# Create cohorts -----------------------------------------------------------------------------------
library(CohortGenerator)
cohortDefinitionSet <- readRDS("extras/CohortDefinitionSet.rds")
negativeControlConceptIds <- c(24134, 73302, 74726, 78162, 79903, 81902, 134441, 139099, 140641, 141663, 141825, 141932, 192279, 193326, 196456, 197032, 197684, 199876, 201606, 256449, 257007, 257012, 260139, 261880, 313459, 313791, 372328, 373478, 374919, 378425, 433736, 435613, 435657, 436073, 437643, 441788, 442077, 443344, 443730, 4002650, 4007453, 4050747, 4077081, 4112853, 4150614, 4153359, 4174977, 4208390, 4242416, 4324765)
negativeControlCohorts <- tibble(
  cohortId = negativeControlConceptIds,
  cohortName = sprintf("Negative control %d", negativeControlConceptIds),
  outcomeConceptId = negativeControlConceptIds
)
cohortTableNames <- getCohortTableNames(cohortTable = cohortTable)
createCohortTables(connectionDetails = connectionDetails,
                   cohortDatabaseSchema = cohortDatabaseSchema,
                   cohortTableNames = cohortTableNames)
generateCohortSet(connectionDetails = connectionDetails,
                  cdmDatabaseSchema = cdmDatabaseSchema,
                  cohortDatabaseSchema = cohortDatabaseSchema,
                  cohortTableNames = cohortTableNames,
                  cohortDefinitionSet = cohortDefinitionSet)
generateNegativeControlOutcomeCohorts(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  negativeControlOutcomeCohortSet = negativeControlCohorts
)

# Fetch CohortMethodData object --------------------------------------------------------------------
library(CohortMethod)
covariateSettings <- createDefaultCovariateSettings(
  excludedCovariateConceptIds = c(1592988, 40228152, 40241331, 43013024, 45775372, 45892847, 1310149),
  addDescendantsToExclude = TRUE
)

cmData <- getDbCohortMethodData(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  targetId = 16329,
  comparatorId = 16330,
  outcomeIds = c(10870, 11024, 11051, 24134, 73302, 74726, 78162, 79903, 81902, 134441, 139099, 140641, 141663, 141825, 141932, 192279, 193326, 196456, 197032, 197684, 199876, 201606, 256449, 257007, 257012, 260139, 261880, 313459, 313791, 372328, 373478, 374919, 378425, 433736, 435613, 435657, 436073, 437643, 441788, 442077, 443344, 443730, 4002650, 4007453, 4050747, 4077081, 4112853, 4150614, 4153359, 4174977, 4208390, 4242416, 4324765),
  exposureDatabaseSchema = cohortDatabaseSchema,
  exposureTable = cohortTable,
  outcomeDatabaseSchema = cohortDatabaseSchema,
  outcomeTable = cohortTable,
  restrictToCommonPeriod = TRUE,
  covariateSettings = covariateSettings,
  maxCohortSize = 10000
)
dir.create(folder, recursive = TRUE)
saveCohortMethodData(cmData, file.path(folder, "cmData.zip"))

# Fit propensity and mediator risk models ----------------------------------------------------------
cmData <- loadCohortMethodData(file.path(folder, "cmData.zip"))

ps <- createPs(cmData,
               control = createControl(noiseLevel = "quiet", 
                                       cvType = "auto", 
                                       seed = 1,
                                       resetCoefficients = TRUE, 
                                       tolerance = 2e-07, 
                                       cvRepetitions = 1, 
                                       fold = 10,
                                       startingVariance = 0.01,
                                       threads = 10))
saveRDS(ps, file.path(folder, "ps.rds"))
ps <- readRDS(file.path(folder, "ps.rds"))
plotPs(ps, showEquiposeLabel = TRUE)

mrs <- createMediatorRiskScore(cohortMethodData = cmData, 
                               mediatorId = 10870,
                               control = createControl(noiseLevel = "quiet", 
                                                       cvType = "auto", 
                                                       seed = 1,
                                                       resetCoefficients = TRUE, 
                                                       tolerance = 2e-07, 
                                                       cvRepetitions = 1, 
                                                       fold = 10,
                                                       startingVariance = 0.01,
                                                       threads = 10))
saveRDS(mrs, file.path(folder, "mrs.rds"))
mrs <- readRDS(file.path(folder, "mrs.rds"))
plotMrs(mrs)
