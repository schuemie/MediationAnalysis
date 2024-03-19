# Implement a real-world example comparing DOACs to warfarin for the outcome of MACE using major bleeds as mediator

# DFetch cohort definitions from WebAPI ------------------------------------------------------------
library(dplyr)
ROhdsiWebApi::authorizeWebApi(baseUrl = Sys.getenv("baseUrl"),
                              authMethod = "windows")
ROhdsiWebApi::insertCohortDefinitionSetInPackage(fileName = "inst/settings/CohortsToCreate.csv",
                                                 baseUrl = Sys.getenv("baseUrl"),
                                                 insertTableSql = FALSE,
                                                 insertCohortCreationR = FALSE,
                                                 generateStats = FALSE,
                                                 packageName = "MediationAnalysis")
cohorts <- tibble(
  cohortIds = c(16329,
                16330,
                10870,
                11024,
                11051),
  cohortNames =c("DOACs", 
                 "Warfarin",
                 "Major bleeding",
                 "MACE",
                 "AMI")
)
cohortDefinitionSet <- ROhdsiWebApi::exportCohortDefinitionSet(baseUrl = Sys.getenv("baseUrl"),
                                                               cohortIds = cohorts$cohortIds)
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
  excludedCovariateConceptIds = c(16329, 16330),
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
  covariateSettings = covariateSettings
)
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




