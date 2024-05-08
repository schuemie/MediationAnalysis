# Code used to prepare for the real-world example, including fetching the cohorts
# from ATLAS, and running the database diagnostics to select the databases to use.

source("RealWorldExample/DatabaseDetails.R")

# Fetch cohort definitions from WebAPI ----------------------------------------
library(dplyr)
ROhdsiWebApi::authorizeWebApi(baseUrl = Sys.getenv("baseUrl"),
                              authMethod = "windows")
cohorts <- tibble(
  cohortId = c(16329,
               16586, 
               16330,
               16853,
               16991,
               2072,
               2087,
               17002),
  cohortName =c("DOACs", 
                "Rivaroxaban",
                "Warfarin",
                "Extracranial Major bleeding IP-ER",
                "Extracranial Major bleeding IP",
                "AMI",
                "Ischemic stroke",
                "Composite")
)

cohortDefinitionSet <- ROhdsiWebApi::exportCohortDefinitionSet(baseUrl = Sys.getenv("baseUrl"),
                                                               cohortIds = cohorts$cohortId)
cohortDefinitionSet <- cohortDefinitionSet %>%
  select(-"cohortName") %>%
  inner_join(cohorts, by = join_by("cohortId"))
saveRDS(cohortDefinitionSet, "RealWorldExample/CohortDefinitionSet.rds")

# Specify target-comparator-mediator-outcome combinations ----------------------
cohortDefinitionSet <- readRDS("RealWorldExample/CohortDefinitionSet.rds") 
targets <- cohortDefinitionSet %>%
  filter(cohortId %in% c(16329, 16586)) %>%
  select(targetId = "cohortId", targetName = "cohortName")
comparator <- cohortDefinitionSet %>%
  filter(cohortId %in% c(16330)) %>%
  select(comparatorId = "cohortId", comparatorName = "cohortName")
mediators <- cohortDefinitionSet %>%
  filter(cohortId %in% c(16853, 16991)) %>%
  select(mediatorId = "cohortId", mediatorName = "cohortName")
outcomes <- cohortDefinitionSet %>%
  filter(cohortId %in% c(2072, 2087, 17002)) %>%
  select(outcomeId = "cohortId", outcomeName = "cohortName")
tcmos <- targets %>%
  cross_join(comparator) %>%
  cross_join(mediators) %>%
  cross_join(outcomes)
saveRDS(tcmos, "RealWorldExample/tcmos.rds") 

# Combine HAS-BLED cohorts into cohort definition set --------------------------
hasBledCohorts <- readr::read_csv("RealWorldExample/HasBled/HasBledCohorts.csv", show_col_types = FALSE)
hasBledCohorts$json <- ""
hasBledCohorts$sql <- ""
for (i in seq_len(nrow(hasBledCohorts))) {
  hasBledCohorts$json[i] <- SqlRender::readSql(file.path("RealWorldExample/HasBled", sprintf("%s.json", hasBledCohorts$cohortName[i])))
  hasBledCohorts$sql[i] <- SqlRender::readSql(file.path("RealWorldExample/HasBled", sprintf("%s.sql", hasBledCohorts$cohortName[i])))
}
saveRDS(hasBledCohorts, "RealWorldExample/HasBledCohortDefinitionSet.rds")

# Database diagnostics ---------------------------------------------------------
library(dplyr)
cohortDefinitionSet <- readRDS("RealWorldExample/CohortDefinitionSet.rds")
ROhdsiWebApi::authorizeWebApi(baseUrl = Sys.getenv("baseUrl"),
                              authMethod = "windows")
targetCohortDefinition <- ROhdsiWebApi::getCohortDefinition(
  cohortId = 16586, 
  baseUrl = Sys.getenv("baseUrl")
)
targetConceptSet <- targetCohortDefinition$expression$ConceptSets[[2]]
targetConceptSet$name
targetConceptIds <- ROhdsiWebApi::resolveConceptSet(
  conceptSetDefinition = targetConceptSet, 
  baseUrl = Sys.getenv("baseUrl")
)
comparatorCohortDefinition <- ROhdsiWebApi::getCohortDefinition(
  cohortId = 16330, 
  baseUrl = Sys.getenv("baseUrl")
)
comparatorConceptSet <- comparatorCohortDefinition$expression$ConceptSets[[1]]
comparatorConceptSet$name
comparatorConceptIds <- ROhdsiWebApi::resolveConceptSet(
  conceptSetDefinition = comparatorConceptSet, 
  baseUrl = Sys.getenv("baseUrl")
)
outcomeCohortDefinition <- ROhdsiWebApi::getCohortDefinition(
  cohortId = 2087, 
  baseUrl = Sys.getenv("baseUrl")
)
outcomeConceptSet <- outcomeCohortDefinition$expression$ConceptSets[[2]]
outcomeConceptSet$name
outcomeConceptIds <- ROhdsiWebApi::resolveConceptSet(
  conceptSetDefinition = outcomeConceptSet, 
  baseUrl = Sys.getenv("baseUrl"))
indicationConceptSet <- targetCohortDefinition$expression$ConceptSets[[1]]
indicationConceptSet$name
indicationConceptIds <- ROhdsiWebApi::resolveConceptSet(
  conceptSetDefinition = indicationConceptSet,
  baseUrl = Sys.getenv("baseUrl")
)
analysisSettings1 <- DbDiagnostics::createDataDiagnosticsSettings(
  analysisId = 1,
  analysisName = "A1",
  minAge = NULL,
  maxAge = NULL,
  genderConceptIds =  NULL,
  raceConceptIds = NULL,
  ethnicityConceptIds = NULL,
  studyStartDate = NULL,
  studyEndDate = NULL,
  requiredDurationDays = 365,
  requiredDomains = c("condition","drug"),
  desiredDomains = NULL,
  requiredVisits = NULL,
  desiredVisits = NULL,
  targetName = targetConceptSet$name,
  targetConceptIds = targetConceptIds,
  comparatorName = comparatorConceptSet$name,
  comparatorConceptIds = comparatorConceptIds,
  indicationName = indicationConceptSet$name,
  indicationConceptIds = indicationConceptIds,
  outcomeName = outcomeConceptSet$name,
  outcomeConceptIds = outcomeConceptIds
)
settingsList <- list(analysisSettings1)
dbProfileConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  user = keyring::key_get("strategusUser"),
  password = keyring::key_get("strategusPassword"),
  connectionString = keyring::key_get("strategusConnectionString")
)
conn <- DatabaseConnector::connect(dbProfileConnectionDetails)
dbDiagnosticResults <- DbDiagnostics::executeDbDiagnostics(
  connectionDetails = dbProfileConnectionDetails,
  resultsDatabaseSchema = "dp_temp",
  resultsTableName = "db_profile_results",
  outputFolder = rootFolder,
  dataDiagnosticsSettingsList = settingsList
)
library(dplyr)
library(tidyr)
dbDiagnosticSummary <- DbDiagnostics::createDataDiagnosticsSummary(dbDiagnosticResults)
CohortGenerator::writeCsv(dbDiagnosticSummary, file.path(rootFolder, "data_diagnostics_summary.csv"))


# Get database descriptions ----------------------------------------------------
library(dplyr)
dbInfo <-list()
for (database in databases) {
  message(sprintf("Fetching metadata from %s", database$databaseId))
  connection <- DatabaseConnector::connect(database$connectionDetails)
  cdmSource <- DatabaseConnector::renderTranslateQuerySql(
    connection = connection,
    sql = "SELECT * FROM @cdm_database_schema.cdm_source;",
    cdm_database_schema = database$cdmDatabaseSchema,
    snakeCaseToCamelCase = TRUE
  )
  dbInfo[[length(dbInfo) + 1]] <- cdmSource
  DatabaseConnector::disconnect(connection)
}
dbInfo <- dbInfo %>%
  bind_rows() %>%
  select("cdmSourceName", "sourceDescription")
saveRDS(dbInfo, "RealWorldExample/DatabaseInfo.rds")
