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


# Fetch cohort definitions from WebAPI ----------------------------------------
library(dplyr)
ROhdsiWebApi::authorizeWebApi(baseUrl = Sys.getenv("baseUrl"),
                              authMethod = "windows")
cohorts <- tibble(
  cohortId = c(16329,
               7835, 
               16330,
               16484,
               2072,
               2087),
  cohortName =c("DOACs", 
                "Rivaroxaban",
                "Warfarin",
                "Major bleeding",
                "AMI",
                "Ischemic stroke")
)
cohortDefinitionSet <- ROhdsiWebApi::exportCohortDefinitionSet(baseUrl = Sys.getenv("baseUrl"),
                                                               cohortIds = cohorts$cohortId)
cohortDefinitionSet <- cohortDefinitionSet %>%
  select(-"cohortName") %>%
  inner_join(cohorts, by = join_by("cohortId"))
saveRDS(cohortDefinitionSet, "RealWorldExample/CohortDefinitionSet.rds")

# Create cohorts ---------------------------------------------------------------
library(CohortGenerator)
library(dplyr)
cohortDefinitionSet <- readRDS("RealWorldExample/CohortDefinitionSet.rds")
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE) %>%
  transmute(cohortId = conceptId,
            cohortName = conceptName,
            outcomeConceptId = conceptId)
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
  negativeControlOutcomeCohortSet = negativeControls
)

# Fetch CohortMethodData object ------------------------------------------------
library(CohortMethod)
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)
covariateSettings <- createDefaultCovariateSettings(
  excludedCovariateConceptIds = c(1592988, 40228152, 40241331, 43013024, 45775372, 45892847, 1310149),
  addDescendantsToExclude = TRUE
)

cmData <- getDbCohortMethodData(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  targetId = 16329,
  comparatorId = 16330,
  outcomeIds = c(2072, 2087, 16484, negativeControls$conceptId),
  exposureDatabaseSchema = cohortDatabaseSchema,
  exposureTable = cohortTable,
  outcomeDatabaseSchema = cohortDatabaseSchema,
  outcomeTable = cohortTable,
  restrictToCommonPeriod = TRUE,
  covariateSettings = covariateSettings
)
dir.create(folder, recursive = TRUE)
saveCohortMethodData(cmData, file.path(folder, "cmData.zip"))

# Fit propensity and mediator risk models --------------------------------------
library(CohortMethod)
library(MediationAnalysis)
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
plotPs(ps, 
       showEquiposeLabel = TRUE,
       showCountsLabel = TRUE,
       fileName = file.path(folder, "PsDistribution.png"))

mrs <- createMediatorRiskScore(cohortMethodData = cmData, 
                               mediatorId = 16484,
                               removeDuplicateSubjects = "keep first",
                               riskWindowStart = 0,
                               startAnchor = "cohort start",
                               riskWindowEnd = 0,
                               endAnchor = "cohort end",
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
plotMrsByExposure(mrs, fileName = file.path(folder, "MrsByExposure.png"))
plotMrsByMediator(mrs, fileName = file.path(folder, "MrsByMediator.png"))

# Compute estimates for negative controls --------------------------------------
library(CohortMethod)
library(MediationAnalysis)
cmData <- loadCohortMethodData(file.path(folder, "cmData.zip"))
ps <- readRDS(file.path(folder, "ps.rds"))
mrs <- readRDS(file.path(folder, "mrs.rds"))

estimates <- list()
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)
for (outcomeId in negativeControls$conceptId) {
  message(sprintf("Fitting model for outcome %d", outcomeId))
  studyPop <- createStudyPopulation(
    cohortMethodData = cmData,
    outcomeId = outcomeId,
    restrictToCommonPeriod = TRUE,
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback = 99999,
    removeDuplicateSubjects = "keep first",
    riskWindowStart = 0,
    startAnchor = "cohort start",
    riskWindowEnd = 0,
    endAnchor = "cohort end"
  )
  model <- fitMediatorModel(
    studyPopulation = studyPop,
    ps = ps,
    mrs = mrs,
    psAdjustment = "matching",
    mrsAdjustment = "model",
    mediatorType = "time-to-event")
  model <- model %>%
    mutate(outcomeId = !!outcomeId)
  estimates[[length(estimates) + 1]] <- model
}
estimates <- bind_rows(estimates)
estimates <- estimates %>%
  mutate(mainSeLogRr = (mainLogUb - mainLogLb) / (2*qnorm(0.975)),
         indirectSeLogRr = (indirectLogUb - indirectLogLb ) / (2*qnorm(0.975)))

EmpiricalCalibration::plotCalibrationEffect(
  logRrNegatives = estimates$mainLogHr,
  seLogRrNegatives = estimates$mainSeLogRr,
  title = "Main effect",
  xLabel = "Hazard ratio",
  showCis = TRUE,
  showExpectedAbsoluteSystematicError = TRUE,
  fileName = file.path(folder, "NegativeControlsMainEffect.png")
)

EmpiricalCalibration::plotCalibrationEffect(
  logRrNegatives = estimates$indirectLogHr,
  seLogRrNegatives = estimates$indirectSeLogRr,
  title = "Indirect effect",
  xLabel = "Hazard ratio",
  showCis = TRUE,
  showExpectedAbsoluteSystematicError = TRUE,
  fileName = file.path(folder, "NegativeControlsIndirectEffect.png")
)

# Compute rough estimate of bleeding incidence rate ----------------------------
library(CohortMethod)
library(dyplr)
cmData <- loadCohortMethodData(file.path(folder, "cmData.zip"))
studyPop <- createStudyPopulation(
  cohortMethodData = cmData,
  outcomeId = 10870,
  restrictToCommonPeriod = TRUE,
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback = 99999,
  removeDuplicateSubjects = "keep first",
  riskWindowStart = 0,
  startAnchor = "cohort start",
  riskWindowEnd = 0,
  endAnchor = "cohort end"
)
studyPop %>%
  summarise(outcomeCount = sum(outcomeCount),
            timeAtRisk = sum(timeAtRisk)) %>%
  mutate(ir = outcomeCount / (timeAtRisk / (1000*365.25)))
# outcomeCount timeAtRisk    ir
# <dbl>      <dbl> <dbl>
#   1         1158   49225009  8.59

# Run CohortDiagnostics --------------------------------------------------------
library(CohortDiagnostics)
cohortDefinitionSet <- readRDS("RealWorldExample/CohortDefinitionSet.rds")
executeDiagnostics(cohortDefinitionSet = cohortDefinitionSet,
                   exportFolder = file.path(folder, "cohortDiagnostics"),
                   databaseId = "CCAE",
                   cohortDatabaseSchema = cohortDatabaseSchema,
                   connectionDetails = connectionDetails,
                   cdmDatabaseSchema = cdmDatabaseSchema,
                   cohortTable = cohortTable)
createMergedResultsFile(dataFolder = file.path(folder, "cohortDiagnostics"),
                        sqliteDbPath = file.path(folder, "MergedCohortDiagnosticsData.sqlite"))
launchDiagnosticsExplorer()

# Create cohort explorer app ---------------------------------------------------
CohortExplorer::createCohortExplorerApp(
  cohortDatabaseSchema = cohortDatabaseSchema,
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortTable = cohortTable,
  cohortDefinitionId = 16484,
  databaseId = "CCAE",
  exportFolder = "d:/temp/ceMajorBleeding")
