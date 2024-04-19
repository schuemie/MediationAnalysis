source("RealWorldExample/DatabaseDetails.R")

# Part 1: Create cohorts --------------------------------------------------
library(CohortGenerator)
library(dplyr)

# database = databases[[1]]
# databases = databases[2:5]
for (database in databases) {
  message(sprintf("Creating cohorts for %s", database$databaseId))
  dir.create(database$outputFolder, showWarnings = FALSE, recursive = TRUE)
  cohortDefinitionSet <- readRDS("RealWorldExample/CohortDefinitionSet.rds")
  negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE) %>%
    transmute(cohortId = conceptId,
              cohortName = conceptName,
              outcomeConceptId = conceptId)
  cohortTableNames <- getCohortTableNames(cohortTable = database$cohortTable)
  createCohortTables(connectionDetails = database$connectionDetails,
                     cohortDatabaseSchema = database$cohortDatabaseSchema,
                     cohortTableNames = cohortTableNames)
  generateCohortSet(connectionDetails = database$connectionDetails,
                    cdmDatabaseSchema = database$cdmDatabaseSchema,
                    cohortDatabaseSchema = database$cohortDatabaseSchema,
                    cohortTableNames = cohortTableNames,
                    cohortDefinitionSet = cohortDefinitionSet)
  generateNegativeControlOutcomeCohorts(
    connectionDetails = database$connectionDetails,
    cdmDatabaseSchema = database$cdmDatabaseSchema,
    cohortDatabaseSchema = database$cohortDatabaseSchema,
    cohortTable = database$cohortTable,
    negativeControlOutcomeCohortSet = negativeControls
  )
  counts <- CohortGenerator::getCohortCounts(
    connectionDetails = database$connectionDetails,
    cohortDatabaseSchema = database$cohortDatabaseSchema,
    cohortTable = database$cohortTable,
    cohortDefinitionSet = cohortDefinitionSet,
    cohortIds = cohortDefinitionSet$cohortId
  )
  counts <- counts %>%
    select("cohortId", "cohortName", "cohortEntries", "cohortSubjects")
  readr::write_csv(counts, file.path(database$outputFolder, "CohortCounts.csv"))
}

# Part 2: Compute diagnostics --------------------------------------------------
library(CohortMethod)
library(MediationAnalysis)

# database = databases[[1]]
for (database in databases) {
  message(sprintf("Computing diagnostics for %s", database$databaseId))
  
  # Fetch CohortMethodData object
  negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)
  covariateSettings <- createDefaultCovariateSettings(
    excludedCovariateConceptIds = c(1592988, 40228152, 40241331, 43013024, 45775372, 45892847, 1310149),
    addDescendantsToExclude = TRUE
  )
  
  cmData <- getDbCohortMethodData(
    connectionDetails = database$connectionDetails,
    cdmDatabaseSchema = database$cdmDatabaseSchema,
    targetId = 16329,
    comparatorId = 16330,
    outcomeIds = c(2072, 2087, 16484, negativeControls$conceptId),
    exposureDatabaseSchema = database$cohortDatabaseSchema,
    exposureTable = database$cohortTable,
    outcomeDatabaseSchema = database$cohortDatabaseSchema,
    outcomeTable = database$cohortTable,
    restrictToCommonPeriod = TRUE,
    covariateSettings = covariateSettings
  )
  saveCohortMethodData(cmData, file.path(database$outputFolder, "cmData.zip"))
  
  # Fit propensity and mediator risk models
  cmData <- loadCohortMethodData(file.path(database$outputFolder, "cmData.zip"))
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
  saveRDS(ps, file.path(database$outputFolder, "ps.rds"))
  ps <- readRDS(file.path(database$outputFolder, "ps.rds"))
  plotPs(ps, 
         showEquiposeLabel = TRUE,
         showCountsLabel = TRUE,
         fileName = file.path(database$outputFolder, "PsDistribution.png"))
  
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
  saveRDS(mrs, file.path(database$outputFolder, "mrs.rds"))
  mrs <- readRDS(file.path(database$outputFolder, "mrs.rds"))
  plotMrsByExposure(mrs, fileName = file.path(database$outputFolder, "MrsByExposure.png"))
  plotMrsByMediator(mrs, fileName = file.path(database$outputFolder, "MrsByMediator.png"))
  
  # Compute estimates for negative controls
  library(CohortMethod)
  library(MediationAnalysis)
  cmData <- loadCohortMethodData(file.path(database$outputFolder, "cmData.zip"))
  ps <- readRDS(file.path(database$outputFolder, "ps.rds"))
  mrs <- readRDS(file.path(database$outputFolder, "mrs.rds"))
  
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
    fileName = file.path(database$outputFolder, "NegativeControlsMainEffect.png")
  )
  
  EmpiricalCalibration::plotCalibrationEffect(
    logRrNegatives = estimates$indirectLogHr,
    seLogRrNegatives = estimates$indirectSeLogRr,
    title = "Indirect effect",
    xLabel = "Hazard ratio",
    showCis = TRUE,
    showExpectedAbsoluteSystematicError = TRUE,
    fileName = file.path(database$outputFolder, "NegativeControlsIndirectEffect.png")
  )
  
  # # Compute rough estimate of bleeding incidence rate ----------------------------
  # library(CohortMethod)
  # library(dyplr)
  # cmData <- loadCohortMethodData(file.path(folder, "cmData.zip"))
  # studyPop <- createStudyPopulation(
  #   cohortMethodData = cmData,
  #   outcomeId = 10870,
  #   restrictToCommonPeriod = TRUE,
  #   removeSubjectsWithPriorOutcome = TRUE,
  #   priorOutcomeLookback = 99999,
  #   removeDuplicateSubjects = "keep first",
  #   riskWindowStart = 0,
  #   startAnchor = "cohort start",
  #   riskWindowEnd = 0,
  #   endAnchor = "cohort end"
  # )
  # studyPop %>%
  #   summarise(outcomeCount = sum(outcomeCount),
  #             timeAtRisk = sum(timeAtRisk)) %>%
  #   mutate(ir = outcomeCount / (timeAtRisk / (1000*365.25)))
  # # outcomeCount timeAtRisk    ir
  # # <dbl>      <dbl> <dbl>
  # #   1         1158   49225009  8.59
  # 
  # # Run CohortDiagnostics --------------------------------------------------------
  # library(CohortDiagnostics)
  # cohortDefinitionSet <- readRDS("RealWorldExample/CohortDefinitionSet.rds")
  # executeDiagnostics(cohortDefinitionSet = cohortDefinitionSet,
  #                    exportFolder = file.path(folder, "cohortDiagnostics"),
  #                    databaseId = "CCAE",
  #                    cohortDatabaseSchema = cohortDatabaseSchema,
  #                    connectionDetails = connectionDetails,
  #                    cdmDatabaseSchema = cdmDatabaseSchema,
  #                    cohortTable = cohortTable)
  # createMergedResultsFile(dataFolder = file.path(folder, "cohortDiagnostics"),
  #                         sqliteDbPath = file.path(folder, "MergedCohortDiagnosticsData.sqlite"))
  # launchDiagnosticsExplorer()
  # 
  # # Create cohort explorer app ---------------------------------------------------
  # CohortExplorer::createCohortExplorerApp(
  #   cohortDatabaseSchema = cohortDatabaseSchema,
  #   connectionDetails = connectionDetails,
  #   cdmDatabaseSchema = cdmDatabaseSchema,
  #   cohortTable = cohortTable,
  #   cohortDefinitionId = 16484,
  #   databaseId = "CCAE",
  #   exportFolder = "d:/temp/ceMajorBleeding")
  
}
