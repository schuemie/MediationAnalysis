source("RealWorldExample/DatabaseDetails.R")

# Part 1: Create cohorts --------------------------------------------------
library(CohortGenerator)
library(dplyr)

# database = databases[[2]]
# databases = databases[2:5]
for (database in databases) {
  message(sprintf("Creating cohorts for %s", database$databaseId))
  dir.create(database$outputFolder, showWarnings = FALSE, recursive = TRUE)
  cohortDefinitionSet <- readRDS("RealWorldExample/CohortDefinitionSet.rds")
  negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE) %>%
    transmute(cohortId = conceptId,
              cohortName = conceptName,
              outcomeConceptId = conceptId)
  hasBledCohorts <- readRDS("RealWorldExample/HasBledCohortDefinitionSet.rds")
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
  generateCohortSet(connectionDetails = database$connectionDetails,
                    cdmDatabaseSchema = database$cdmDatabaseSchema,
                    cohortDatabaseSchema = database$cohortDatabaseSchema,
                    cohortTableNames = cohortTableNames,
                    cohortDefinitionSet = hasBledCohorts)
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
library(tidyr)
source("RealWorldExample/HasBled/HasBledCovariateBuilder.R")
tcmos <- readRDS("RealWorldExample/tcmos.rds") 
tcs <- tcmos %>%
  distinct(targetId, targetName, comparatorId, comparatorName)
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)

# database = databases[[2]]
# i = 1
for (database in databases) {
  message(sprintf("Computing diagnostics for %s", database$databaseId))
  for (i in seq_len(nrow(tcs))) {
    tc <- tcs[i, ]
    
    cmDataFileName <- file.path(database$outputFolder, sprintf("cmData_t%d_c%s.zip", tc$targetId, tc$comparatorId))
    
    if (!file.exists(cmDataFileName)) {
      defaultCovariateSettings <- createDefaultCovariateSettings(
        excludedCovariateConceptIds = c(1592988, 40228152, 40241331, 43013024, 45775372, 45892847, 1310149),
        addDescendantsToExclude = TRUE
      )
      hasBledCovariateSettings <- createHasBledCovariateSettings()
      covariateSettings <- list(defaultCovariateSettings, hasBledCovariateSettings) 
      cmData <- getDbCohortMethodData(
        connectionDetails = database$connectionDetails,
        cdmDatabaseSchema = database$cdmDatabaseSchema,
        targetId = tc$targetId,
        comparatorId = tc$comparatorId,
        outcomeIds = c(unique(tcmos$outcomeId), unique(tcmos$mediatorId), negativeControls$conceptId),
        exposureDatabaseSchema = database$cohortDatabaseSchema,
        exposureTable = database$cohortTable,
        outcomeDatabaseSchema = database$cohortDatabaseSchema,
        outcomeTable = database$cohortTable,
        restrictToCommonPeriod = TRUE,
        studyStartDate = "20101101",
        studyEndDate = "20221231",
        covariateSettings = covariateSettings
      )
      # cmData$outcomes %>% distinct(outcomeId) %>% pull()
      # Censor at switch:
      cohorts <- cmData$cohorts %>%
        arrange(personSeqId, cohortStartDate) %>%
        collect()
      idx <- which(duplicated(cohorts$personSeqId))
      cohorts$daysToCohortEnd[idx - 1] <- pmin(cohorts$daysToCohortEnd[idx - 1],
                                               as.numeric(cohorts$cohortStartDate[idx] - cohorts$cohortStartDate[idx - 1]) - 1)
      cohorts <- cohorts %>%
        filter(daysToCohortEnd > 0)
      counts <- CohortMethod:::getCounts(cohorts, "Censor at switch")
      cmData$cohorts <- cohorts
      attr(cmData, "metaData")$attrition <- bind_rows(attr(cmData, "metaData")$attrition,
                                                      counts)
      saveCohortMethodData(cmData, cmDataFileName)
    } 
    cmData <- loadCohortMethodData(cmDataFileName)   
    psFileName <- file.path(database$outputFolder, sprintf("ps_t%d_c%s.rds", tc$targetId, tc$comparatorId))
    if (!file.exists(psFileName)) {
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
      plotPs(ps, 
             showEquiposeLabel = TRUE,
             showCountsLabel = TRUE,
             fileName = file.path(database$outputFolder, sprintf("ps_t%d_c%s.png", tc$targetId, tc$comparatorId)))
      equipoise <- bind_cols(tc, data.frame(equipoise = computeEquipoise(ps)))
      saveRDS(equipoise, file.path(database$outputFolder, sprintf("equipoise_t%d_c%s.rds", tc$targetId, tc$comparatorId)))
      saveRDS(ps, psFileName)
    } else {
      ps <- readRDS(psFileName)
    }
    balanceFileName <- file.path(database$outputFolder, sprintf("balance_t%d_c%s.rds", tc$targetId, tc$comparatorId))
    if (!file.exists(balanceFileName)) {
      studyPop <- createStudyPopulation(
        cohortMethodData = cmData,
        population = ps,
        restrictToCommonPeriod = TRUE,
        removeDuplicateSubjects = "keep first",
        riskWindowStart = 0,
        startAnchor = "cohort start",
        riskWindowEnd = 0,
        endAnchor = "cohort end"
      )
      strataPop <- matchOnPs(studyPop, maxRatio = 100)
      balance <- computeCovariateBalance(strataPop, cmData)
      plotCovariateBalanceScatterPlot(
        balance = balance,
        threshold = 0.1,
        showCovariateCountLabel = TRUE,
        showMaxLabel = TRUE,
        fileName = file.path(database$outputFolder, sprintf("balanceScatterplot_t%d_c%s.png", tc$targetId, tc$comparatorId))
      )
      saveRDS(balance, balanceFileName)
    }
    mdrrFileName <- file.path(database$outputFolder, sprintf("mdrr_t%d_c%s.csv", tc$targetId, tc$comparatorId))
    if (!file.exists(mdrrFileName)) {
      outcomes <- tcmos %>%
        inner_join(tcs, by = join_by(targetId, targetName, comparatorId, comparatorName)) %>%
        distinct(outcomeId, outcomeName)
      mdrrs <- list()
      for (j in seq_len(nrow(outcomes))) {
        outcome <- outcomes[j, ]
        studyPop <- createStudyPopulation(
          cohortMethodData = cmData,
          population = ps,
          outcomeId = outcome$outcomeId,
          restrictToCommonPeriod = TRUE,
          removeSubjectsWithPriorOutcome = TRUE,
          priorOutcomeLookback = 365,
          removeDuplicateSubjects = "keep first",
          riskWindowStart = 0,
          startAnchor = "cohort start",
          riskWindowEnd = 0,
          endAnchor = "cohort end"
        )
        strataPop <- matchOnPs(studyPop, maxRatio = 100)
        mdrr <- computeMdrr(strataPop)
        followup <- studyPop %>%
          group_by(treatment) %>%
          summarise(p10FollowUp = quantile(survivalTime, 0.1),
                    p25FollowUp = quantile(survivalTime, 0.25),
                    medianFollowup = median(survivalTime),
                    p75FollowUp = quantile(survivalTime, 0.75),
                    p90FollowUp = quantile(survivalTime, 0.9)) %>%
          mutate(treatment = if_else(treatment == 1, "Target", "Comparator")) %>%
          pivot_wider(names_from = treatment,
                      values_from = c("p10FollowUp", "p25FollowUp", "medianFollowup", "p75FollowUp", "p90FollowUp"))
        mdrr <- bind_cols(tc, outcome, mdrr, followup)
        mdrrs[[j]] <- mdrr
      }
      mdrrs <- bind_rows(mdrrs)
      readr::write_csv(mdrrs, mdrrFileName)
    }
    mediators <- tcmos %>%
      inner_join(tc, by = join_by(targetId, targetName, comparatorId, comparatorName)) %>%
      distinct(mediatorId, mediatorName)
    # j = 1
    for (j in seq_len(nrow(mediators))) {
      mediator <- mediators[j, ]
      mrsFileName <- file.path(database$outputFolder, sprintf("mrs_t%d_c%s_m%d.png", tc$targetId, tc$comparatorId, mediator$mediatorId))
      if (!file.exists(mrsFileName)) {
        mrs <- createMediatorRiskScore(cohortMethodData = cmData, 
                                       mediatorId = mediator$mediatorId,
                                       removeDuplicateSubjects = "keep first",
                                       removeSubjectsWithPriorOutcome = TRUE,
                                       priorOutcomeLookback = 30,
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
        plotMrsByExposure(mrs, fileName = file.path(database$outputFolder, sprintf("mrsByExposure_t%d_c%s_m%d.png", tc$targetId, tc$comparatorId, mediator$mediatorId)))
        plotMrsByMediator(mrs, fileName = file.path(database$outputFolder, sprintf("mrsByMediator_t%d_c%s_m%d.png", tc$targetId, tc$comparatorId, mediator$mediatorId)))
        saveRDS(mrs, mrsFileName)
      } else {
        mrs <- readRDS(mrsFileName)
      }
    
      ncsFileName <- file.path(database$outputFolder, sprintf("ncs_t%d_c%s_m%d.rds", tc$targetId, tc$comparatorId, mediator$mediatorId))

      if (!file.exists(ncsFileName)) {
        estimates <- list()
        negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)
        for (outcomeId in negativeControls$conceptId) {
          message(sprintf("Fitting model for outcome %d", outcomeId))
          studyPop <- createStudyPopulation(
            cohortMethodData = cmData,
            outcomeId = outcomeId,
            restrictToCommonPeriod = TRUE,
            removeSubjectsWithPriorOutcome = TRUE,
            priorOutcomeLookback = 365,
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
                 indirectSeLogRr = (indirectLogUb - indirectLogLb ) / (2*qnorm(0.975)),
                 directSeLogRr = (directLogUb - directLogLb ) / (2*qnorm(0.975)))
        nullMain <- EmpiricalCalibration::fitMcmcNull(
          logRr = estimates$mainLogHr,
          seLogRr = estimates$mainSeLogRr
        )
        EmpiricalCalibration::plotCalibrationEffect(
          logRrNegatives = estimates$mainLogHr,
          seLogRrNegatives = estimates$mainSeLogRr,
          null = nullMain,
          title = "Main effect",
          xLabel = "Hazard ratio",
          showCis = TRUE,
          showExpectedAbsoluteSystematicError = TRUE,
          fileName = file.path(database$outputFolder, sprintf("ncsMainEffect_t%d_c%s_m%d.png", tc$targetId, tc$comparatorId, mediator$mediatorId))
        )
        nullDirect <- EmpiricalCalibration::fitMcmcNull(
          logRr = estimates$directLogHr,
          seLogRr = estimates$directSeLogRr
        )
        EmpiricalCalibration::plotCalibrationEffect(
          logRrNegatives = estimates$directLogHr,
          seLogRrNegatives = estimates$directSeLogRr,
          null = nullDirect,
          title = "Direct effect",
          xLabel = "Hazard ratio",
          showCis = TRUE,
          showExpectedAbsoluteSystematicError = TRUE,
          fileName = file.path(database$outputFolder, sprintf("ncsDirectEffect_t%d_c%s_m%d.png", tc$targetId, tc$comparatorId, mediator$mediatorId))
        )      
        nullIndirect <- EmpiricalCalibration::fitMcmcNull(
          logRr = estimates$indirectLogHr,
          seLogRr = estimates$indirectSeLogRr
        )
        EmpiricalCalibration::plotCalibrationEffect(
          logRrNegatives = estimates$indirectLogHr,
          seLogRrNegatives = estimates$indirectSeLogRr,
          null = nullIndirect,
          title = "Indirect effect",
          xLabel = "Hazard ratio",
          showCis = TRUE,
          showExpectedAbsoluteSystematicError = TRUE,
          fileName = file.path(database$outputFolder, sprintf("ncsIndirectEffect_t%d_c%s_m%d.png", tc$targetId, tc$comparatorId, mediator$mediatorId))
        )
        ease <- bind_cols(tc, mediator) %>%
          mutate(easeMain = EmpiricalCalibration::computeExpectedAbsoluteSystematicError(nullMain)$ease,
                 easeDirect = EmpiricalCalibration::computeExpectedAbsoluteSystematicError(nullDirect)$ease,
                 easeIndirect = EmpiricalCalibration::computeExpectedAbsoluteSystematicError(nullIndirect)$ease) 
        readr::write_csv(ease, file.path(database$outputFolder, sprintf("ease_t%d_c%s_m%d.png", tc$targetId, tc$comparatorId, mediator$mediatorId)))
        saveRDS(estimates, ncsFileName)
      }
    }
  }
}
    