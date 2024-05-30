source("RealWorldExample/DatabaseDetails.R")

# Part 1: Create cohorts --------------------------------------------------
library(CohortGenerator)
library(dplyr)

# database = databases[[1]]
# databases = databases[c(5,4)]
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

# database = databases[[1]]
# i = 2
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
      hasBledCovariateSettings <- createHasBledCovariateSettings(
        hasBledCohortDatabaseSchema = database$cohortDatabaseSchema,
        hasBledCohortTable = database$cohortTable)
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
    # Magical code to fix mysterious error (cannot get a slot ("slots"))
    attr(cmData, "metaData")$call <- NULL
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
      mrsFileName <- file.path(database$outputFolder, sprintf("mrs_t%d_c%s_m%d.rds", tc$targetId, tc$comparatorId, mediator$mediatorId))
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
      
      # outcomes <- tcmos %>%
      #   filter(targetId == tc$targetId,
      #          comparatorId == tc$comparatorId,
      #          mediatorId == mediator$mediatorId)
      # # outcomeId = outcomes$outcomeId[1]
      # for (outcomeId in outcomes$outcomeId) {
      #   table1FileName <- file.path(database$outputFolder, sprintf("table1_t%d_c%s_m%d_o%d.rds", tc$targetId, tc$comparatorId, mediator$mediatorId, outcomeId))
      #   if (!file.exists(table1FileName)) {
      #     studyPop <- createStudyPopulation(
      #       cohortMethodData = cmData,
      #       outcomeId = outcomeId,
      #       restrictToCommonPeriod = TRUE,
      #       removeSubjectsWithPriorOutcome = TRUE,
      #       priorOutcomeLookback = 365,
      #       removeDuplicateSubjects = "keep first",
      #       riskWindowStart = 0,
      #       startAnchor = "cohort start",
      #       riskWindowEnd = 0,
      #       endAnchor = "cohort end"
      #     )
      #     studyPop <- studyPop %>%
      #       inner_join(mrs %>%
      #                    select("rowId"),
      #                  by = join_by(rowId)) %>%
      #       inner_join(ps %>%
      #                    select("rowId", "propensityScore"),
      #                  by = join_by(rowId))
      #     strataPop <- matchOnPs(studyPop, maxRatio = 100)
      #     covariateData1 <- FeatureExtraction::filterByRowId(cmData,
      #                                                        strataPop %>%
      #                                                          filter(treatment == 1) %>%
      #                                                          pull(rowId))
      #     covariateData2 <- FeatureExtraction::filterByRowId(cmData,
      #                                                        strataPop %>%
      #                                                          filter(treatment == 0) %>%
      #                                                          pull(rowId))
      #     covariateData1 <- FeatureExtraction::aggregateCovariates(covariateData1)
      #     covariateData2 <- FeatureExtraction::aggregateCovariates(covariateData2)
      #     specs <- FeatureExtraction::getDefaultTable1Specifications() %>%
      #       bind_rows(tibble(label = "HAS-BLED",
      #                        analysisId = 999,
      #                        covariateIds = "999"))
      #     
      #     table1 <- FeatureExtraction::createTable1(covariateData1 = covariateData1,
      #                                               covariateData2 = covariateData2,
      #                                               specifications = specs)
      #     saveRDS(table1, table1FileName)
      #   }
      # }
      
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
        readr::write_csv(ease, file.path(database$outputFolder, sprintf("ease_t%d_c%s_m%d.csv", tc$targetId, tc$comparatorId, mediator$mediatorId)))
        saveRDS(estimates, ncsFileName)
      }
    }
  }
}

# Fix filenames ---------------------------------------------------------------
for (database in databases) {
  message(sprintf("Fixing filenames for %s", database$databaseId))
  toFix <- list.files(database$outputFolder, "^mrs_t.*png$", full.names = TRUE)
  if (length(toFix) > 0) {
    newNames <- gsub(".png$", ".rds", toFix)
    file.rename(toFix, newNames)
  }
}
for (database in databases) {
  message(sprintf("Fixing filenames for %s", database$databaseId))
  toFix <- list.files(database$outputFolder, "^ease_t.*png$", full.names = TRUE)
  if (length(toFix) > 0) {
    newNames <- gsub(".png", ".csv", toFix)
    file.rename(toFix, newNames)
  }
}

# Back up old cmData files before refetching:
for (database in databases) {
  message(sprintf("Fixing filenames for %s", database$databaseId))
  toFix <- list.files(database$outputFolder, "^cmData_t.*zip$", full.names = TRUE)
  if (length(toFix) > 0) {
    newNames <- paste0(toFix, ".bak")
    file.rename(toFix, newNames)
  }
}

# Part 3: Generate results for outcomes of interest ------------------------------------
library(CohortMethod)
library(MediationAnalysis)
library(tidyr)
source("RealWorldExample/HasBled/HasBledCovariateBuilder.R")
tcmos <- readRDS("RealWorldExample/tcmos.rds") 
tcs <- tcmos %>%
  distinct(targetId, targetName, comparatorId, comparatorName)
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)
table1Specs <- getDefaultCmTable1Specifications() %>%
  bind_rows(tibble(label = "Mean score", analysisId = 999, covariateId = "999"))
# database = databases[[1]]
# i = 1
for (database in databases) {
  message(sprintf("Computing results for %s", database$databaseId))
  for (i in seq_len(nrow(tcs))) {
    tc <- tcs[i, ]
    cmDataFileName <- file.path(database$outputFolder, sprintf("cmData_t%d_c%s.zip", tc$targetId, tc$comparatorId))
    cmData <- loadCohortMethodData(cmDataFileName)   
    # Magical code to fix mysterious error (cannot get a slot ("slots"))
    attr(cmData, "metaData")$call <- NULL
    psFileName <- file.path(database$outputFolder, sprintf("ps_t%d_c%s.rds", tc$targetId, tc$comparatorId))
    ps <- readRDS(psFileName)
    mediators <- tcmos %>%
      inner_join(tc, by = join_by(targetId, targetName, comparatorId, comparatorName)) %>%
      distinct(mediatorId, mediatorName)
    # j = 1
    for (j in seq_len(nrow(mediators))) {
      mediator <- mediators[j, ]
      mrsFileName <- file.path(database$outputFolder, sprintf("mrs_t%d_c%s_m%d.rds", tc$targetId, tc$comparatorId, mediator$mediatorId))
      mrs <- readRDS(mrsFileName)
      ncsFileName <- file.path(database$outputFolder, sprintf("ncs_t%d_c%s_m%d.rds", tc$targetId, tc$comparatorId, mediator$mediatorId))
      ncEstimates <- readRDS(ncsFileName)
      hoisFileName <- file.path(database$outputFolder, sprintf("hois_t%d_c%s_m%d.rds", tc$targetId, tc$comparatorId, mediator$mediatorId))
      # if (!file.exists(hoisFileName)) {
      if (TRUE) {
        outcomes <- tcmos %>%
          filter(targetId == tc$targetId,
                 comparatorId == tc$comparatorId,
                 mediatorId == mediator$mediatorId)
        estimates <- list()
        # k = 1
        for (k in seq_len(nrow(outcomes))) {
          outcome = outcomes[k, ]
          studyPop <- createStudyPopulation(
            cohortMethodData = cmData,
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
          table1FileName <- file.path(rootFolder, sprintf("Table1_t%d_c%s_m%d_o%d_%s.csv", tc$targetId, tc$comparatorId, mediator$mediatorId, outcome$outcomeId, database$databaseId))
          if (!file.exists(table1FileName)) {
            strataPop <- matchOnPs(inner_join(studyPop, select(ps, "rowId", "propensityScore"), by = join_by(rowId)), maxRatio = 100)
            balance <- computeCovariateBalance(
              population = strataPop,
              cohortMethodData = cmData,
              covariateFilter = table1Specs
            )
            cohorts <- cmData$cohorts %>%
              collect()
            table1 <- createCmTable1(balance, table1Specs)
            table1 <- rbind(
              table1[1, ],
              c("", 
                paste0("n=", format(sum(cohorts$treatment), big.mark = ",")), 
                paste0("n=", format(sum(!cohorts$treatment), big.mark = ",")), 
                "", 
                paste0("n=", format(sum(strataPop$treatment), big.mark = ",")), 
                paste0("n=", format(sum(!strataPop$treatment), big.mark = ",")), 
                ""),
              table1[2:nrow(table1), ]
            )
            table1[nrow(table1), c(2, 3, 5, 6)] <- format(as.numeric(table1[nrow(table1), c(2, 3, 5, 6)]) / 100, digits = 2)
            table1[nrow(table1), 1] <- "HAS-BLED"
            readr::write_csv(table1, table1FileName)
          }
          set.seed(123) # Make Bootstrap reproducible
          model <- fitMediatorModel(
            studyPopulation = studyPop,
            ps = ps,
            mrs = mrs,
            psAdjustment = "matching",
            mrsAdjustment = "model",
            mediatorType = "time-to-event")
          model <- model %>%
            mutate(outcomeId = outcome$outcomeId,
                   outcomeName = outcome$outcomeName)
          estimates[[k]] <- model
        }
        estimates <- bind_rows(estimates)
        set.seed(123) # Make MCMC reproducible
        nullMain <- EmpiricalCalibration::fitMcmcNull(
          logRr = ncEstimates$mainLogHr,
          seLogRr = ncEstimates$mainSeLogRr
        )
        nullDirect <- EmpiricalCalibration::fitMcmcNull(
          logRr = ncEstimates$directLogHr,
          seLogRr = ncEstimates$directSeLogRr
        )
        nullIndirect <- EmpiricalCalibration::fitMcmcNull(
          logRr = ncEstimates$indirectLogHr,
          seLogRr = ncEstimates$indirectSeLogRr
        )
        estimates <- estimates %>%
          mutate(mainSeLogRr = (mainLogUb - mainLogLb) / (2*qnorm(0.975)),
                 indirectSeLogRr = (indirectLogUb - indirectLogLb ) / (2*qnorm(0.975)),
                 directSeLogRr = (directLogUb - directLogLb ) / (2*qnorm(0.975)))
        mainCalibratedP <- EmpiricalCalibration::calibrateP(nullMain, estimates$mainLogHr, estimates$mainSeLogRr)
        indirectCalibratedP <- EmpiricalCalibration::calibrateP(nullIndirect, estimates$indirectLogHr, estimates$indirectSeLogRr)
        directCalibratedP <- EmpiricalCalibration::calibrateP(nullDirect, estimates$directLogHr, estimates$directSeLogRr)
        mainCalibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(estimates$mainLogHr, estimates$mainSeLogRr, EmpiricalCalibration::convertNullToErrorModel(nullMain))
        indirectCalibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(estimates$indirectLogHr, estimates$indirectSeLogRr, EmpiricalCalibration::convertNullToErrorModel(nullIndirect))
        directCalibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(estimates$directLogHr, estimates$directSeLogRr, EmpiricalCalibration::convertNullToErrorModel(nullDirect))
        estimates <- estimates %>%
          mutate(mainCalibratedHr = exp(mainCalibratedCi$logRr),
                 mainCalibratedLb = exp(mainCalibratedCi$logLb95Rr),
                 mainCalibratedUb = exp(mainCalibratedCi$logUb95Rr),
                 mainCalibratedP = !!mainCalibratedP$p,
                 indirectCalibratedHr = exp(indirectCalibratedCi$logRr),
                 indirectCalibratedLb = exp(indirectCalibratedCi$logLb95Rr),
                 indirectCalibratedUb = exp(indirectCalibratedCi$logUb95Rr),
                 indirectCalibratedP =!!indirectCalibratedP$p,
                 directCalibratedHr = exp(directCalibratedCi$logRr),
                 directCalibratedLb = exp(directCalibratedCi$logLb95Rr),
                 directCalibratedUb = exp(directCalibratedCi$logUb95Rr),
                 directCalibratedP = !!directCalibratedP$p)
        saveRDS(estimates, hoisFileName)
      }
    }
  }
}

# Part 4: Generate forest plots ------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
tcmos <- readRDS("RealWorldExample/tcmos.rds") 
tcms <- tcmos %>%
  distinct(targetId, 
           targetName,
           comparatorId,
           comparatorName,
           mediatorId,
           mediatorName)
# database = databases[[1]]
estimates <- list()
for (database in databases) {
  for (i in seq_len(nrow(tcms))) {
    tcm <- tcms[i, ]
    hoisFileName <- file.path(database$outputFolder, sprintf("hois_t%d_c%s_m%d.rds", tcm$targetId, tcm$comparatorId, tcm$mediatorId))
    temp <- readRDS(hoisFileName)
    temp <- temp %>%
      bind_cols(tcm) %>%
      mutate(database = database$databaseId)
    estimates[[length(estimates) + 1]] <- temp
  }
}
estimates <- bind_rows(estimates)

# Forest plot
# i = 1
for (i in seq_len(nrow(tcmos))) {
  tcmo <- tcmos[i, ]
  plotFileName <- file.path(rootFolder, sprintf("hr_t%d_c%d_m%d_o%d.png",
                                                tcmo$targetId,
                                                tcmo$comparatorId,
                                                tcmo$mediatorId,
                                                tcmo$outcomeId))
  data <- estimates %>%
    inner_join(tcmo, by = join_by(outcomeId, outcomeName, targetId, targetName, comparatorId, comparatorName, mediatorId, mediatorName))
  data <- bind_rows(
    data %>%
      transmute(estimand = "Main effect",
                hr = mainCalibratedHr,
                lb = mainCalibratedLb,
                ub = mainCalibratedUb,
                database = database),
    data %>%
      transmute(estimand = "Direct effect",
                hr = directCalibratedHr,
                lb = directCalibratedLb,
                ub = directCalibratedUb,
                database = database),
    data %>%
      transmute(estimand = "Indirect effect",
                hr = indirectCalibratedHr,
                lb = indirectCalibratedLb,
                ub = indirectCalibratedUb,
                database = database)
  )
  data$database <- sprintf("%s    %0.2f (%0.2f-%0.2f)",
                           data$database,
                           data$hr,
                           data$lb,
                           data$ub)
  data$estimand <- factor(data$estimand, levels = c("Main effect", "Direct effect", "Indirect effect"))
  data$database <- factor(data$database, levels = rev(sort(unique(data$database))))
  title <- sprintf("%s vs %s for %s\nMediator: %s",
                   tcmo$targetName,
                   tcmo$comparatorName,
                   tcmo$outcomeName,
                   tcmo$mediatorName)
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  plot <- ggplot(data, aes(x = hr, y = database)) +
    geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
    geom_vline(xintercept = 1, size = 0.5) +
    geom_point(size = 3, shape = 18) +
    geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.15) +
    scale_x_log10("Hazard Ratio", breaks = breaks) +
    coord_cartesian(xlim = c(0.4, 2.5)) +
    facet_grid(estimand~., scales = "free_y") +
    ggtitle(title) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_line(colour = "white"),
      plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines"),
      plot.title = element_text(size = 12, hjust = 0.5)
    )
  plot
  ggsave(plotFileName, plot, width = 5, height = 3.5, dpi = 200)
}

# Table of counts
table <- estimates %>%
  select(Target = "targetName",
         Comparator ="comparatorName",
         Mediator = "mediatorName",
         Outcome = "outcomeName",
         Database = "database",
         "targetSubjects",
         "comparatorSubjects",
         "targetMediators",
         "comparatorMediators",
         "targetOutcomes",
         "comparatorOutcomes",
         "targetMediatorOutcomes",
         "comparatorMediatorOutcomes")
readr::write_csv(table, file.path(rootFolder, "Counts.csv"))

# Tables for AHA abstract:
filteredEstimates <- estimates %>%
  filter(mediatorId == 17003, targetId == 16586)

hrTable <- filteredEstimates %>%
  mutate(hrMain = sprintf("%0.2f (%0.2f - %0.2f)", exp(mainLogHr), exp(mainLogLb), exp(mainLogUb)),
         hrIndirect = sprintf("%0.2f (%0.2f - %0.2f)", exp(indirectLogHr), exp(indirectLogLb), exp(indirectLogUb))) %>%
  select("database", "outcomeName", "hrMain", "hrIndirect") %>%
  pivot_wider(names_from = "outcomeName", values_from = c("hrMain", "hrIndirect"), names_sort = TRUE) %>%
  arrange("database")
hrTable


countTable <- filteredEstimates %>%
  group_by(outcomeName) %>%
  summarise(outcomes = format(sum(targetOutcomes + comparatorOutcomes), big.mark = ",", trim = TRUE),
            mediators = format(sum(targetMediators + comparatorMediators), big.mark = ",", trim = TRUE),
            mediatorOutcomes = format(sum(targetMediatorOutcomes + comparatorMediatorOutcomes), big.mark = ",", trim = TRUE)) %>%
  t()
countTable <- cbind(countTable, rownames(countTable))
colnames(countTable) <- c(paste0("hrMain_", countTable[1, seq_len(ncol(countTable) -1)]), "database")
countTable <- as_tibble(countTable)


combinedTable <- countTable %>% 
  rename
