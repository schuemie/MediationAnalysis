# Copyright 2024 Observational Health Data Sciences and Informatics
#
# This file is part of MediationAnalysis
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
# library(dplyr)

#' Run a set of simulations
#'
#' @param folder                 Folder to write (intermediate) results.
#' @param simulationSettingsList A list of simulation settings objects as created
#'                               by `createSimulationSettings()`.
#' @param modelSettingsList      A list of model fitting settings as created by
#'                               `createModelSettings()`.
#' @param nSimulations           Number of times to repeat each simulation
#' @param maxCores               Maximum number of CPU cores to use at the same 
#'                               time.
#'
#' @return
#' Does not return anything. Writes results to the folder.
#' 
#' @export
runSetOfSimulations <- function(folder, 
                                simulationSettingsList, 
                                modelSettingsList,
                                nSimulations = 1000,
                                maxCores = 4) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  ParallelLogger::addDefaultFileLogger(file.path(folder, "log.txt"))
  ParallelLogger::addDefaultErrorReportLogger(file.path(folder, "errorReportR.txt"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER", silent = TRUE))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER", silent = TRUE), add = TRUE)
  
  cluster <- ParallelLogger::makeCluster(maxCores)
  on.exit(ParallelLogger::stopCluster(cluster))
  allResults <- list()
  # simulationSettings = simulationSettingsList[[1]]
  # modelSettings = modelSettingsList[[1]]
  for (i in seq_along(simulationSettingsList)) {
    simulationSettings <- simulationSettingsList[[i]]
    message(sprintf("Running scenario %d of %d", i, length(simulationSettingsList)))
    fileName <- file.path(folder, 
                          sprintf("result_%s.rds", 
                                  digest::digest(list(simulationSettings))))
    if (file.exists(fileName)) {
      results <- readRDS(fileName)
    } else {
      estimates <- ParallelLogger::clusterApply(cluster, 
                                                seq_len(nSimulations), 
                                                runOneSimulation, 
                                                simulationSettings = simulationSettings,
                                                modelSettingsList = modelSettingsList)
      results <- list()
      for (j in seq_along(modelSettingsList)) {
        modelSettings <- modelSettingsList[[j]]
        modelEstimates <- lapply(estimates, function(x) x[[j]]) %>%
          bind_rows()
        results[[j]] <- evaluateSingleResult(simulationSettings, modelSettings, modelEstimates)
      }
      results <- bind_rows(results)
      saveRDS(results, fileName)
    }
    allResults[[i]] <- results
  }
  allResults <- bind_rows(allResults)
  readr::write_csv(allResults, file.path(folder, "Results.csv"))
}

runOneSimulation <- function(seed, simulationSettings, modelSettingsList) {
  set.seed(seed)
  data <- simulateData(simulationSettings)
  estimates <- list()
  for (i in seq_along(modelSettingsList)) {
    modelSettings <- modelSettingsList[[i]]
    estimates[[i]] <- fitModel(data, modelSettings)
  }
  return(estimates)
}

evaluateSingleResult <- function(simulationSettings, modelSettings, estimates) {
  nonEstimableIdx <- is.na(estimates$mainLogDiff)
  estimates <- estimates[!nonEstimableIdx, ]
  hasIndirectEffect <- simulationSettings$mA != 0 & simulationSettings$yM != 0
  results <- tibble(
    coverageDirectEffect = mean(simulationSettings$yA >= estimates$directLogLb & 
                                  simulationSettings$yA <= estimates$directLogUb),
    coverageMediatorEffect = mean(simulationSettings$yM >= estimates$mediatorLogLb & 
                                    simulationSettings$yM <= estimates$mediatorLogUb),
    covarageMainEffect = mean(log(estimates$hrMain) >= estimates$mainLogLb &
                                log(estimates$hrMain) <= estimates$mainLogUb),
    covarageIndirectEffect = mean(log(estimates$hrIndirect) >= estimates$mainLogLbDiff &
                                    log(estimates$hrIndirect) <= estimates$mainLogUbDiff),
    biasDirectEffect = mean(simulationSettings$yA - estimates$directLogHr),
    biasMediatorEffect = mean(simulationSettings$yA - estimates$mainLogHr),
    biasMainEffect = mean(log(estimates$hrMain) - estimates$mainLogHr),
    biasIndirectEffect = mean(log(estimates$hrIndirect) - estimates$mainLogDiff),
    mseDirectEffect = mean((simulationSettings$yA - estimates$directLogHr)^2),
    msesMediatorEffect = mean((simulationSettings$yA - estimates$mainLogHr)^2),
    mseMainEffect = mean((log(estimates$hrMain) - estimates$mainLogHr)^2),
    mseIndirectEffect = mean((log(estimates$hrIndirect) - estimates$mainLogDiff)^2),
    indirectType1Error = if_else(hasIndirectEffect, NA, mean(estimates$mainLogLbDiff > 0 | estimates$mainLogUbDiff < 0)),
    indirectType2Error = if_else(hasIndirectEffect, mean(estimates$mainLogLbDiff <= 0 & estimates$mainLogUbDiff >= 0), NA),
    nonEstimableFraction = mean(nonEstimableIdx)
  )
  simSettingsForOutput <- simulationSettings
  class(simSettingsForOutput) <- "list"
  if (is(simSettingsForOutput, "SimulationSettings")) {
    simSettingsForOutput$aX <- paste(simSettingsForOutput$aX, collapse = ", ")
    simSettingsForOutput$mX <- paste(simSettingsForOutput$mX, collapse = ", ")
    simSettingsForOutput$yX <- paste(simSettingsForOutput$yX, collapse = ", ")
  }
  results <- results %>% 
    bind_cols(as_tibble(modelSettings)) %>%
    bind_cols(as_tibble(simSettingsForOutput))
  return(results)
}
