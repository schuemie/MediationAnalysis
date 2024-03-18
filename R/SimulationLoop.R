# Copyright 2023 Observational Health Data Sciences and Informatics
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
  nScenarios <- length(simulationSettingsList) * length(modelSettingsList)
  i <- 1
  # ParallelLogger::clusterRequire(cluster, "MediationAnalysis")
  # simulationSettings = simulationSettingsList[[1]]
  # modelSettings = modelSettingsList[[1]]
  for (simulationSettings in simulationSettingsList) {
    for (modelSettings in modelSettingsList) {
      message(sprintf("Running scenario %d of %d", i, nScenarios))
      fileName <- file.path(folder, 
                            sprintf("result_%s.rds", 
                                    digest::digest(list(modelSettings, simulationSettings))))
      if (file.exists(fileName)) {
        summaryResults <- readRDS(fileName)
      } else {
        results <- ParallelLogger::clusterApply(cluster, 
                                                seq_len(nSimulations), 
                                                MediationAnalysis:::runOneSimulation, 
                                                simulationSettings = simulationSettings,
                                                modelSettings = modelSettings)
        results <- results %>%
          bind_rows()
        hasIndirectEffect <- simulationSettings$mA != 0 & simulationSettings$yM != 0
        summaryResults <- tibble(
          coverageDirectEffect = mean(simulationSettings$yA >= results$directLogLb & 
                                      simulationSettings$yA <= results$directLogUb),
          coverageMediatorEffect = mean(simulationSettings$yM >= results$mediatorLogLb & 
                                          simulationSettings$yM <= results$mediatorLogUb),
          covarageMainEffect = mean(log(results$hrMain) >= results$mainLogLb &
                                          log(results$hrMain) <= results$mainLogUb),
          covarageIndirectEffect = mean(log(results$hrIndirect) >= results$mainLogLbDiff &
                                          log(results$hrIndirect) <= results$mainLogUbDiff),
          biasDirectEffect = mean(simulationSettings$yA - results$directLogHr),
          biasMediatorEffect = mean(simulationSettings$yA - results$mainLogHr),
          biasMainEffect = mean(log(results$hrMain) - results$mainLogHr),
          biasIndirectEffect = mean(log(results$hrIndirect) - results$mainLogDiff),
          mseDirectEffect = mean((simulationSettings$yA - results$directLogHr)^2),
          msesMediatorEffect = mean((simulationSettings$yA - results$mainLogHr)^2),
          mseMainEffect = mean((log(results$hrMain) - results$mainLogHr)^2),
          mseIndirectEffect = mean((log(results$hrIndirect) - results$mainLogDiff)^2),
          indirectType1Error = if_else(hasIndirectEffect, NA, mean(results$mainLogLbDiff > 0 | results$mainLogUbDiff < 0)),
          indirectType2Error = if_else(hasIndirectEffect, mean(results$mainLogLbDiff <= 0 & results$mainLogUbDiff >= 0), NA)
        )
        simSettingsForOutput <- simulationSettings
        simSettingsForOutput$aX <- paste(simSettingsForOutput$aX, collapse = ", ")
        simSettingsForOutput$mX <- paste(simSettingsForOutput$mX, collapse = ", ")
        simSettingsForOutput$yX <- paste(simSettingsForOutput$yX, collapse = ", ")
        summaryResults <- summaryResults %>% 
          bind_cols(as_tibble(modelSettings)) %>%
          bind_cols(as_tibble(simSettingsForOutput))
        saveRDS(summaryResults, fileName)
      }
      allResults[[i]] <- summaryResults
      i <- i + 1
    }
  }
  allResults <- bind_rows(allResults)
  readr::write_csv(allResults, file.path(folder, "Results.csv"))
}

runOneSimulation <- function(seed, simulationSettings, modelSettings) {
  set.seed(seed)
  data <- simulateData(simulationSettings)
  estimates <- fitModel(data, modelSettings)
  return(estimates)
}
