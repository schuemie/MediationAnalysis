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
    dir.create(folder)
  }
  cluster <- ParallelLogger::makeCluster(maxCores)
  # ParallelLogger::clusterRequire(cluster, "MediationAnalysis")
  for (simulationSettings in simulationSettingsList) {
    for (modelSettings in modelSettingsList) {
      fileName <- file.path(folder, 
                            sprintf("result_%s.rds", 
                                    digest::digest(list(modelSettings, simulationSettings))))
      if (!file.exists(fileName)) {
        results <- ParallelLogger::clusterApply(cluster, 
                                                seq_len(nSimulations), 
                                                runOneSimulation, 
                                                simulationSettings = simulationSettings,
                                                modelSettings = modelSettings)
        results <- results %>%
          bind_rows()
        summaryResults <- tibble(
          coverageMainEffect = mean(simulationSettings$yA >= results$mainLogLb & 
                                      simulationSettings$yA <= results$mainLogUb),
          coverageMediatorEffect = mean(simulationSettings$yM >= results$mediatorLogLb & 
                                          simulationSettings$yM <= results$mediatorLogUb),
          coverageMainEffectNoM = mean(simulationSettings$yA >= results$mainLogLb & 
                                         simulationSettings$yA <= results$mainLogUb),
          meanMainEffect = mean(results$mainLogHr),
          meanMediatorEffect = mean(results$mediatorLogHr),
          meanMainEffectNoM = mean(results$mainLogHrNoM),
        )
        simulationSettings$aX <- paste(simulationSettings$aX, collapse = ", ")
        simulationSettings$mX <- paste(simulationSettings$mX, collapse = ", ")
        simulationSettings$yX <- paste(simulationSettings$yX, collapse = ", ")
        summaryResults <- summaryResults %>% 
          bind_cols(as_tibble(modelSettings)) %>%
          bind_cols(as_tibble(simulationSettings))
        saveRDS(summaryResults, fileName)
      }
    }
  }
  ParallelLogger::stopCluster(cluster)
}

runOneSimulation <- function(seed, simulationSettings, modelSettings) {
  set.seed(seed)
  data <- simulateData(simulationSettings)
  estimates <- fitModel(data, modelSettings)
  return(estimates)
}
