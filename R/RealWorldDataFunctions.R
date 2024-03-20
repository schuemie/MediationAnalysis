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

#' Create the mediator risk score
#'
#' @param cohortMethodData An object of type `CohortMethodData` as generated using `CohortMethod::getDbCohortMethodData()`.
#' @param mediatorId       The outcome ID corresponding to the mediator.
#' @param prior            The prior used to fit the model. See `Cyclops::createPrior()` for details.
#' @param control	         The control object used to control the cross-validation used to determine the hyperparameters 
#'                         of the prior (if applicable). See `Cyclops::createControl()` for details.
#' @return
#' A tibble with the mediator risk score (MRS).
#' 
#' @export
createMediatorRiskScore <- function(cohortMethodData, 
                                    mediatorId,
                                    removeDuplicateSubjects = "keep all",
                                    riskWindowStart = 0,
                                    startAnchor = "cohort start",
                                    riskWindowEnd = 0,
                                    endAnchor = "cohort end",
                                    prior = createPrior("laplace", exclude = c(0), useCrossValidation = TRUE),
                                    control = createControl(noiseLevel = "silent", 
                                                            cvType = "auto", 
                                                            seed = 1,
                                                            resetCoefficients = TRUE, 
                                                            tolerance = 2e-07, 
                                                            cvRepetitions = 1, 
                                                            fold = 10,
                                                            startingVariance = 0.01)) {
  covariateData <- FeatureExtraction::tidyCovariateData(cohortMethodData)
  studyPop <- CohortMethod::createStudyPopulation(
    cohortMethodData = cohortMethodData,
    outcomeId = mediatorId,
    removeDuplicateSubjects = removeDuplicateSubjects,
    removeSubjectsWithPriorOutcome = FALSE,
    riskWindowStart = riskWindowStart,
    startAnchor = startAnchor,
    riskWindowEnd = riskWindowEnd,
    endAnchor = endAnchor
  )
  covariateData$outcomes <- studyPop %>%
    mutate(y = outcomeCount > 0) %>%
    select("rowId", time = "survivalTime", "y")
  cyclopsData <- Cyclops::convertToCyclopsData(
    outcomes = covariateData$outcomes,
    covariates = covariateData$covariates,
    modelType = "pr"
  )
  fit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control) 
  covariateData$outcomes <- cohortMethodData$cohorts %>%
    mutate(time = 1)
  mrs <- predict(fit, 
                 newOutcomes = covariateData$outcomes,
                 newCovariates = covariateData$covariates)
  
  # daysToMediator <- cohortMethodData$outcomes %>%
  #   filter(.data$outcomeId == mediatorId, .data$daysToEvent >= 0) %>%
  #   group_by(.data$rowId) %>%
  #   summarise(daysToMediator = min(.data$daysToEvent, na.rm = TRUE)) %>%
  #   collect()
  mrs <- covariateData$outcomes %>%
    select("rowId", "treatment") %>%
    collect() %>%
    mutate(mrs = mrs) %>%
    left_join(studyPop %>%
                filter(.data$outcomeCount > 0) %>%
                select("rowId", daysToMediator = "survivalTime"),
              by = join_by("rowId"))
  return(mrs)
}

#' Plot the mediator risk score
#'
#' @param mrs             A tibble containing the MRS as created by `createMrs()`.
#' @param targetLabel     A label to us for the target cohort.
#' @param comparatorLabel A label to us for the comparator cohort.
#'
#' @return
#' A ggplot object. Use the `ggplot2::ggsave()` function to save to file in a different format.
#' 
#' @export
plotMrs <- function(mrs,
                    targetLabel = "Target",
                    comparatorLabel = "Comparator",
                    showFraction = TRUE) {
  mrs <- mrs %>%
    mutate(label = if_else(treatment == 1, targetLabel, comparatorLabel))
  mrs$label <- factor(mrs$label, levels = c(targetLabel, comparatorLabel))
  plot <- ggplot2::ggplot(mrs, ggplot2::aes(x = mrs)) +
    ggplot2::geom_density(ggplot2::aes(color = label, group = label, fill = label)) +
    ggplot2::scale_fill_manual(values = c(
      rgb(0.8, 0, 0, alpha = 0.5),
      rgb(0, 0, 0.8, alpha = 0.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      rgb(0.8, 0, 0, alpha = 0.5),
      rgb(0, 0, 0.8, alpha = 0.5)
    )) +
    ggplot2::scale_x_log10("Mediator risk score") +
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "top")
  if (showFraction) {
    labelData <- data.frame(text = sprintf("%0.2f%% have the mediator", 100*mean(!is.na(mrs$daysToMediator))))
    y <- max(density(log10(mrs$mrs[mrs$treatment == 1]))$y, density(log10(mrs$mrs[mrs$treatment == 0]))$y)
    # Weird bug in ggplot2: log transform is not applied to geom_label x coordinates:
    plot <- plot + ggplot2::geom_label(x = log10(max(mrs$mrs)), y = y, hjust = "right", vjust = "top", alpha = 0.8, ggplot2::aes(label = text), data = labelData, size = 3.5)
  }
  # plot
  return(plot)
}

fitMediatorModel <- function(studyPopulation, 
                             ps, 
                             mrs,  
                             psAdjustment = "matching", 
                             mrsAdjustment = "model", 
                             mediatorType = "time-to-event") {
  
}


