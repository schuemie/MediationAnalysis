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
#' @param removeDuplicateSubjects          Remove subjects that are in both the target and comparator
#'                                         cohort? See details for allowed values.
#' @param riskWindowStart                  The start of the risk window (in days) relative to the `startAnchor`.
#' @param startAnchor                      The anchor point for the start of the risk window. Can be `"cohort start"`
#'                                         or `"cohort end"`.
#' @param riskWindowEnd                    The end of the risk window (in days) relative to the `endAnchor`.
#' @param endAnchor                        The anchor point for the end of the risk window. Can be `"cohort start"`
#'                                         or `"cohort end"`.
#' @param prior            The prior used to fit the model. See `Cyclops::createPrior()` for details.
#' @param control	         The control object used to control the cross-validation used to determine the hyperparameters 
#'                         of the prior (if applicable). See `Cyclops::createControl()` for details.
#'                         
#' @details
#' The `removeduplicateSubjects` argument can have one of the following values:
#'
#' - `"keep all"`: Do not remove subjects that appear in both target and comparator cohort
#' - `"keep first"`: When a subjects appear in both target and comparator cohort, only keep whichever cohort is first in time. If both cohorts start simultaneous, the person is removed from the analysis.
#' - `"remove all"`: Remove subjects that appear in both target and comparator cohort completely from the analysis."
#' 
#' @return
#' A tibble with the mediator risk score (MRS).
#' 
#' @export
createMediatorRiskScore <- function(cohortMethodData, 
                                    mediatorId,
                                    removeDuplicateSubjects = "keep all",
                                    censorAtNewRiskWindow = FALSE,
                                    removeSubjectsWithPriorOutcome = TRUE,
                                    priorOutcomeLookback = 999999,
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
    censorAtNewRiskWindow = censorAtNewRiskWindow,
    removeSubjectsWithPriorOutcome = removeSubjectsWithPriorOutcome,
    priorOutcomeLookback = priorOutcomeLookback,
    riskWindowStart = riskWindowStart,
    startAnchor = startAnchor,
    riskWindowEnd = riskWindowEnd,
    endAnchor = endAnchor
  )
  covariateData$outcomes <- studyPop %>%
    mutate(y = .data$outcomeCount > 0) %>%
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

#' Plot the mediator risk score by exposure
#'
#' @param mrs             A tibble containing the MRS as created by `createMrs()`.
#' @param targetLabel     A label to us for the target cohort.
#' @param comparatorLabel A label to us for the comparator cohort.
#' @param showFraction    Add a label to the plot showing what fraction of the population has the 
#'                        mediator during the time-at-risk?
#' @param title             Optional: the main title for the plot.
#' @param fileName        Name of the file where the plot should be saved, for example 'plot.png'.
#'                        see the function [ggplot2::ggsave()] for supported file formats.
#'
#' @return
#' A ggplot object. Use the `ggplot2::ggsave()` function to save to file in a different format.
#' 
#' @export
plotMrsByExposure <- function(mrs,
                              targetLabel = "Target",
                              comparatorLabel = "Comparator",
                              showFraction = TRUE,
                              title = NULL,
                              fileName = NULL) {
  mrs <- mrs %>%
    mutate(label = if_else(.data$treatment == 1, targetLabel, comparatorLabel))
  mrs$label <- factor(mrs$label, levels = c(targetLabel, comparatorLabel))
  plot <- ggplot2::ggplot(mrs, ggplot2::aes(x = .data$mrs)) +
    ggplot2::geom_density(ggplot2::aes(color = .data$label, group = .data$label, fill = .data$label)) +
    ggplot2::scale_fill_manual(values = c(
      rgb(0.8, 0, 0, alpha = 0.5),
      rgb(0, 0, 0.8, alpha = 0.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      rgb(0.8, 0, 0, alpha = 0.5),
      rgb(0, 0, 0.8, alpha = 0.5)
    )) +
    ggplot2::scale_x_log10("Mediator risk score", labels = scales::comma) +
    ggplot2::scale_y_continuous("Density") + 
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "top")
  if (showFraction) {
    labelData <- data.frame(text = sprintf("%0.2f%% have the mediator", 100*mean(!is.na(mrs$daysToMediator))))
    y <- max(density(log10(mrs$mrs[mrs$treatment == 1]))$y, density(log10(mrs$mrs[mrs$treatment == 0]))$y)
    # Weird bug in ggplot2: log transform is not applied to geom_label x coordinates:
    plot <- plot + ggplot2::geom_label(x = log10(max(mrs$mrs)), 
                                       y = y,
                                       hjust = "right", 
                                       vjust = "top", 
                                       alpha = 0.8, 
                                       ggplot2::aes(label = .data$text), 
                                       data = labelData, 
                                       size = 3.5)
  }
  # plot
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 5, height = 3.5, dpi = 400)
  }
  return(plot)
}

#' Plot the mediator risk score by mediator status
#'
#' @param mrs             A tibble containing the MRS as created by `createMrs()`.
#' @param targetLabel     A label to us for the target cohort.
#' @param comparatorLabel A label to us for the comparator cohort.
#' @param showFraction    Add a label to the plot showing what fraction of the population has the 
#'                        mediator during the time-at-risk?
#' @param title             Optional: the main title for the plot.
#' @param fileName        Name of the file where the plot should be saved, for example 'plot.png'.
#'                        see the function [ggplot2::ggsave()] for supported file formats.
#'
#' @return
#' A ggplot object. Use the `ggplot2::ggsave()` function to save to file in a different format.
#' 
#' @export
plotMrsByMediator <- function(mrs,
                              showFraction = TRUE,
                              title = NULL,
                              fileName = NULL) {
  mediatorLabel <- "With mediator"
  withoutMediatorLabel <- "Without mediator"
  
  mrs <- mrs %>%
    mutate(label = if_else(is.na(.data$daysToMediator), withoutMediatorLabel, mediatorLabel))
  mrs$label <- factor(mrs$label, levels = c(mediatorLabel, withoutMediatorLabel))
  plot <- ggplot2::ggplot(mrs, ggplot2::aes(x = .data$mrs)) +
    ggplot2::geom_density(ggplot2::aes(color = .data$label, group = .data$label, fill = .data$label)) +
    ggplot2::scale_fill_manual(values = c(
      rgb(0.4, 0.09, 0.30, alpha = 0.5),
      rgb(0.80, 0.70, 0.15, alpha = 0.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      rgb(0.4, 0.09, 0.30, alpha = 0.5),
      rgb(0.96, 0.82, 0.26, alpha = 0.5)
    )) +
    ggplot2::scale_x_log10("Mediator risk score", labels = scales::comma) +
    ggplot2::scale_y_continuous("Density") + 
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "top")
  plot
  if (showFraction) {
    labelData <- data.frame(text = sprintf("%0.2f%% have the mediator", 100*mean(!is.na(mrs$daysToMediator))))
    y <- max(density(log10(mrs$mrs[is.na(mrs$daysToMediator)]))$y, density(log10(mrs$mrs[!is.na(mrs$daysToMediator)]))$y)
    # Weird bug in ggplot2: log transform is not applied to geom_label x coordinates:
    plot <- plot + ggplot2::geom_label(x = log10(max(mrs$mrs)), 
                                       y = y,
                                       hjust = "right", 
                                       vjust = "top", 
                                       alpha = 0.8, 
                                       ggplot2::aes(label = .data$text), 
                                       data = labelData, 
                                       size = 3.5)
  }
  plot
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 5, height = 3.5, dpi = 400)
  }
  return(plot)
}

#' Fit mediator model using real-world data
#'
#' @param studyPopulation The study population as created using the `CohortMethod::createStudyPopulation()` 
#'                        function, with `outcomeId` specifying the outcome of interest.
#' @param ps              The propensity score object as created using the `CohortMethod::createPs()` function.
#' @param mrs             The mediator risk score object as created using the `createMrs()` function.
#' @param psAdjustment    Type of PS adjustment. Can be "matching", "model", 
#'                        "2D matching", or "none".
#' @param mrsAdjustment   Type of MRS adjustment. Can be "model", "2D matching", 
#'                        or "none".
#' @param mediatorType    Type of mediator in the model. Can be "time-to-event" or
#'                        "binary".
#'
#' @return
#' Returns a data frame with one row, and columns for the various model estimands.
#' 
#' @export
fitMediatorModel <- function(studyPopulation, 
                             ps, 
                             mrs,  
                             psAdjustment = "matching",
                             mrsAdjustment = "model",
                             mediatorType = "time-to-event",
                             bootstrapSettings = createBootstrapSettings()) {
  settings <- createModelsettings(
    psAdjustment = psAdjustment,
    mrsAdjustment = mrsAdjustment,
    mediatorType = mediatorType,
    bootstrapSettings = bootstrapSettings
  )
  data <- studyPopulation %>%
    inner_join(ps %>%
                 select("rowId", "propensityScore"),
               by = join_by("rowId")) %>%
    inner_join(mrs %>%
                 select("rowId", "mrs", "daysToMediator"),
               by = join_by("rowId")) %>%
    transmute(a = .data$treatment == 1,
              tStart = 0,
              tEnd = .data$survivalTime,
              m = !is.na(.data$daysToMediator) & .data$daysToMediator < .data$survivalTime,
              y = .data$outcomeCount > 0,
              tM = .data$daysToMediator,
              pA = .data$propensityScore,
              hMstar = .data$mrs,
              hrIndirect = NA,
              hrMain = NA,
              personSeqId = .data$rowId)
  model <- fitModel(data, settings)
  return(model)
}
