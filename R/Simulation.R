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

#' Create simulation settings
#'
#' @param n          Number of subject to simulate.
#' @param nX         Number of covariates to simulate.
#' @param hCensor    Hazard for random censoring.
#' @param aIntercept Intercept for probability of treatment (log scale).  
#' @param aX         Coefficients for probability of treatment for the covariates.
#'                   (log scale).
#' @param mIntercept Intercept for hazard of mediator (log scale)
#' @param mX         Coefficients for hazard of mediator for the covariates.
#'                   (log scale).
#' @param mA         Coefficient for hazard of mediator for the treatment. (log
#'                   scale).
#' @param yIntercept Intercept for hazard of outcome. (log scale)
#' @param yX         Coefficients for hazard of outcome for the covariates. (log
#'                   scale).
#' @param yA         Coefficient for the hazard of the outcome for treatment. 
#'                   (log scale).
#' @param yM         Coefficeint for the hazard of the outcome for the mediator.
#'                   (log scale).
#'
#' @return A simulation settings object.
#' 
#' @export
createSimulationSettings <- function(n = 1000,
                                     nX = 8,
                                     hCensor = 0.1,
                                     aIntercept = log(0.5),
                                     aX = c( 1, -1, 1, -1, 0, 0, 0, 0),
                                     mIntercept = log(0.1),
                                     mX = c(0, 0, 0.5, -1, 1.5, -1, 0, 0),
                                     mA = log(2),
                                     yIntercept = log(0.05),
                                     yX = c(0.5, 0, 0.5, 0, -1, 0, 0, 0),
                                     yA = log(0.5),
                                     yM = log(1.5)) {
  settings <- list()
  for (name in names(formals("createSimulationSettings"))) {
    settings[[name]] <- get(name)
  }
  return(settings)
}

logistic <- function(x) {
  return(1/(1+exp(-x)))
}


#' Simulate data
#'
#' @param settings A simulation settings object as created using `createSimulationSettings()`.
#'
#' @return
#' A tibble with one row per subject and the following columns:
#' 
#' - a: Was the patient treated (i.e. in the target group)?
#' - tStart: Start of follow-up time (currently always 0).
#' - tEnd: End of follow-up time (time to either outcome or random censor)
#' - m: Did the patient experience the mediator during follow-up?
#' - y: Did the patient experience the outcome at the end of follow-up?
#' - tM: Time to mediator (truncated to t).
#' - tY: Time to outcome (truncated to t).
#' - pA: Probability of treatment.
#' - hM_star: The hazard of the mediator, based on the counterfactual model (a
#'   model without the exposure status).
#' - x1 ... xn: The covariates.   
#' 
#' @export
simulateData <- function(settings) {
  x <- matrix(runif(settings$n * settings$nX), ncol = settings$nX)
  colnames(x) <- paste0("x", seq_len(settings$nX))
  pA <- logistic(settings$aIntercept + x %*% settings$aX)[, 1]
  a <- runif(settings$n) < pA
  hM <- exp(settings$mIntercept + x %*% settings$mX + a * settings$mA)[, 1]
  hY_M <- exp(settings$yIntercept + x %*% settings$yX + a * settings$yA)[, 1]
  hYM <- exp(settings$yIntercept + x %*% settings$yX + a * settings$yA + 1 * settings$yM)[, 1]
  tCensor <- rexp(settings$n, settings$hCensor)
  tM <- rexp(settings$n, hM)
  tY <- rexp(settings$n, hY_M)
  idx <- tM < tY
  tY[idx] <- tM[idx] + rexp(sum(idx), hYM[idx])
  t <- pmin(tCensor, tY)
  hMstar <- exp(settings$mIntercept + x %*% settings$mX)[, 1]
  hrIndirect <- computeTrueIndirectEffect(hCensor = settings$hCensor, 
                                          hM = hM, 
                                          hMstar = hMstar, 
                                          hY_M = hY_M, 
                                          hYM = hYM, 
                                          a = a)
  data <- tibble(
    a = a,
    tStart = 0,
    tEnd = t,
    m = tM < t,
    y = tY == t,
    tM = pmin(tM, t),
    tY = pmin(tY, t),
    pA = pA,
    hMstar = hMstar,
    hrIndirect = hrIndirect
  ) %>%
    bind_cols(x)
  return(data)
}

computeTrueIndirectEffect <- function(hCensor, hM, hMstar, hY_M, hYM, a) {
  # Only exposed have (indirect) effect:
  idx <- a == 1
  hM <- hM[idx]
  hYM <- hYM[idx]
  hY_M <- hY_M[idx]
  hMstar <- hMstar[idx]
  nearInfinity <- 999
  hFull <- cumHazard(nearInfinity, hCensor, hM, hY_M, hYM) - cumHazard(0, hCensor, hM, hY_M, hYM)
  hNoMediation <- cumHazard(nearInfinity, hCensor, hMstar, hY_M, hYM) - cumHazard(0, hCensor, hMstar, hY_M, hYM)
  return(mean(hFull / hNoMediation))
}



cumHazard <- function(t, hCensor, hM, hY_M, hYM) {
  # Closed-form integral of hazard over time
  part1 <- ((hYM - hY_M) * exp(t*(-(hM + hCensor)))) / (hM + hCensor)
  part2 <- hYM * exp(-hCensor*t) / hCensor
  return(part1 - part2)
  
  # Test code
  # hazard <- function(t, hCensor, hM, hY_M, hYM) {
  #   return((exp(-hM*t)*hY_M + (1-exp(-hM*t))*hYM) * exp(-hCensor*t))
  # }
  # 
  # hCensor = settings$hCensor
  # hM = hM[1]
  # hYM = hYM[1]
  # hY_M = hY_M[1]
  # hMstar = hMstar[1]
  # integrate(hazard, 0, 2, hCensor = hCensor, hM = hM, hY_M = hY_M, hYM = hYM)
}