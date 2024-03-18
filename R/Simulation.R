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
  # hrIndirect <- computeTrueIndirectEffect(tM, t, hYM, hY_M, n)
  hrIndirect <- computeTrueIndirectEffect(settings$hCensor, hM, hY_M, hYM, settings$n)
  
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

# computeTrueIndirectEffect <- function(tM, t, hYM, hY_M, n) {
#   tM <- pmin(t, tM)
#   totalH <- sum(tM*hY_M + (t-tM) * hYM)
#   hrIndirect <- totalH / sum(hY_M*t)
#   return(hrIndirect)
# }

computeTrueIndirectEffect <- function(hCensor, hM, hY_M, hYM, n) {
  sample <- sapply(seq_len(1000),
                   computeTrueIndirectEffectOnce,
                   hCensor = hCensor,
                   hM = hM,
                   hY_M = hY_M,
                   hYM = hYM,
                   n = n)
  return(median(sample))
}

computeTrueIndirectEffectOnce <- function(dummy, hCensor, hM, hY_M, hYM, n) {
  tCensor <- rexp(n, hCensor)
  tM <- rexp(n, hM)
  tY <- rexp(n, hY_M)
  idx <- tM < tY
  tY[idx] <- tM[idx] + rexp(sum(idx), hYM[idx])
  t <- pmin(tCensor, tY)
  tM <- pmin(t, tM)
  totalH <- sum(tM*hY_M + (t-tM) * hYM)
  return(totalH / sum(hY_M*t))
}



# computeTrueIndirectEffect <- function(hM, hY_M, hYM, hC, n) {
#   # Well this was stupid: cannot use discrete time for continuous model
#   # p(M,t) = 1-Prod_i=1^t(1-hM)
#   # h*Y(t) = p(M,t)*hYm + (1-p(M,t))*hY_m
#   # p(Y,t) = 1-Prod_i=1^t(1-h*Y(i))
#   # p(c,t) = 1-Prod_i=1^t(1-hC)
#   # p(o,t) = 1-((1-p(Y,t)) * (1-p(c,t))
#   # hr(t) = h*Y(t) / hY_m
#   # hr = Sum_t=1^inf(hr(t)*p(o,t)) / Sum_t=1^inf(p(o,t)) 
#   sumHrW <- rep(0, n)
#   sumW <- rep(0, n)
#   p_M <- rep(1, n)
#   p_Y <- rep(1, n)
#   p_C <- rep(1, n)
#   for (t in seq_len(100)) {
#     p_M <- p_M * (1-hM)
#     pM <- 1-p_M  
#     hStarY <- pM * hYM + p_M * hY_M
#     p_Y <- p_Y * (1-hStarY)
#     p_C <- p_C * (1-hC)
#     pO <- p_Y * p_C
#     hr <- hStarY / hY_M
#     sumHrW <- sumHrW + hr * pO
#     sumW <- pO
#   }
#   hrIndirect <- sumHrW / sumW
#   
#   
# }
