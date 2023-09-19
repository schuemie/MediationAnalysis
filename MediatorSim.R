# A simulation of a time-to-event mediator for a time-to-event outcome
# Not using Cyclops because of https://github.com/OHDSI/Cyclops/issues/69

library(dplyr)
library(survival)
library(splines)
library(ParallelLogger)
library(CohortMethod)

# x: covariates, a: treatment, m: mediator, y: outcome
# n: count, p: probability, h: hazard, t: time

# Simulation parameters ----------------------------------------------
settings = list(
  n = 1000,
  nX = 8,
  
  # Censor model parameters
  hCensor = 0.1,
  
  # Exposure model parameters
  aIntercept = log(0.5),
  aX = c( 1, -1, 1, -1, 0, 0, 0, 0),
  
  # Mediator model parameters
  mIntercept = log(0.1),
  mX = c(0, 0, 0.5, -1, 1.5, -1, 0, 0),
  mA = log(2),
  
  # Outcome hazard model parameters
  yIntercept = log(0.05),
  yX = c(0.5, 0, 0.5, 0, -1, 0, 0, 0),
  yA = log(0.5),
  yM = log(1.5)
)

# Simulation ----------------------------------------------------------
# seed = 0
runOneSimulation <- function(seed, settings, psAdjustment = "matching", mAdjustment = "mrs") {
  set.seed(seed)
  logistic <- function(x) {
    return(1/(1+exp(-x)))
  }
  
  x <- matrix(runif(settings$n * settings$nX), ncol = settings$nX)
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
  
  # Survival analysis ---------------------------------------------------
  # Perfect oracle: propensity score (PS) and mediator risk score (MRS) known:
  ps <- pA
  mrs <- hM
  
  # Sanity checks:
  # plot(ps, mrs)
  # writeLines(sprintf("Percent with the outcome: %%%0.1f", 100 * mean(tY < tCensor)))
  # writeLines(sprintf("Percent with the mediator: %%%0.1f", 100 * mean(tM < t)))
  # writeLines(sprintf("Percent with the outcome after the mediator: %%%0.1f", 100 * mean(tM < t & t == tY)))
  # writeLines(sprintf("Percent time after the mediator: %%%0.1f", 100 * sum(pmax(0, t - tM)) / sum(t)))
  
  data <- tibble(
    a = a,
    tStart = 0,
    tEnd = t,
    m = tM < t,
    y = tY == t,
    ps = ps,
    mrs = mrs,
    tM = tM
  )
  
  f <- formula(Surv(tStart, tEnd, y) ~ a)
  
  if (mAdjustment == "mrs") {
    f <- update(f, ~ . + ns(mrs, 5)) 
  } else if (mAdjustment == "covariates") {
    colnames(x) <- paste0("x", seq_len(settings$nX))
    data <- bind_cols(data, x)
    newF <- as.formula(paste("~ . + ", paste(colnames(x), collapse = " + ")))
    f <- update(f, newF) 
  } else if (mAdjustment == "none") {
    # Do nothing
  } else {
    stop("Unknown mAdjustment: ", mAdjustment) 
  }
  if (psAdjustment == "matching") {
    data <- data %>%
      mutate(propensityScore = ps,
             treatment = a,
             rowId = row_number()) %>%
      matchOnPs(maxRatio = 1) %>%
      select(-"propensityScore", -"treatment", -"rowId", -"stratumId")
  } else if (psAdjustment == "model") {
    f <- update(f, ~ . + ns(ps, 5))
  } else if (psAdjustment == "none") {
    # Do nothing
  } else {
    stop("Unknown psAdjustment: ", psAdjustment) 
  }
  
  # Split time before / after mediator
  data1 <- data %>%
    filter(!m)
  data2 <- data %>%
    filter(m) %>%
    mutate(tEnd = tM,
           m = FALSE,
           y = FALSE)
  data3 <- data %>%
    filter(m) %>%
    mutate(tStart = tM)
  data <- bind_rows(data1, data2, data3)
  
  # With mediator:
  fit <- coxph(update(f, ~ . + m), data = data)
  ci <- confint(fit)
  coverageA1 <- settings$yA >= ci["aTRUE",1] & settings$yA <= ci["aTRUE",2]
  coverageM1 <- settings$yM >= ci["mTRUE",1] & settings$yM <= ci["mTRUE",2]
  e <- coef(fit)
  eYa1 <- e["aTRUE"]
  eYm1 <- e["mTRUE"]
  
  # Without mediator:
  fit <- coxph(f, data = data)
  ci <- confint(fit)
  coverageA2 <- settings$yA >= ci["aTRUE",1] & settings$yA <= ci["aTRUE",2]
  e <- coef(fit)
  eYa2 <- e["aTRUE"]
  
  tibble(coverageA1 = coverageA1,
         coverageM1 = coverageM1,
         eYa1 = eYa1,
         eYm1 = eYm1,
         coverageA2 = coverageA2,
         eYa2) %>%
    return()
}

cluster <- makeCluster(20)
clusterRequire(cluster, "survival")
clusterRequire(cluster, "splines")
clusterRequire(cluster, "dplyr")
clusterRequire(cluster, "CohortMethod")

# CoverageA1: Coverage of the CI for the main effect when including the mediator in the model
# CoverageM1: Coverage of the CI for the mediator effect when including the mediator in the model
# eYa1: Mean estimate of the log main effect when including the mediator in the model
# eYm1: Mean estimate of the log mediator effect when including the mediator in the model
# CoverageA2: Coverage of the CI for the main effect when *not* including the mediator in the model
# eYa2: Mean estimate of the log main effect when *not* including the mediator in the model

clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "matching", mAdjustment = "mrs") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.7160000  0.9410000 -0.4699055  0.3905300  0.7080000 -0.4673302 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "model", mAdjustment = "mrs") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.4820000  0.9390000 -0.3963065  0.3979469  0.4790000 -0.3940503 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "none", mAdjustment = "mrs") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.4130000  0.9360000 -0.3666600  0.3943466  0.4080000 -0.3645749 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "matching", mAdjustment = "covariates") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9380000  0.9390000 -0.7041261  0.4001168  0.9290000 -0.6356787 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "model", mAdjustment = "covariates") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9260000  0.9350000 -0.7022992  0.4034742  0.9090000 -0.6333219 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "none", mAdjustment = "covariates") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9240000  0.9380000 -0.6996240  0.4026938  0.9070000 -0.6307743 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "matching", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9390000  0.8320000 -0.6502500  0.2489746  0.9110000 -0.6084248 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "model", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9340000  0.7150000 -0.6511412  0.2242072  0.8980000 -0.6137354 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "none", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.8300000  0.8000000 -0.5721684  0.2605790  0.7090000 -0.5216629 

# Without confounding:
settingsNoConfounding <- settings 
settingsNoConfounding$aX <- rep(0, settings$nX)
settingsNoConfounding$mX <- rep(0, settings$nX)
settingsNoConfounding$yX <- rep(0, settings$nX)
clusterApply(cluster, 1:1000, runOneSimulation, settings = settingsNoConfounding, psAdjustment = "matching", mAdjustment = "mrs") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9390000  0.9430000 -0.7011284  0.4114851  0.9250000 -0.6243665 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settingsNoConfounding, psAdjustment = "none", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9380000  0.9490000 -0.7017348  0.4071058  0.9010000 -0.6260749 

settings$yA
settings$yM

stopCluster(cluster)
