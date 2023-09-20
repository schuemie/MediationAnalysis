# A simulation of a time-to-event mediator for a time-to-event outcome
# Not using Cyclops because of https://github.com/OHDSI/Cyclops/issues/69

library(dplyr)
library(survival)
library(splines)
library(ParallelLogger)
library(CohortMethod)
library(MatchIt)

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
runOneSimulation <- function(seed, settings, psAdjustment = "matching", mAdjustment = "mrs", mType = "time-to-event") {
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
  # mrs <- hM # This won't work. Need counterfactual MRS:
  mrs <- exp(settings$mIntercept + x %*% settings$mX)[, 1]
  
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
  } else if (mAdjustment == "2D matching") {
    if (psAdjustment != "2D matching") {
      stop("If mAdjustment = '2D matching' then psAdjustment should also be '2D matching'")
    }
  }  else if (mAdjustment == "none") {
    # Do nothing
  } else {
    stop("Unknown mAdjustment: ", mAdjustment) 
  }
  if (psAdjustment == "matching") {
    data <- data %>%
      mutate(propensityScore = ps,
             treatment = a,
             rowId = row_number()) %>%
      matchOnPs(maxRatio = 100) %>%
      select(-"propensityScore", -"treatment", -"rowId")
    f <- update(f, ~ . + strata(stratumId))
  } else if (psAdjustment == "model") {
    f <- update(f, ~ . + ns(ps, 5))
  } else if (psAdjustment == "covariates") {
    if (mAdjustment != "covariates") {
      stop("If mAdjustment = 'covariates' then psAdjustment should also be 'covariates'")
    }
  } else if (psAdjustment == "2D matching") {
    matchit <- matchit(a ~ ps + mrs, 
                       method = "nearest", 
                       distance ="euclidean", 
                       caliper = c("ps" = 0.02, "mrs" = 0.02),
                       std.caliper = FALSE,
                       ratio = 100,
                       data = data)
    data <- match.data(matchit, data = data, subclass = "stratumId")
    f <- update(f, ~ . + strata(stratumId))
  } else if (psAdjustment == "none") {
    # Do nothing
  } else {
    stop("Unknown psAdjustment: ", psAdjustment) 
  }
  
  if (mType == "time-to-event") {
    # Split time before / after mediator
    data <- bind_rows(
      data %>%
        filter(!m),
      data %>%
        filter(m) %>%
        mutate(tEnd = tM,
               m = FALSE,
               y = FALSE),
      data %>%
        filter(m) %>%
        mutate(tStart = tM)
    )
  } else if (mType == "binary") {
    # Do nothing 
  } else {
    stop("Unknown mType: ", mType) 
  }
  # Cleanup: remove (near) zero-length intervals:
  data <- data %>%
    filter(tEnd - tStart > 0.0001)
  
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
clusterRequire(cluster, "MatchIt")

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
# 0.9350000  0.9400000 -0.7137519  0.4166025  0.9270000 -0.6449663 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "model", mAdjustment = "mrs") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9290000  0.9410000 -0.6900137  0.3992755  0.8980000 -0.6220887 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "none", mAdjustment = "mrs") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.8240000  0.9340000 -0.5693873  0.3949595  0.6590000 -0.5029238 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "covariates", mAdjustment = "covariates") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9240000  0.9380000 -0.6996478  0.4026710  0.9070000 -0.6308008 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "matching", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9260000  0.8710000 -0.6706329  0.2447659  0.9190000 -0.6291749 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "model", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9290000  0.7140000 -0.6523017  0.2247884  0.8970000 -0.6147374 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "none", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.8300000  0.8000000 -0.5721839  0.2605555  0.7090000 -0.5216825 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "2D matching", mAdjustment = "2D matching") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9390000  0.9350000 -0.6928786  0.3937783  0.9220000 -0.6273094 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "covariates", mAdjustment = "covariates", mType = "binary") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.5560000  0.0000000 -0.4498662 -0.8682919  0.9070000 -0.6307832 

# Without confounding:
settingsNoConfounding <- settings 
settingsNoConfounding$aX <- rep(0, settings$nX)
settingsNoConfounding$aX[1] <- 1
settingsNoConfounding$mX <- rep(0, settings$nX)
settingsNoConfounding$mX[2] <- 1
settingsNoConfounding$yX <- rep(0, settings$nX)
settingsNoConfounding$yX[3] <- 1
clusterApply(cluster, 1:1000, runOneSimulation, settings = settingsNoConfounding, psAdjustment = "matching", mAdjustment = "mrs") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9380000  0.9410000 -0.7067062  0.4068102  0.9090000 -0.6313603 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settingsNoConfounding, psAdjustment = "none", mAdjustment = "none") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.9420000  0.9360000 -0.6761984  0.3950845  0.8440000 -0.6056510 
clusterApply(cluster, 1:1000, runOneSimulation, settings = settings, psAdjustment = "none", mAdjustment = "none", mType = "binary") %>%
  bind_rows() %>%
  colMeans()
# coverageA1 coverageM1       eYa1       eYm1 coverageA2       eYa2 
# 0.2040000  0.0000000 -0.3233349 -0.8556581  0.7090000 -0.5216636 

settings$yA
settings$yM

stopCluster(cluster)


settings$aX <- rnorm(settings$nX)
settings$mX <- rnorm(settings$nX)
settings$yX <- rnorm(settings$nX)

