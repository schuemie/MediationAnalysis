# A simulation of a time-to-event mediator for a time-to-event outcome

library(dplyr)
library(survival)
library(Cyclops)
library(splines)

# x: covariates, a: treatment, m: mediator, y: outcome
# n: count, p: probability, h: hazard, t: time

# Simulation parameters ----------------------------------------------
n <- 10000 # Number of subjects
nX <- 10 # Number of covariates

# Censor model parameters
hCensor <- 0.1 # Hazard of censoring

# Exposure model parameters
aIntercept <- log(0.5)
aX <- c( 1, -1, 1, -1, 0, 0, 0, 0, 0, 0)

# Mediator model parameters
mIntercept <- log(0.1)
mX <- c(0, 0, 0.5, -1, 1.5, -1, 0, 0, 0 ,0)
mA <- log(2)

# Outcome hazard model parameters
yIntercept <- log(0.05)
yX <- c(0.5, 0, 0.5, 0, -1, 0, 0, 0, 0, 0)
yA <- log(2)
yM <- log(1.5)

# Simulation ----------------------------------------------------------
logistic <- function(x) {
  return(1/(1+exp(-x)))
}

x <- matrix(runif(n * nX) < 0.5, ncol = nX)
pA <- logistic(aIntercept + x %*% aX)[, 1]
a <- runif(n) < pA
hM <- exp(mIntercept + x %*% mX + a * mA)[, 1]
hY_M <- exp(yIntercept + x %*% yX + a * yA)[, 1]
hYM <- exp(yIntercept + x %*% yX + a * yA + 1 * yM)[, 1]
tCensor <- rexp(n, hCensor)
tM <- rexp(n, hM)
tY <- rexp(n, hY_M)
idx <- tM < tY
tY[idx] <- tM[idx] + rexp(sum(idx), hYM[idx])

# Survival analysis ---------------------------------------------------
# Perfect oracle: propensity score and mediator risk score known:
ps <- pA
mrs <- hM

t <- pmin(tCensor, tY)
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

# Adjusted:
cyclopsData <- createCyclopsData(Surv(tStart, tEnd, y) ~ a + m + ps + mrs, 
                                 data = data,
                                 modelType = "cox")
fit <- fitCyclopsModel(cyclopsData)
exp(confint(fit, parm = "aTRUE"))
exp(confint(fit, parm = "mTRUE"))

# Unadjusted:
cyclopsData <- createCyclopsData(Surv(tStart, tEnd, y) ~ a + m, 
                                 data = data,
                                 modelType = "cox")
fit <- fitCyclopsModel(cyclopsData)
exp(confint(fit, parm = "aTRUE"))
exp(confint(fit, parm = "mTRUE"))
