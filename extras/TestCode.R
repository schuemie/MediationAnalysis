# Evaluating using Cyclops instead of coxph for outcome model fitting:

# Using Cyclops ----------------------------------------------------
system.time({
  bigData <- survival::survSplit(Surv(tStart, tEnd, y) ~ a + m + mrs + stratumId,
                                 data = data,
                                 cut = unique(
                                   data %>% filter(y == 1) %>% pull(tM)),
                                 episode = "id2") %>%
    mutate(tEnd = tEnd - tStart,
           tStart = 0)
  cd <- Cyclops::createCyclopsData(Surv(tStart, tEnd, y) ~ a + ns(mrs, 5) + m + strata(stratumId) + strata(id2), data = bigData, modelType = "cox")
  fit2 <- Cyclops::fitCyclopsModel(cd)
  coef(fit2)
})
# user  system elapsed 
# 1.449   0.022   1.479 

# Using coxph ---------------------------------------------------------
system.time({
  fit1 <- coxph(update(f, ~ . + m), data = data)
  fit1
})
# user  system elapsed 
# 0.008   0.001   0.009 

# Using agreg.fit ----------------------------------------------------
x <- cbind(data$a, ns(data$mrs, 5), data$m)
y <- Surv(data$tStart, data$tEnd, data$y)
system.time({
  fit <- agreg.fit(x, y, data$stratumId, control = coxph.control(), method = "efron", rownames = seq_len(nrow(x)),  init = rep(0, ncol(x)))
  fit$coefficients
})
# user  system elapsed 
# 0.001   0.000   0.001 

# Debug coverage issue for indirect effect ---------------------------------------------------------
library(MediationAnalysis)
folder <- tempfile()

dir.create(folder)
simulationSettings <- createAbstractSimulationSettings(
  n = 1000,
  confoundingAySd = 0.1,
  confoundingmYSd = 0.1,
  confoundingAymSd = 0.1,
  aIntercept = log(0.1),
  mIntercept = log(0.01),
  yIntercept = log(0.01),
  mA = log(0.8),
  yA = log(0.9),
  yM = log(1.0)
)
# modelSettings <- createModelsettings(psAdjustment = "none",
#                                      mrsAdjustment = "none")
modelSettings1 <- createModelsettings()
modelSettings2 <- createModelsettings(bootstrapSettings = createBootstrapSettings(bootstrapType = "reduced bias-corrected"))

runSetOfSimulations(folder = folder, 
                    simulationSettingsList = list(simulationSettings), 
                    modelSettingsList = list(modelSettings1, modelSettings2),
                    nSimulations = 1000,
                    maxCores = 10) 
results <- read.csv(file.path(folder, "Results.csv"))
results
unlink(folder, recursive = TRUE)

seed <- seed + 1
set.seed(seed)
data <- simulateData(simulationSettings)
# sum(data$m)
# sum(data$y)
# sum(data$m & data$y)
model <- fitModel(data, modelSettings)
writeLines(sprintf("Mediated proportion: %0.2f (%0.2f - %0.2f)", model$mediatedProportion, model$mediatedProportionLb, model$mediatedProportionUb))

settings <- modelSettings
librarlibrarlibrary(dplyr)
sampling <- "strata" # strata or person
bootstrapType <- "percentile" # percentile or pivoted
sprintf("CI: %0.2f-%0.2f, true HR: %0.2f", exp(model$indirectLogLb), exp(model$indirectLogUb), model$trueIndirectHr)

indirectLogHr
fit3<- coxph(Surv(tStart, tEnd, y) ~ a + ns(log(mrs), 5) + strata(stratumId) + m + m*a, data = data)
fit3

# Investigate invividual components of sample --------------------------------------
source("R/ModelFitting.R")
data <- simulateData(simulationSettings)
settings <- modelSettings
library(dplyr)

bootstrap <- sapply(seq_len(10000), singleBootstrapSample, x = x, y = y, stratumIds = stratumIds, uniqueStratumIds =uniqueStratumIds)  

# Main effect:
hist(bootstrap[1, ], breaks = 100)
median(bootstrap[1, ])
log(data$hrMain[1])

# Direct effect:
hist(bootstrap[2, ], breaks = 100)
median(bootstrap[2, ])
simulationSettings$yA

# Indirect effect:
hist(bootstrap[1, ] - bootstrap[2, ], breaks = 100)
median(bootstrap[1, ] - bootstrap[2, ])
log(data$hrIndirect[1])



