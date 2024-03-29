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
  confoundingAySd = 1,
  confoundingmYSd = 0.1,
  confoundingAymSd = 1,
  aIntercept = log(0.5),
  mIntercept = log(0.01),
  yIntercept = log(0.01),
  mA = log(2),
  yA = log(0.5),
  yM = log(2)
)
modelSettings <- createModelsettings()
runSetOfSimulations(folder = folder, 
                    simulationSettingsList = list(simulationSettings), 
                    modelSettingsList = list(modelSettings),
                    nSimulations = 10,
                    maxCores = 10) 
results <- read.csv(file.path(folder, "Results.csv"))
results
unlink(folder, recursive = TRUE)
