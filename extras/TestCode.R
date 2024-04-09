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
  confoundingAySd = 0.1,
  confoundingmYSd = 0.1,
  confoundingAymSd = 0.1,
  aIntercept = log(0.5),
  mIntercept = log(0.01),
  yIntercept = log(0.1),
  mA = log(2),
  yA = log(0.5),
  yM = log(2)
)
modelSettings <- createModelsettings()
runSetOfSimulations(folder = folder, 
                    simulationSettingsList = list(simulationSettings), 
                    modelSettingsList = list(modelSettings),
                    nSimulations = 100,
                    maxCores = 10) 
results <- read.csv(file.path(folder, "Results.csv"))
results
unlink(folder, recursive = TRUE)


data <- simulateData(simulationSettings)
sum(data$m)
sum(data$y)
sum(data$m & data$y)
model <- fitModel(data, modelSettings)
settings <- modelSettings
sprintf("CI: %0.2f-%0.2f, true HR: %0.2f", exp(model$indirectLogLb), exp(model$indirectLogUb), model$trueIndirectHr)

indirectLogHr
fit3<- coxph(Surv(tStart, tEnd, y) ~ a + ns(log(mrs), 5) + strata(stratumId) + m + m*a, data = data)
fit3

# Investigate invividual components of sample --------------------------------------
singleBootstrapSample <- function(dummy, x, y, stratumIds, uniqueStratumIds) {
  if (is.null(stratumIds)) {
    idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
  } else {
    if (sampling == "strata") {
      # sampledStratumIds <- sample(uniqueStratumIds, size = length(uniqueStratumIds), replace = TRUE)
      sampledStratumIds <- sample(uniqueStratumIds$stratumId, 
                                  size = nrow(uniqueStratumIds), 
                                  prob = uniqueStratumIds$weight,
                                  replace = TRUE)
      idx <- inner_join(tibble(stratumId = sampledStratumIds), 
                        stratumIds, 
                        by = join_by("stratumId"), 
                        relationship = "many-to-many") %>%
        pull("idx")
      stratumIds <- stratumIds$stratumId[idx]
    } else {
      idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
      # idx <- c(
      #   sample(which(x[,7] == 1), sum(x[, 7] == 1), replace = TRUE),
      #   sample(which(x[,7] == 0), sum(x[, 7] == 0), replace = TRUE)
      # )
      # idx <- c(
      #   sample(which(y[, 3] == 1), sum(y[, 3] == 1), replace = TRUE),
      #   sample(which(y[, 3] == 0), sum(y[, 3] == 0), replace = TRUE)
      # )
      stratumIds <- stratumIds[idx]
    }
  }
  x <- x[idx, ]
  y <- y[idx, ]
  control <- coxph.control()
  tryCatch({
    suppressWarnings({
      fit1 <- agreg.fit(x, y, stratumIds, control = control, method = "efron", rownames = seq_along(idx), init = rep(0,ncol(x)))
      fit2 <- agreg.fit(x[, -ncol(x), drop = FALSE], y, stratumIds, control = control, method = "efron", rownames = seq_along(idx),  init = rep(0, ncol(x)-1))
    })
    return(tibble(direct = fit1$coefficients[1], main = fit2$coefficients[1]))
  },
  error = function(e) {
    return(NA)
  })
}

bootstrap <- lapply(seq_len(10000), singleBootstrapSample, x = x, y = y, stratumIds = stratumIds, uniqueStratumIds =uniqueStratumIds)  
bootstrap <- bind_rows(bootstrap)
hist(bootstrap$main, breaks = 100)
median(bootstrap$main)
log(data$hrMain[1])

hist(bootstrap$direct, breaks = 100)
median(bootstrap$direct)
simulationSettings$yA



