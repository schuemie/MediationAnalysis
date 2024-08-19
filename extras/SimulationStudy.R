# Code to execute the full simulation study

library(MediationAnalysis)

# Pilot simulation study -------------------------------------------------------
folder <- "PilotSimulation"
ssList <- list()
for (yA in log(c(0.5))) {
  for (mA in log(c(2.0))) {
    for (yM in log(c(2.0))) {
      for (confounding in c(0.1, 1)) {
        for (aIntercept in log(c(0.5))) {
          for (mIntercept in log(c(0.01, 0.1))) {
            for (yIntercept in log(c(0.01, 0.1))) {
              ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings(
                aIntercept = aIntercept,
                confoundingAySd = confounding,
                confoundingmYSd = confounding,
                confoundingAymSd = confounding,
                mIntercept = mIntercept,
                yIntercept = yIntercept,
                yA = yA,
                yM = yM,
                mA = mA
              )
            }              
          }              
        }           
      }
    }
  }
}

msList <- list()
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "strata"))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "strata",
                                                                                                bootstrapType = 'reduced bias-corrected'))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "strata",
                                                                                                bootstrapType = 'bias-corrected'))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "strata",
                                                                                                bootstrapType = 'pivoted'))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "person"))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "person",
                                                                                                bootstrapType = 'reduced bias-corrected'))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "person",
                                                                                                bootstrapType = 'bias-corrected'))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(sampling = "person",
                                                                                                bootstrapType = 'pivoted'))
print(length(ssList))
print(length(msList))
runSetOfSimulations(folder = folder, 
                    simulationSettingsList = ssList, 
                    modelSettingsList = msList,
                    nSimulations = 1000,
                    maxCores = 25) 

results <- readr::read_csv(file.path(folder, "Results.csv"))
library(dplyr)
results |>
  group_by(sampling, bootstrapType) |>
  summarise(coverageIndirectEffect = mean((0.95-coverageIndirectEffect)^2), 
            coverageMediatedProportion = mean((0.95-coverageMediatedProportion)^2)) |>
  arrange(coverageIndirectEffect)
prepareForShinyApp(folder)

# Full simulation study --------------------------------------------------------
folder <- "Simulation"
ssList <- list()
for (yA in log(c(0.5, 1.0, 2.0))) {
  for (mA in log(c(0.5, 1.0, 2.0))) {
    for (yM in log(c(0.5, 1.0, 2.0))) {
      for (confounding in c(0.1, 1)) {
        for (aIntercept in log(c(0.1, 0.5, 0.9))) {
          for (mIntercept in log(c(0.01, 0.1))) {
            for (yIntercept in log(c(0.01, 0.1))) {
              ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings(
                aIntercept = aIntercept,
                confoundingAySd = confounding,
                confoundingmYSd = confounding,
                confoundingAymSd = confounding,
                mIntercept = mIntercept,
                yIntercept = yIntercept,
                yA = yA,
                yM = yM,
                mA = mA
              )
            }              
          }              
        }           
      }
    }
  }
}
print(length(ssList))

msList <- list()
msList[[1]] <- createModelsettings(ps = "fit",
                                   mrs = "fit",
                                   psAdjustment = "matching",
                                   mrsAdjustment = "model",
                                   mediatorType = "time-to-event",
                                   bootstrapSettings = createBootstrapSettings(sampling = "person",
                                                                               bootstrapType = 'reduced bias-corrected'))

runSetOfSimulations(folder = folder, 
                    simulationSettingsList = ssList, 
                    modelSettingsList = msList,
                    nSimulations = 1000,
                    maxCores = 25) 

prepareForShinyApp(folder)

