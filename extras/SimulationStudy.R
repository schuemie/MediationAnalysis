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
print(length(ssList))
msList <- list()
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event")
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(bootstrapType = 'reduced bias-corrected'))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(bootstrapType = 'bias-corrected'))
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "model",
                                                    mediatorType = "time-to-event",
                                                    bootstrapSettings = createBootstrapSettings(bootstrapType = 'pivoted'))
runSetOfSimulations(folder = folder, 
                    simulationSettingsList = ssList, 
                    modelSettingsList = msList,
                    nSimulations = 1000,
                    maxCores = 25) 

results <- readr::read_csv(file.path(folder, "Results.csv"))
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
                                   mediatorType = "time-to-event")

runSetOfSimulations(folder = folder, 
                    simulationSettingsList = ssList, 
                    modelSettingsList = msList,
                    nSimulations = 1000,
                    maxCores = 25) 

prepareForShinyApp(folder)

