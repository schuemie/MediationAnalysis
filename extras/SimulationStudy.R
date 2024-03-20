# Code to execute the full simulation study

library(MediationAnalysis)

folder <- "Simulation"
ssList <- list()
for (yA in log(c(0.5, 1.0, 2.0))) {
  for (mA in log(c(0.5, 1.0, 2.0))) {
    for (yM in log(c(0.5, 1.0, 2.0))) {
      for (confoundingAySd in c(0.5, 2)) {
        for (confoundingmYSd in c(0.5, 2)) {
          for (confoundingAymSd in c(0.5, 2)) {
            for (aIntercept in log(c(0.1, 0.5, 0.9))) {
              for (mIntercept in log(c(0.01, 0.1))) {
                for (yIntercept in log(c(0.01, 0.1))) {
                  ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings(
                    aIntercept = aIntercept,
                    confoundingAySd = confoundingAySd,
                    confoundingmYSd = confoundingmYSd,
                    confoundingAymSd = confoundingAymSd,
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
                                                    mediatorType = "binary")
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    psAdjustment = "matching",
                                                    mrsAdjustment = "none",
                                                    mediatorType = "time-to-event")

runSetOfSimulations(folder = folder, 
                    simulationSettingsList = ssList, 
                    modelSettingsList = msList,
                    nSimulations = 1000,
                    maxCores = 25) 

