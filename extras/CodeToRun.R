library(MediationAnalysis)

folder <- "Simulation"
ssList <- list()
msList <- list()

ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings()
ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings(mA = log(0.5))
ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings(mA = log(0.5), confoundingAymSd = 2)
ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings(yA = log(2))
ssList[[length(ssList) + 1]] <- createAbstractSimulationSettings(yA = log(2), mA = log(0.5))

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
                    nSimulations = 100,
                    maxCores = 10) 


