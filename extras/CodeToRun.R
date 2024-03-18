library(MediationAnalysis)

folder <- "Simulation"
ssList <- list()
msList <- list()

ssList[[length(ssList) + 1]] <- createSimulationSettings()
ssList[[length(ssList) + 1]] <- createSimulationSettings(mA = log(4))
ssList[[length(ssList) + 1]] <- createSimulationSettings(mA = log(0.5))

msList[[length(msList) + 1]] <- createModelsettings() 
# msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
#                                                     mrs = "fit",
#                                                     mediatorType = "binary") 

runSetOfSimulations(folder = folder, 
                    simulationSettingsList = ssList, 
                    modelSettingsList = msList,
                    nSimulations = 1000,
                    maxCores = 10) 


