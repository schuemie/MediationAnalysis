library(MediationAnalysis)

folder <- "d:/temp/MediationAnalysis"
ssList <- list()
msList <- list()

ssList[[length(ssList) + 1]] <- createSimulationSettings()
msList[[length(msList) + 1]] <- createModelsettings() 
msList[[length(msList) + 1]] <- createModelsettings(ps = "fit",
                                                    mrs = "fit") 

runSetOfSimulations(folder = folder, 
                    simulationSettingsList = ssList, 
                    modelSettingsList = msList,
                    nSimulations = 1000,
                    maxCores = 20) 
