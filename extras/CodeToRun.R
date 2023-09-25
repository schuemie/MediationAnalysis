library(MediationAnalysis)

folder <- "c:/temp/MediationAnalysis"
simulationSettingsList <- list(createSimulationSettings())
modelSettingsList <- list(createModelsettings())

runSetOfSimulations(folder, 
                    simulationSettingsList, 
                    modelSettingsList,
                    nSimulations = 10,
                    maxCores = 4) 