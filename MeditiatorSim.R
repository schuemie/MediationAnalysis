# A simulation of a time-to-event mediator for a time-to-event outcome

# Simulation parameters ----------------------------------------------
n <- 100000 # Number of subjects
pExposure <- 0.5
hCensor <- 0.1 # Hazard of censoring

# Outcome hazard model parameters
betaIntercept <- log(0.05)
betaExposure <- log(1.1)
betaMediator <- log(1.0) # Simulating no effect of mediator on outcome

# Mediator hazard model parameters
gammaIntercept <- log(0.05)
gammaExposure <- log(2.0)

# Simulation ----------------------------------------------------------
exposure <- runif(n) < pExposure 
tCensor <- rexp(n, hCensor)
hMediator <- exp(gammaIntercept + gammaExposure*exposure)
tMediator <- rexp(n, hMediator)
hOutcomeBeforeMediator <- exp(betaIntercept + betaExposure*exposure)
hOutcomeAfterMediator <- exp(betaIntercept + betaExposure*exposure + betaMediator)
tOutcome <- rexp(n, hOutcomeBeforeMediator)
idx <- tOutcome > tMediator
tOutcome[idx] <- tMediator[idx] + rexp(sum(idx), hOutcomeAfterMediator)

# Survival analysis ---------------------------------------------------
time <- pmin(tCensor, tOutcome)
outcome <- tOutcome < tCensor
mediator <- tMediator < time

library(Cyclops)
# Using Poisson regression censoring at outcome to implement survival analysis:
cyclopsData <- createCyclopsData(outcome ~ exposure, offset = log(time), modelType = "pr")
fit <- fitCyclopsModel(cyclopsData)
exp(coef(fit))
# (Intercept) exposureTRUE 
# 0.049726     1.084155 

cyclopsDataDirect <- createCyclopsData(outcome ~ exposure + mediator, offset = log(time), modelType = "pr")
fitDirect <- fitCyclopsModel(cyclopsDataDirect)
exp(coef(fitDirect))
# (Intercept) exposureTRUE mediatorTRUE 
# 0.06753598   1.29442749   0.39577760  # Would not have expected exposure coefficient to be different from above, since outcome does not depend on mediator

cyclopsDataMediator <- createCyclopsData(mediator ~ exposure, offset = log(time), modelType = "pr")
fitMediator <- fitCyclopsModel(cyclopsDataMediator)
exp(coef(fitMediator))
# (Intercept) exposureTRUE 
# 0.03741103   1.62996980 
