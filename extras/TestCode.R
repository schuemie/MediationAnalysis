# Evaluating using Cyclops instead of coxph for outcome model fitting:

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

system.time({
  fit1 <- coxph(update(f, ~ . + m), data = data)
  fit1
})
# user  system elapsed 
# 0.008   0.001   0.009 
