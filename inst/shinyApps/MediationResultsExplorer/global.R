
results <- readRDS("data/simulation.rds")
results$type <- as.factor(results$type)
results$metric <- as.factor(results$metric)
types <- unique(results$type)
metrics <- unique(results$metric)
simParams <- colnames(results)[!(colnames(results) %in% c("type", "metric", "value"))]


