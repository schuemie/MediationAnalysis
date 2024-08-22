library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggh4x)

# Plots from real-world example ------------------------------------------------
source("RealWorldExample/DatabaseDetails.R")

target <- 16329 # DOACs
comparator <- 16330 # Warfarin
mediator <- 16853 # Extracranial Major bleeding IP-ER
outcome1 <- 2072 # AMI                              
outcome2 <- 2087 # Ischemic stroke                  

database = databases[[1]]
for (database in databases) {
  
  # Combined PS and MRS plot
  psFileName <- file.path(database$outputFolder, sprintf("ps_t%d_c%s.rds", target, comparator))
  ps <- readRDS(psFileName)
  mrsFileName <- file.path(database$outputFolder, sprintf("mrs_t%d_c%s_m%d.rds", target, comparator, mediator))  
  mrs <- readRDS(mrsFileName)
  psPlot <- CohortMethod::plotPs(ps,
                                 scale = "propensity",
                                 targetLabel = "DOACs",
                                 comparatorLabel = "Warfarin")
  psPlot <- psPlot + 
    # scale_fill_manual(values = c(alpha("#336B91", 0.7), alpha("#EB6622", 0.7))) +
    scale_fill_manual(values = c(alpha("#EB6622", 0.6), alpha("#336B91", 0.6))) +
    scale_color_manual(values = c(alpha("#EB6622", 0.7), alpha("#336B91", 0.7))) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())
  
  mrsPlot <- MediationAnalysis::plotMrsByMediator(mrs, showFraction = FALSE)
  mrsPlot <- mrsPlot +    
    scale_fill_manual(values = c(alpha("#6d904f", 0.6), alpha("#e5ae38", 0.6))) + 
    scale_color_manual(values = c(alpha("#6d904f", 0.7), alpha("#e5ae38", 0.7))) + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) 
  plot <- grid.arrange(psPlot, mrsPlot, ncol = 2)
  ggsave(file.path(rootFolder, sprintf("PsMrs_%s.png", database$databaseId)), plot, width = 6, height = 3, dpi = 300)
  ggsave(file.path(rootFolder, sprintf("PsMrs_%s.pdf", database$databaseId)), plot, width = 6, height = 3, dpi = 300)
  
  # Balance scatter plot
  balanceFileName <- file.path(database$outputFolder, sprintf("balance_t%d_c%s.rds", target, comparator))
  balance <- readRDS(balanceFileName)
  plot <- CohortMethod::plotCovariateBalanceScatterPlot(balance,
                                                        showCovariateCountLabel = TRUE,
                                                        threshold = 0.1,
                                                        title = "")
  ggsave(file.path(rootFolder, sprintf("Balance_%s.png", database$databaseId)), plot, width = 3.5, height = 3.5, dpi = 300)
  ggsave(file.path(rootFolder, sprintf("Balance_%s.pdf", database$databaseId)), plot, width = 3.5, height = 3.5, dpi = 300)
  
  # Negative controls
  ncsFileName <- file.path(database$outputFolder, sprintf("ncs_t%d_c%s_m%d.rds", target, comparator, mediator))
  ncs <- readRDS(ncsFileName)  
  negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)
  ncs <- ncs |>
    inner_join(negativeControls |>
                 select(outcomeId = conceptId,
                        outcomeName = conceptName),
               by = join_by(outcomeId))
  vizData <- bind_rows(
    ncs |>
      select(logRr = mainLogHr,
             lb = mainLogLb,
             ub = mainLogUb,
             outcomeName) |>
      mutate(label = "Main effect"),
    ncs |>
      select(logRr = indirectLogHr,
             lb = indirectLogLb,
             ub = indirectLogUb,
             outcomeName) |>
      mutate(label = "Indirect effect")
  )
  vizData$outcomeName <- factor(vizData$outcomeName, levels = sort(unique(vizData$outcomeName), decreasing = TRUE))
  vizData$label <- factor(vizData$label, levels = sort(unique(vizData$label), decreasing = TRUE))
  
  plot <- ggplot(vizData, aes(x = logRr, y = outcomeName)) +
    geom_vline(xintercept = 0) + 
    geom_point(shape = 16, color = "#336B91", alpha = 0.6) + 
    geom_errorbarh(aes(xmin = lb, xmax = ub), color = "#336B91", alpha = 0.6) +
    coord_cartesian(xlim = c(log(0.25), log(4))) +
    scale_x_continuous("Hazard ratio", breaks = log(c(0.5, 1, 2)), labels = c(0.5, 1, 2)) +
    facet_grid(~label) + 
    theme(axis.title.y = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(file.path(rootFolder, sprintf("Ncs_%s.png", database$databaseId)), plot, width = 7, height = 8, dpi = 300)
  ggsave(file.path(rootFolder, sprintf("Ncs_%s.pdf", database$databaseId)), plot, width = 7, height = 8, dpi = 300)
  
  # Percent negative controls signifcant
  ncs |>
    summarise(
      coverageMain = mean(mainLogLb <= 0 & mainLogUb >= 0, na.rm = TRUE),
      coverageIdirect = mean(indirectLogLb <= 0 & indirectLogUb >= 0, na.rm = TRUE)
    )
}



# Forest plot
library(ggplot2)
library(dplyr)
library(tidyr)
# database = databases[[1]]
estimates <- list()
for (database in databases) {
  hoisFileName <- file.path(database$outputFolder, sprintf("hois_t%d_c%s_m%d.rds", target, comparator, mediator))
  temp <- readRDS(hoisFileName)
  temp <- temp |>
    mutate(database = database$databaseId)
  estimates[[length(estimates) + 1]] <- temp
}
estimates <- bind_rows(estimates)
estimates <- estimates |>
  filter(outcomeId %in% c(outcome1, outcome2)) |>
  mutate(outcome = case_when(outcomeId == outcome1 ~ "Acute\nmyocardial infarction",
                             TRUE ~ "Ischemic stroke"),
         database = case_when(database == "OptumDod" ~ "Clinformatics",
                              database == "OptumEhr" ~ "Optum EHR",
                              TRUE ~ database))
estimates$database <- factor(estimates$database, levels = sort(unique(estimates$database), decreasing = TRUE))

vizData1 <- estimates |>
    transmute(estimand = "Main effect",
              hr = exp(mainLogHr),
              lb = exp(mainLogLb),
              ub = exp(mainLogUb),
              database = database,
              outcome)
vizData2 <- estimates |>
    transmute(estimand = "Indirect effect",
              hr = exp(indirectLogHr),
              lb = exp(indirectLogLb),
              ub = exp(indirectLogUb),
              database = database,
              outcome)
vizData3 <- estimates |>
    transmute(estimand = "Mediated proportion",
              hr = mediatedProportion,
              lb = mediatedProportionLb,
              ub = mediatedProportionUb,
              database = database,
              outcome)
breaks <- c(0.5, 1, 2)
plot1 <- ggplot(vizData1, aes(x = hr, y = database)) +
  geom_vline(xintercept = 1, size = 0.5) +
  geom_point(size = 3, shape = 18, color = "#336B91") +
  geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.15, color = "#336B91") +
  scale_x_log10("Hazard ratio", breaks = breaks) +
  coord_cartesian(xlim = c(0.4, 2.5)) +
  facet_grid(outcome~estimand) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines"),
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.text.y.right = element_blank()
  )
plot1

breaks <- c(0.9, 1, 1.1)
plot2 <- ggplot(vizData2, aes(x = hr, y = database)) +
  geom_vline(xintercept = 1, size = 0.5) +
  geom_point(size = 3, shape = 18, color = "#336B91") +
  geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.15, color = "#336B91") +
  scale_x_log10("Hazard ratio", breaks = breaks) +
  coord_cartesian(xlim = c(0.8, 1.2)) +
  facet_grid(outcome~estimand) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines"),
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.text.y.right = element_blank()
  )
plot2

breaks <- c(-0.5, 0, 0.5, 1, 1.5)
plot3 <- ggplot(vizData3, aes(x = hr, y = database)) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(size = 3, shape = 18, color = "#336B91") +
  geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.15, color = "#336B91") +
  scale_x_continuous("Proportion", breaks = breaks) +
  coord_cartesian(xlim = c(-1, 2)) +
  facet_grid(outcome~estimand) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines"),
    plot.title = element_text(size = 12, hjust = 0.5)
  )
plot3

plot <- grid.arrange(plot1, plot2, plot3, ncol = 3, widths = c(1, 0.8, 0.8))
ggsave(file.path(rootFolder, "ForestForMethods.png"), plot = plot, width = 7, height = 3.2, dpi = 300)
ggsave(file.path(rootFolder, "ForestForMethods.pdf"), plot = plot, width = 7, height = 3.2, dpi = 300)


# Table with counts

table <- estimates |>
  select("database",
         "outcome", 
         "targetSubjects",
         "comparatorSubjects",
         "targetMediators",
         "comparatorMediators",
         "targetOutcomes",
         "comparatorOutcomes",
         "targetMediatorOutcomes",
         "comparatorMediatorOutcomes")
readr::write_csv(table, file.path(rootFolder, "CountsForMethods.csv"))

# Plots from simulation studies ------------------------------------------------

# Table showing results of pilot simulations:
folder <- "PilotSimulation"
results <- readr::read_csv(file.path(folder, "Results.csv"))
# results |>
#   group_by(sampling, bootstrapType) |>
#   summarise(mseCoverageIndirectEffect = mean((0.95-coverageIndirectEffect)^2), 
#             mseCoverageMediatedProportion = mean((0.95-coverageMediatedProportion)^2),
#             mseIndirectEffect = mean(mseIndirectEffect),
#             mseMediatedProportion = mean(mseMediatedProportion)) |>
#   arrange(mseCoverageIndirectEffect)
table <- results |>
  filter(sampling != "weighted strata") |>
  group_by(sampling, bootstrapType) |>
  summarise(mseCoverageIndirectEffect = mean((0.95-coverageIndirectEffect)^2), 
            mseCoverageMediatedProportion = mean((0.95-coverageMediatedProportion)^2),
            .groups = "drop") 
readr::write_csv(table, file.path(folder, "SummaryTable.csv"))

# Violin plot showing pilot simulation results:
folder <- "PilotSimulation"
results <- readRDS("inst/shinyApps/MediationResultsExplorer/data/pilotSimulation.rds")
vizData <- results |>
  filter(metric %in% c("Coverage indirect effect",
                       "Coverage mediated proportion",
                       "Bias indirect effect",
                       "Bias mediated proportion"),
         sampling %in% c("strata", "person")) |>
  mutate(sampling = stringr::str_to_sentence(sampling),
         bootstrapType = gsub("^Reduced ", "Reduced\n", stringr::str_to_sentence(bootstrapType)),
         metricType = gsub(" .*$", "" ,metric),
         metricEstimand = gsub("^[^ ]* ", "", metric)) |>
  mutate(metricEstimand = case_when(
    metricEstimand == "indirect effect" ~ "Indirect\neffect",
    metricEstimand == "mediated proportion" ~ "Mediated\nproportion"
  )) |>
  filter(value < 5 & value > -5)

reference <- vizData |>
  distinct(sampling, bootstrapType, metricType, metricEstimand) |>
  mutate(y = case_when(metricType == "Bias" ~ 0,
                       metricType == "Coverage" ~ 0.95))

plot <- ggplot(vizData, aes(x = 0, y = value)) +
  geom_hline(aes(yintercept = y), data = reference, linetype = "dashed") +
  geom_violin(scale = "width", color = "#336B91", fill = "#336B91", alpha = 0.6) +
  expand_limits(y = 0) + 
  facet_nested(metricType + metricEstimand ~ sampling + bootstrapType, scales = "free_y") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
ggsave(file.path(folder, "PilotSimulations.png"), plot, width = 8.5, height = 6, dpi = 300)
ggsave(file.path(folder, "PilotSimulations.pdf"), plot, width = 8.5, height = 6, dpi = 300)


# Violin plot showing main simulation results:
folder <- "Simulation"
results <- readRDS("inst/shinyApps/MediationResultsExplorer/data/simulation.rds")
vizData <- results |>
  filter(metric %in% c("Coverage main effect",
                       "Coverage indirect effect",
                       "Coverage mediated proportion",
                       "Bias main effect",
                       "Bias indirect effect",
                       "Bias mediated proportion")) |>
  mutate(confounding = sprintf("Confounding = %0.1f", Confounding),
         outcomePrevalence = sprintf("Outcome baseline rate = %0.2f", `Baseline outcome prevalence`),
         metricType = gsub(" .*$", "" ,metric),
         metricEstimand = gsub("^[^ ]* ", "", metric)) |>
  mutate(metricEstimand = case_when(
    metricEstimand == "indirect effect" ~ "Indirect\neffect",
    metricEstimand == "main effect" ~ "Main\neffect",
    metricEstimand == "mediated proportion" ~ "Mediated\nproportion"
  )) |>
  filter(value < 5 & value > -5)
vizData$metricEstimand <- factor(vizData$metricEstimand, levels = c("Main\neffect", "Indirect\neffect", "Mediated\nproportion"))

reference <- vizData |>
  distinct(confounding, outcomePrevalence, metricType, metricEstimand) |>
  mutate(y = case_when(metricType == "Bias" ~ 0,
                       metricType == "Coverage" ~ 0.95))

plot <- ggplot(vizData, aes(x = 0, y = value)) +
  geom_hline(aes(yintercept = y), data = reference, linetype = "dashed") +
  geom_violin(scale = "width", color = "#336B91", fill = "#336B91", alpha = 0.6) +
  expand_limits(y = 0) + 
  facet_nested(metricType + metricEstimand ~ outcomePrevalence + confounding, scales = "free_y") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
plot
ggsave(file.path(folder, "MainSimulations.png"), plot, width = 7, height = 6, dpi = 300)
ggsave(file.path(folder, "MainSimulations.pdf"), plot, width = 7, height = 6, dpi = 300)

