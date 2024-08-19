
library(dplyr)
library(ggplot2)
library(gridExtra)

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
                                 targetLabel = "DOACs",
                                 comparatorLabel = "Warfarin")
  psPlot <- psPlot + 
    # scale_fill_manual(values = c(alpha("#336B91", 0.7), alpha("#EB6622", 0.7))) +
    scale_fill_manual(values = c(alpha("#336B91", 0.6), alpha("#EB6622", 0.6))) +
    scale_color_manual(values = c(alpha("#336B91", 0.7), alpha("#EB6622", 0.7))) +
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
}

# Forest plot
# See last part of Execute.R

