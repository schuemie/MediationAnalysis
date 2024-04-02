library(shiny)
library(ggplot2)

singleRankVector <- function(row) {
  return(data.frame(type = row$type,
                    rank = row$start:row$end,
                    weight = 1 / (1 + row$end - row$start),
                    row.names = NULL))
}

computeRankVectors <- function(subgroup, descending = TRUE) {
  totalLength <- nrow(subgroup)
  valueCounts <- aggregate(rep(1, totalLength) ~ value, data = subgroup, FUN = sum) 
  if (nrow(valueCounts) == 1) {
    # One unique value. Distribute weights evenly over all ranks
    valueCounts$start <- 1
    valueCounts$end <- totalLength
  } else if (nrow(valueCounts) == totalLength) {
    # All values are unique. Simple ranking without consideration of ties
    rankVectors <- data.frame(type = subgroup$type,
                              rank = order(subgroup$value, decreasing = descending),
                              weight = 1,
                              row.names = NULL)
    return(rankVectors)
  } else {
    # Handle ties
    if (descending) {
      valueCounts <- valueCounts[order(-valueCounts$value), ]
    } else {
      valueCounts <- valueCounts[order(valueCounts$value), ]
    }
    valueCounts$end <- cumsum(valueCounts[, 2])
    valueCounts$start <-  c(1, valueCounts$end[1:(nrow(valueCounts) - 1)] + 1 )
    valueCounts[, 2] <- NULL
  }
  valueCounts <- merge(subgroup[, c("type", "value")], valueCounts) 
  rankVectors <- lapply(split(valueCounts, 1:totalLength), singleRankVector)
  rankVectors <- do.call("rbind", rankVectors)
  rankVectors$start <- NULL
  rankVectors$end <- NULL
  return(rankVectors)
}

shinyServer(function(input, output, session) {
  
  pivotData <- function(simParam, subset, dropIfUnique = TRUE) {
    if (dropIfUnique && length(unique(subset[[simParam]])) == 1) {
      return(NULL)
    } else {
      temp <- subset
      maxValue <- max(subset[[simParam]])
      temp$parameterValue <- subset[[simParam]]
      temp$jitter <- temp$parameterValue
      # temp$jitter <- temp$parameterValue + runif(nrow(subset), -0.02 * maxValue, 0.02 * maxValue)
      temp$simParam <- simParam
      temp[simParams] <- NULL
      return(temp)
    }
  }
  
  filteredResults <- reactive({
    subset <- results
    subset <- subset[subset$metric %in% input$metric, ]
    subset <- subset[subset$type %in% input$type, ]
    for (simParam in simParams) {
        subset <- subset[as.character(subset[[simParam]]) %in% input[[paste0(simParam, "Main")]], ]
    }
    return(subset)
  })
  
  filteredPivotedResults <- reactive({
    subset <- filteredResults()
    vizData <- lapply(simParams, pivotData, subset = subset)
    vizData <- do.call(rbind, vizData)
    return(vizData)
  })
  
  filteredViolinPivotedResults <- reactive({
    subset <- filteredResults()
    vizData <- pivotData(input$simParamRadioButton, subset, dropIfUnique = FALSE)
    return(vizData)
  })
  
  getReferenceValues <- function(metrics) {
    ref <- data.frame()
    idx <- grepl("Bias", metrics)
    if (any(idx)) {
      ref <- rbind(ref, data.frame(value = 0,
                                   metric = metrics[idx])) 
    }
    idx <- grepl("Coverage", metrics)
    if (any(idx)) {
      ref <- rbind(ref, data.frame(value = 0.95,
                                   metric = metrics[idx])) 
    }
    idx <- grepl("MSE", metrics)
    if (any(idx)) {
      ref <- rbind(ref, data.frame(value = 0,
                                   metric = metrics[idx])) 
    }
    idx <- grepl("Non Estimable", metrics)
    if (any(idx)) {
      ref <- rbind(ref, data.frame(value = 0,
                                   metric = metrics[idx])) 
    }
    return(ref)
  }
  
  output$mainPlot  <- renderPlot({
    subset <- filteredPivotedResults()
    if (nrow(subset) == 0) {
      return(NULL)
    } else {
      ref <- getReferenceValues(subset$metric)
      subset$type <- gsub(" ", "\n", subset$type)
      plot <- ggplot2::ggplot(subset, ggplot2::aes(x = jitter, y = value, group = type, color = type)) 
      if (nrow(ref) > 0) {
        plot <- plot + geom_hline(aes(yintercept = value), data = ref, linetype = "dashed")
      }
      plot <- plot + ggplot2::geom_point(alpha = 0.4) +
        ggplot2::scale_fill_manual(values = c("#69AED5", "#FBC511", "#EB6622", "#11A08A", "#336B91")) +
        ggplot2::scale_color_manual(values = c("#69AED5", "#FBC511", "#EB6622", "#11A08A", "#336B91")) +
        ggplot2::facet_grid(metric~simParam, scales = "free", switch = "both") +
        ggplot2::theme(legend.position = "top",
                       legend.title = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       strip.placement = "outside",
                       strip.background = ggplot2::element_blank())
      return(plot)
    }
  },
  res = 125,
  height = 800)
  
  output$mainViolinPlot  <- renderPlot({
    subset <- filteredViolinPivotedResults()
    if (nrow(subset) == 0) {
      return(NULL)
    } else {
      ref <- getReferenceValues(subset$metric)
      subset$type <- gsub(" ", "\n", subset$type)
      # subset <- subset[subset$simParam == input$simParamRadioButton, ]
      plot <- ggplot(subset, aes(x = factor(parameterValue), y = value, fill = type)) 
      if (nrow(ref) > 0) {
        plot <- plot + geom_hline(aes(yintercept = value), data = ref, linetype = "dashed")
      }
      plot <- plot + geom_violin(position = position_dodge(0.9), scale = "width", alpha = 0.4) +
        ggplot2::scale_fill_manual(values = c("#69AED5", "#FBC511", "#EB6622", "#11A08A", "#336B91")) +
        ggplot2::scale_color_manual(values = c("#69AED5", "#FBC511", "#EB6622", "#11A08A", "#336B91")) +
        facet_grid(metric~., scales = "free", switch = "both") +
        ggplot2::expand_limits(y = 0) + 
        theme(legend.position = "top",
              legend.title = element_blank(),
              axis.title = element_blank(),
              strip.placement = "outside",
              strip.background = element_blank())
      return(plot)
    }
  },
  res = 125,
  height = 800)
  
  output$hoverInfoPlot <- renderUI({
    subset <- filteredPivotedResults()
    if (nrow(subset) == 0) {
      return(NULL)
    } 
    hover <- input$plotHoverMainPlot
    point <- nearPoints(subset, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px - 361, "px; top:", top_px - 150, "px; width:250px;")
    
    # Unpivot:
    unpivotedRow <- results[results[, point$simParam] == point$parameterValue & 
                              results$type == point$type & 
                              results$metric == point$metric & 
                              results$value == point$value, ]
    unpivotedRow <- unpivotedRow[1, ]
    allMetrics <- merge(results, unpivotedRow[, c(simParams, "type")])
    
    lines <- sprintf("<b> Type: </b>%s", point$type)
    lines <- c(lines, "")
    lines <- c(lines, sprintf("<b> %s: </b>%s", simParams, unpivotedRow[, simParams]))
    lines <- c(lines, "")
    lines <- c(lines, sprintf("<b> %s: </b>%.2f", allMetrics$metric, allMetrics$value))
    
    div(
      style = "position: relative; width: 0; height: 0",
      wellPanel(
        style = style,
        p(HTML(paste(lines, collapse = "<br/>")))))
  })
  
  output$mainCaption <- renderUI({
    subset <- filteredPivotedResults()
    if (nrow(subset) == 0) {
      return(NULL)
    }
    count <- sum(subset$type == subset$type[1] & subset$metric == subset$metric[1] & subset$simParam == subset$simParam[1]) 
    HTML(sprintf("<strong>Figure S1.2. </strong>Each dot represents one of the %s selected simulation scenarios. The y-axes represent the various metrics 
    as estimated over 1,000 iterations per scenario, and the x-axes represent the various simulation parameters. Color indicates the various tested
                 meta-analysis algorithms. Hover over a data point to reveal more details.", count))
  })
  
  output$mainViolinCaption <- renderUI({
    subset <- filteredPivotedResults()
    if (nrow(subset) == 0) {
      return(NULL)
    }
    count <- sum(subset$type == subset$type[1] & subset$metric == subset$metric[1] & subset$simParam == subset$simParam[1]) 
    HTML(sprintf("<strong>Figure S1.1. </strong>Violin plots showing the performance accross the %s selected simulation scenarios. The y-axes represent the various metrics 
    as estimated over 1,000 iterations per scenario, and the x-axes represent %s. Color indicates the various tested
                 meta-analysis algorithms.", count, input$simParamRadioButton))
  })
  
  output$rankPlot  <- renderPlot({
    subset <- filteredResults()
    subset <- subset[!grepl("Bias", subset$metric), ]
    
    processMetric <- function(metricSubset) {
      metric <- metricSubset$metric[1]
      descending <- grepl("Precision", metric)
      if (grepl("Coverage", metric)) {
        metricSubset$value <- abs(0.95 - metricSubset$value)
      }
      subgroups <- split(metricSubset, apply(metricSubset[, c(simParams, "metric")],1,paste,collapse = " "))
      names(subgroups) <- NULL
      metricSubset <- lapply(subgroups, computeRankVectors, descending = descending)
      metricSubset <- do.call(rbind, metricSubset)  
      results <- aggregate(weight ~ type + rank, data = metricSubset, sum)
      metricSubset$metric <- metric
      return(metricSubset)
    }
    
    rankedSubset <- lapply(split(subset, subset$metric, drop = TRUE), processMetric)
    rankedSubset <- do.call(rbind, rankedSubset)
    rankedSubset$type <- gsub(" ", "\n", rankedSubset$type)
    plot <- ggplot2::ggplot(rankedSubset, ggplot2::aes(x = rank, y = weight)) +
      ggplot2::geom_col(color = rgb(0, 0, 0.8, alpha = 0), fill = rgb(0, 0, 0.8), alpha = 0.6) +
      ggplot2::scale_x_continuous("Rank (lower is better)", breaks = min(rankedSubset$rank):max(rankedSubset$rank)) +
      ggplot2::scale_y_continuous("Count") +
      ggplot2::facet_grid(type~metric) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank())
    
    return(plot)
  },
  res = 110,
  height = 800)
  
  output$rankCaption <- renderUI({
    subset <- filteredPivotedResults()
    if (nrow(subset) == 0) {
      return(NULL)
    }
    text <- "<strong>Figure S1.3. </strong>Histograms of algorithm ranks. Each bar represents the number of simulation scenarios where the algorithm on the 
    right achieved that rank on the metric at the top, compared to the other selected algorithms."
    if (any(grepl("coverage", subset$metric))) {
      text <- paste(text, "For coverage, algorithms were ranked by absolute difference between the estimated coverage and 95 percent.") 
    }
    HTML(text)
  })
})
