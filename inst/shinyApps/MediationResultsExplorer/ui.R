library(shiny)
library(DT)
source("widgets.R")

shinyUI(
  fluidPage(style = "width:1500px;",
            titlePanel("Simulation Results"),
            tabsetPanel(id = "mainTabsetPanel",
                        tabPanel("About",
                                 HTML("<p>TODO</p>"),
                                 HTML("<p>For review use only")
                        ),
                        tabPanel("Simulations",       
                                 fluidRow(
                                   column(3,
                                          checkboxGroupInput("type", "Type", choices = types, selected = types[!grepl("Pooled", types)]),
                                          lapply(simParams, createSimParamWidget, results = results, suffix = "Main"),
                                          checkboxGroupInput("metric", "Metric", choices = metrics, selected = metrics[grepl("indirect effect", metrics)])
                                          ),  
                                   column(9,
                                          tabsetPanel(id = "EffectsTabsetPanel",
                                                      type = "pills",
                                                      tabPanel("Violin plots",
                                                               plotOutput("mainViolinPlot", height = 800),
                                                               div(align = "center",
                                                                   radioButtons("simParamRadioButton", "", choices = simParams, selected = simParams[1], inline = TRUE)
                                                               ),
                                                               uiOutput("mainViolinCaption")
                                                      ),tabPanel("Scatter plots",
                                                               uiOutput("hoverInfoPlot"),
                                                               plotOutput("mainPlot", height = 800, hover = hoverOpts("plotHoverMainPlot", delay = 100, delayType = "debounce")),
                                                               uiOutput("mainCaption")
                                                      ),
                                                      tabPanel("Rankings",
                                                               plotOutput("rankPlot", height = 800),
                                                               uiOutput("rankCaption")
                                                      )
                                          )
                                   )
                                 )
                        )
            )
  )
)
