#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(dplyr)

source("plotsAndTables.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  cache <- new.env()
  
  # some basic reactives from input -----
  exposureId <- reactive(
    exposure %>%
      filter(.data$exposureName == input$exposure) %>%
      pull(.data$exposureId)
  )
  
  baseExposureId <- reactive(
    exposure %>%
      filter(.data$exposureName == input$exposure) %>%
      pull(.data$baseExposureId)
  )
  
  analysisId <- reactive(
    analysis %>%
      filter(.data$method %in% input$method, .data$timeAtRisk %in% input$timeAtRisk) %>%
      pull(.data$analysisId)
  )
  
  # type 1 error ----
  type1Plot <- reactive({
    
  })
  
  output$type1Plot <- renderPlot({
    return(type1Plot())
  })
  
  # estimation metrics ----
  filterMSEs <- reactive({
    subsetMSEs <- mses %>%
      filter(.data$databaseId == input$database,
             .data$exposureId == exposureId,
             .data$method %in% input$method, 
             .data$analysisId %in% analysisId)
    if (input$trueRr != 'Any'){
      subsetMSEs <- subsetMSEs %>%
        filter(.data$effectSize == as.numeric(input$trueRr))
    }
    return(subsetMSEs)
  })
  
  output$mseCoverageTable <- renderDataTable({
    data <- filterMSEs()
    if (nrow(data) == 0) {
      return(data)
    }
    data <- data %>% inner_join(analysis, by = c('method', 'analysisId')) %>%
      select(method, description, timeAtRisk, effectSize, maxsprtMse, adjustedMse, maxsprtCoverage, adjustedCovearge)
    
    selection = list(mode = "single", target = "row")
    options = list(pageLength = 10, 
                   searching = TRUE, 
                   lengthChange = TRUE,
                   ordering = TRUE,
                   digits = 3,
                   columnDefs = list(list(width = '20%', targets = 1)))
    
    table <- DT::datatable(data, 
                           selection = selection, 
                           options = options, 
                           colnames = c("Design", "Description", "Time-at-risk", 
                                        "Effect size", "MSE (MaxSPRT)", "MSE (Bayesian)", 
                                        "Coverage (MaxSPRT)", "Coverage (Bayesian)"),
                           rownames = FALSE, 
                           escape = FALSE) 
    
    colors <- c("#b7d3e6", "#b7d3e6", "#f2b4a9", "#f2b4a9")
    mins <- c(0, 0, 0, 0)
    maxs <- c(max(data[, 5:6]), max(data[, 5:6]), 1, 1)
    for (i in 1:length(colors)) {
      table <- DT::formatStyle(table = table, 
                               columns = i + 4,
                               background = styleColorBar(c(mins[i], maxs[i]), colors[i]),
                               backgroundSize = '98% 88%',
                               backgroundRepeat = 'no-repeat',
                               backgroundPosition = 'center')
    }
    return(table)
  })
  
  # database information -----
  output$databaseInfoTable <- renderDataTable({
    
    table <- database %>%
      select(.data$databaseId, .data$databaseName, .data$description, .data$vocabularyVersion)
    options = list(pageLength = 25,
                   searching = TRUE,
                   lengthChange = FALSE,
                   ordering = TRUE,
                   paging = FALSE,
                   columnDefs = list(list(width = '30%', targets = 1),
                                     list(width = '60%', targets = 2))
    )
    table <- datatable(table,
                       options = options,
                       colnames = c("ID", "Name", "Description", "Vocabulary"),
                       rownames = FALSE,
                       class = "stripe compact")
    return(table)
  })
  
  # vaccine information table ----
  output$exposureInfoTable <- renderDataTable({
    
    table <- exposure %>%
      select(.data$exposureName, .data$shot, 
             .data$startDate, .data$endDate,
             .data$historyStartDate, .data$historyEndDate)
    options = list(pageLength = 25,
                   searching = TRUE,
                   lengthChange = FALSE,
                   ordering = TRUE,
                   paging = FALSE,
                   columnDefs = list(list(width = '30%', targets = 0))
    )
    table <- datatable(table,
                       options = options,
                       colnames = c("Name", "Shot", "Start Date", "End Date",
                                    "History Start Date", "History End Date"),
                       rownames = FALSE,
                       class = "stripe compact")
    return(table)
  })
  
})
