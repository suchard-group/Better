#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

titleName <- "BETTER: Bayesian Evaluation of Time-To-Event and Reliability (for vaccine surveillance)"
methods <- c("HistoricalComparator", "SCCS")
timeAtRisks <- c("1-28", "1-42")

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    
  # style and title
  style = "width:1500px;",
  titlePanel(title = div(img(height = 50, src="OHDSI_toc_header.png"),
                         HTML('&emsp;'), titleName, HTML('&emsp;')),
             windowTitle = titleName),
  tags$head(tags$style(type = "text/css", "
             #loadmessage {
                                 position: fixed;
                                 top: 0px;
                                 left: 0px;
                                 width: 100%;
                                 padding: 5px 0px 5px 0px;
                                 text-align: center;
                                 font-weight: bold;
                                 font-size: 100%;
                                 color: #000000;
                                 background-color: #ADD8E6;
                                 z-index: 105;
                                 }
                                 ")),
  
    tabsetPanel(
      id = "mainTabsetPanel",
      tabPanel("About",
               includeMarkdown("md/about.md")
      ),
      tabPanel(
        "Estimation metrics",
        fluidRow(
          column(
            2,
            style = "background-color:#e8e8e8;",
            selectInput("exposure", label = "Vaccine:", choices = exposure$exposureName),
            selectInput("database", label = "Database:", choices = database$databaseId),
            selectInput("trueRr", label = "True effect size:", choices = trueRrs),
            checkboxGroupInput("method", label = "Design:", choices = methods, selected = methods),
            checkboxGroupInput("timeAtRisk", label = "Time at risk:", choices = timeAtRisks, selected = timeAtRisks),
          ),
          column(
            10,
            dataTableOutput("mseCoverageTable"),
            div(strong("Table:"), 
                "Estimation-oriented metrics. Mean-squared errors (MSEs) and coverage rate of 95% CIs.")
          )
        )
      ),
      tabPanel("Database information",
               dataTableOutput("databaseInfoTable"),
               div(strong("Table:"),"Information about each database.")
      ),
      tabPanel("Vaccine information",
               dataTableOutput("exposureInfoTable"),
               div(strong("Table:"),"Information about each vaccine exposure.")
      )
    )
  )
)
