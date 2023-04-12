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
sensitivity_levels <- c(.25, .5)

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
        "Testing metrics",
        fluidRow(
          column(
            2,
            style = "background-color:#e8e8e8;",
            selectInput("exposureTest", label = "Vaccine:", choices = exposure$exposureName),
            selectInput("databaseTest", label = "Database:", choices = database$databaseId),
            radioButtons("methodTest", label = "Design:", choices = methods, selected = methods[1]),
            radioButtons("timeAtRiskTest", label = "Time at risk:", choices = timeAtRisks, selected = timeAtRisks[1]),
          ),
          column(
            10,
            tabsetPanel(type = "pills",
                        tabPanel("Type 1 error",
                                 plotOutput("type1Plot"),
                                 div(strong("Plot:"), 
                                     "Empirical Type 1 error rate over analysis periods. Type 1 error rates are measured by fraction of H0 rejected over all negative control outcomes.")
                        ),
                        tabPanel("Statistical power",
                                 selectizeInput("analysis", "Design variant:", choices = c("choose" = "", unique(analysis$description[analysis$method == methods[1]]))),
                                 plotOutput("powerPlot"),
                                 div(strong("Plot:"), 
                                     "Statistical power over analysis periods. Power is measured by fraction of positive control outcomes with H0 rejected, stratified by effect sizes.")
                                 ),
                        tabPanel("Time-to-signal",
                                 selectInput("sensitivity", 'Sensitivity:', choices = sensitivity_levels, selected = sensitivity_levels[2]),
                                 plotOutput("ttsPlot"),
                                 div(strong("Plot:"), 
                                     "Timeliness. Measured by number of analyses (in months) needed to reach a desired sensitivity level for detecting true positive signals.")
                                 )
                        )
            )
          )
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
