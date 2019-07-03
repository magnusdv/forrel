#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

defaultPedigrees = c("Nuclear family (1 son)" = "nucPed-1s",
                     "Nuclear family (1 daughter)" = "nucPed-1d",
                     "Unrelated" = "unrelated",
                     "Custom (from .ped file)" = "pedfile")

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  includeCSS('styles.css'),
  fluidRow(
    column(12, h3('Exclusion Power'))
  ),
  sidebarLayout(
    sidebarPanel(
      verticalLayout(
        p(strong('Claim pedigree')),
        plotOutput('pedClaimPlot'),
        p(strong('True pedigree')),
        plotOutput('pedTruePlot'))
      ),
    mainPanel(
      tabsetPanel(
        tabPanel('Pedigrees',
                 fluidRow(
                   column(6,
                          selectInput("pedClaim", "Claim pedigree", defaultPedigrees, width = '100%'),
                          HTML('<p class="text-center">- or -</p>'),
                          fileInput('pedClaimFile', ".ped file",  width = '100%')),
                   column(6,
                          selectInput("pedTrue", "True pedigree", defaultPedigrees, width = '100%'),
                          HTML('<p class="text-center">- or -</p>'),
                          fileInput('pedTrueFile', ".ped file",  width = '100%'))
                 )),
        tabPanel('Genetic data',
                 fluidRow(
                   column(4,
                          p(strong('Load frequency database file')),
                          p(advancedTableFileLoaderInput('frequencyDbFile', 'Select file'))),
                   column(4,
                          p(strong('Load reference file(s)')),
                          p(advancedTableFileLoaderInput('referenceFiles', 'Select file'))),
                          # helpText('Loading multiple reference files is supported')),
                   column(4,
                          checkboxGroupInput("ids", "Individuals available for genotyping"))
                 )),
        tabPanel('Settings',
                 fluidRow(
                   column(6,
                          p(strong('Mutation model'))),
                   column(6,
                          p(strong('Sex-linked markers')))
                 ),
                 fluidRow(
                   column(12, h4('Inbreeding parameters'))),
                 fluidRow(
                   column(6,
                          p(strong('Claim pedigree'))),
                   column(6,
                          p(strong('True pedigree')))
                 )),
        tabPanel('Results',
                 fluidRow(
                   column(12,
                          actionButton('computeButton', 'Compute exclusion power'))
                 ),
                 fluidRow(
                   column(12,
                          tableOutput('exclusionPowerResults'))
                 ))
      ),
      fluidRow(
        column(8, p(" ")),
        column(4,
               wellPanel(
                 p(strong('Save workspace')),
                 downloadButton('saveWorkspaceButton', 'Download'),
                 p(fileInput('loadWorkspace', 'Load workspace'))
               ))
        )
    )
  )
))
