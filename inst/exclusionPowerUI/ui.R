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
                     "Nuclear family (2 sons)" = "nucPed-2s",
                     "Nuclear family (2 daughters)" = "nucPed-2d",
                     "Nuclear family (1 son, 1 daughter)" = "NucPed-1s-1d",
                     "Unrelated" = "unrelated",
                     "Custom (from .ped file)" = "pedfile")

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  includeCSS('styles.css'),
  fluidRow(
    column(12, titlePanel('Exclusion Power'))
  ),
  fluidRow(
    column(12, h3('Step 1. Choose pedigrees'))
  ),
  fluidRow(
    column(6,
           selectInput("pedClaim", "Claim pedigree", defaultPedigrees, width = '100%'),
           HTML('<p class="text-center">- or -</p>'),
           fileInput('pedClaimFile', ".ped file",  width = '100%'),
           plotOutput('pedClaimPlot'),
           checkboxGroupInput("ids", "Individuals available for genotyping")),
    column(6,
           selectInput("pedTrue", "True pedigree", defaultPedigrees, width = '100%'),
           HTML('<p class="text-center">- or -</p>'),
           fileInput('pedTrueFile', ".ped file",  width = '100%'),
           plotOutput('pedTruePlot'))
  ),
  fluidRow(
    column(12,
           h3('Step 2. Provide available genetic data'))
  ),
  fluidRow(
    column(4,
           fileInput('frequencyDbFile', 'Allele frequency database', width = '100%'),
           textOutput('frequencyDbSummary')),
    column(4,
           fileInput('referenceFiles', 'Reference files', width = '100%'))
  )
))
