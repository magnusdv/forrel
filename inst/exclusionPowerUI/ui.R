library(shiny)

#' UI definition of the shiny application
#'
#' Evaluating this file should yield a shiny interface definition object as
#' returned by [shiny::shinyUI()].
#'
#' This file mimics an HTML file describing the layout of the UI. The contents
#' are either static text labels refering to input widgets or references to
#' output widgets defined in the server.R file.
#'
#' In what follows are some comments explaining design decisions or
#' non-idiomatic uses of shiny functions.
#'
#' @seealso [shiny::shinyServer()]
#'
#' @references For a gentle introduction to shiny see
#'   http://shiny.rstudio.com/tutorial/
#' @references For detailed API documentation see
#'   http://shiny.rstudio.com/reference/shiny/
#'
#'
#' @author Elias Hernandis <eliashernandis@gmail.com>


shinyUI(fluidPage(
  # include some CSS needed by the application and associated modules
  includeCSS('styles.css'),
  fluidRow(
    column(12, h3('Exclusion Power'))
  ),
  sidebarLayout(
    mainPanel(
      tabsetPanel(
        ### Pedigree data
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

        ### Genetic data: frequencies, known genotypes and individuals available for genotyping
        tabPanel('Database',
                 fluidRow(
                   column(4,
                          wellPanel(
                            p(strong('Load frequency database file')),
                            p(advancedTableFileLoaderInput('frequencyDbFile', 'Select file')),
                            HTML('<p class="text-center">- or -</p>'),
                            fileInput('familiasFrequencyFile', 'Load Familias DNA data')
                          ),
                          radioButtons('frequencySource',
                                       'Use frequency data from',
                                       choices = c('Frequency DB file' = 'file',
                                                   'Familias DNA data' = 'fam')),
                                                   #'Pedigree file' = 'ped')),

                          # Help message describing the loaded frequency DB file
                          uiOutput('frequencyDbDescription'),
                          tabularDataPreviewInput('frequencies',
                                                  label = 'View loaded frequency data')),
                   column(4,
                          wellPanel(
                            p(strong('Load reference file')),
                            p(advancedTableFileLoaderInput('referenceFile', 'Select file')),
                            HTML('<p class="text-center">- or -</p>'),
                            fileInput('familiasRefereceFile', 'Load Familias case data')
                          ),
                          radioButtons('referenceSource',
                                       'Use known genotype data from',
                                       choices = c('Do not use reference data' = 'none',
                                                   'Reference file' = 'file',
                                                   'Familias case data' = 'fam')),
                                                   #'Pedigree file' = 'ped')),

                          # Help message describing the loaded reference file
                          uiOutput('referenceFileDescription'),
                          tabularDataPreviewInput('references',
                                                  label = 'View loaded reference data')),
                   column(4,
                          checkboxGroupInput('ids', 'Individuals available for genotyping'))
                 )),

        ### Settings: mutation model, inbreeding parameters and sex-linked markers
        tabPanel('Parameters',
                 p(strong('Marker settings')),
                 fluidRow(
                   column(12,
                          settingsTableInput('markerSettingsTable'))
                 )),

        ### Results: exclusion power calculation
        tabPanel('Results',
                 fluidRow(
                   column(12,
                          actionButton('computeButton', 'Compute exclusion power'))
                 ),
                 fluidRow(
                   column(12,
                          tableOutput('exclusionPowerResults'))
                 ))
      )

      # ### Save/restore workspace (visible across all tabs)
      # fluidRow(
      #   column(8, p(" ")),
      #   column(4,
      #          wellPanel(
      #            p(strong('Save workspace')),
      #            downloadButton('saveWorkspaceButton', 'Download'),
      #            p(fileInput('loadWorkspace', 'Load workspace'))
      #          ))
      #   )
    ),
    ### Pedigree plots, visible across all tabs
    sidebarPanel(
      verticalLayout(
        HTML('<h5 class="text-center">Claim pedigree</h5>'),
        plotOutput('pedClaimPlot'),
        HTML('<h5 class="text-center">True pedigree</h5>'),
        plotOutput('pedTruePlot'))
    )
  )
))
