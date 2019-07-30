library(shiny)
library(forrel)

#' Server definition of the shiny application
#'
#' Evaluating this file should yield a shiny server definition object as
#' returned by [shiny::shinyServer()].
#'
#' This file mostly contains an enumeration of outputs that are used by the
#' `ui.R` file. Keep in mind that the order of definition of this enumeration is
#' not important and that all output definitions are evaluated whenever an input
#' they depend on changes per the principles of reactive programming. An
#' important exception to this is the calculation of the exclusion power, which
#' is triggered by a button in the last tab of the UI. This design decision is
#' justified by the facts that the exclusion power cannot be calculated until
#' all parameters are set and because it may take a long time to compute.
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

shinyServer(function(input, output, session) {

  # obtain the claim pedigree
  claimPedigree <- reactive({
    pedigreeFromUI(input$pedClaim, pedfile = input$pedClaimFile)
  })

  pedPlotColors = reactive({
    list(red = input$ids)
  })

  pedPlotShaded = reactive({
    unique(references()[,1])
  })

  # render the pedigree plot when the user chooses a Claim pedigree or updates
  # individuals available for genotyping
  output$pedClaimPlot <- renderPlot({
    if (!is.pedList(claimPedigree())) {
      plot(claimPedigree(),
           col = list(red = intersect(input$ids, labels(claimPedigree()))),
           shaded = intersect(pedPlotShaded(), labels(claimPedigree())))
    } else {
      plot.arg.list = lapply(claimPedigree(), function(x) {
        list(x = x,
             col = list(red = intersect(input$ids, labels(x))),
             shaded = intersect(pedPlotShaded(), labels(x)))
      })
      pedtools::plotPedList(plot.arg.list,
                            frames = FALSE)
    }
  })

  # obtain the true pedigree
  truePedigree <- reactive({
    if (input$pedTrue == 'unrelated') {
      unrelatedPedFromClaimPed(claimPedigree(), input$ids)
    } else {
      pedigreeFromUI(input$pedTrue, pedfile = input$pedTrueFile)
    }
  })

  # render the pedigree plot when the user chooses a True pedigree
  output$pedTruePlot <- renderPlot({
    if (!is.pedList(truePedigree())) {
      plot(truePedigree(),
           col = list(red = intersect(input$ids, labels(truePedigree()))),
           shaded = intersect(pedPlotShaded(), labels(truePedigree())))
    } else {
      plot.arg.list = lapply(truePedigree(), function(x) {
        list(x = x,
             col = list(red = intersect(input$ids, labels(x))),
             shaded = intersect(pedPlotShaded(), labels(x)))
      })
      pedtools::plotPedList(plot.arg.list,
                            frames = FALSE)
    }
  })

  # update list of persons available for genotyping when the user chooses a pedigree
  observe({
    if (!is.null(claimPedigree())) {
      updateCheckboxGroupInput(session, "ids", choices = getMembers(claimPedigree()))
    }
  })

  # choose "Custom .ped file" in the claim pedigree dropdown whenever the user
  # selects a .ped file using the file input.
  observe({
    if (!is.null(input$pedClaimFile)) {
      updateSelectInput(session, 'pedClaim', selected = 'pedfile')
    }
  })

  # choose "Custom .ped file" in the true pedigree dropdown whenever the user
  # selects a .ped file using the file input.
  observe({
    if (!is.null(input$pedTrueFile)) {
      updateSelectInput(session, 'pedTrue', selected = 'pedfile')
    }
  })

  # load frequency file
  frequencyDB = callModule(advancedTableFileLoader, 'frequencyDbFile', id = 'frequencyDbFile',
                           columnHeaders = TRUE,
                           rowHeaders = TRUE)

  # load reference file
  references <- callModule(advancedTableFileLoader, 'referenceFile', id = 'referenceFile',
                           columnHeaders = TRUE)

  # update list of sex-linked markers when new allele denomination data becomes available
  observe({
    updateCheckboxGroupInput(session, 'sexLinkedMarkers',
                             choices = getMarkerNames(computedClaimPed()),
                             selected = input$sexLinkedMarkers)
  })

  # data provided thorugh tabular files is attached to the claim pedigree here
  computedClaimPed = reactive({
    ped = claimPedigree()

    switch (input$frequencySource,
      'file' = {
        if (isTruthy(frequencyDB())) {
          ped = attachAlleleFrequenciesToPedigree(ped, df = frequencyDB(), Xchrom = input$sexLinkedMarkers)
        }
      }
    )

    switch (input$referenceSource,
      'file' = {
        if (isTruthy(references()) && isTruthy(frequencyDB())) {
          ped = attachGenotypeToPedigree(ped, df = references())
        }
      }
    )

    ped
  })

  # compose the description for the frequency file
  output$frequencyDbDescription = renderUI({
    if (!isTruthy(frequencyDB())) {
      p('No frequency database loaded.')
    } else {
      p(sprintf('Allele frequency loaded for %d markers', length(computedClaimPed()$markerdata)))
    }
  })

  # compose the description for the reference file
  output$referenceFileDescription = renderUI({
    if (!isTruthy(references())) {
      p('No known genotype data loaded.')
    } else {
      p(sprintf('Known genotypes loaded for %s',
                paste(unique(references()[,1]), collapse = ', ')))
    }
  })

  callModule(tabularDataPreview, 'frequencies',
             id = 'frequencies',
             df = frequencyDB(),
             title = 'Frequency data used for calculation')

  callModule(tabularDataPreview, 'references',
             id = 'references',
             df = references(),
             title = 'Reference data used for calculation')

  # compute exclusion power
  output$exclusionPowerResults <- renderTable({
    if (input$computeButton < 1) return(NULL)

    isolate({
      withProgress({
        markerNames = getMarkerNames(computedClaimPed())
        Nmarkers = length(markerNames)
        exclusionProbabilities = vector('numeric', Nmarkers)

        # compute exclusion power
        for (i in 1:Nmarkers) {
          markerName = markerNames[i]

          EP = exclusionPower(ped_claim = computedClaimPed(),
                              ped_true = truePedigree(),
                              ids = input$ids,
                              markerindex = i,
                              plot = FALSE)
          exclusionProbabilities[i] = EP

          incProgress(1/Nmarkers)
        }

        # TODO: consider building this row by row
        data.frame("Marker" = markerNames, "Exclusion probability" = exclusionProbabilities)
      }, message = 'Calculating exclusion power of all markers...')
    })
  }, digits = 4)
})
