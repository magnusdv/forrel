#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(forrel)

pedigreeFromUI = function(pedigreeID, pedfile = NULL) {
  if (pedigreeID == 'nucPed-1s') {
    return(nuclearPed(1, father = "Father", mother = "Mother", children = c("Son")))
  } else if (pedigreeID == 'nucPed-1d') {
    return (nuclearPed(1, sex = 2, father = "Father", mother = "Mother", children = c("Daughter")))
  } else if (pedigreeID == 'pedfile') {
    if (is.null(pedfile)) {
      return()
    }

    return(as.ped(read.table(pedfile$datapath)))
  } else if (pedigreeID == "unrelated") {
    return(list(singleton(id = 1), singleton(id = 3, sex = 2)))
    }
}

shinyServer(function(input, output, session) {

  # obtain the claim pedigree
  claimPedigree <- reactive({
    pedigreeFromUI(input$pedClaim, pedfile = input$pedClaimFile)
  })

  # render the pedigree plot when the user chooses a Claim pedigree or updates
  # individuals available for genotyping
  output$pedClaimPlot <- renderPlot({
    colors = ifelse(labels(claimPedigree()) %in% input$ids, 2, 1)
    plot(claimPedigree(), col = colors)
  })

  # update list of persons available for genotyping when the user chooses a pedigree
  observe({
    if (!is.null(claimPedigree())) {
      updateCheckboxGroupInput(session, "ids", choices = claimPedigree()$ID)
    }

    if (!is.null(input$pedClaimFile)) {
      updateSelectInput(session, 'pedClaimFile', selected = 'pedfile')
    }
  })

  # obtain the true pedigree
  truePedigree <- reactive({
    pedigreeFromUI(input$pedTrue, pedfile = input$pedTrueFile)
  })

  # render the pedigree plot when the user chooses a True pedigree
  output$pedTruePlot <- renderPlot({
    colors = ifelse(labels(truePedigree()) %in% input$ids, 2, 1)
    plot(truePedigree(), col = colors)
  })

  # load frequency file
  frequencyDB = callModule(advancedTableFileLoader, 'frequencyDbFile', id = 'frequencyDbFile')

  # load reference file(s)
  references <- callModule(advancedTableFileLoader, 'referenceFiles', id = 'referenceFiles')

  # compute exclusion power
  output$exclusionPowerResults <- renderTable({
    if (input$computeButton < 1) return(NULL);

    isolate({
      withProgress({
        Nmarkers = ncol(frequencyDB())
        markerNames = colnames(frequencyDB())
        exclusionProbabilities = vector('numeric', Nmarkers)

        for (i in 1:Nmarkers) {
          marker = markerNames[i]
          alleles = rownames(frequencyDB()[!is.na(frequencyDB()[marker]),])
          frequencies = frequencyDB()[!is.na(frequencyDB()[marker]),marker]
          EP = exclusionPower(ped_claim = claimPedigree(),
                              ped_true = truePedigree(),
                              ids = input$ids,
                              alleles = alleles,
                              afreq = frequencies,
                              plot = FALSE)
          exclusionProbabilities[i] = EP

          incProgress(1/Nmarkers)
        }

        data.frame("Marker" = markerNames, "Exclusion probability" = exclusionProbabilities)
      }, message = 'Calculating exclusion power of all markers...')
    })
  }, digits = 4)
})
