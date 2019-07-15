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

unrelatedPedFromClaimPed = function(claimPed, ids) {
  if (length(ids) == 0) return(list())

  sexes = getSex(claimPed, ids)
  tuples = cbind(ids, sexes)
  return(apply(tuples, 1, function(row) {pedtools::singleton(row[1], sex = row[2])}))
}

pedigreeFromUI = function(pedigreeID, pedfile = NULL, ids = NULL, sids = NULL) {
  if (pedigreeID == 'nucPed-1s') {
    return(nuclearPed(1, father = "Father", mother = "Mother", children = c("Son")))
  } else if (pedigreeID == 'nucPed-1d') {
    return (nuclearPed(1, sex = 2, father = "Father", mother = "Mother", children = c("Daughter")))
  } else if (pedigreeID == 'pedfile') {
    if (is.null(pedfile)) {
      return()
    }
    return(as.ped(read.table(pedfile$datapath)))
  }
}

shinyServer(function(input, output, session) {

  # obtain the claim pedigree
  claimPedigree <- reactive({
    pedigreeFromUI(input$pedClaim, pedfile = input$pedClaimFile, input$ids)
  })

  # render the pedigree plot when the user chooses a Claim pedigree or updates
  # individuals available for genotyping
  output$pedClaimPlot <- renderPlot({
    if (input$pedClaim != "unrelated") {
      colors = ifelse(labels(claimPedigree()) %in% input$ids, 2, 1)
      plot(claimPedigree(), col = colors)
    } else {
      pedtools::plotPedList(claimPedigree(), frames = FALSE)
    }
  })

  # update list of persons available for genotyping when the user chooses a pedigree
  observe({
    if (!is.null(claimPedigree())) {
      updateCheckboxGroupInput(session, "ids", choices = claimPedigree()$ID)
    }
  })

  # update selected entry in dropdowns when user uploads pedigree file
  observe({
    if (!is.null(input$pedTrueFile)) {
      updateSelectInput(session, 'pedTrue', selected = 'pedfile')
    }
  })

  observe({
    if (!is.null(input$pedClaimFile)) {
      updateSelectInput(session, 'pedClaim', selected = 'pedfile')
    }
  })

  # obtain the true pedigree
  truePedigree <- reactive({
    if (input$pedTrue == 'unrelated') {
      unrelatedPedFromClaimPed(claimPedigree(), input$ids)
    } else {
      pedigreeFromUI(input$pedTrue, pedfile = input$pedTrueFile, input$ids)
    }
  })

  # render the pedigree plot when the user chooses a True pedigree
  output$pedTruePlot <- renderPlot({
    if (input$pedTrue != "unrelated") {
      colors = ifelse(labels(truePedigree()) %in% input$ids, 2, 1)
      plot(truePedigree(), col = colors)
    } else {
      pedtools::plotPedList(truePedigree(), frames = FALSE)
    }
  })

  # load frequency file
  frequencyDB = callModule(advancedTableFileLoader, 'frequencyDbFile', id = 'frequencyDbFile',
                           columnHeaders = TRUE,
                           rowHeaders = TRUE)

  output$frequencyDbDescription = renderUI({
    if (!isTruthy(frequencyDB())) {
      p('No frequency database loaded.')
    } else {
      p(sprintf('Allele frequency loaded for %d markers', ncol(frequencyDB())))
    }
  })

  # update list of sex-linked markers when a new reference file is loaded
  # TODO:
  # this will probably need to be changed if we support loading of allele
  # denomination data from other sources such as pedigrees
  observe({
    if (isTruthy(frequencyDB())) {
      updateCheckboxGroupInput(session, 'sexLinkedMarkers', choices = colnames(frequencyDB()))
    }
  })

  # load reference file
  references <- callModule(advancedTableFileLoader, 'referenceFile', id = 'referenceFile',
                           columnHeaders = TRUE)

  output$referenceFileDescription = renderUI({
    if (!isTruthy(references())) {
      p('No known genotype data loaded.')
    } else {
      p(sprintf('Known genotypes loaded for %s', paste(unique(references()[,1]), sep = ', ')))
    }
  })

  # compute exclusion power
  output$exclusionPowerResults <- renderTable({
    if (input$computeButton < 1) return(NULL)

    isolate({
      withProgress({
        Nmarkers = ncol(frequencyDB())
        markerNames = colnames(frequencyDB())
        exclusionProbabilities = vector('numeric', Nmarkers)

        # compute exclusion power
        for (i in 1:Nmarkers) {
          markerName = markerNames[i]

          alleles = rownames(frequencyDB()[!is.na(frequencyDB()[markerName]),])
          frequencies = frequencyDB()[!is.na(frequencyDB()[markerName]),markerName]

          # load relevant genotype data for this marker
          knownGen = NULL
          if (isTruthy(references())) {
            ref = references()
            # TODO: turn this into functional code
            relevantRef = ref[grep(markerName, ref[,2]),]
            if (nrow(relevantRef) > 0) {
              knownGen = vector("list", nrow(relevantRef))
              for (i in 1:nrow(relevantRef)) {
                knownGen[[i]] = c(as.character(relevantRef[i,1]), relevantRef[i,3], relevantRef[i,4])
              }
            } else {
              knownGen = NULL
            }
          }

          EP = exclusionPower(ped_claim = claimPedigree(),
                              ped_true = truePedigree(),
                              ids = input$ids,
                              known_genotypes = knownGen,
                              alleles = alleles,
                              afreq = frequencies,
                              Xchrom = markerName %in% input$sexLinkedMarkers,
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
