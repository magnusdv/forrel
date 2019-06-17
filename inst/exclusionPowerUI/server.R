#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

pedigreeFromUI = function(pedigreeID, pedfile = NULL) {
  if (pedigreeID == 'nucPed-1s') {
    return(nuclearPed(1, father = "Father", mother = "Mother", children = c("Son")))
  } else if (pedigreeID == 'nucPed-1d') {
    return (nuclearPed(1, sex = 2, father = "Father", mother = "Mother", children = c("Daughter")))
  } else if (pedigreeID == 'pedfile') {
    if (is.null(pedfile)) return();

    return(as.ped(read.table(pedfile$datapath)))
  } else if (pedigreeID == "unrelated") {
    return(list(singleton(id = 1), singleton(id = 3, sex = 2)))
    }
}

shinyServer(function(input, output, session) {
  # update list of persons available for genotyping when the user chooses a pedigree
  observe({
    pedigree = pedigreeFromUI(input$pedClaim, input$pedClaimFile)
    if (!is.null(pedigree)) {
      updateCheckboxGroupInput(session, "ids", choices = pedigree$ID)

      # render the pedigree plot when the user chooses a Claim pedigree
      output$pedClaimPlot <- renderPlot({
        plot(pedigree)
      })
    }

    if (!is.null(input$pedClaimFile)) {
      updateSelectInput(session, 'pedClaimFile', selected = 'pedfile')
    }
  })

  # render the pedigree plot when the user chooses a True pedigree
  output$pedTruePlot <- renderPlot({
    pedigree = pedigreeFromUI(input$pedTrue, input$pedTrueFile)
    plot(pedigree)
  })

  # load frequency file
  frequencyDB <- reactiveVal()
  output$frequencyDbSummary <- renderText({
    if (is.null(input$frequencyDbFile)) {
      return("Load a frequency database file")
    } else {
      frequencyDB(read.table(input$frequencyDbFile$datapath))
      return(sprintf("Loaded frequency database with %d markers", ncol(frequencyDB())))
    }
  })
})
