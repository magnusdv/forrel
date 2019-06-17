#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

pedigreeFromUI = function(pedigreeID) {
  if (pedigreeID == 'nucPed-1s') {
    return(nuclearPed(1, father = "Father", mother = "Mother", children = c("Son")))
  } else if (pedigreeID == 'nucPed-1d') {
    return (nuclearPed(1, sex = 2, father = "Father", mother = "Mother", children = c("Daughter")))
  } else if (pedigreeID == "unrelated") {
    return(list(singleton(id = 1), singleton(id = 3, sex = 2)))
    }
}

shinyServer(function(input, output, session) {
  observe({
    pedigree = pedigreeFromUI(input$pedClaim)
    updateCheckboxGroupInput(session, "ids", choices = pedigree$ID)
    output$pedClaimPlot <- renderPlot({
      plot(pedigree)
    })
  })

  output$pedTruePlot <- renderPlot({
    pedigree = pedigreeFromUI(input$pedTrue)
    plot(pedigree)
  })

})
