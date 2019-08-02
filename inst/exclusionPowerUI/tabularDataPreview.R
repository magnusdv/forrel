tabularDataPreviewInput = function(id, label = 'Preview') {
  # create a namespace to avoid name collisions with other modules / other instances of this one
  ns <- NS(id)

  tagList(
    actionLink(ns('modalLink'), label = label)
  )
}

tabularDataPreview = function(input, output, session, id = 'namespace', title = NULL, df = NULL) {
  ns <- NS(id)

  inputModal <- function() {
    modalDialog(
      fluidPage(
        fluidRow(
          column(12,
                 tableOutput(ns('tablePreview')))
        )
      ),
      title = title,
      footer = tagList(
        actionButton(ns('ok'), 'OK')
      ))
  }

  output$tablePreview = renderTable(df, rownames = TRUE, decimals = 6)

  # show modal when action link is clicked
  observeEvent(input$modalLink, {
    showModal(inputModal())
  })

  # dismiss modal when the user clicks OK
  observeEvent(input$ok, {
    removeModal()
  })
}
