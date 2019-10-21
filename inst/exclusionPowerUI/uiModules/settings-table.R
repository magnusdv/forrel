settingsTableInput = function(inputId) {
  # create a namespace to avoid name collisions with other modules / other
  # instances of this one
  ns <- NS(inputId)

  tagList(
    shiny::singleton(tags$head(
      tags$script(src = 'settingsTable/js/binding.js'),
      tags$link(rel = 'stylesheet',
                type = 'text/css',
                href = 'settingsTable/css/style.css')
    )),
    tags$div(
      tags$table(
        # the table content is build on the client through the shiny binding
        class = 'table table-bordered table-condensed'
      ),
      class = 'table-responsive settings-table',
      id = ns('settingsTable')
    )
  )
}

hashToJSONlike = function(hash) {
  values = names(hash)
  labels = unname(hash)

  as.data.frame(list('value' = values, 'label' = labels))
}

STlabel = function(name) {
  list('name' = name,
       'type' = 'label')
}

STcheckbox = function(name) {
  list('name' = name,
       'type' = 'checkbox')
}

STdropdown = function(name, options) {
  list('name' = name,
       'type' = 'dropdown',
       'options' = hashToJSONlike(options))
}

STradio = function(name, options) {
  list('name' = name,
       'type' = 'radio',
       'options' = hashToJSONlike(options))
}

settingsTable = function(input, output, session, id = 'namespace', fields) {
  session$sendInputMessage(id, list('fields' = fields, 'data' = list()))

  ns <- NS(id)

  return(reactive({
    if (isTruthy(input$settingsTable))
      jsonlite::fromJSON(input$settingsTable)
    else
      NULL
  }))
}

updateSettingsTableInput = function(session, id = 'namespace', fields, data) {
  ns = NS(id)

  session$sendInputMessage(ns('settingsTable'),
                           jsonlite::toJSON(list('fields' = fields, 'data' = data),
                                            auto_unbox = F))
}
