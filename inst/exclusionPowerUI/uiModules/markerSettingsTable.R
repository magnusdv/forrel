
#' Given a pedigree, it returns a `data.frame` with the relevant marker
#' metadata: sex-linkedness (as chromosome number), mutation model...
#'
#' The purpose of this function is to construct an R object that can be easily
#' serialised as JSON.
#'
#' @param ped a `ped` object, or a list of such (in which case the marker
#'   settings will be deduced from only the first pedigree).
#' @return a `data.frame`
pedToMarkerSettingsState = function(ped) {
  if (is.pedList(ped)) return(pedToMarkerSettingsState(ped[[1]]))

  lapply(ped$markerdata, function(m) {
    list('markerName' = attr(m, 'name'),
         'chrom' = if (is.na(attr(m, 'chrom'))) '0' else attr(m, 'chrom'),
         'includeInCalculation' = '1',
         'mutationModel' = if (is.null(attr(m, 'mutmod'))) 'None' else attr(m, 'mutmod'),
         'rate' = attr(m, 'rate'))
  })
}

#' Initialises a marker settings table UI-widget.
#'
#' @param inputId a unique identifier for the input. Must be unique per
#'   shiny-app.
#' @return an HTML skeleton
#'
#' @note This function only initialises the widget, to have it display some
#'   data, [updateMarkerSettingsTable()] must be called with an appropriate
#'   initial state, and as long as updates are needed.
markerSettingsTable = function(inputId) {
  tagList(
    shiny::singleton(tags$head(
      tags$script(src = 'markerSettingsTableComponent/js/binding.js'),
      tags$link(rel = 'stylesheet',
                type = 'text/css',
                href = 'markerSettingsTableComponent/css/style.css')
    )),

    tags$table(
      tags$thead(
        tags$tr(
          tags$th('Marker'),
          tags$th('Sex-linked?'),
          tags$th('Mutation model'),
          tags$th('Include in calculation?')
        )
      ),
      tags$tbody(), # table body is built on the client, through the shiny binding
      class = 'table table-bordered table-condensed marker-settings-table',
      id = inputId
    )
  )
}

# send an input message to the client-side widget with new marker data when the
# user updates another part of the input.
updateMarkerSettingsTable = function(session, inputId, ped) {
  session$sendInputMessage(inputId, pedToMarkerSettingsState(ped))
}

#' Returns the list of markers included in the calculation, according to the
#' marker settings table.
#'
#' @param markerSettings a `data.frame` returned by the markerSettings module
#' @return a list of marker names (strings) included in the calculation
getIncludedMarkers = function(markerSettings) {
  if (!is.data.frame(markerSettings)) return(list());

  markerSettings[markerSettings$includeInCalculation == TRUE, ]$markerName
}
