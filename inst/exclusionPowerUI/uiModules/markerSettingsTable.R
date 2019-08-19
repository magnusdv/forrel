# mutationModels = tagList(tags$option('None', value = 'None'),
#                          tags$option('Equal', value = 'Equal'))

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

# buildTableRow = function(markerInfo) {
#   tagList(
#     tags$tr(
#       tags$td(markerInfo$markerName),
#       tags$td(
#         tags$input(type = 'checkbox',
#                    value = if (markerInfo$chrom == 23) '1' else '0',
#                    class = 'sexLinked-checkbox')),
#       tags$td(
#         tags$select(mutationModels,
#                     class = 'mutationModel-select')),
#       tags$td(
#         tags$input(type = 'checkbox',
#                    value = markerInfo$includeInCalculation,
#                    class = 'includeInCalculation-checkbox')),
#       'data-marker' = markerInfo$markerName
#     )
#   )
# }

markerSettingsTable = function(inputId, initialState) {
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
      tags$tbody(),
      class = 'table table-bordered table-condensed marker-settings-table',
      id = inputId
    )
  )
}

updateMarkerSettingsTable = function(session, inputId, ped) {
  session$sendInputMessage(inputId, pedToMarkerSettingsState(ped))
}

getIncludedMarkers = function(markerSettings) {
  if (!is.data.frame(markerSettings)) return(list());

  markerSettings[markerSettings$includeInCalculation == TRUE, ]$markerName
}
