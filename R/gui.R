#' GUI launcher
#'
#' Launches the Shiny-based graphical user interface for the exclusion power functionality.
#'
#' @export

lauchGUI = function() {
  shiny::runApp(system.file('exclusionPowerUI', package='forrel'))
}
