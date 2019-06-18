#' GUI launcher
#'
#' Launches the Shiny-based graphical user interface for the exclusion power functionality.
#'
#' @export

launchGUI = function() {
  shiny::runApp(system.file('exclusionPowerUI', package='forrel'))
}
