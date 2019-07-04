#' Advanced Table File Loader Input
#'
#' UI function for the Advanced Table File Loader shiny module.
#'
#' The Advanced Table File Loader shiny module makes it easy for users to load
#' tabular data from delimited files by allowing them to customise loading
#' parameters such as column separator, header presence and more. When invoked,
#' it shows an interactive modal with a preview of the table that is being
#' loaded that is updated when the user modifies loading settings. The module
#' returns a dataframe from the data loaded using the user-specified parameters.
#'
#' This UI function must be included inside the `ui.R` file or `ui()` function of
#' the *shiny* interface. See examples. Also, the server part of this module,
#' [advancedTableFileLoader()] must be included in the server section of the
#' shiny application.
#'
#' @param id a stirng that is unique among the widgets used in the interface
#' @param label the string shown as a link to open the modal dialog for loading
#'   the file
#'
#' @author Elias Hernandis
#'
#' @seealso [advancedTableFileLoader()] for the server part of this module
#'
#' @example \dontrun{ # Include the following in the relvant part of the UI
#' code: advancedTableFileLoaderInput('referenceFiles', 'Select file')
#'
#' # Include the call to the server function in the server part of the code }
#'
advancedTableFileLoaderInput <- function(id, label = 'Table loader') {
  # create a namespace to avoid name collisions with other modules / other instances of this one
  ns <- NS(id)

  tagList(
    actionLink(ns('modalLink'), label = label)
  )
}

#' Advanced Table File Loader
#'
#' Server function for the Advanced Table File Loader shiny module.
#'
#' The Advanced Table File Loader shiny module makes it easy for users to load
#' tabular data from delimited files by allowing them to customise loading
#' parameters such as column separator, header presence and more. When invoked,
#' it shows an interactive modal with a preview of the table that is being
#' loaded that is updated when the user modifies loading settings. The module
#' returns a dataframe from the data loaded using the user-specified parameters.
#'
#' This server function must be included inside the `server.R` file or `server()` function of
#' the *shiny* interface. See examples. Also, the UI part of this module,
#' [advancedTableFileLoaderInput()] must be included in the UI section of the
#' shiny application.
#'
#' @param input the input param of the *shiny* server function
#' @param output the output param of the *shiny* server function
#' @param session the session param of the *shiny* server function
#' @param id a string that is unique among the widget names of the UI. Must be the same as the one passed as `id` to the UI function [advancedTableFileLoaderInput()].
#' @param columnHeaders whether the option to treat the first row as headers should be enabled by default. Defaults to `FALSE`. See [read.table()].
#' @param rowHeaders whether the option to treat the first column as row headers should be enabled by default. When set to TRUE, [read.table()] is called with `row.headers = 1`.
#' @param sep the default value for the column separator. Defaults to a comma (`,`).
#' @param quote the default value for the quoting character. Defaults to a simple double quote (`"`).
#' @param na.strings the default value for the string indicating that a value should be treated as NA. Defaults to `NA`.
#' @param dec the default value for the decimal separator character. Defaults to a dot (`.`).
#' @param transpose whether the loaded table should be transposed. Defaults to `FALSE`.
#'
#' @author Elias Hernandis
#'
#' @seealso [advancedTableFileLoader()] for the server part of this module
#'
#' @example \dontrun{ # Include the following in the relvant part of the UI
#' code: advancedTableFileLoaderInput('referenceFiles', 'Select file')
#'
#' # Include the call to the server function in the server part of the code }
#'
advancedTableFileLoader <- function(input, output, session, id = 'namespace',
                                    columnHeaders = FALSE,
                                    rowHeaders = FALSE,
                                    sep = ',',
                                    quote = '"',
                                    na.strings = 'NA',
                                    dec = '.',
                                    transpose = FALSE) {
  ns <- NS(id)

  inputModal <- function() {
    modalDialog(
      fluidPage(
        fluidRow(
          column(12,
                 fileInput(ns('inputFile'), 'File', width = '100%'))
        ),
        fluidRow(
          column(12,
                 helpText('Select file, then adjust parameters until the table preview looks right. Table should contain one marker per column and one allele per row.'))
        ),
        fluidRow(
          column(4,
                 radioButtons(ns('sep'), 'Column separator', selected = sep, choices = c(', (Comma)' = ',', '; (Semicolon)' = ';', '<TAB>' = '\t', '<SPACE>' = ' '))),
          column(4,
                 radioButtons(ns('quote'), 'Cell quoting', selected = quote, choices = c('" (Double quotes)' = '"', "' (Single quotes)" = "'"))),
          column(4,
                 p(strong('Headers')),
                 checkboxInput(ns('columnHeaders'), 'Column headers', value = columnHeaders),
                 checkboxInput(ns('rowHeaders'), 'Row headers', value = rowHeaders))
        ),
        fluidRow(
          column(4,
                 radioButtons(ns('dec'), 'Decimal separator', selected = dec, choices = c('. (Dot)' = '.', ', (Comma)' = ','))),
          column(4,
                 textInput(ns('na.strings'), 'NA strings', value = na.strings),
                 helpText('String used in place of a value when none applies.')),
          column(4,
                 p(strong('Transpose')),
                 checkboxInput(ns('transpose'), 'Exchange rows and columns', value = transpose))
        ),
        fluidRow(
          column(12,
                 tableOutput(ns('tablePreview'))),
          column(12,
                 helpText('Only first rows are shown in table preview.'))
        )
      ),
      title = 'Load a table file',
      footer = tagList(
        modalButton('Cancel'), # no id required as it always dismisses the modal
        actionButton(ns('ok'), 'OK')
      )
    )
  }

  # load file when one is chose or the options change
  dataframe <- reactive({
    # don't do anything until the user chooses a file
    req(input$inputFile)

    table = read.table(input$inputFile$datapath,
                       header = input$columnHeaders,
                       sep = input$sep,
                       quote = input$quote,
                       row.names = if (input$rowHeaders) 1 else NULL,
                       dec = input$dec,
                       na.strings = input$na.strings)

    if (input$transpose) t(table) else table
  })

  # display file preview
  output$tablePreview <- renderTable(head(dataframe()), rownames = TRUE, decimals = 4)

  # show modal when action link is clicked
  observeEvent(input$modalLink, {
    showModal(inputModal())
  })

  # dismiss modal when the user clicks OK
  observeEvent(input$ok, {
    removeModal()
  })

  return(dataframe)
}
