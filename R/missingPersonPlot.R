#' Missing person plot
#'
#' Plot the competing hypotheses of a family reunion case. A plot with two
#' panels is generated. The left panel shows a pedigree in which the _person of
#' interest_ (POI) is identical to the _missing person_ (MP). The right panel
#' shows the situation where these two are unrelated.
#'
#' A standard family reunification case involves the following ingredients:
#'
#' * A family in which a single member ("MP") is missing.
#'
#' * Some of the family members have been genotyped
#'
#' * A person of interest ("POI") may be the missing person.
#'
#' After genotyping of POI, the genetic evidence is typically assessed by
#' computing the likelihood ratio of the following hypotheses:
#'
#' * H1: POI is MP
#'
#' * H2: POI is unrelated to the family
#'
#' The goal of this function is to illustrate the above hypotheses, using
#' labels, colours and shading to visualise the different aspects of the situation.
#'
#' This function cannot handle cases with more complicated hypotheses (e.g.
#' multiple missing persons, or where H2 specifies a different relationship).
#' However, as it is basically a wrapper of [pedtools::plotPedList()], an
#' interested user should be able to extend the source code to such cases
#' without too much trouble.
#'
#' @param reference A [pedtools::ped()] object
#' @param missing The ID label of the missing pedigree member
#' @param id.labels A character vector with labels for the pedigree members. See
#'   [pedtools::plot.ped()] for
#' @param MP.label The label of the missing member. Default: "MP"
#' @param POI.label The label of the person of interest. Default: "POI"
#' @param marker Optional vector of marker indices to be included in the plot
#' @param shaded A vector of ID labels indicating who should appear with shaded
#'   symbols. By default, all typed members.
#' @param POI.col The plot colour of POI. Default: red
#' @param POI.shaded A logical: If TRUE (default), the POI is plotted with a
#'   shaded symbol
#' @param POI.height A numeric controlling the vertical placement of the POI
#'   singleton (in the right panel)
#' @param width A positive number controlling the width of the plot. More
#'   specifically this number is the relative width of the reference pedigree,
#'   compared to a singleton. Default: 4
#' @param newdev A logical: If TRUE the plot is created in a new plot window
#' @param frametitles A character of length 2, with subtitles for the two frames
#' @param ... Extra parameters passed on to [pedtools::plotPedList()]
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", children = c("b1", "b2"))
#'
#' # Default plot
#' missingPersonPlot(x, missing = "b2")
#'
#' # A bit nicer using various options
#' missingPersonPlot(x, missing = "b2", MP.label = "Missing", id.label = NULL,
#'                   shaded = "b1", POI.shaded = TRUE,
#'                   width = 2,      # adjust internal spacing (see above)
#'                   dev.width = 7,  # device width (see ?plotPedList())
#'                   dev.height = 3, # device height (see ?plotPedList())
#'                   fmar = 0.02,    # adjust frame margin (see ?plotPedList())
#'                   cex = 1.5,      # larger symbols and label font (see ?par())
#'                   cex.main = 1.3  # larger frame titles (see ?par())
#'                   )
#'
#' @importFrom stats setNames
#' @export
missingPersonPlot = function(reference, missing, id.labels = labels(reference),
                             MP.label = "MP", POI.label = "POI", marker = NULL,
                             shaded = "typed", POI.col = "red", POI.shaded = FALSE,
                             POI.height = 8, width = 4, newdev = TRUE,
                             frametitles = c(expression(H[1] * ": POI is MP"),
                                           expression(H[2] * ": POI unrelated")),
                             ...) {

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")
  if(length(POI.label) != 1 || POI.label == "")
    stop2("`POI.label` must be a non-empty character string")
  if(!is_number(width, minimum = 1))
    stop2("`width` must be a number larger than 1")

  nInd = pedsize(reference)
  if(identical(id.labels, "num")) {
    id.labels = labels(reference)
    names(id.labels) = as.character(1:nInd)
  }

  # Shading
  if(identical(shaded, "typed"))
    shaded = typedMembers(reference)

  ### Hypothesis 1: Related
  ped_related = reference

  # Labels
  mp_poi = if(MP.label != "") paste(MP.label, "=", POI.label) else POI.label
  mp_label = setNames(missing, mp_poi)
  if (is.null(id.labels) || identical(id.labels, ""))
    labs = mp_label
  else
    labs = c(mp_label, id.labels[id.labels != missing])

  # Colour MP = POI red
  col1 = list(red = missing)

  # Shading
  shaded1 = shaded
  if(POI.shaded)
    shaded1 = c(shaded1, missing)

  # Build plot 1
  plot1 = list(ped_related, id.labels = labs, col = col1, shaded = shaded1)


  ### Hypothesis 2: Unrelated
  # First part

  # MP label: Either full vector or
  names(labs)[match(missing, labs)] = MP.label

  plot2 = list(reference, id.labels = labs, shaded = shaded)

  # Second part: POI singleton
  s = singleton(id = missing, sex = getSex(reference, missing))

  if(!is.null(marker))
    s = transferMarkers(from = reference, to = s)

  # Relabel POI (above MP label was used to enable transfer/sex)
  s = relabel(s, POI.label)
  plot3 = list(s, col = POI.col, shaded = if(POI.shaded) POI.label,
               margins = c(POI.height,0,0,0))

  ### Plot
  plotPedList(list(plot1, plot2, plot3),
              widths = c(width, width, 1),
              frames = list(1, 2:3),
              frametitles = frametitles,
              marker = marker,
              skip.empty.genotypes = TRUE,
              newdev = newdev, ...)
}
