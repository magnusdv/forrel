#' Missing person plot
#'
#' Visualises the competing hypotheses of a family reunion case. A plot with two
#' panels is generated. The left panel shows a pedigree in which the _person of
#' interest_ (POI) is identical to the _missing person_ (MP). The right panel
#' shows the situation where these two are unrelated. See Details for further
#' explanations.
#'
#' A standard family reunification case involves the following ingredients:
#'
#' * A reference family with a single missing person ("MP").
#'
#' * Some of the family members have been genotyped
#'
#' * A person of interest ("POI") is to be matched against the reference family
#'
#' After genotyping of POI, the genetic evidence is typically assessed by
#' computing the likelihood ratio of the following hypotheses:
#'
#' * H1: POI is MP
#'
#' * H2: POI is unrelated to the family
#'
#' The goal of this function is to illustrate the above hypotheses, using
#' labels, colours and shading to visualise the different aspects of the
#' situation.
#'
#' This function cannot handle cases with more complicated hypotheses (e.g.
#' multiple missing persons, or where H2 specifies a different relationship).
#' However, as it is basically a wrapper of [pedtools::plotPedList()], an
#' interested user should be able to extend the source code to such cases
#' without too much trouble.
#'
#' @param reference A [pedtools::ped()] object.
#' @param missing The ID label of the missing pedigree member.
#' @param labs A character vector with labels for the pedigree members. See
#'   [pedtools::plot.ped()].
#' @param marker Optional vector of marker indices to be included in the plot.
#' @param hatched A vector of ID labels indicating who should appear with
#'   hatched symbols in the plot. By default, all typed members.
#' @param MP.label,POI.label Custom labels of the missing person and the POI.
#'   Default: "MP" and "POI".
#' @param MP.col,POI.col Fill colours for MP and POI.
#' @param POI.sex The sex of POI. This defaults to that of the missing person,
#'   but may be set explicitly. This is particularly useful when the missing
#'   person has unknown sex.
#' @param POI.hatched Deprecated (ignored).
#' @param titles A character of length 2, with subtitles for the two frames.
#' @param width A positive number controlling the width of the plot. More
#'   specifically this number is the relative width of the reference pedigree,
#'   compared to a singleton.
#' @param cex Expansion factor for pedigree symbols and font size.
#' @param ... Extra parameters passed on to [pedtools::plotPedList()].
#'
#' @return None
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", children = c("b1", "b2"))
#'
#' # Default plot
#' missingPersonPlot(x, missing = "b2")
#'
#' # Open in separate window; explore various options
#' missingPersonPlot(x,
#'                   missing = "b2",
#'                   hatched = "b1",
#'                   deceased = c("fa", "mo"),
#'                   cex = 1.5,      # larger symbols and labels (see ?par())
#'                   cex.main = 1.3, # larger frame titles (see ?par())
#'                   dev.width = 7,  # device width (see ?plotPedList())
#'                   dev.height = 3  # device height (see ?plotPedList())
#'                   )
#'
#' @importFrom stats setNames
#' @export
missingPersonPlot = function(reference, missing, labs = labels(reference),
                             marker = NULL, hatched = typedMembers(reference),
                             MP.label = "MP", POI.label = "POI",
                             MP.col = "#FF9999", POI.col = "lightgreen",
                             POI.sex = getSex(reference, missing),
                             POI.hatched = NULL,
                             titles = c(expression(H[1] * ": POI = MP"),
                                        expression(H[2] * ": POI unrelated")),
                             width = NULL, cex = 1.2, ...) {

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")
  if(length(POI.label) != 1 || POI.label == "")
    stop2("`POI.label` must be a non-empty character string")

  nInd = pedsize(reference)
  if(identical(labs, "num"))
    labs = setNames(labels(reference), as.character(1:nInd))

  ### Hypothesis 1: Related
  ped_related = reference

  # Ensure MP has same sex as POI (relevant if MP has unknown sex)
  ped_related = setSex(ped_related, ids = missing, sex = POI.sex)

  # Labels
  mp_poi = if(MP.label != "") sprintf("%s=%s", POI.label, MP.label) else POI.label
  mp_label = setNames(missing, mp_poi)
  if (is.null(labs) || identical(labs, ""))
    labs1 = mp_label
  else
    labs1 = c(mp_label, labs[labs != missing])

  # Shading: Remove MP/POI
  hatched = .mysetdiff(hatched, missing)

  # Colour of MP/POI in first ped
  fill1 = setNames(POI.col, missing)

  # Build plot 1
  plot1 = list(ped_related, labs = labs1, fill = fill1, hatched = hatched)


  ### Hypothesis 2: Unrelated

  # First part
  if(missing %in% labs)  # fix MP label
    names(labs)[match(missing, labs)] = MP.label

  # Fill color MP
  fill2 = setNames(MP.col, missing)

  plot2 = list(reference, labs = labs, hatched = hatched, fill = fill2)

  # Second part: POI singleton
  s = singleton(id = missing, sex = POI.sex)

  if(!is.null(marker))
    s = transferMarkers(from = reference, to = s)

  fill3 = setNames(POI.col, POI.label)

  # Relabel POI (above MP label was used to enable transfer/sex)
  s = relabel(s, POI.label)
  plot3 = list(s, fill = fill3)

  widths = if(!is.null(width)) c(width, width, 1) else NULL

  ### Plot
  plotPedList(list(plot1, plot2, plot3),
              widths = widths,
              groups = list(1, 2:3),
              titles = titles,
              marker = marker,
              cex = cex, ...)
}
