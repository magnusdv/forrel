#' Exclusion/inclusion power plots
#'
#' This function offers three different visualisations of exclusion/inclusion
#' powers, particularly for missing person cases.
#'
#' The plot types are as follows:
#'
#' `type = 1`: x = Exclusion power; y = Inclusion power
#'
#' `type = 2`: x = Exclusion odds ratio; y = Inclusion odds ratio
#'
#' `type = 3`: x = Expected number of exclusions; y = average log(LR)
#'
#' For each `EPresult` object in `ep`, and the corresponding element of `ip`,
#' the relevant data is extracted from each, producing a single point the final
#' plot.
#'
#' In the most general case `ep` (and similarly for `ip`) can be a list of lists
#' of `EPresult` objects. To simplify the discussion we refer to the inner lists
#' as "groups". A group may consist of a single point, or several (typically
#' many simulations of the same situation). Points within the same group are
#' always drawn with the same color and shape.
#'
#' When plotting several groups, two sets of points are drawn:
#'
#' * Major points: Group means.
#'
#' * Minor points: Individual points in groups with more than one element.
#'
#' @param ep,ip Lists of equal length, with outputs from one or more runs of
#'   [missingPersonEP()] and [missingPersonIP()] respectively. Alternatively,
#'   `ep` can be a single output from [MPPsims()], in which case `ip` should be
#'   NULL. See Examples.
#' @param type Plot type; either 1, 2 or 3.
#' @param ellipse A logical. If TRUE, data ellipsis are drawn for each inner
#'   list of `ep`/`ip` containing more than 1 element.
#' @param col A color vector, recycle to match the top level length of `ep`.
#' @param labs A character of the same length as `ep`. If NULL, the names of
#'   `ep` are used, if present.
#' @param alpha Transparency for minor points (see Details).
#' @param shape Either "circle", "square", "diamond", "triangleUp" or
#'   "triangleDown", determining the shapes of both minor and major points.
#' @param size Point size.
#' @param hline,vline Single numbers indicating positions for
#'   horizontal/vertical "threshold" lines. When `type = 1`, these determine the
#'   shaded vertical regions, by default starting at 0.95.
#' @param xlim,ylim Axis limits; automatically chosen if NULL.
#' @param xlab,ylab Axis labels; automatically chosen if NULL.
#' @param legendOrder A permutation of `1,2,...,L`, where `L=length(ep)`.
#'
#' @return A `ggplot2` plot object.
#' @seealso [missingPersonEP()], [missingPersonEP()], [MPPsims()]
#'
#' @examples
#'
#' ref = nuclearPed(father = "fa", mother = "mo", child = "MP")
#' ref = setMarkers(ref, marker(ref, alleles = 1:5))
#'
#' # Alternatives for genotyping
#' sel = list("fa", c("fa", "mo"))
#'
#' # Simulate power for each selection
#' simData = MPPsims(ref, selections = sel, nProfiles = 10,
#'                   lrSims = 10, thresholdIP = 2)
#'
#' # Power plot 1: EP vs IP
#' powerPlot(simData, type = 1)
#'
#' # Change shape, and modify legend order
#' powerPlot(simData, type = 1, shape = "square", legendOrder = c(1,3,2))
#'
#' # Zoom in, and adjust the blue strips
#' powerPlot(simData, type = 1, xlim = c(0.5, 1), ylim = c(0.5, 1),
#'           hline = 0.8, vline = 0.9)
#'
#' # Power plot 3: Expected number of exclusions vs E[LR]
#' powerPlot(simData, type = 3)
#'
#' # With horizontal/vertical lines
#' powerPlot(simData, type = 3, hline = log10(2), vline = 1)
#'
#' # Plot 4: Illustrating the general inequality ELR > 1/(1-EP)
#' powerPlot(simData, type = 4)
#'
#' @importFrom stats aggregate
#' @export
powerPlot = function(ep, ip, type = 1, ellipse = FALSE, col = NULL, labs = NULL, alpha = 1,
                     shape = "circle", size = 2, hline = NULL, vline = NULL, xlim = NULL, ylim = NULL,
                     xlab = NULL, ylab = NULL, legendOrder = NULL) {
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop2("Package `ggplot2` is not installed. Please install this and try again.")

  if(inherits(ep, "MPPsim")) {
    ip = lapply(ep, "[[", "ip")
    ep = lapply(ep, "[[", "ep")
  }

  if(is.null(labs) && is.list(ep) && !isEP(ep))
    labs = names(ep)

  if(isEP(ep)) ep = list(ep)
  if(isIP(ip)) ip = list(ip)


  ### Ensure each element of ep (ip) is a *list* of EPresult objects.
  for(i in seq_along(ep)) {
    if(isEP(ep[[i]])) ep[[i]] = list(ep[[i]])
  }
  for(i in seq_along(ip)) {
    if(isIP(ip[[i]])) ip[[i]] = list(ip[[i]])
  }

  L = lengths(ep)
  if(!all(L == lengths(ip)))
    stop2("Arguments `ep` and `ip` are incompatible")

  ### Group labels legend ordering
  if(is.null(labs))
    labs = seq_along(L)
  if(is.null(legendOrder))
    legendOrder = seq_along(labs)
  group = factor(rep(labs, times = L), levels = labs[legendOrder])

  ### COlors
  if(is.null(col)) {
    col = c("lightgreen", "firebrick1", "deepskyblue", "#FFFF33", "gray70", "#F781BF", "cyan", "wheat", "#FF7F00")

    # Use white for "Baseline", if present
    if(!is.na(baseIdx <- match("Baseline", labs)))
      col = append(col, "white", after = baseIdx - 1)
  }

  col = setNames(rep_len(col, length(labs)), labs)

  ### Point shapes of minor and major points

  majorShapes = c(circle=21, square=22, diamond=23, triangleUp=24, triangleDown=25)
  minorShapes = c(circle=1, square=0, diamond=5, triangleUp=2, triangleDown=6)
  shapeIdx = pmatch(shape, names(majorShapes), duplicates.ok = TRUE)
  if(anyNA(shapeIdx))
    stop2("Unrecognized `shape` entry: ", shape[is.na(shapeIdx)],
          "\nMust match uniquely to one of: ", choices)
  shapeIdx = rep_len(shapeIdx, length(labs))

  # Major points: Use as group aestethic
  shapeMapMajor = setNames(majorShapes[shapeIdx], labs)

  # Minor points: Add column; map explicitly
  shapeMapMinor = setNames(minorShapes[shapeIdx], labs)

  ### Extract numbers and build data frame with plotting data
  if(type == 1) {
    epnum = vapply(unlist(ep, recursive = FALSE), function(a) a$EPtotal, FUN.VALUE = 1)
    ipnum = tryCatch(vapply(unlist(ip, recursive = F), function(a) a$IP, FUN.VALUE = 1),
                     error = function(e) stop2("Missing IP data. Are you sure you wanted power plot type 1?"))


    alldata = data.frame(ep = epnum, ip = ipnum, group = group)
    if(is.null(xlab)) xlab = "Exclusion power"
    if(is.null(ylab)) ylab = "Inclusion power"
    if(is.null(xlim)) xlim = c(0, 1)
    if(is.null(ylim)) ylim = c(0, 1)
  }
  else if(type == 2) {
    epnum = vapply(unlist(ep, recursive = F), function(a) a$EPtotal, FUN.VALUE = 1)
    ipnum = tryCatch(vapply(unlist(ip, recursive = F), function(a) a$IP, FUN.VALUE = 1),
                     error = function(e) stop2("Missing IP data. Are you sure you wanted power plot type 2?"))

    # Odds ratios: Assumes first entry is baseline!
    ep.odds = epnum/(1 - epnum)
    ip.odds = ipnum/(1 - ipnum)
    ep.OR = ep.odds/ep.odds[1]
    ip.OR = ip.odds/ip.odds[1]

    alldata = data.frame(ep = ep.OR, ip = ip.OR, group = group)

    if(is.null(xlab)) xlab = "Exclusion odds ratio"
    if(is.null(ylab)) ylab = "Inclusion odds ratio"
    if(is.null(xlim)) xlim = c(0, max(ep.OR))
    if(is.null(ylim)) ylim = c(0, max(ip.OR))
  }
  else if(type == 3) {
    epnum = vapply(unlist(ep, recursive = F), function(a) a$expectedMismatch, FUN.VALUE = 1)
    ipnum = vapply(unlist(ip, recursive = F), function(a) a$meanLogLR, FUN.VALUE = 1)

    alldata = data.frame(ep = epnum, ip = ipnum, group = group)

    if(is.null(xlab)) xlab = "Expected number of exclusions"
    if(is.null(ylab)) ylab = expression(log[10]~LR)
    if(is.null(xlim)) xlim = c(min(0, epnum), max(1, epnum))
    if(is.null(ylim)) ylim = c(min(0, ipnum), max(1, ipnum))
  }
  else if(type == 4) {
    epnum = vapply(unlist(ep, recursive = F), function(a) a$EPtotal, FUN.VALUE = 1)
    ipnum = vapply(unlist(ip, recursive = F), function(a) a$meanLR, FUN.VALUE = 1)

    alldata = data.frame(ep = epnum, ip = ipnum, group = group)

    if(is.null(xlab)) xlab = "Exclusion power"
    if(is.null(ylab)) ylab = "LR"
    if(is.null(xlim)) xlim = c(min(0, epnum), max(1, epnum))
    if(is.null(ylim)) ylim = c(min(0, ipnum), max(1, ipnum))
  }


  ### Major data points: Mean points of each group
  major = aggregate(. ~ group, alldata, mean)

  ### Minor data points: Individual entries in groups with more than one
  minor = subset(alldata, table(group)[group] > 1)

  ### Plot
  p = ggplot2::ggplot() +
    ggplot2::aes(x = ep, y = ip) +
    ggplot2::geom_point(data = minor, ggplot2::aes(colour = group), size = size,
                        shape = shapeMapMinor[minor$group], alpha = alpha) +
    ggplot2::geom_point(data = major, ggplot2::aes(fill = group, shape = group),
                        size = 2*size, colour = "black", stroke = 1.5) +
    ggplot2::labs(x = xlab, y = ylab, fill = NULL, colour = NULL) +
    ggplot2::scale_colour_manual(values = col) +
    ggplot2::scale_fill_manual(values = col) +
    ggplot2::scale_shape_manual(values = shapeMapMajor) +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::guides(colour = FALSE,
                    fill = ggplot2::guide_legend(title = "", reverse = TRUE),
                    shape = ggplot2::guide_legend(title = "", reverse = TRUE)
                    ) +
    ggplot2::coord_cartesian(clip = 'off') +
    ggplot2::theme_bw(base_size = 14)

  if(ellipse) {
    p = p + ggplot2::stat_ellipse(data = minor, na.rm = TRUE)
  }

  if(type == 1) {
    xmin = xlim[1]; if(is.na(xmin) || xmin > 0) xmin = -Inf
    ymin = ylim[1]; if(is.na(ymin) || ymin > 0) ymin = -Inf
    ethresh = if(is.null(vline)) 0.95 else vline
    ithresh = if(is.null(hline)) 0.95 else hline
    np = length(p$layers)
    p = p +
      ggplot2::annotate("rect", xmin = ethresh, xmax = 1, ymin = ymin, ymax = 1,  fill = "lightblue", alpha = 0.3) +
      ggplot2::annotate("rect", xmin = xmin, xmax = 1, ymin = ithresh, ymax = 1,  fill = "lightblue", alpha = 0.3)
    p$layers[] = p$layers[c(np + 1:2, 1:np)]
  }
  if(type == 2) {
    # Nothing to do
  }
  if(type == 3) {
    # Nothing to do
  }
  if(type == 4) {
    epvec = seq(xlim[1], xlim[2], length = 100)[-100]
    asympt = data.frame(ep = epvec, ip = 1/(1 - epvec))
    np = length(p$layers)
    p = p +
      ggplot2::geom_line(data = asympt, ggplot2::aes(colour = NULL, fill = NULL))
      #ggplot2::annotate("text", asympt$ep[1], asympt$ip[1], label = "ELR == frac(1, 1 - EP)",
      #                  hjust = -0.05, vjust = -0.5, parse = TRUE)
      p$layers[] = p$layers[c(np + 1, 1:np)]
    p = suppressMessages(p + ggplot2::scale_y_log10())
  }

  if(type > 1) {
    if(!is.null(hline)) {
      np = length(p$layers)
      p = p + ggplot2::geom_hline(yintercept = hline, linetype = 2, size = 1)
      p$layers[] = p$layers[c(np + 1, 1:np)]
    }
    if(!is.null(vline)) {
      np = length(p$layers)
      p = p + ggplot2::geom_vline(xintercept = vline, linetype = 2, size = 1)
      p$layers[] = p$layers[c(np + 1, 1:np)]
    }
  }

  p
}

