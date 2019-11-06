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
#' For each `mpEP` object in `ep`, and the corresponding element of `ip`, the
#' relevant data is extracted from each, producing a single point the final plot.
#'
#' In the most general case `ep`
#' (and similarly for `ip`) can be a list of lists of `mpEP` objects. To
#' simplify the discussion we refer to the inner lists as "groups". A group
#' may consist of a single point, or several (typically many simulations of the same
#' situation). Points within the same group are always drawn with the same
#' color and shape.
#'
#' When plotting several groups, two sets of points are drawn:
#'
#' * Major points: Group means.
#'
#' * Minor points: Individual points in groups with more than one element.
#'
#'
#' @param ep Exclusion power data, typically in the form of one or several
#'   `mpEP` objects (output from [missingPersonEP()]. See Details and Examples.
#' @param ip Inclusion power data, typically in the form of one or several
#'   `mpIP` objects (output from [missingPersonIP()]. See Details and Examples.
#'   The list structure of `ip` must be identical to that of `ep`.
#' @param type Plot type; either 1, 2 or 3.
#' @param LRthresh A single positive number; the LR threshold for "inclusion".
#' @param ellipse A logical. If TRUE, data ellipsis are drawn for each inner
#'   list of `ep`/`ip` containing more than 1 element.
#' @param col A color vector, recycle to match the top level length of `ep`.
#' @param labs A character of the same length as `ep`. If NULL, the names of
#'   `ep` are used, if present.
#' @param alpha Transparency for minor points (see Details).
#' @param shape,size Plotting parameters
#' @param xlim,ylim Axis limits; automatically chosen if NULL.
#' @param xlab,ylab Axis labels; automatically chosen if NULL.
#'
#' @return
#' @export
#'
#' @examples
powerPlot = function(ep, ip, type = 1, LRthresh = 1e4, ellipse = TRUE, col = NULL,
                                  labs = NULL, alpha = 1, shape = 1, size = 2,
                                  xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL) {
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop2("Package `ggplot2` is not installed. Please install this and try again.")

  if(is.null(labs))
    labs = names(ep)

  if(isEP(ep)) ep = list(ep)
  if(isIP(ip)) ip = list(ip)

  # Ensure each element of ep (ip) is a *list* of mpEP objects.
  for(i in seq_along(ep)) {
    if(isEP(ep[[i]])) ep[[i]] = list(ep[[i]])
  }
  for(i in seq_along(ip)) {
    if(isIP(ip[[i]])) ip[[i]] = list(ip[[i]])
  }

  L = lengths(ep)
  if(!all(L == lengths(ip)))
    stop2("Arguments `ep` and `ip` are incompatible")

  # Group labels and colors
  if(is.null(labs))
    labs = seq_along(L)
  group = as.factor(rep(labs, times = L))

  if(is.null(col)) {
    Set1 = c("#377EB8", "#4DAF4A", "#E41A1C", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    col = Set1[seq_along(L)]
  }

  # Extract numbers and build data frame with plotting data
  if(type == 1) {
    epnum = vapply(unlist(ep, recursive = F), function(a) a$EPtotal, FUN.VALUE = 1)
    ipnum = vapply(unlist(ip, recursive = F), function(a) a$IP, FUN.VALUE = 1)
    alldata = data.frame(ep = epnum, ip = ipnum, group = group)
    if(is.null(xlab)) xlab = "Exclusion power"
    if(is.null(ylab)) ylab = "Inclusion power"
    if(is.null(xlim)) xlim = c(0, 1)
    if(is.null(ylim)) ylim = c(0, 1)
  }
  else if(type == 2) {
    epnum = vapply(unlist(ep, recursive = F), function(a) a$EPtotal, FUN.VALUE = 1)
    ipnum = vapply(unlist(ip, recursive = F), function(a) a$IP, FUN.VALUE = 1)

    # Odds ratios: Assumes first entry is baseline!
    ep.odds = epnum/(1 - epnum)
    ip.odds = ipnum/(1 - ipnum)
    ep.OR = ep.odds/ep.odds[1]
    ip.OR = ip.odds/ip.odds[1]

    alldata = data.frame(ep = ep.OR, ip = ip.OR, group = group)

    if(is.null(xlab)) xlab = "Exclusion odds ratio"
    if(is.null(ylab)) ylab = "Inclusion odds ratio"
    if(is.null(xlim)) xlim = c(0, NA)
    if(is.null(ylim)) ylim = c(0, NA)
  }
  else if(type == 3) {
    epnum = vapply(unlist(ep, recursive = F), function(a) a$expectedMismatch, FUN.VALUE = 1)
    ipnum = vapply(unlist(ip, recursive = F), function(a) a$meanLogLR, FUN.VALUE = 1)

    alldata = data.frame(ep = epnum, ip = ipnum, group = group)

    if(is.null(xlab)) xlab = "Expected number of exclusions"
    if(is.null(ylab)) ylab = expression("Average "~log[10](LR))
    if(is.null(xlim)) xlim = c(0, if(max(epnum) < 1) 1 else NA)
    if(is.null(ylim)) ylim = c(0, if(max(ipnum) < 1) 1 else NA)
  }

  # Major data points: Mean points of each group
  major = aggregate(. ~ group, alldata, mean)

  # Minor data points: Individual entries in groups with more than one
  minor = subset(alldata, table(group)[group] > 1)

  # Plot
  p = ggplot(NULL, aes(x = ep, y = ip, color = group, fill = group)) +
    geom_point(data = minor, size = size, shape = shape, alpha = alpha) +
    geom_point(data = major, size = 2*size, shape = 21, col = 1, stroke = 1.5) +
    labs(x = xlab, y = ylab, fill = NULL, color = NULL) +
    guides(color = F) +
    scale_color_manual(limits = labs, values = col) +
    scale_fill_manual(limits = labs, values = col) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    coord_cartesian(clip = 'off') +
    theme_bw(base_size = 14)

  if(ellipse) {
    p = p + stat_ellipse(data = minor, na.rm = TRUE)
  }

  np = length(p$layers)

  if(type == 1) {
    xmin = xlim[1]; if(is.na(xmin) || xmin > 0) xmin = -Inf
    ymin = ylim[1]; if(is.na(ymin) || ymin > 0) ymin = -Inf
    p = p +
      annotate("rect", xmin = 0.95, xmax = 1, ymin = ymin, ymax = 1,  fill = "lightblue", alpha = 0.3) +
      annotate("rect", xmin = xmin, xmax = 1, ymin = 0.95, ymax = 1,  fill = "lightblue", alpha = 0.3)
    p$layers[] = p$layers[c(np + 1:2, 1:np)]
  }
  if(type == 2) {
    p = p +
      geom_hline(yintercept = 1, linetype = 2, size = 1) +
      geom_vline(xintercept = 1, linetype = 2, size = 1)
    p$layers[] = p$layers[c(np + 1:2, 1:np)]
  }
  if(type == 3) {
    p = p +
      geom_hline(yintercept = log10(LRthresh), linetype = 2, size = 1) +
      geom_vline(xintercept = 1, linetype = 2, size = 1)
    p$layers[] = p$layers[c(np + 1:2, 1:np)]
  }
  p
}

