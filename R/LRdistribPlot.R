#' Plot LR distributions under two hypotheses
#'
#' Simulates the LR H1/H2 under each hypothesis and plots the two log10 LR
#' distributions, including their density overlap.
#'
#' @param numeratorPed,denominatorPed Pedigrees describing H1 and H2.
#' @param ids Individuals to simulate.
#' @param markers Marker names or indices to include, or a named list of
#'   frequency vectors defining new markers. By default all attached markers.
#' @param nsim Number of simulations under each hypothesis.
#' @param seed Integer seed for the random number generator.
#' @param data Precomputed data, either output from `LRdistribPlot(...,
#'   returnData = TRUE)` or a list of two `LRpowerResult` objects, with H1 true
#'   first and H2 true second.
#' @param returnData If TRUE, return the simulated log10 LRs instead of a plot.
#' @param title Plot title.
#' @param bw Density bandwidth on the log10 LR scale. By default it is estimated
#'   from the pooled simulations. Increase `bw` for smoother curves and decrease
#'   it to show more detail.
#' @param col Two colours for the H1-true and H2-true distributions.
#' @param verbose A logical.
#'
#' @return A `ggplot` object, or if `returnData = TRUE`, a data frame with
#'   columns `hypothesis`, `sim` and `log10LR`.
#'
#' @examples
#' if(requireNamespace("ggplot2", quietly = TRUE)) {
#'
#' db = NorwegianFrequencies[1:15]
#'
#' # Example 1: Sibs vs unrelated (increase nsim for better results)
#' ids = c("A", "B")
#' H1 = nuclearPed(children = ids)
#' H2 = singletons(ids)
#' LRdistribPlot(H1, H2, ids = ids, markers = db, nsim = 10, seed = 123)
#'
#' # Example 2: Full sibs vs half sibs
#' ids = c("A", "B")
#' H1 = nuclearPed(children = ids)
#' H2 = halfSibPed() |> relabel(old = 4:5, new = ids)
#' LRdistribPlot(H1, H2, ids = ids, markers = db, nsim = 10, seed = 123,
#'               title = "LR distributions for H1: Full sibs, H2: Half sibs")
#'
#' # Example 3: Full sibs vs half sibs, including shared parent
#' ids = c("A", "B", "C")
#' H1 = nuclearPed(fa = ids[1], children = ids[2:3])
#' H2 = halfSibPed() |> relabel(old = c(2,4:5), new = ids)
#' LRdistribPlot(H1, H2, ids = ids, markers = db, nsim = 10, seed = 123,
#'               title = "Full vs. half sibs, when parent is available")
#' }
#'
#' @export
LRdistribPlot = function(numeratorPed = NULL, denominatorPed = NULL, ids = NULL,
                         markers = NULL, nsim = 500, seed = NULL, data = NULL,
                         returnData = FALSE, title = "LR distributions", bw = NULL,
                         col = c("#E69F00", "#0072B2"), verbose = TRUE) {
  if(is.null(data)) {
    # Simulate under each hyp
    r1 = LRpower(numeratorPed, denominatorPed, truePed = numeratorPed,
                 ids = ids, markers = markers, source = "numerator",
                 nsim = nsim, seed = seed, verbose = verbose)
    r2 = LRpower(numeratorPed, denominatorPed, truePed = denominatorPed,
                 ids = ids, markers = markers, source = "numerator",
                 nsim = nsim, verbose = verbose)
    data = list(r1, r2)
  }

  # Convert LRpower output to plotting data
  if(is.list(data) && !is.data.frame(data)) {
    if(length(data) != 2 || !all(vapply(data, inherits, logical(1), "LRpowerResult")))
      stop2("`data` must contain two `LRpowerResult` objects")

    n = lengths(lapply(data, `[[`, "log10LRperSim"))
    nmark = length(data[[1]]$params$markers)
    data = data.frame(
      hypothesis = rep(c("H1", "H2"), n),
      sim = sequence(n),
      log10LR = c(data[[1]]$log10LRperSim, data[[2]]$log10LRperSim)
    )
    attr(data, "nMarkers") = nmark
  }

  if(returnData) return(data)

  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop2("Package `ggplot2` is needed for plotting")
  if(!all(c("hypothesis", "log10LR") %in% names(data)))
    stop2("Invalid `data` input")

  # Estimate finite densities on the same grid
  x1 = data$log10LR[data$hypothesis == "H1"]
  x2 = data$log10LR[data$hypothesis == "H2"]
  if(any(!is.finite(c(x1, x2))))
    stop2("Cannot plot distributions containing LR = 0 or Inf")

  bw = bw %||% stats::bw.nrd0(c(x1, x2))
  xr = range(x1, x2) + c(-3, 3) * bw
  d1 = stats::density(x1, bw = bw, from = xr[1], to = xr[2], n = 1024)
  d2 = stats::density(x2, bw = bw, from = xr[1], to = xr[2], n = 1024)

  d = data.frame(x = d1$x, H1 = d1$y, H2 = d2$y)
  d$overlap = pmin(d$H1, d$H2)
  ov = sum(d$overlap) * diff(d$x)[1]

  nsim = table(data$hypothesis)
  nmark = attr(data, "nMarkers")

  subtitParts = c(
    sprintf("%s sims per hypothesis", if(nsim[1] == nsim[2]) nsim[1] else toString(nsim)),
    sprintf("%d markers", attr(data, "nMarkers")),
    sprintf("Distribution overlap: %.1f%%", 100 * ov)
  )
  subtitle = paste(subtitParts, collapse = " | ")

  # Plot distributions and their overlap
  x = H1 = H2 = overlap = NULL
  p = ggplot2::ggplot(d, ggplot2::aes(x)) +
    ggplot2::geom_area(ggplot2::aes(y = H1, fill = "H1 true"), alpha = 0.5) +
    ggplot2::geom_area(ggplot2::aes(y = H2, fill = "H2 true"), alpha = 0.5) +
    ggplot2::geom_area(ggplot2::aes(y = overlap), fill = "grey50", alpha = 0.75) +
    ggplot2::scale_fill_manual(values = c("H1 true" = col[1], "H2 true" = col[2]),
                               name = NULL) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = expression(log[10](LR)), y = "Density") +
    ggplot2::theme_classic()

  p
}
