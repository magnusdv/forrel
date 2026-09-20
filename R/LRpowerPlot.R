#' Plot LR distributions under two hypotheses
#'
#' Simulates the LR comparing two pedigree hypotheses and plots the two log10 LR
#' distributions, including their density overlap.
#'
#' @param numeratorPed,denominatorPed Pedigrees describing H1 and H2. If `denominatorPed`
#'   is NULL, an 'unrelated' hypothesis will be created as `singletons(ids)`.
#' @param ids Individuals to simulate.
#' @param markers Marker names or indices to include, or a named list of frequency vectors
#'   defining new markers. By default all attached markers.
#' @param nsim Number of simulations under each hypothesis.
#' @param seed Integer seed for the random number generator.
#' @param threshold An LR threshold. If given, the plot includes exceedance probabilities.
#'   Default: 10000.
#' @param data Precomputed data, either output from `LRpowerPlot(..., returnData = TRUE)`
#'   or a list of two `LRpowerResult` objects, with H1 true first and H2 true second.
#' @param returnData If TRUE, return the simulated log10 LRs instead of a plot.
#' @param title Plot title.
#' @param bw Density bandwidth on the log10 LR scale. By default it is estimated from the
#'   pooled simulations. Increase `bw` for smoother curves and decrease it to show more
#'   detail.
#' @param col Two colours for the H1-true and H2-true distributions.
#' @param verbose A logical.
#'
#' @return A `ggplot` object, or if `returnData = TRUE`, a data frame with columns
#'   `hypothesis`, `sim` and `log10LR`.
#'
#' @seealso [LRpower()]
#'
#' @examples
#' if(requireNamespace("ggplot2", quietly = TRUE)) {
#'
#' db = NorwegianFrequencies[1:10]
#'
#' ### Example 1: Sibs vs unrelated (increase nsim!)
#' ids = c("A", "B")
#' H1 = nuclearPed(children = ids)
#'
#' LRpowerPlot(H1, ids = ids, markers = db, nsim = 10, seed = 123)
#'
#'
#' ### Example 2: Full sibs vs half sibs
#' ids = c("A", "B")
#' H1 = nuclearPed(children = ids)
#' H2 = halfSibPed() |> relabel(old = 4:5, new = ids)
#'
#' LRpowerPlot(H1, H2, ids = ids, markers = db, nsim = 10, seed = 123,
#'             title = "H1: Full sibs, H2: Half sibs")
#'
#'
#' ### Example 3: Full sibs vs half sibs, including shared parent
#' ids = c("A", "B", "C")
#' H1 = nuclearPed(fa = ids[1], children = ids[2:3])
#' H2 = halfSibPed() |> relabel(old = c(2,4:5), new = ids)
#'
#' LRpowerPlot(H1, H2, ids = ids, markers = db, nsim = 10, seed = 123,
#'             title = "Full vs. half sibs, when parent is available")
#'
#'
#' # Example 4: Paternity case (requires mutation modelling!)
#' H1 = nuclearPed() |>
#'   setMarkers(locusAttributes = db) |>
#'   setMutmod(model = "equal", rate = 0.01)
#'
#' LRpowerPlot(H1, ids = c(1,3), nsim = 10, seed = 123)
#'
#' # Alternative syntax: With returnData = TRUE
#' dat = LRpowerPlot(H1, ids = c(1,3), nsim = 10, returnData = TRUE)
#'
#' LRpowerPlot(data = dat, threshold = 1e6, col = 2:3)
#'
#' }
#'
#' @export
LRpowerPlot = function(numeratorPed = NULL, denominatorPed = NULL, ids = NULL,
                         markers = NULL, nsim = 500, seed = NULL, threshold = 1e4,
                         data = NULL, returnData = FALSE, title = NULL,
                         bw = NULL, col = c("#E69F00", "#0072B2"), verbose = TRUE) {
  if(is.null(data)) {

    if(is.list(ids))
      stop2("`ids` must be a vector")

    if(is.null(denominatorPed)) {
      if(verbose)
        message(paste("Creating H2: Unrelated singletons", toString(ids)))
      denominatorPed = singletons(ids)
    }

    # Simulate under each hyp
    if(verbose) message("\n--- H1 ---")
    r1 = LRpower(numeratorPed, denominatorPed, truePed = numeratorPed,
                 ids = ids, markers = markers, source = "numerator",
                 nsim = nsim, seed = seed, verbose = verbose)
    if(verbose) message("\n--- H2 ---")
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

  # Estimate both densities on the same grid
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

  caption = sprintf("%d sims under each hypothesis  |  %d markers  |  Density overlap: %.1f%%",
                    nsim[1], nmark, 100 * ov)

  subtitle = NULL
  if(!is.null(threshold)) {
    thr = log10(threshold)
    ep1 = mean(x1 >= thr)
    ep2 = mean(x2 >= thr)
    subtitle = sprintf("Exceedance (LR \u2265 %g): %.0f%% under H1; %.0f%% under H2",
                       threshold, 100 * ep1, 100 * ep2)
  }
  # Plot distributions and their overlap
  x = H1 = H2 = overlap = NULL
  p = ggplot2::ggplot(d, ggplot2::aes(x)) +
    ggplot2::geom_area(ggplot2::aes(y = H1, fill = "H1 true"), alpha = 0.5) +
    ggplot2::geom_area(ggplot2::aes(y = H2, fill = "H2 true"), alpha = 0.5) +
    ggplot2::geom_area(ggplot2::aes(y = overlap), fill = "grey50", alpha = 0.75) +
    ggplot2::scale_fill_manual(values = c("H1 true" = col[1], "H2 true" = col[2]),
                               name = NULL) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = expression(log[10](LR)), y = "Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      plot.caption = ggplot2::element_text(size = 10, hjust = 0)
    )

  if(!is.null(threshold))
    p = p + ggplot2::geom_vline(xintercept = log10(threshold), linetype = 2, col = "red")

  p
}
