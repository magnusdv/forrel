#' Check pedigree data for relationship errors
#'
#' This function provides a convenient way to check for pedigree errors, given
#' the available marker data. The function calls [ibdEstimate()] to estimate IBD
#' coefficients for all pairs of typed pedigree members, and computes the
#' likelihood ratio (LR) comparing each estimate to the coefficients implied by
#' the pedigree. By default, the estimates are shown in a colour-coded plot
#' where unlikely relationships are easy to spot.
#'
#' By default, inbred individuals are excluded from the analysis, since pairwise
#' relationships involving inbred individuals have undefined kappa coefficients
#' (and therefore no position in the triangle). In some cases it may still be
#' informative to include their estimates; set `excludeInbred = FALSE` to
#' achieve this.
#'
#' @param x A `ped` object or a list of such.
#' @param excludeInbred A logical, by default TRUE, indicating if inbred
#'   individuals should be excluded from the analysis.
#' @param plot A logical (default: TRUE). If TRUE, a plot is produced, showing
#'   the IBD estimates in the IBD triangle.
#' @param plotType Either "base" (default), "ggplot2", "plotly" or "none".
#'   Abbreviations are allowed.
#' @param labels A logical (default: FALSE). If TRUE, labels are included in the
#'   IBD triangle plot.
#' @param LRthreshold A positive number (default: 1000). IBD estimates whose LR
#'   exceed this, when compared to the coefficients implied by the pedigree, are
#'   encircled in the plot.
#' @param ... Further parameters passed on to [ribd::ibdTriangle()].
#'
#' @return A data frame containing both the estimated and pedigree-based IBD
#'   coefficients for each pair of typed individuals. The last column contains
#'   the likelihood ratio comparing the estimated coefficients to the
#'   pedigree-based ones.
#'
#' @seealso [ibdEstimate()]
#'
#' @examples
#'
#' ### Example with realistic data
#'
#' x = avuncularPed() |>
#'   profileSim(markers = NorwegianFrequencies, seed = 1729)
#'
#' checkPairwise(x)
#'
#' ### Create an error: sample swap 1 <-> 3
#' als = getAlleles(x)
#' als[c(1,3), ] = als[c(3,1), ]
#' y = setAlleles(x, alleles = als)
#'
#' checkPairwise(y)
#'
#' \donttest{
#' # Combined plot of pedigree and IBD estimates
#' dev.new(height = 5, width = 8, noRStudioGD = TRUE)
#' layout(rbind(1:2), widths = 2:3)
#' plot(y, margins = c(4,2,4,2), title = "Swapped 1 - 3")
#' checkPairwise(y, labels = TRUE)
#' }
#'
#' @importFrom ribd inbreeding kappaIBD ibdTriangle showInTriangle
#' @importFrom graphics legend points
#' @export
checkPairwise = function(x, excludeInbred = TRUE, plot = TRUE,
                         plotType = c("base", "ggplot2", "plotly", "none"),
                         labels = FALSE, LRthreshold = 1000, ...) {
  includeIds = typedMembers(x)
  plotType = match.arg(plotType)
  if(isFALSE(plot)) {
    message("Argument `plot` is replaced with `plotType`. For `plot = F` use `plotType = 'none'`")
    plotType = "none"
  }

  if(excludeInbred) {
    inbr = names(which(inbreeding(x) > 0))
    if(length(inbr))
      message("Excluding inbred individuals: ", inbr)
    includeIds = setdiff(includeIds, inbr)
  }

  if(length(includeIds) < 2) {
    message("No relationships to check")
    return(invisible())
  }

  # Estimated coefficients
  kEst = ibdEstimate(x, ids = includeIds, verbose = FALSE)

  # Pedigree coefficients
  kTrue = kappaIBD(x, ids = includeIds, simplify = FALSE)

  # Merge (to ensure same pairing)
  kMerge = merge(kEst, kTrue, by = 1:2)
  k0 = kMerge$k0
  k2 = kMerge$k2
  kappa0 = kMerge$kappa0
  kappa2 = kMerge$kappa2

  # LR comparing estimate against pedigree claim
  kMerge$LR = vapply(1:nrow(kMerge), function(i) {
    ids = kMerge[i, 1:2]
    loglik1 = .IBDlikelihood(x, ids = ids, kappa = c(k0[i], k2[i]), log = TRUE)
    loglik2 = .IBDlikelihood(x, ids = ids, kappa = c(kappa0[i], kappa2[i]), log = TRUE)
    exp(loglik1 - loglik2)
  }, FUN.VALUE = 1)

  # Relationship according to kappa (pedigree)
  kStr = paste(kappa0, kappa2, sep = "-")
  pedrel = factor(kStr, levels = unique(kStr[order(kappa0, kappa2)]))
  levels(pedrel)[levels(pedrel) == "0-0"] = "Parent-offspring"
  levels(pedrel)[levels(pedrel) == "0.25-0.25"] = "Full siblings"
  levels(pedrel)[levels(pedrel) == "0.5-0"] = "Half/Uncle/Grand"
  levels(pedrel)[levels(pedrel) == "0.75-0"] = "First cousins"
  levels(pedrel)[levels(pedrel) == "1-0"] = "Unrelated"
  levels(pedrel)[levels(pedrel) == "NA-NA"] = "NA (inbred)"
  kMerge$pedrel = pedrel

  # Test for large LRs
  kMerge$err = err = kMerge$LR > LRthreshold

  if(any(err)) {
    errDat = kMerge[err, , drop = FALSE]
    errDat$labs = labs = paste(errDat$id1, "-", errDat$id2)
  }
  else
    errDat = NULL

  # Plot --------------------------------------------------------------------

  if(plotType == "base") {
    cols = pchs = as.integer(pedrel) + 1
    nlev = nlevels(pedrel)

    # Legend specifications
    legcol = legpch = seq_len(nlev) + 1
    legcex = rep(1, nlev)
    legtxt = levels(pedrel)

    if(any(err, na.rm = T)) {
      legcol = c(legcol, NA, 1)
      legpch = c(legpch, NA, 1)
      legcex = c(legcex, NA, 3)
      legtxt = c(legtxt, NA, sprintf("LR > %d", LRthreshold))
    }

    ribd::showInTriangle(kMerge[1:6], plotType = "base", col = cols, pch = pchs,
                         labels = labels, ...)
    points(errDat$k0, errDat$k2, pch = 1, lwd = 1.8, cex = 3)

    legend("topright", title = " According to pedigree", title.adj = 0,
           bg = "whitesmoke", legend = legtxt, col = legcol, pch = legpch,
           pt.cex = legcex, lty = NA, lwd = 2)
  }

  if(plotType == "ggplot2") {
    if(!requireNamespace("ggplot2", quietly = TRUE))
        stop2("Package `ggplot2` must be installed for this option to work")

    p = ribd::ibdTriangle(plotType = "ggplot2", ...) +
      ggplot2::geom_point(data = kMerge, ggplot2::aes(k0, k2, color = pedrel, shape = pedrel),
                          size = 2, stroke = 1.5) +
      ggplot2::labs(color = "According to pedigree", shape = "According to pedigree") +
      ggplot2::scale_shape_manual(values = 1+seq_len(nlevels(pedrel))) +
      ggplot2::scale_color_manual(values = 1+seq_len(nlevels(pedrel))) +
      ggplot2::theme(legend.position = c(1,1),
                     legend.justification = c(1, 1.1),
                     legend.margin = ggplot2::margin(3,6,3,6),
                     legend.background = ggplot2::element_rect(fill = "whitesmoke"))

    if(!is.null(errDat)) {
      if(!requireNamespace("ggrepel", quietly = TRUE))
        stop2("Package `ggrepel` must be installed for this option to work")

      p = p +
        ggplot2::geom_point(data = errDat, ggplot2::aes(k0, k2, size = I(8)),
                          shape = 1, size = 8, stroke = 1, col = 1) +
        ggrepel::geom_text_repel(ggplot2::aes(k0, k2, label = labs, color = pedrel),
                                 data = errDat, size = 4, max.overlaps = Inf,
                                 box.padding = 1, show.legend = FALSE)
    }
    print(p)
  }

  if(plotType == "plotly")
    stop2("Not implemented yet")

  kMerge
}
