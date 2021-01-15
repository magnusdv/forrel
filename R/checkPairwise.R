#' Check pedigree data for relationship errors
#'
#' This function provides a convenient way to check for pedigree errors, given
#' the available marker data. The function calls [ibdEstimate()] to estimate IBD
#' coefficients for all pairs of typed pedigree members, and computes the
#' likelihood ratio (LR) comparing each estimate to the coefficients implied by
#' the pedigree. By default, the estimates are shown in a colour-coded plot
#' where unlikely relationships are easy to spot.
#'
#' @param x A `ped` object or a list of such.
#' @param plot A logical (default: TRUE). If TRUE, a plot is produced, showing
#'   the IBD estimates in the IBD triangle.
#' @param labels A logical (default: FALSE). If TRUE, labels are included in the
#'   IBD triangle plot.
#' @param LRthreshold A positive number (default: 1000). IBD estimates whose LR
#'   exceed this, when compared to the coefficients implied by the pedigree, are
#'   encircled in the plot.
#' @param ... further parameters passed on to [ribd::ibdTriangle()].
#'
#' @return A data frame containing both the estimated and pedigree-based IBD
#'   coefficients for each pair of typed individuals. The last column contains
#'   the likelihood ratio comparing the estimated coefficients to the
#'   pedigree-based ones.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [ibdEstimate()]
#'
#' @examples
#'
#' ### Example with realistic data
#'
#' x = addSon(nuclearPed(nch = 2), parent = 4)
#' x = setMarkers(x, locus = NorwegianFrequencies)
#' x = profileSim(x, N = 1, seed = 1729)[[1]]
#'
#' checkPairwise(x)
#'
#' ### Create sample swap between 1 and 3
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
#' plot(x, margins = c(4,2,4,2))
#' checkPairwise(x, labels = TRUE)
#' }
#'
#' @importFrom ribd kappaIBD ibdTriangle showInTriangle
#' @export
checkPairwise = function(x, plot = TRUE, labels = FALSE, LRthreshold = 1000, ...) {
  # Estimated coefficients
  kEst = ibdEstimate(x)

  # Pedigree coefficients
  if(is.ped(x)) {
    kTrue = kappaIBD(x)
  } else { #TODO: fix in ribd
    kTrue = do.call(rbind, lapply(x, kappaIBD))
    nPed = length(x)
    for(i in seq_len(nPed - 1)) for(j in seq(i+1, nPed))
      for(a in labels(x[[i]]))
        kTrue = rbind(kTrue, data.frame(id1 = a, id2 = labels(x[[j]]),
                                        kappa0 = 1, kappa1 = 0, kappa2 = 0,
                                        stringsAsFactors = FALSE))
  }

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


  if(plot) {
    # Factor defining colours
    kStr = paste(kappa0, kappa2, sep = "-")
    kFac = factor(kStr, levels = unique(kStr[order(kappa0, kappa2)]))
    levels(kFac)[levels(kFac) == "0-0"] = "Parent-offspring"
    levels(kFac)[levels(kFac) == "0.25-0.25"] = "Full siblings"
    levels(kFac)[levels(kFac) == "0.5-0"] = "Half/Uncle/Grand"
    levels(kFac)[levels(kFac) == "0.75-0"] = "First cousins"
    levels(kFac)[levels(kFac) == "1-0"] = "Unrelated"
    levels(kFac)[levels(kFac) == "NA-NA"] = "NA (inbred)"

    cols = pchs = as.integer(kFac) + 1
    nlev = nlevels(kFac)

    # Legend specifications
    legcol = legpch = seq_len(nlev) + 1
    legcex = rep(1, nlev)
    legtxt = levels(kFac)

    # Test for large LRs
    err = kMerge$LR > LRthreshold

    if(any(err, na.rm = T)) {
      legcol = c(legcol, NA, 1)
      legpch = c(legpch, NA, 1)
      legcex = c(legcex, NA, 3)
      legtxt = c(legtxt, NA, sprintf("LR > %d", LRthreshold))
    }

    ribd::ibdTriangle(...)
    ribd::showInTriangle(kMerge[1:6], col = cols, pch = pchs, labels = labels, new = FALSE)
    points(k0[err], k2[err], pch = 1, lwd = 2, cex = 3)

    legend("topright", title = " According to pedigree", title.adj = 0,
           bg = "whitesmoke", legend = legtxt, col = legcol, pch = legpch,
           pt.cex = legcex, lty = NA, lwd = 2)
  }

  kMerge
}
