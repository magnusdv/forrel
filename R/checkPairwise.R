#' Check pedigree for relationship errors
#'
#' This function provides a convenient way to check for pedigree errors, given
#' the available marker data. The function calls [IBDestimate()] to estimate IBD
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
#' @param LRthreshold A positive number (default: 1000). Estimates whose LR
#'   (compared to the pedigree) exceeds
#'
#' @return A data frame containing both the estimated and pedigree-based IBD
#'   coefficients for each pair of typed individuals.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [IBDestimate()], [IBDtriangle()], [showInTriangle()]
#'
#' @examples
#'
#' ### Example with realistic data
#' x = addSon(nuclearPed(nch = 2), parent = 4)
#' x = setMarkers(x, locus = NorwegianFrequencies)
#' x = profileSim(x, N = 1, numCores = 1, seed = 1729)[[1]]
#'
#' checkPairwise(x)
#'
#' ### Create sample swap between 3 and 7
#' als = getAlleles(x)
#' als[c(2,5), ] = als[c(5,2), ]
#' y = setAlleles(x, alleles = als)
#'
#' checkPairwise(y)
#'
#' \donttest{
#' # Nice plot of pedigree + triangle
#' dev.new(height = 5, width = 8, noRStudioGD = TRUE)
#' layout(rbind(1:2), widths = 2:3)
#' plot(y, margins = c(4,2,4,2))
#' checkPairwise(y)
#' }
#'
#' @importFrom ribd kappaIBD
#' @export
checkPairwise = function(x, plot = TRUE, labels = FALSE, LRthreshold = 1000) {
  kEst = IBDestimate(x)
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
  kMerge = merge(kEst, kTrue, by = 1:2, suffixes = c(".est", ".ped")) # Ensure same pairing

  # LR comparing estimate against pedigree claim
  k0est = kMerge$kappa0.est
  k2est = kMerge$kappa2.est
  k0ped = kMerge$kappa0.ped
  k2ped = kMerge$kappa2.ped
  kMerge$LR = vapply(1:nrow(kMerge), function(i) {
    ids = kMerge[i, 1:2]
    loglik1 = .IBDlikelihood(x, ids = ids, kappa = c(k0est[i], k2est[i]), log = TRUE)
    loglik2 = .IBDlikelihood(x, ids = ids, kappa = c(k0ped[i], k2ped[i]), log = TRUE)
    exp(loglik1 - loglik2)
  }, FUN.VALUE = 1)


  if(plot) {
    # Factor defining colours
    kStr = paste(k0ped, k2ped, sep = "-")
    kFac = factor(kStr, levels = unique(kStr[order(k0ped, k2ped)]))
    levels(kFac)[levels(kFac) == "0-0"] = "Parent-offspring"
    levels(kFac)[levels(kFac) == "0.25-0.25"] = "Full siblings"
    levels(kFac)[levels(kFac) == "0.5-0"] = "Half/Uncle/Grand"
    levels(kFac)[levels(kFac) == "0.75-0"] = "First cousins"
    levels(kFac)[levels(kFac) == "1-0"] = "Unrelated"

    cols = pchs = as.integer(kFac) + 1
    nlev = nlevels(kFac)

    # Legend specifications
    legcol = legpch = seq_len(nlev) + 1
    legcex = rep(1, nlev)
    legtxt = levels(kFac)

    # Test for large LRs
    err = kMerge$LR > LRthreshold

    if(any(err)) {
      legcol = c(legcol, NA, 1)
      legpch = c(legpch, NA, 1)
      legcex = c(legcex, NA, 3)
      legtxt = c(legtxt, NA, sprintf("LR > %d", LRthreshold))
    }

    IBDtriangle()
    showInTriangle(kEst, col = cols, pch = pchs, labels = labels, new = FALSE)
    points(k0est[err], k2est[err], pch = 1, lwd = 2, cex = 3)

    legend("topright", title = " According to pedigree", title.adj = 0,
           bg = "whitesmoke", legend = legtxt, col = legcol, pch = legpch,
           pt.cex = legcex, lty = NA, lwd = 2)
  }

  kMerge
}
