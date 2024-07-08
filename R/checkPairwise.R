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
#' @param ids A vector of ID labels; the individuals to include in the check.
#'   Default: All typed members of `x`.
#' @param excludeInbred A logical, by default TRUE, indicating if inbred
#'   individuals should be excluded from the analysis.
#' @param plotType Either "base" (default), "ggplot2", "plotly" or "none".
#'   Abbreviations are allowed.
#' @param labels A logical (default: FALSE). If TRUE, labels are included in the
#'   IBD triangle plot.
#' @param LRthreshold A positive number (default: 1000). IBD estimates whose LR
#'   exceed this when compared to the coefficients implied by the pedigree, are
#'   flagged as potential errors in the output table and encircled in the plot.
#' @param pvalThreshold A positive number, or NULL (default). If given, this is
#'   used instead of `LRthreshold` to identify potential errors. Ignored if
#'   `nsim = 0`.
#' @param nsim A nonnegative number; the number of simulations used to estimate
#'   p-values. If 0 (default), this step is skipped.
#' @param plot Deprecated. To suppress the triangle plot, use `plotType =
#'   "none"`
#' @param verbose A logical.
#' @param ... Further parameters passed on to [ribd::ibdTriangle()].
#'
#' @return If `plotType` is "none" or "base": A data frame containing both the
#'   estimated and pedigree-based IBD coefficients for each pair of typed
#'   individuals. The last columns contains LRs and test results comparing the
#'   estimated coefficients to the pedigree-based ones.
#'
#'   If `plotType` is "ggplot2" or "plotly", only the plot objects are returned.
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
#' # Using p-values instead of LR (increase nsim!)
#' nsim = 10 # increase!
#' checkPairwise(y, nsim = nsim, pvalThreshold = 0.05)
#'
#' \donttest{
#' # Combined plot of pedigree and IBD estimates
#' dev.new(height = 5, width = 8, noRStudioGD = TRUE)
#' layout(rbind(1:2), widths = 2:3)
#' plot(y, margins = 2, title = "Swapped 1 - 3")
#' checkPairwise(y, labels = TRUE)
#' }
#'
#' @importFrom ribd inbreeding kappaIBD ibdTriangle showInTriangle
#' @importFrom graphics legend points
#' @importFrom grDevices palette
#' @export
checkPairwise = function(x, ids = typedMembers(x), excludeInbred = TRUE,
                         plotType = c("base", "ggplot2", "plotly", "none"),
                         labels = FALSE, LRthreshold = 1000, pvalThreshold = NULL,
                         nsim = 0, plot = TRUE, verbose = TRUE, ...) {

  includeIds = .myintersect(ids, typedMembers(x))
  if(length(includeIds) < 2) {
    if(verbose) message("No relationships to check: Less than 2 typed individuals included")
    return(invisible())
  }

  plotType = match.arg(plotType)
  if(isFALSE(plot)) {
    message("Argument `plot` is replaced with `plotType`. Use `plotType = 'none'` to suppress plotting")
    plotType = "none"
  }

  # If pedlist, allow duplicated names among non-included indivs
  if(is.pedList(x) && anyDuplicated(labs <- labels(x, unlist = TRUE))) {
    dups = unique.default(labs[duplicated.default(labs)])
    if(any(dups %in% includeIds))
      stop2("Duplicated ID label: ", intersect(dups, includeIds))

    # Other dups: Give unique labels
    for(i in seq_along(x)) {
      cmp = x[[i]]
      old = intersect(cmp$ID, dups)
      x[[i]] = relabel(cmp, old = old, new = sprintf("%s (Fam %d)", old, i))
    }
  }

  if(excludeInbred) {
    inbr = names(which(inbreeding(x, includeIds) > 0))
    if(length(inbr) && verbose)
      message("Excluding inbred individuals: ", toString(inbr))
    includeIds = setdiff(includeIds, inbr)
  }

  if(length(includeIds) < 2) {
    if(verbose) message("No relationships to check: Less than 2 typed individuals included")
    return(invisible())
  }

  # Estimated coefficients
  kEst = ibdEstimate(x, ids = includeIds, verbose = FALSE)

  # Pedigree coefficients
  kTrue = kappaIBD(x, ids = includeIds, inbredAction = as.integer(verbose), simplify = FALSE)

  # Merge (to ensure same pairing)
  kMerge = merge(kEst, kTrue, by = 1:2)
  k0 = kMerge$k0
  k2 = kMerge$k2
  kappa0 = kMerge$kappa0
  kappa1 = kMerge$kappa1 # needed in GLR
  kappa2 = kMerge$kappa2
  NR = nrow(kMerge)

  # LR comparing estimate against pedigree claim
  logLR = vapply(1:NR, function(i) {
    if(is.na(kappa0[i]) || is.na(kappa2[i]))
      return(NA_real_)
    ids = kMerge[i, 1:2]
    llFun = ibdLoglikFUN(x, ids = ids, input = "kappa02")
    loglik1 = llFun(c(k0[i], k2[i]))
    loglik2 = llFun(c(kappa0[i], kappa2[i]))
    loglik1 - loglik2
  }, FUN.VALUE = numeric(1))

  kMerge$LR = exp(logLR)

  # Relationship according to kappa (pedigree)
  kStr = paste(kappa0, kappa2, sep = "-")
  relLabs = c("0-0"       = "Parent-offspring",
              "0.25-0.25" = "Full siblings",
              "0.5-0"     = "Half/Uncle/Grand",
              "0.75-0"    = "First cousins",
              "1-0"       = "Unrelated",
              "NA-NA"     = "Other")
  pedrel = as.character(relLabs[kStr])
  pedrel[is.na(pedrel)] = "Other"
  kMerge$pedrel = pedrel = factor(pedrel, levels = .myintersect(relLabs, pedrel))

  # GLR
  kMerge$GLR = GLR = 1/kMerge$LR

  # Empirical p-value
  pval = rep(NA_real_, NR)
  if(nsim > 0) {
    ecdfList = list()
    db = getFreqDatabase(x)

    for(i in 1:NR) {
      kap = c(kappa0[i], kappa1[i], kappa2[i])
      if(anyNA(kap))
        next
      kapStr = paste(kap, collapse = "-")

      # If ecdf already computed, get it - otherwise simulate
      if(kapStr %in% names(ecdfList))
        cdf = ecdfList[[kapStr]]
      else {
        if(verbose) cat("Simulating null distribution for GLR at kappa =", kapStr, "\n")
        cdf = ecdfList[[kapStr]] = ecdfGLR(kappa = kap, nsim = nsim, freqList = db)
      }
      pval[i] = cdf(log(GLR[i]))
    }
  }
  kMerge$pval = pval


  # Test for errors
  if(isNumber(pvalThreshold, minimum = 0, maximum = 1) && nsim > 0) {
    kMerge$err = err = !is.na(kMerge$pval) & kMerge$pval < pvalThreshold
    errtxt = sprintf("pval < %g", pvalThreshold)
  }
  else {
    kMerge$err = err = !is.na(kMerge$LR) & kMerge$LR > LRthreshold
    errtxt = sprintf("LR > %g", LRthreshold)
  }

  if(plotType == "none")
    return(kMerge)

  # Plot --------------------------------------------------------------------

  if(any(err)) {
    errDat = kMerge[err, , drop = FALSE]
    errDat$labs = labs = paste(errDat$id1, "-", errDat$id2)
  }
  else
    errDat = NULL

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
      legtxt = c(legtxt, NA, errtxt)
    }

    ribd::showInTriangle(kMerge[1:6], plotType = "base", col = cols, pch = pchs,
                         labels = labels, ...)
    points(errDat$k0, errDat$k2, pch = 1, lwd = 1.8, cex = 3)

    legend("topright", title = " According to pedigree", title.adj = 0,
           bg = "whitesmoke", legend = legtxt, col = legcol, pch = legpch,
           pt.cex = legcex, lty = NA, lwd = 2)

    return(kMerge)
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
    return(p)
  }

  if(plotType == "plotly") {
    if(!requireNamespace("plotly", quietly = TRUE))
        stop2("Package `plotly` must be installed for this option to work")

    dat = kMerge
    dat$idx = seq_len(nrow(kMerge))
    dat$labs = paste0("ID1: ", kMerge$id1, "<br>", "ID2: ", kMerge$id2)

    symbs = c("Parent-offspring" = "triangle-up-open",
              "Full siblings" = "cross-thin-open",
              "Half/Uncle/Grand" = "x-thin-open",
              "First cousins" = "diamond-open",
              "Unrelated" = "triangle-down-open",
              "Other" = "asterisk-open")

    cols = c(palette()[c(2,3,4,7,5)],"#C8A2C8")
    names(cols) = names(symbs)

    p = ribd::ibdTriangle(plotType = "plotly", ...)
    for(r in rev(levels(pedrel))) {
      datr = dat[dat$pedrel == r, , drop = FALSE]
      p = p |> plotly::add_markers(data = datr, x = ~k0, y = ~k2, customdata = ~idx,
                                   symbol = I(symbs[r]), color = I(cols[r]),
                  marker = list(size = 12,line = list(width = if(r == "Other") 1 else 2)),
                  text= ~labs, hoverinfo = "text", name = r)
    }
    if(any(dat$err)) {
      p = p |> plotly::add_markers(data = errDat, x = ~k0, y = ~k2,
                                   name = errtxt,
                                   symbol = I("circle-open"), color = I("black"),
                                   marker = list(size = 20,line = list(width = 1)))
    }

    p = p |> plotly::layout(
      legend = list(title = list(text = 'According to pedigree'),
                    x = 1, y = 1, xanchor = 'right', yanchor = 'top',
                    bgcolor = "whitesmoke",
                    traceorder = "reversed")
      )
    return(p)
  }

  invisible(NULL)
}


# Empiric cdf of GLR
#' @importFrom stats ecdf
ecdfGLR = function(kappa = NULL, nsim = 1000, freqList = NULL, log = TRUE) {

  # Insert default allele labels (1,2,3,...) where needed
  freqList = lapply(freqList, function(fr)
    if(is.null(names(fr))) .setnames(fr, seq_along(fr)) else fr)

  # Simulate from kappa; return only allele indices (most efficient)
  sims = profileSimParametric(kappa = kappa, N = nsim, freqList = freqList,
                              returnValue = "internal")

  # GLR for each sim
  glr = lapply(sims, function(s) {
    numer = .ibdLoglikFromAlleles(s, freqList, kappa = kappa)
    denom = .ibdEstimFromAlleles(s, freqList, param = "kappa", start = kappa, returnValue = "loglik")
    numer - denom
  })

  glrvec = unlist(glr, use.names = FALSE)
  if(!log)
    glrvec = exp(glrvec)

  ecdf(glrvec)
}
