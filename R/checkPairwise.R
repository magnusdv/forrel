#' Check pedigree data for relationship errors
#'
#' The `checkPairwise()` function provides a convenient way to check for
#' pedigree errors, given the available marker data. The function calls
#' [ibdEstimate()] to estimate IBD coefficients for all pairs of typed pedigree
#' members, and uses the estimates to test for potential errors. By default, the
#' results are shown in a colour-coded plot (based on [ribd::ibdTriangle()])
#' where unlikely relationships are easy to spot.
#'
#' To identify potential pedigree errors, the function calculates the
#' *generalised likelihood ratio* (GLR) of each pairwise relationship.
#' This compares the likelihood of the estimated coefficients with that of the
#' coefficients implied by the pedigree. By default, relationships whose GLR
#' exceed 1000 are flagged as errors and shown with a circle in the plot.
#' Alternatively, if arguments `nsim` and `pvalThreshold` are supplied, the
#' p-value of each score is estimated by simulation, and used as threshold for
#' calling errors.
#'
#' By default, inbred individuals are excluded from the analysis, since pairwise
#' relationships involving inbred individuals have undefined kappa coefficients
#' (and therefore no position in the triangle). In some cases it may still be
#' informative to include their estimates; set `excludeInbred = FALSE` to
#' enforce this.
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
#' @param GLRthreshold A positive number, by default 1000. Threshold for the
#'   generalised likelihood ratio (see Details). Scores exceeding this are
#'   flagged as potential errors in the output table and encircled in the plot.
#' @param pvalThreshold A positive number, or NULL (default). If given, this is
#'   used instead of `GLRthreshold` to identify potential errors. Ignored if
#'   `nsim = 0`.
#' @param nsim A nonnegative number; the number of simulations used to estimate
#'   p-values. If 0 (default), this step is skipped.
#' @param seed An integer seed for the random number generator (optional, and
#'   only relevant if `nsim > 0`).
#' @param cpRes A data frame: the output from `checkPairwise()`.
#' @param errtxt A character string to use for the error legend.
#' @param plot Deprecated. To suppress the triangle plot, use `plotType =
#'   "none"`
#' @param verbose A logical.
#' @param ... Further parameters passed on to [ribd::ibdTriangle()].
#'
#' @return If `plotType` is "none" or "base": A data frame containing both the
#'   estimated and pedigree-based IBD coefficients for each pair of typed
#'   individuals. The last columns (`GLR`, `pval` and `err`) contain test
#'   results using the GLR scores to identify potential pedigree errors.
#'
#'   If `plotType` is "ggplot2" or "plotly", the plot objects are returned.
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
#' # Using p-values instead of GLR
#' nsim = 10 # increase!
#' checkPairwise(y, nsim = nsim, pvalThreshold = 0.05)
#'
#' # Plot can be done separately
#' res = checkPairwise(y, nsim = nsim, pvalThreshold = 0.05, plotType = "none")
#' plotCP(res, plotType = "base", errtxt = "Not good!")
#'
#' \donttest{
#' # Combined plot of pedigree and check results
#' dev.new(height = 5, width = 8, noRStudioGD = TRUE)
#' layout(rbind(1:2), widths = 2:3)
#' plot(y, margins = 2, title = "Swapped 1 - 3")
#' plotCP(res, labels = TRUE)
#' }
#'
#' @importFrom ribd inbreeding kappaIBD
#' @importFrom verbalisr verbalise
#' @export
checkPairwise = function(x, ids = typedMembers(x), excludeInbred = TRUE,
                         plotType = c("base", "ggplot2", "plotly", "none"),
                         GLRthreshold = 1000, pvalThreshold = NULL, nsim = 0,
                         seed = NULL, plot = TRUE, verbose = TRUE, ...) {

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
  cpRes = merge(kEst, kTrue, by = 1:2)
  k0 = cpRes$k0
  k2 = cpRes$k2
  kappa0 = cpRes$kappa0
  kappa1 = cpRes$kappa1 # needed in GLR
  kappa2 = cpRes$kappa2
  NR = nrow(cpRes)

  pedrel = character(NR)
  for(i in 1:NR) {
    pedrel[i] = verbalise(x, ids = cpRes[i, 1:2]) |>
      format(cap = TRUE, simplify = TRUE, abbreviate = TRUE, collapse = " & ")
  }
  cpRes$pedrel = pedrel

  # Relationship group (for legend)
  kStr = paste(kappa0, kappa2, sep = "-")
  relLabs = c("0-0"       = "Parent-offspring",
              "0.25-0.25" = "Full siblings",
              "0.5-0"     = "Half/Uncle/Grand",
              "0.75-0"    = "First cousins/etc",
              "NA-NA"     = "Other",
              "1-0"       = "Unrelated")
  relgroup = as.character(relLabs[kStr])
  relgroup[is.na(relgroup)] = "Other"
  relgroup = factor(relgroup, levels = .myintersect(relLabs, relgroup))
  cpRes$relgroup = relgroup

  # GLR: compare estimate against pedigree claim
  logGLR = vapply(1:NR, function(i) {
    khat = c(k0[i], k2[i])
    kped = c(kappa0[i], kappa2[i])
    if(anyNA(kped))
      return(NA_real_)
    llFun = ibdLoglikFUN(x, ids = cpRes[i, 1:2], input = "kappa02")
    llFun(khat) - llFun(kped)
  }, FUN.VALUE = numeric(1))

  cpRes$GLR = exp(logGLR)

  # Empirical p-values
  pval = rep(NA_real_, NR)

  if(nsim > 0) {
    ecdfList = list()
    db = getFreqDatabase(x)

    for(i in which(!is.na(logGLR))) { # only rows with GLR result
      ks = kStr[i]

      # If ecdf not already computed: simulate
      if(!ks %in% names(ecdfList)) {
        if(verbose)
          cat("Simulating null distribution for GLR at kappa =", ks, "\n")
        kap = c(kappa0[i], kappa1[i], kappa2[i])
        ecdfList[[ks]] = ecdfGLR(kap, nsim = nsim, freqList = db, log = TRUE, seed = seed)
      }
      cdf = ecdfList[[ks]]
      pval[i] = 1 - cdf(logGLR[i])
    }
  }
  cpRes$pval = pval

  # Flag potential pedigree errors
  if(isNumber(pvalThreshold, minimum = 0, maximum = 1) && nsim > 0) {
    cpRes$err = err = !is.na(pval) & pval < pvalThreshold
    errtxt = sprintf("pval < %g", pvalThreshold)
  }
  else {
    cpRes$err = err = !is.na(cpRes$GLR) & cpRes$GLR >= GLRthreshold
    errtxt = sprintf("GLR > %g", GLRthreshold)
  }

  if(plotType != "none")
    plotCP(cpRes, plotType = plotType, errtxt = errtxt, seed = seed, ...)
  else
    return(cpRes)
}

# Plot methods ----------------------------------------------------------------

#' @rdname checkPairwise
#' @importFrom ribd ibdTriangle showInTriangle
#' @importFrom graphics legend points
#' @importFrom grDevices palette
#' @export
plotCP = function(cpRes, plotType = c("base", "ggplot2", "plotly"),
                  labels = FALSE, errtxt = "Potential error",
                  seed = NULL, ...) {
  err = cpRes$err
  relgroup = cpRes$relgroup
  k0 = cpRes$k0
  k2 = cpRes$k2
  plotType = match.arg(plotType)

  errDat = NULL
  if(any(err)) {
    errDat = cpRes[err, , drop = FALSE]
    errDat$labs = labs = paste(errDat$id1, "-", errDat$id2)
  }

  # Copied from checkPairwise()
  relLabs = c("0-0"       = "Parent-offspring",
              "0.25-0.25" = "Full siblings",
              "0.5-0"     = "Half/Uncle/Grand",
              "0.75-0"    = "First cousins/etc",
              "NA-NA"     = "Other",
              "1-0"       = "Unrelated")
  ALLCOLS = .setnames(2:7, relLabs)
  ALLSHAPES = .setnames(c(2,3,4,5,8,6), relLabs)

  if(plotType == "base") {
    cols = ALLCOLS[as.character(relgroup)]
    pchs = ALLSHAPES[as.character(relgroup)]

    # Legend specifications
    legtxt = levels(relgroup)
    legcol = ALLCOLS[levels(relgroup)]
    legpch = ALLSHAPES[levels(relgroup)]
    legcex = rep(1, length(legcol))

    if(any(err, na.rm = T)) {
      legtxt = c(legtxt, NA, errtxt)
      legcol = c(legcol, NA, 1)
      legpch = c(legpch, NA, 1)
      legcex = c(legcex, NA, 3)
    }

    ribd::showInTriangle(cpRes[1:6], plotType = "base", col = cols, pch = pchs,
                         labels = labels, ...)
    points(errDat$k0, errDat$k2, pch = 1, lwd = 1.8, cex = 3)

    legend("topright", title = " According to pedigree", title.adj = 0,
           bg = "whitesmoke", legend = legtxt, col = legcol, pch = legpch,
           pt.cex = legcex, lty = NA, lwd = 2)

    return(cpRes)
  }

  if(plotType == "ggplot2") {
    if(!requireNamespace("ggplot2", quietly = TRUE))
        stop2("Package `ggplot2` must be installed for this option to work")

    p = ribd::ibdTriangle(plotType = "ggplot2", ...) +
      ggplot2::geom_point(data = cpRes, ggplot2::aes(k0, k2, color = relgroup, shape = relgroup),
                          size = 2, stroke = 1.5) +
      ggplot2::labs(color = "According to pedigree", shape = "According to pedigree") +
      ggplot2::scale_shape_manual(values = ALLSHAPES) +
      ggplot2::scale_color_manual(values = ALLCOLS) +
      ggplot2::theme(legend.position = c(1,1),
                     legend.justification = c(1, 1.1),
                     legend.margin = ggplot2::margin(3,6,3,6),
                     legend.background = ggplot2::element_rect(fill = "whitesmoke"),
                     legend.text = ggplot2::element_text(size = 10),
                     legend.title = ggplot2::element_text(size = 12)
                     )

    if(isTRUE(labels) || !is.null(errDat)) {
      if(!requireNamespace("ggrepel", quietly = TRUE))
        stop2("Package `ggrepel` must be installed for this option to work")
    }

    # TODO: ad hoc!
    if(isTRUE(labels)) {
      cpRes$labs = paste(cpRes$id1, "-", cpRes$id2)

      p = p + ggrepel::geom_text_repel(data = cpRes,
        ggplot2::aes(k0, k2, label = labs, color = relgroup), min.segment.length = 1.5, seed = seed, ...,
        size = 4, max.overlaps = Inf, box.padding = 1, show.legend = FALSE)
    }

    if(!is.null(errDat)) {
      p = p +
        ggplot2::geom_point(data = errDat, ggplot2::aes(k0, k2, size = "big"),
                          shape = 1, stroke = 1, col = 1) +
        ggplot2::scale_size_manual(values = c(big = 8), labels = errtxt, name = NULL) +
        ggplot2::guides(shape = ggplot2::guide_legend(order = 1),
                        colour = ggplot2::guide_legend(order = 1))
        if(!isTRUE(labels) && !is.null(labels)) # if TRUE, all labs anyway
          p = p + ggrepel::geom_text_repel(ggplot2::aes(k0, k2, label = labs, color = relgroup),
                                 data = errDat, size = 4, max.overlaps = Inf,
                                 box.padding = 1, show.legend = FALSE)
    }

    return(p)
  }

  if(plotType == "plotly") {
    if(!requireNamespace("plotly", quietly = TRUE))
        stop2("Package `plotly` must be installed for this option to work")

    dat = cpRes
    dat$idx = seq_len(nrow(cpRes))
    dat$labs = paste0("ID1: ", cpRes$id1, "<br>", "ID2: ", cpRes$id2)

    symbs = c("Parent-offspring" = "triangle-up-open",
              "Full siblings" = "cross-thin-open",
              "Half/Uncle/Grand" = "x-thin-open",
              "First cousins/etc" = "diamond-open",
              "Other" = "asterisk-open",
              "Unrelated" = "triangle-down-open")

    cols = c(palette()[2:7]) #c(2,3,4,7,5)],"#C8A2C8")
    names(cols) = names(symbs)

    p = ribd::ibdTriangle(plotType = "plotly", ...)

    for(r in rev(levels(relgroup))) {
      datr = dat[dat$relgroup == r, , drop = FALSE]
      p = p |> plotly::add_markers(data = datr, x = ~k0, y = ~k2, customdata = ~idx,
                                   symbol = I(symbs[r]), color = I(cols[r]),
                  marker = list(size = 12,line = list(width = if(r == "Other") 1 else 2)),
                  text= ~labs, hoverinfo = "text", name = r, legendrank = 2)
    }
    if(any(dat$err)) {

      # Invisible spacer trace
      p = p |> plotly::add_markers(x = 0, y = 0, name = " ", hoverinfo = "none",
                             marker = list(size = 0, color = 'rgba(0,0,0,0)'),
                             showlegend = TRUE, legendrank = 1.5)

      # Error circles
      p = p |> plotly::add_markers(data = errDat, x = ~k0, y = ~k2,
                                   name = errtxt,
                                   symbol = I("circle-open"), color = I("black"),
                                   marker = list(size = 22, line = list(width = 1)),
                                   hoverinfo = "none", legendrank = 1)
    }

    p = p |> plotly::layout(
      legend = list(title = list(text = 'According to pedigree'),
                    x = 1, y = 1, xanchor = 'right', yanchor = 'top',
                    bgcolor = "whitesmoke",
                    traceorder = "reversed")
      )
    return(p)
  }
}



# Empiric cdf of GLR
#' @importFrom stats ecdf
ecdfGLR = function(kappa = NULL, nsim = 1000, freqList = NULL, log = TRUE, seed = NULL) {

  # Insert default allele labels (1,2,3,...) where needed
  freqList = lapply(freqList, function(fr)
    if(is.null(names(fr))) .setnames(fr, seq_along(fr)) else fr)

  # Simulate from kappa; return only allele indices (most efficient)
  sims = profileSimParametric(kappa = kappa, N = nsim, freqList = freqList,
                              returnValue = "internal", seed = seed)

  # GLR for each sim
  glr = lapply(sims, function(s) {
    numer = .ibdEstimFromAlleles(s, freqList, param = "kappa", start = kappa, returnValue = "loglik")
    denom = .ibdLoglikFromAlleles(s, freqList, kappa = kappa)
    numer - denom
  })

  glrvec = unlist(glr, use.names = FALSE)
  if(!log)
    glrvec = exp(glrvec)

  ecdf(glrvec)
}
