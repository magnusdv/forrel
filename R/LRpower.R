#' Power simulation for kinship LR
#'
#' This function uses simulations to estimate the likelihood ratio (LR)
#' distribution in a given kinship testing scenario. In the most general
#' setting, three pedigrees are involved: the two pedigrees being compared, and
#' the true relationship (which may differ from the other two). A subset of
#' individuals are available for genotyping. Some individuals may already be
#' genotyped; all simulations are then conditional on these.
#'
#' @inheritParams exclusionPower
#' @param numeratorPed,denominatorPed `ped` objects (or lists of such),
#'   describing the two relationships under comparison.
#' @param truePed A `ped` object (or a list of such), describing the true
#'   relationship. By default equal to `numeratorPed`.
#' @param ids Individuals available for genotyping.
#' @param source Either "true" (default), "numerator" or "denominator",
#'   indicating which pedigree is used as source for marker data.
#' @param nsim A positive integer: the number of simulations.
#' @param threshold A numeric vector with one or more positive numbers used as
#'   LR thresholds.
#' @param disableMutations Not implemented yet.
#' @param alleles,afreq,Xchrom If these are given, they are used (together with
#'   `knownGenotypes`) to create a marker object on the fly.
#' @param seed An integer seed for the random number generator (optional).
#'
#' @return A `LRpowerResult` object, which is essentially a list with the
#'   following entries:
#'
#'   * `LRperSim`: A numeric vector of length `nsim` containing the total LR for
#'   each simulation.
#'
#'   * `meanLRperMarker`: The mean LR per marker, over all simulations.
#'
#'   * `meanLR`: The mean total LR over all simulations.
#'
#'   * `meanLogLR`: The mean total `log10(LR)` over all simulations.
#'
#'   * `IP`: A named numeric of the same length as `threshold`. For each element
#'   of `threshold`, the fraction of simulations resulting in a LR exceeding the
#'   given number.
#'
#'   * `time`: The total computation time.
#'
#'   * `params`: A list containing the input parameters `missing`, `markers`,
#'   `nsim`, `threshold` and `disableMutations`
#'
#' @examples
#'
#' # Paternity LR of siblings
#' ids = c("A", "B")
#' truth = nuclearPed(children = ids)
#' claim = nuclearPed(fa = "A", mo = "NN", children = "B")
#' unrel = singletons(ids)
#'
#' # Simulation parameters
#' nsim = 10   # increase!
#' thresh = 1
#'
#' # Simulation 1:
#' als = 1:5
#' afr = runif(5)
#' afr = afr/sum(afr)
#'
#' pow1 = LRpower(claim, unrel, truth, ids = ids, nsim = nsim,
#'                threshold = thresh, alleles = als, afreq = afr,
#'                seed = 123)
#' pow1
#'
#' # Simulation 2: Same, but using an attached marker
#' truth = addMarker(truth, alleles = als, afreq = afr)
#'
#' pow2 = LRpower(claim, unrel, truth, ids = ids, nsim = nsim,
#'                threshold = thresh, markers = 1, seed = 123)
#'
#' stopifnot(identical(pow1$LRperSim, pow2$LRperSim))
#'
#' \donttest{
#' # True pedigree has inbred founders
#' truth2 = setFounderInbreeding(truth, value = 0.5)
#'
#' pow3 = LRpower(claim, unrel, truth2, ids = ids, nsim = nsim,
#'                threshold = thresh, markers = 1, seed = 123) # plot = TRUE
#' pow3
#' }
#'
#' @export
LRpower = function(numeratorPed, denominatorPed, truePed = numeratorPed, ids, markers = NULL,
                   source = "true", nsim = 1, threshold = NULL, disableMutations = NA,
                   alleles = NULL, afreq = NULL, Xchrom = FALSE, knownGenotypes = NULL,
                   plot = FALSE, plotMarkers = NULL, seed = NULL, verbose = TRUE) {
  st = Sys.time()
  if(is.list(ids)) {
    ids = lapply(ids, as.character)
    allids = unique.default(unlist(ids))
  }
  else {
    allids = as.character(ids)
  }

  ### Input option 1: Single marker with attributes given directly
  if(!is.null(alleles) || !is.null(afreq)) {

    # Convert knownGenotypes to allele matrix
    am = NULL
    if(!is.null(knownGenotypes)) {

      # Check format
      if(!is.list(knownGenotypes) || !all(lengths(knownGenotypes) == 3))
        stop2("`knownGenotypes` must be a list of vectors of the form `c(ID, allele1, allele2)`")

      am = do.call(rbind, knownGenotypes)
      rownames(am) = as.character(am[, 1])
      am = am[, -1, drop = FALSE]
    }

    # If `alleles` is single integer, convert to sequence
    if(isNumber(alleles))
      alleles = seq_len(alleles)

    # Create and attach locus to both pedigrees
    locus = list(alleles = alleles, afreq = afreq, chrom = if (Xchrom) 23 else NA)
    truePed = setMarkers(truePed, alleleMatrix = am, locusAttributes = locus)

    markers = 1
    typed = typedMembers(truePed)
    disableMutations = FALSE # don't do anything
  }
  else {
    source = match.arg(source, c("true", "numerator", "denominator"))
    sourcePed = switch(source, true = truePed, numerator = numeratorPed, denominator = denominatorPed,
                       stop2("`source` must be either 'true', 'numerator' or 'denominator': ", source))

    nmTot = nMarkers(sourcePed)
    if(nmTot == 0)
      stop2("No markers attached to the source pedigree ('", source, "')")

    # If `markers` not given, use all attached. Use names if present.
    if(is.null(markers)) {
      if(verbose)
        message("Using all ", nmTot, " attached markers")
      markers = name(sourcePed, 1:nmTot)
      if(anyNA(markers))
        markers = 1:nmTot
    }

    # Check for already typed members. TODO: Support for partially typed members
    typed = typedMembers(sourcePed)
    if(length(bad <- intersect(allids, typed)))
      stop2("Individual is already genotyped: ", toString(bad))

    # Select markers from source and transfer to truePed (if neccessary)
    truePed = switch(source,
       true = selectMarkers(truePed, markers),
       numerator = transferMarkers(from = selectMarkers(numeratorPed, markers),
                                   to = truePed),
       denominator = transferMarkers(from = selectMarkers(denominatorPed, markers),
                                     to = truePed))
  }

  # Plot
  if (isTRUE(plot) || plot == "plotOnly") {
    tp = selectMarkers(truePed, NULL)
    if(identical(tp, selectMarkers(numeratorPed, NULL))) {
      peds = list(numeratorPed, denominatorPed)
      frms = c("Numerator (True)", "Denominator")
    }
    else if(identical(tp, selectMarkers(denominatorPed, NULL))) {
      peds = list(numeratorPed, denominatorPed)
      frms = c("Numerator", "Denominator (True)")
    }
    else {
      peds = list(numeratorPed, denominatorPed, truePed)
      frms = c("Numerator", "Denominator", "True")
    }

    plotPedList(peds, newdev = TRUE, titles = frms,
                hatched = typedMembers,
                fill = list(green = allids),
                marker = match(plotMarkers, markers))

    if (plot == "plotOnly")
      return()
  }

  # Set seed once
  set.seed(seed)

  # Simulate nsim complete profiles from truePed
  if(verbose)
    message(sprintf("Simulating %d profile%s from the true pedigree ...\n",
                    nsim, pluralise(nsim)),
            appendLF = FALSE)

  allsims = profileSim(truePed, N = nsim, ids = allids, simplify1 = FALSE, verbose = FALSE)

  if(verbose)
    message("done")

  # List of input parameters
  params = list(markers = markers, nsim = nsim, threshold = threshold, seed = seed,
                disableMutations = disableMutations, typed = typed, allids = allids)

  if(is.list(ids))
    res = lapply(ids, function(idvec) lrPowerCompute(allsims, numeratorPed, denominatorPed,
                                                     idvec, params = params, verbose = verbose))
  else
    res = lrPowerCompute(allsims, numeratorPed, denominatorPed, ids,
                         params = params, verbose = verbose)

  # Timing
  if(verbose)
    message("Total time used: ", ftime(st))

  # Return results
  res
}


lrPowerCompute = function(sims, numeratorPed, denominatorPed, ids, params, verbose = TRUE) {

  if(verbose)
    message(sprintf("Computing LR distribution for individual%s %s ... ",
                    pluralise(length(ids)), toString(ids)), appendLF = FALSE)

  markers = params$markers
  threshold = params$threshold

  # Target individuals: Pre-typed and simulated
  targetsIds = c(params$typed, ids)

  # Trans
  # Compute LR distribution using the ids individuals only
  lrs = vapply(sims, function(s) {
    numerSim = transferMarkers(from = s, to = numeratorPed, ids = targetsIds)
    denomSim = transferMarkers(from = s, to = denominatorPed, ids = targetsIds)
    #print(numerSim);print(denomSim)
    lr = kinshipLR(list(numerSim, denomSim), ref = 2)
    lr$LRperMarker[,1]
  }, FUN.VALUE = numeric(length(markers)))

  # Ensure matrix
  if(length(markers) == 1)
    lrs = rbind(lrs, deparse.level = 0)

  rownames(lrs) = markers

  if(verbose)
    message("done")

  # Results
  LRperSim = apply(lrs, 2, prod)
  meanLRperMarker = apply(lrs, 1, mean)
  meanLR = mean(LRperSim)
  meanLogLR = mean(log10(LRperSim))
  IP = sapply(threshold, function(thr) mean(LRperSim >= thr))
  names(IP) = threshold

  params$ids = ids
  structure(list(LRperSim = LRperSim, meanLRperMarker = meanLRperMarker,
                 meanLR = meanLR, meanLogLR = meanLogLR, IP = IP,
                 params = params), class = "LRpowerResult")

}

#' @export
print.LRpowerResult = function(x, ...) {
  cat("Mean LR:", round(x$meanLR, 3), "\n")
  cat("Mean log10(LR):", round(x$meanLogLR, 3), "\n")
  ip = x$IP
  cat("Estimated power:", if(!length(ip)) NA, "\n")
  for(i in seq_along(ip))
    cat(sprintf("  P(LR >= %s) = %.3g\n", names(ip)[i], ip[i]))
}


