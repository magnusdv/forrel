#' Power simulation for kinship LR
#'
#' This functions uses simulations to estimate the distribution of LR in a given
#' kinship testing scenario, given a certain subset of individuals who are
#' available for testing. Some individuals individuals may already be genotyped;
#' all simulations are then conditional on these. In the most general setting,
#' three pedigrees may be involved: The two pedigrees being compared, and the
#' true relationship (which may differ from the other two).
#'
#' @param numeratorPed,denominatorPed `ped` objects (or lists of such),
#'   describing the two relationships under comparison.
#' @param truePed A `ped` object (or a list of such), describing the true
#'   relationship. By default equal to `numeratorPed`.
#' @param ids Individuals available for genotyping.
#' @param source Either "true" (default), "numerator" or "denominator",
#'   indicating which pedigree is used as source for marker data.
#' @param nsim A positive integer: the number of simulations.
#' @param threshold A numeric vector with one or more positive numbers used as
#'   LR tresholds.
#' @param disableMutations Not implemented yet.
#' @param seed A numeric seed for the random number generator (optional)
#'
#' @inheritParams exclusionPower
#' @examples
#'
#' # Paternity LR of siblings
#' claim = nuclearPed(fa = "A", mo = "NN", children = "B")
#' unrel = list(singleton("A"), singleton("B"))
#' truth = nuclearPed(children = c("A", "B"))
#'
#' # Simulation parameters
#' nsim = 10
#' thresh = 1
#' ids = c("A", "B")
#'
#' # Simulation 1:
#' als = 1:5
#' afr = runif(5)
#' afr = afr/sum(afr)
#'
#' pow1 = LRpower(claim, unrel, truth, ids = ids, nsim = nsim, threshold = thresh,
#'                alleles = als, afreq = afr, seed = 123)
#' pow1
#'
#' # Simulation 2: Same, but using an attached marker
#' m = marker(truth, alleles = als, afreq = afr)
#' truth = setMarkers(truth, m)
#'
#' pow2 = LRpower(claim, unrel, truth, ids = ids, nsim = nsim, threshold = thresh,
#'                markers = 1, seed = 123)
#' stopifnot(identical(pow1$LRperSim, pow2$LRperSim))
#'
#' # Founder inbreeding in true pedigree
#' founderInbreeding(truth, founders(truth)) = 0.5
#' pow3 = LRpower(claim, unrel, truth, ids = ids, nsim = nsim, threshold = thresh,
#'                markers = 1, seed = 123, plot = TRUE)
#' pow3
#'
#'
#' @export
LRpower = function(numeratorPed, denominatorPed, truePed = numeratorPed, ids, markers = NULL,
                   source = "true", nsim = 1, threshold = NULL, disableMutations = NA,
                   alleles = NULL, afreq = NULL, Xchrom = FALSE, knownGenotypes = NULL,
                   plot = FALSE, plotMarkers = NULL, seed = NULL, verbose = TRUE) {
  st = Sys.time()
  ids = as.character(ids)
  source = match.arg(source, c("true", "numerator", "denominator"))

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
    if(length(alleles) == 1 && is.numeric(alleles))
      alleles = seq_len(alleles)

    # Create and attach locus to both pedigrees
    locus = list(alleles = alleles, afreq = afreq, chrom = if (Xchrom) 23 else NA)
    truePed = setMarkers(truePed, alleleMatrix = am, locusAttributes = locus)

    markers = 1
    typed = typedMembers(truePed)
    hasMut = FALSE
    disableMutations = FALSE # don't do anything
  }
  else {
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
    if(length(bad <- intersect(ids, typed)))
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

    plotPedList(peds, newdev = TRUE, frametitles = frms,
                shaded = function(p) c(ids, typedMembers(p)),
                col = list(red = ids),
                marker = match(plotMarkers, markers))

    if (plot == "plotOnly")
      return()
  }

  # Set seed once
  set.seed(seed)

  # Simulate nsim complete profiles from truePed
  if(verbose)
    message("Simulating ", nsim, " profiles from true pedigree...", appendLF = F)

  allsims = profileSim(truePed, N = nsim, ids = ids)

  if(verbose)
    message("done\nComputing likelihood ratios...", appendLF = F)

  # Compute the exclusion power of each marker
  lrs = vapply(allsims, function(s) {
    numerSim = transferMarkers(from = s, to = numeratorPed)
    denomSim = transferMarkers(from = s, to = denominatorPed)

    lr = LR(list(numerSim, denomSim), ref = 2)
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

  # Timing
  if(verbose)
    message("Total time used: ", format(Sys.time() - st, digits = 3))

  # Lits of input parameters
  params = list(ids = ids, markers = markers,
                nsim = nsim, threshold = threshold, seed = seed,
                disableMutations = disableMutations)

  structure(list(LRperSim = LRperSim, meanLRperMarker = meanLRperMarker,
                 meanLR = meanLR, meanLogLR = meanLogLR, IP = IP,
                 params = params), class = "LRpowerResult")
}

#' @export
print.LRpowerResult = function(x, ...) {
  cat("Mean total LR:", round(x$meanLR, 3), "\n")
  cat("Mean total log10(LR):", round(x$meanLogLR, 3), "\n")
  ip = x$IP
  cat("Estimated inclusion powers:", if(!length(ip)) NA, "\n")
  for(i in seq_along(ip))
    cat(sprintf("  P(LR > %s) = %.3g\n", names(ip)[i], ip[i]))
}


