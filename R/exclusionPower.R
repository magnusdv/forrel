#' Power of exclusion
#'
#' Computes the power (of a single marker, or for a collection of markers) of
#' excluding a claimed relationship, given the true relationship.
#'
#' This function implements the formula for exclusion power as defined and
#' discussed in (Egeland et al., 2014).
#'
#' @param claimPed A `ped` object (or a list of such), describing the claimed
#'   relationship. If a list, the sets of ID labels must be disjoint, that is,
#'   all ID labels must be unique.
#' @param truePed A `ped` object (or a list of such), describing the true
#'   relationship. ID labels must be consistent with `claimPed`.
#' @param ids Individuals available for genotyping.
#' @param markers A vector indicating the names or indices of markers attached
#'   to the source pedigree. If NULL (default), then all markers attached to the
#'   source pedigree are used. If `alleles` or `afreq` is non-NULL, then this
#'   parameter is ignored.
#' @param source Either "claim" (default) or "true", deciding which pedigree is
#'   used as source for marker data.
#' @param disableMutations This parameter determines how mutation models are
#'   treated. Possible values are as follows:
#'
#'   * `NA` (the default): Mutations are disabled only for those markers whose
#'   known genotypes are compatible with both `claimPed` and `truePed`. This is
#'   determined by temporarily removing all mutation models and checking which
#'   markers have nonzero likelihood in both alternatives.
#'
#'   * `TRUE`: Mutations are disabled for all markers.
#'
#'   * `FALSE`: No action is done to disable mutations.
#'
#'   * A vector containing the names or indices of those markers for which
#'   mutations should be disabled.
#' @param exactMaxL A positive integer, or `Inf` (default). Exact EPs are
#'   calculated for markers whose number of alleles is less or equal to
#'   `exactMaxL`; remaining markers are handled by simulation.
#' @param nsim A positive integer; the number of simulations used for markers
#'   whose number of alleles exceeds `exactMaxL`.
#' @param seed A numeric seed for the random number generator (optional).
#' @param alleles,afreq,Xchrom If these are given, they are used (together with
#'   `knownGenotypes`) to create a marker object on the fly.
#' @param knownGenotypes A list of triplets `(a, b, c)`, indicating that
#'   individual `a` has genotype `b/c`. Ignored unless `alleles` or `afreq` is
#'   non-NULL.
#' @param plot Either a logical or the character "plotOnly". If the latter, a
#'   plot is drawn, but no further computations are done.
#' @param plotMarkers A vector of marker names or indices whose genotypes are to
#'   be included in the plot.
#' @param verbose A logical.
#'
#' @return If `plot = "plotOnly"`, the function returns NULL after producing the
#'   plot.
#'
#'   Otherwise, the function returns an `EPresult` object, which is essentially
#'   a list with the following entries:
#'
#'   * `EPperMarker`: A numeric vector containing the exclusion power of each
#'   marker. If the known genotypes of a marker are incompatible with the true
#'   pedigree, the corresponding entry is `NA`.
#'
#'   * `EPtotal`: The total exclusion power, computed as `1 - prod(1 -
#'   EPperMarker, na.rm = TRUE)`.
#'
#'   * `expectedMismatch`: The expected number of markers giving exclusion,
#'   computed as `sum(EPperMarker, na.rm = TRUE)`.
#'
#'   * `distribMismatch`: The probability distribution of the number of markers
#'   giving exclusion. This is given as a numeric vector of length `n+1`, where
#'   `n` is the number of nonzero elements of `EPperMarker`. The vector has
#'   names `0:n`.
#'
#'   * `time`: The total computation time.
#'
#'   * `params`: A list containing the (processed) parameters `ids`, `markers`
#'   and `disableMutations`.
#'
#'
#' @author Magnus Dehli Vigeland
#'
#' @references T. Egeland, N. Pinto and M.D. Vigeland, *A general approach to
#'   power calculation for relationship testing.* Forensic Science
#'   International: Genetics 9 (2014): 186-190.
#'   \doi{10.1016/j.fsigen.2013.05.001}
#'
#' @examples
#'
#' ############################################
#' ### A standard case paternity case:
#' ### Compute the power of exclusion when the claimed father is in fact
#' ### unrelated to the child.
#' ############################################
#'
#' # Claim: Individual 1 is the father of individual 3
#' claim = nuclearPed(nch = 1, sex = 2)
#'
#' # Truth: 1 and 3 are unrelated
#' true = list(singleton(id = 1), singleton(id = 3, sex = 2))
#'
#' # Attach 3 markers
#' m1 = marker(claim, alleles = 1:2)
#' m2 = marker(claim, alleles = 1:3)
#' m3 = marker(claim, alleles = 1:2, chrom = "X")
#' claim = setMarkers(claim, list(m1, m2, m3))
#'
#' # Compute EP when father and child is available for genotyping
#' exclusionPower(claim, true, ids = c(1,3))
#'
#' # Suppose child is genotyped
#' genotype(claim, marker = 1, id = 3) = c(1, 1)
#' genotype(claim, marker = 2, id = 3) = c(1, 1)
#' genotype(claim, marker = 3, id = 3) = c(1, 2)
#'
#' # Compute EP when father is available
#' exclusionPower(claim, true, ids = 1)
#'
#'
#' ############################################
#' ### Two females claim to be mother and daughter, but are in reality sisters.
#' ### We compute the power of various markers to reject the claim.
#' ############################################
#'
#' mother_daughter = nuclearPed(1, sex = 2)
#' sisters = relabel(nuclearPed(2, sex = c(2, 2)), c(101, 102, 2, 3))
#' ids = 2:3
#'
#' # SNP with MAF = 0.1:
#' PE1 = exclusionPower(claimPed = mother_daughter, truePed = sisters,
#'                      ids = ids, alleles = 2, afreq = c(0.9, 0.1))
#'
#' # Tetra-allelic marker with one major allele:
#' PE2 = exclusionPower(claimPed = mother_daughter, truePed = sisters,
#'                      ids = ids, alleles = 4, afreq = c(0.7, 0.1, 0.1, 0.1))
#'
#' stopifnot(all.equal(c(PE1$EPtotal, PE2$EPtotal), c(0.00405, 0.03090)))
#'
#' ### How does the power change if the true pedigree is inbred?
#' sisters_LOOP = addParents(sisters, 101, father = 201, mother = 202)
#' sisters_LOOP = addParents(sisters_LOOP, 102, father = 201, mother = 203)
#'
#'
#' # SNP with MAF = 0.1:
#' PE3 = exclusionPower(claimPed = mother_daughter, truePed = sisters_LOOP,
#'                      ids = ids, alleles = 2, afreq = c(0.9, 0.1))
#'
#' stopifnot(all.equal(PE3$EPtotal, 0.00765))
#'
#' @importFrom pedprobr likelihood
#' @export
exclusionPower = function(claimPed, truePed, ids, markers = NULL, source = "claim",
                          disableMutations = NA, exactMaxL = Inf, nsim = 1000, seed = NULL,
                          alleles = NULL, afreq = NULL, knownGenotypes = NULL, Xchrom = FALSE,
                          plot = FALSE, plotMarkers = NULL, verbose = TRUE) {

  st = Sys.time()

  # Ensure `ids` is a list of character vectors
  if(!is.list(ids))
    ids = list(ids)
  ids = lapply(ids, as.character)
  idslabs = sapply(ids, paste, collapse = "-")
  allids = unique.default(unlist(ids))

  # Number of `ids` vectors
  NI = length(ids)

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
    claimPed = setMarkers(claimPed, alleleMatrix = am, locusAttributes = locus)
    truePed = setMarkers(truePed, alleleMatrix = am, locusAttributes = locus)

    markers = 1
    typed = typedMembers(truePed)
    hasMut = FALSE
    disableMutations = FALSE # don't do anything
  }
  else {
    sourcePed = switch(source, claim = claimPed, true = truePed,
                       stop2("`source` must be either 'claim' or 'true': ", source))

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

    # Transfer marker data to the other pedigree
    switch(source,
           claim = {
             claimPed = selectMarkers(claimPed, markers)
             truePed = transferMarkers(from = claimPed, to = truePed)
           },
           true = {
             truePed = selectMarkers(truePed, markers)
             claimPed = transferMarkers(from = truePed, to = claimPed)
           })
  }

  # Plot
  if (isTRUE(plot) || plot == "plotOnly") {
    plotPedList(list(claimPed, truePed),
                newdev = TRUE,
                frametitles = c("Claim", "True"),
                shaded = function(p) c(allids, typedMembers(p)),
                col = list(red = allids),
                marker = match(plotMarkers, markers))

    if (plot == "plotOnly")
      return()
  }

  nMark = length(markers)

  ### Mutation disabling

  # Which of the markers allow mutations?
  hasMut = allowsMutations(claimPed)

  if(!is.null(disableMutations) && !isFALSE(disableMutations)) {

    if(isTRUE(disableMutations))
      disableMutations = which(hasMut)
    else if(identical(disableMutations, NA)) {
      # Disable mutations for markers consistent in both true and claim
      cons = consistentMarkers(claimPed, hasMut) & consistentMarkers(truePed, hasMut)

      disableTF = logical(nMark)
      disableTF[hasMut] = cons
      disableMutations = which(disableTF)
    }
    else {# if numeric or character: Use original pedigree, i.e. before selection
      midx_orig = whichMarkers(sourcePed, markers)
      dis_idx = whichMarkers(sourcePed, disableMutations)
      disableMutations = match(dis_idx, midx_orig) # = index in the after-selection peds!
    }

    # Disable
    if(length(disableMutations)) {
      mutmod(claimPed, disableMutations) = NULL
      mutmod(truePed, disableMutations) = NULL
    }
  }

  ### Baseline likelihoods
  trueBase = likelihood(truePed, seq_len(nMark))
  claimBase = likelihood(claimPed, seq_len(nMark))

  ### Compute the exclusion power of each marker.
  # The result is rectangular, with nMark columns and NI rows
  ep = vapply(seq_len(nMark), function(i) {

    m = markers[i]
    if(verbose)
      message("Marker ", m, ". ", appendLF = FALSE)

    # Disable mutations?
    if(verbose) {
      action = if(i %in% disableMutations) "Mutation model disabled. "
      else if(hasMut[i]) "Mutation model applied. "
      else "No mutation model. "
      message(action, appendLF = FALSE)
    }

    # If impossible in true, return NA
    if(trueBase[i] == 0) {
      if(verbose)
        message("INCOMPATIBLE WITH TRUE PEDIGREE! EP = NA")
      return(rep(NA_real_, NI))
    }

    # If impossible in claim, return 1
    if(claimBase[i] == 0) {
      if(verbose)
        message("INCOMPATIBLE WITH CLAIMED PEDIGREE! EP = 1")
      return(rep(1, NI))
    }

    # Exact of simulation?
    exactCalc = nAlleles(truePed, marker = i) <= exactMaxL

    # Compute/estimate EP value for each `ids` entry
    if(exactCalc) {
      if(verbose)
        message("Computing exact EP.")
      this.ep = unlist(lapply(ids, function(idvec) {
        .EPsingleMarker(claimPed, truePed, idvec, marker = i, verbose = FALSE)
      }))
    }
    else {
      if(verbose)
        message("Estimating EP by simulation.")
      trueSims = markerSim(truePed, ids = allids, N = nsim, partialmarker = i,
                           seed = seed, verbose = FALSE)

      this.ep = unlist(lapply(ids, function(idvec) {
        claimSims = transferMarkers(trueSims, claimPed, ids = c(typed, idvec))
        liks = likelihood(claimSims, markers = 1:nsim)
        mean(liks == 0)
      }))
    }

    if(verbose)
      if(NI == 1) cat(round(this.ep, 3), "\n") else print(setNames(round(this.ep, 3), idslabs))

    this.ep
  }, FUN.VALUE = numeric(NI))

  # List of parameters
  params = list(markers = markers, disableMutations = disableMutations,
                exactMaxL = exactMaxL, nsim = nsim, seed = seed, allids = allids, typed = typed)

  # Summarise over all markers
  if(length(ids) == 1) { # in this case `ep` is just a vector
    res = summariseEP(ep)
    res$params = c(params, list(ids = ids[[1]]))
  }
  else {
    res = lapply(seq_along(ids), function(i) {
      resi = summariseEP(ep[i, ])
      resi$params = c(params, list(ids = ids[[i]]))
      resi
    })
  }

  # Timing
  time = Sys.time() - st
  if(verbose)
    message("Total time used: ", format(time, digits = 3))

  res
}


# Summarise a vector of single-marker exclusion powers
summariseEP = function(epvec) {
  # Total EP
  tot = 1 - prod(1 - epvec, na.rm = TRUE)

  # Expected number of exclusions
  expMis = sum(epvec, na.rm = TRUE)

  # Distribution of number of mismatches
  # This is a sum of different Bernoulli variables, i.e., Poisson binomial.
  n.nonz = sum(epvec > 0, na.rm = TRUE)
  if(n.nonz == 0)
    distrib = c(`0` = 1)
  else {
    distrib = setNames(rep(NA_real_, n.nonz + 1), 0:n.nonz)
    if (requireNamespace("poibin", quietly = TRUE))
      distrib[] = poibin::dpoibin(kk = 0:n.nonz, pp = epvec[!is.na(epvec) & epvec > 0])
    else
      warning("Package `poibin` not found. Cannot compute the distribution of exclusion counts without this; returning NA's")
  }

  # Ad hoc fix: If all epvec are NA, the other results should also be NA (not 0)
  if(all(is.na(epvec))) {
    tot = expMis = distrib = NA
  }

  structure(list(EPperMarker = epvec, EPtotal = tot, expectedMismatch = expMis, distribMismatch = distrib),
            class = "EPresult")
}

#' @export
print.EPresult = function(x, ...) {
  incons = is.na(x$EPperMarker)
  incons_names = if(any(incons)) sprintf("(%s)", toString(x$params$markers[incons])) else ""

  poss = !incons & x$EPperMarker > 0
  poss_names = if(any(poss)) sprintf("(%s)", toString(x$params$markers[poss])) else ""
  nPoss = if(all(incons)) NA else sum(poss)

  expMis = round(x$expectedMismatch, 3)

  cat("Exclusion power with individuals ", toString(x$params$ids), ".\n", sep = "")
  cat("Total EP:", round(x$EPtotal, 3), "\n")
  cat("Markers inconsistent with pedigree:", sum(incons), incons_names, "\n")
  cat("Markers with potential mismatch:", nPoss, poss_names, "\n")
  cat("Expected number of mismatches:", expMis, "\n")
}

  ###############################
  ### Computations start here ###
  ###############################

.EPsingleMarker = function(claimPed, truePed, ids, marker, verbose = TRUE) {

  # Step 1: Genotype combinations incompatible with "claim"
  # Step 2: Probability of incomp under `true` pedigree

  if (is.ped(claimPed))
    claimPed = list(claimPed)
  if (is.ped(truePed))
    truePed = list(truePed)

  claim = selectMarkers(claimPed, marker)
  true = selectMarkers(truePed, marker)

  # Extract ped component of each id, and check that all were found
  compsClaim = getComponent(claim, ids, checkUnique = TRUE)
  if(anyNA(compsClaim))
    stop2("Individuals not found in `claimPed`: ", ids[is.na(compsClaim)])

  # Check that all `ids` are in `true`
  compsTrue = getComponent(true, ids, checkUnique = TRUE)
  if(anyNA(compsTrue))
    stop2("Individuals not found in `truePed`: ", ids[is.na(compsTrue)])

  names(compsClaim) = names(compsTrue) = ids

  ### Step 1
  # Incompatible combinations in each component
  I.g.list = lapply(seq_along(claim), function(i) {
    if(!any(compsClaim == i))
      return(NULL)

    omd = oneMarkerDistribution(claim[[i]], ids = ids[compsClaim == i],
                                partialmarker = 1, verbose = FALSE, eliminate = 0)
    omd == 0
    })

  # Remove NULLs
  I.g.list = I.g.list[lengths(I.g.list) > 0]

  # outer product OR
  outerOR = function(x,y) outer(x, y, FUN = `|`)
  I.g = Reduce(outerOR, I.g.list)

  # If no incompatibilities, return 0
  if(!any(I.g))
    return(0)

  # Extract the TRUE positions (= incompatible combinations)
  incomp = which(I.g, arr.ind = TRUE, useNames = FALSE)
  colnames(incomp) = ids

  if(verbose) {
    cat("Impossible genotype combinations of `ids` given the claimed pedigree:\n")
    inc = sapply(seq_along(ids), function(i) dimnames(I.g)[[i]][incomp[, i]])
    inc = as.data.frame(inc)
    names(inc) = ids
    print(inc)
  }


  ### Step 2: Probability of incomp under `true` pedigree

  # Grid to be used as input to oneMarkerDistribution
  # Recall: Entries now refers to rows in allGenotypes(n)
  incomp.grid = incomp

  # If X: Modify male columns.  TODO: would be nice to clean this up a bit.
  m = getMarkers(claim[[1]], 1)[[1]]
  if(isXmarker(m)) {
    sx = sapply(seq_along(ids), function(i) getSex(true[[compsTrue[i]]], ids[i]))
    nall = nAlleles(true[[1]], 1)
    if(nall > 1) {
      homoz = c(1, cumsum(nall:2) + 1)
      incomp.grid[, sx == 1] = homoz[incomp.grid[, sx == 1]]
    }
  }

  # In the true ped: Sum probs of the incompat combinations
  p.g.list = lapply(seq_along(true), function(i) {
    if(!any(compsTrue == i))
      return(NULL)

    ids.i = ids[compsTrue == i]
    grid.i = unique.matrix(incomp.grid[, ids.i, drop = FALSE])

    omd = oneMarkerDistribution(true[[i]], ids.i, partialmarker = 1,
                          grid.subset = grid.i, verbose = FALSE,
                          eliminate = 1)
    omd
  })

  # Remove NULLs
  p.g.list = p.g.list[lengths(p.g.list) > 0]

  # Reduce to array
  p.g = Reduce("%o%", p.g.list)

  # The entries corresponding to exclusions
  excl = p.g[I.g]

  if(verbose) {
    cat("\nProbabilities under `true` ped:\n")
    inc = cbind(inc, prob = excl)
    print(inc)
  }

  sum(excl)
}

