#' Likelihood ratios for kinship testing
#'
#' This function computes likelihood ratios (LRs) for a list of pedigrees. One
#' of the pedigrees (the last one, by default) is designated as 'reference', to
#' be used in the denominator in all LR calculations. To ensure that all
#' pedigrees use the same data set, one of the pedigrees may be chosen as
#' 'source', from which data is transferred to all the other pedigrees.
#'
#' By default, all markers are assumed to be unlinked. To accommodate linkage, a
#' genetic map may be supplied with the argument `linkageMap`. This requires the
#' software MERLIN to be installed.
#'
#' @param ... Pedigree alternatives. Each argument should be either a single
#'   `ped` object or a list of such. The pedigrees may be named; otherwise they
#'   are assigned names "H1", "H2", ... automatically.
#'
#'   It is also possible to pass a single `list` containing all the pedigrees.
#' @param ref An index or name indicating which of the input pedigrees should be
#'   used as "reference pedigree", i.e., used in the denominator of each LR. If
#'   NULL (the default), the last pedigree is used as reference.
#' @param source An index or name designating one of the input pedigrees as
#'   source for marker data. If given, marker data is transferred from this to
#'   all the other pedigrees (replacing any existing markers). The default
#'   action (`source = NULL`) is as follows: If all pedigree have attached
#'   markers, no transfers are done. If exactly one of the pedigrees have
#'   attached markers, these are transferred to the others. all other cases give
#'   an error.
#' @param markers A vector of marker names or indices indicating which markers
#'   should be included. If NULL (the default) all markers are used.
#' @param linkageMap If this is non-NULL, the markers are interpreted as being
#'   linked, and likelihoods will be computed by an external call to MERLIN.
#'
#'   The supplied object should be either:
#'
#'   * a data frame, whose first three columns must be (i) chromosome (ii)
#'   marker name (iii) centiMorgan position, or
#'
#'   * a map object created with `ibdsim2::uniformMap()` or
#'   `ibdsim2::loadMap()`. This will internally be applied to the attached
#'   markers to produce a suitable data frame as above.
#'
#' @param keepMerlin Either NULL (default) or the path to an existing folder. If
#'   given, MERLIN files are stored here, typically for debugging purposes.
#' @param verbose A logical.
#'
#' @seealso [LRpower()], [pedprobr::likelihoodMerlin()]
#'
#' @return A `LRresult` object, which is essentially a list with entries
#'
#'   * `LRtotal` : A vector of length `L`, where `L` is the number of input
#'   pedigrees. The i'th entry is the total LR (i.e., the product over all
#'   markers) comparing pedigree `i` to the reference pedigree. The entry
#'   corresponding to the reference will always be 1.
#'
#'   * `LRperMarker` : A numerical matrix, where the i'th column contains the
#'   marker-wise LR values comparing pedigree `i` to the reference. The product
#'   of all entries in a column should equal the corresponding entry in
#'   `LRtotal`.
#'
#'   * `likelihoodsPerMarker` : A numerical matrix of the same dimensions as
#'   `LRperMarker`, but where the entries are likelihood of each pedigree for
#'   each marker.
#'
#'   * `time` : Elapsed time
#'
#' @author Magnus Dehli Vigeland and Thore Egeland
#'
#' @examples
#'
#' ### Example 1: Full vs half sibs
#'
#' # Simulate 5 markers for a pair of full sibs
#' ids = c("A", "B")
#' sibs = nuclearPed(children = ids)
#' sibs = simpleSim(sibs, N = 5, alleles = 1:4, ids = ids, seed = 123)
#'
#' # Create two alternative hypotheses
#' halfsibs = relabel(halfSibPed(), old = 4:5, new = ids)
#' unrel = list(singleton("A"), singleton("B"))
#'
#' # Compute LRs. By default, the last ped is used as reference
#' kinshipLR(sibs, halfsibs, unrel)
#'
#' # Input pedigrees can be named, reflected in the output
#' kinshipLR(S = sibs, H = halfsibs, U = unrel)
#'
#' # Select non-default reference (by index or name)
#' kinshipLR(S = sibs, H = halfsibs, U = unrel, ref = "H")
#'
#' # Alternative syntax: List input
#' peds = list(S = sibs, H = halfsibs, U = unrel)
#' kinshipLR(peds, ref = "H", source = "S", verbose = TRUE)
#'
#' # Detailed results
#' res = kinshipLR(peds)
#' res$LRperMarker
#' res$likelihoodsPerMarker
#'
#'
#' ### Example 2: Separating grandparent/halfsib/uncle-nephew
#' \donttest{
#' # Requires ibdsim2 and MERLIN
#' if(requireNamespace("ibdsim2", quietly = TRUE) && pedprobr::checkMerlin()) {
#'
#'   # Define pedigrees
#'   ids = c("A", "B")
#'   H = relabel(halfSibPed(),   old = c(4,5), new = ids)
#'   U = relabel(avuncularPed(), old = c(3,6), new = ids)
#'   G = relabel(linearPed(2),   old = c(1,5), new = ids)
#'
#'   # Attach FORCE panel of SNPs to G
#'   G = setSNPs(G, FORCE[1:10, ])  # use all for better results
#'
#'   # Simulate recombination pattern in G
#'   map = ibdsim2::loadMap("decode19", uniform = TRUE)   # unif for speed
#'   ibd = ibdsim2::ibdsim(G, N = 1, ids = ids, map = map)
#'
#'   # Simulate genotypes conditional on pattern
#'   G = ibdsim2::profileSimIBD(G, ibdpattern = ibd)
#'
#'   # Compute LR (genotypes are automatically transferred to H and U)
#'   kinshipLR(H, U, G, linkageMap = map)
#' }}
#'
#' @export
kinshipLR = function(..., ref = NULL, source = NULL, markers = NULL, linkageMap = NULL,
                     keepMerlin = NULL, verbose = FALSE) {
  st = proc.time()

  x = list(...)
  if(length(x) == 1)
    x = x[[1]]

  if(length(x) < 2 || is.ped(x))
    stop2("The input must contain at least two pedigrees")

  # Check that all ... are pedigrees. NB: Abbrev args will end up there!
  if(!all(sapply(x, function(ped) is.ped(ped) || is.pedList(ped)))) {
    msg = NULL

    if(!is.null(nms <- names(x))) {
      if(any(err <- startsWith(nms, "verb")))
        msg = sprintf("Did you accidentally write `%s` instead of `verbose`?", nms[err][1])
      else if(any(err <- startsWith(nms, "mark")))
        msg = sprintf("Did you accidentally write `%s` instead of `markers`?", nms[err][1])

      if(!is.null(msg))
        msg = paste("\n\nThis function disallows abbreviated argument names.", msg)
    }
    stop2("The input is not a list of pedigrees", msg)
  }

  # Check `ref`
  if(is.null(ref))
    ref = length(x)
  if(isNumber(ref, 1, length(x)))
    refIdx = as.integer(ref)
  else if(is.character(ref) && ref %in% names(x))
    refIdx = match(ref, names(x))
  else
    stop2("Invalid value for `ref`: ", ref)
  if(verbose)
    cat("Reference pedigree: ", ref, "\n")

  if(is.null(source)) {
    # Identify peds without marker data
    empty = !sapply(x, hasMarkers)

    if(all(empty))
      stop2("None of the pedigrees has attached marker data")

    if(sum(!empty) == 1)
      source = which(!empty)
    else if(any(empty))
      stop2("Data transfer needed, but more than one possible source. Please specify `source`")
  }

  # If source given, transfer to all
  if(!is.null(source)) {
    if(length(source) != 1)
      stop2("`source` must have length 1: ", source)
    if(verbose)
      cat("Source pedigree: ", source, "\n")
    srcPed = x[[source]]
    if(is.null(srcPed))
      stop2("Unknown source pedigree: ", source)
    if(nMarkers(srcPed) == 0)
      stop2("The source pedigree has no attached markers")
    if(!is.null(markers)) {
      srcPed = selectMarkers(srcPed, markers)
      markers = NULL # important!
    }
    if(!is.null(linkageMap))
      srcPed = lumpAlleles(srcPed, always = TRUE, verbose = verbose)

    x = lapply(x, transferMarkers, from = srcPed)
  }

  # By default, use all markers
  if(is.null(markers)) {
    nm = vapply(x, nMarkers, FUN.VALUE = 1)
    if(!all(nm == nm[1]))
      stop2("When `markers = NULL`, all pedigrees must have the same number of attached markers: ", nm)
    markers = seq_len(nm[1])
  }
  else {
    x = lapply(x, selectMarkers, markers)
  }

  if(verbose)
    cat("Number of markers: ", length(markers), "\n")

  # Fix hypothesis names
  hypnames = names(x)
  if(is.null(hypnames))
    hypnames = paste0("H", seq_along(x))
  else {
    unnamed = hypnames == "" | is.na(hypnames)
    hypnames[unnamed] = paste0("H", which(unnamed))
  }
  if(dup <- anyDuplicated(hypnames))
    stop2("Duplicated hypothesis name: ", hypnames[dup])

  names(x) = hypnames


  # Linked markers: MERLIN ---------------------------------------------
  if(hasLinkedMarkers(x[[1]]) && is.null(linkageMap))
    stop2("Linked markers detected, but no `linkageMap` provided")

  if(!is.null(linkageMap)) {
    if(verbose)
      cat("\nLinkage map detected - preparing MERLIN or MINX\n")

    if(!checkMerlin("merlin", version = FALSE, error = FALSE))
      stop2("Kinship analysis with linked markers requires MERLIN to be installed and available in the search path")

    # If ibdsim recombination map, apply it to the attached markers
    if(inherits(linkageMap, "genomeMap") || inherits(linkageMap, "chromMap")) {
      if(!requireNamespace("ibdsim2", quietly = TRUE))
        stop2("The supplied `linkageMap` requires the `ibdsim2` package, which seems unavailable. Install `ibdsim2` and try again.")
      if(verbose)
        cat("Converting attached marker positions from MB to CM\n")

      physMap = getMap(x[[1]])
      linkageMap = ibdsim2:::convertMap(physMap, genomeMap = linkageMap)
    }

    # Lump all peds, if not already done (a bit of a hack)
    #if(is.null(source))
    #  x = lapply(x, lumpAlleles, verbose = verbose)

    if(!is.null(keepMerlin)) {
      if(!dir.exists(keepMerlin))
        stop2("`keepMerlin` must point to an existing folder: ", keepMerlin)
      cleanup = FALSE
      dir = keepMerlin
    }
    else {
      cleanup = TRUE
      dir = tempdir()
    }

    merlinFUN = function(x, verb)
      likelihoodMerlin(x, markers = markers, linkageMap = linkageMap, perChrom = TRUE, logbase = exp(1),
                       checkpath = FALSE, cleanup = cleanup, dir = dir, verbose = verb)

    lnLikList = lapply(seq_along(x), function(i) merlinFUN(x[[i]], verbose && (i == 1)))

    lnLikChrom = do.call(cbind, lnLikList)
    rownames(lnLikChrom) = names(lnLikList[[1]]) # chrom labels

    lnLik = colSums(lnLikChrom)
    lnDiff = lnLik - lnLik[[refIdx]]
    LRtotal = signif(exp(lnDiff), 3)
    names(LRtotal) = paste0(hypnames, ":", hypnames[refIdx])

    if(any(LRtotal == Inf & lnDiff < Inf))
      warning("Overflow!\nOutput entries `lnLik` and `lnLikChrom` may still be informative",
              immediate. = TRUE, call. = FALSE)
    if(any(LRtotal == 0 & lnDiff > -Inf))
      warning("Underflow!\nOutput entries `lnLik` and `lnLikChrom` may still be informative",
              immediate. = TRUE, call. = FALSE)

    lnDiffChrom = lnLikChrom - lnLikChrom[,refIdx]
    LRchrom = signif(exp(lnDiffChrom), 3)

    res = list(LRtotal = LRtotal, lnLik = lnLik,
               LRchrom = LRchrom, lnLikChrom = lnLikChrom,
               time = proc.time() - st)
    return(structure(res, class = "LRresult"))
  }

  # Break all loops (NB: rapply() doesn't work here, since is.list(ped) = TRUE)
  breaklp = function(a) {
    if(is.pedList(a))
      return(lapply(a, breaklp))
    if(is.singleton(a))
      return(a)
    breakLoops(a, verbose = FALSE)
  }

  x_loopfree = lapply(x, breaklp)

  # compute likelihoods
  liks = lapply(x_loopfree, function(xx) likelihood(xx, markers = markers))
  likelihoodsPerMarker = do.call(cbind, liks)

  # LR per marker and total
  LRperMarker = do.call(cbind, lapply(1:length(x), function(j) liks[[j]]/liks[[refIdx]]))

  # Create names for LR comparisons
  colnames(LRperMarker) = paste0(hypnames, ":", hypnames[refIdx])

  # Total LR
  LRtotal = apply(LRperMarker, 2, prod)

  # Use marker names in output
  markernames = name(x[[1]], markers)
  if(any(NAnames <- is.na(markernames)))
    markernames[NAnames] = sprintf("<%d>", which(NAnames))
  rownames(likelihoodsPerMarker) = rownames(LRperMarker) = markernames

  # Output `LRresult` object
  structure(list(
    LRtotal = LRtotal,
    LRperMarker = LRperMarker,
    likelihoodsPerMarker = likelihoodsPerMarker,
    time = proc.time() - st),
    class = "LRresult")
}

#' @export
print.LRresult = function(x, ...) {
  print(x$LRtotal)
}

