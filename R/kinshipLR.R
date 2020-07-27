#' Likelihood ratios for kinship testing
#'
#' This function computes likelihood ratios (LRs) for a list of pedigrees. One
#' of the pedigrees (the last one, by default) is designated as 'reference', to
#' be used in the denominator in all LR calculations. To ensure that all
#' pedigrees use the same data set, one of the pedigrees may be chosen as
#' 'source', from which data is transferred to all the other pedigrees.
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
#' @param verbose A logical.
#'
#' @seealso [LRpower()]
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
#' @export
kinshipLR = function(..., ref = NULL, source = NULL, markers = NULL, verbose = FALSE) {
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
    message("Reference pedigree: ", ref)

  if(is.null(source)) {
    # Identity peds without marker data
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
      message("Source pedigree: ", source)
    srcPed = x[[source]]
    if(is.null(srcPed))
      stop2("Unknown source pedigree: ", source)
    if(nMarkers(srcPed) == 0)
      stop2("The source pedigree has no attached markers")
    if(!is.null(markers)) {
      srcPed = selectMarkers(srcPed, markers)
      markers = NULL # important!
    }
    x = lapply(x, transferMarkers, from = srcPed)
  }

  # By default, use all markers
  if(is.null(markers)) {
    nm = vapply(x, nMarkers, FUN.VALUE = 1)
    if(!all(nm == nm[1]))
      stop2("When `markers = NULL`, all pedigrees must have the same number of attached markers: ", nm)
    markers = seq_len(nm[1])
  }

  if(verbose)
    message("Number of markers: ", length(markers))

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

  # Break all loops (NB: would use rapply, but doesnt work since is.list(ped) = TRUE
  breaklp = function(a) {
    if(is.pedList(a))
      return(lapply(a, breaklp))
    if(is.singleton(a))
      return(a)
    breakLoops(a, verbose = FALSE)
  }

  x_loopfree = lapply(x, breaklp)

  # compute likelihoods
  liks = lapply(x_loopfree, function(xx) likelihood(xx, markers))
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
