#' Likelihood ratios for kinship testing
#'
#' This function computes likelihood ratios (LRs) for a given a list of
#' pedigrees with attached markers. The user must indicate which of the
#' pedigrees is the 'reference', which will be used in the denominator in each
#' LR.
#'
#' @param ... A list of pedigree alternatives. Each alternative should be either
#'   a single `ped` object or a list of such.
#' @param ref An index or name indicating which of the input pedigrees should be
#'   used as "reference pedigree", i.e., used in the denominator of each LR. If
#'   NULL (the default), the last pedigree is used as reference.
#' @param source An index or name designating one of the input pedigrees as
#'   source for marker data. If given, marker data is transferred from this to
#'   all the other pedigrees (replacing any existing markers). The default
#'   action (`source = NULL`) is as follows: If all pedigrees already have
#'   attached markers, no transfers are done. If some pedigrees are empty,
#'   marker data is transferred to those, using the first nonempty pedigree as
#'   source.
#' @param markers A vector of marker names or indices indicating which markers
#'   should be included. If NULL (the default) all markers
#'   are used.
#' @param verbose A logical.
#'
#' @seealso [LRpower()], [pedtools::transferMarkers()]
#'
#' @return A `LRresult`object, which is essentially a list with entries
#'
#'   * `LRtotal` : Total likelihood ratios
#'
#'   * `LRperMarker` : Likelihood ratios for each marker
#'
#'   * `likelihoodsPerMarker` : Likelihoods for each marker
#'
#'   * `time` : Elapsed time
#'
#' @author Magnus Dehli Vigeland and Thore Egeland
#'
#' @examples
#'
#' # Simulate 5 markers for a pair of full sibs
#' sibs = nuclearPed(children = c("A", "B"))
#' sibs = simpleSim(sibs, N = 5, alleles = 1:4, ids = c("A", "B"), seed = 123)
#'
#' # Create two alternative hypotheses
#' halfsibs = relabel(halfSibPed(), old = 4:5, new = c("A", "B"))
#' unrel = list(singleton("A"), singleton("B"))
#'
#' # Compute LR with 'unrel' as reference (by default, the last ped is used)
#' res = kinshipLR(sibs, halfsibs, unrel)
#' res
#'
#' # Detailed results
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

  if(!all(sapply(x, function(ped) is.ped(ped) || is.pedList(ped))))
    stop2("The input is not a list of pedigrees")

  # Check `ref`
  if(is.null(ref))
    refIdx = length(x)
  else if(is_number(ref, 1, length(x)))
    refIdx = as.integer(ref)
  else if(is.character(ref) && ref %in% names(x))
    refIdx = match(ref, names(x))
  else
    stop2("Invalid value for `ref`: ", ref)
  if(verbose)
    message("Using pedigree ", refIdx, " as reference")

  # If explicit source given, transfer to all
  if(!is.null(source)) {
    if(verbose)
      message("Transfering marker data from pedigree ", source, " to all others")
    srcPed = x[[source]]
    if(is.null(srcPed))
      stop2("Unknown source pedigree: ", source)
    if(nMarkers(srcPed) == 0)
      stop2("The source pedigree has no attached markers")
    if(!is.null(markers))
      srcPed = selectMarkers(srcPed, markers)
    x = lapply(x, transferMarkers, from = srcPed)
  }
  else {
    # Any empty pedigrees?
    empty = !sapply(x, hasMarkers)

    if(all(empty))
      stop2("None of the pedigrees has attached marker data")

    if(any(empty)) {
      source = which(!empty)[1] # use first non-empty as source
      srcPed = x[[source]]
      if(!is.null(markers))
        srcPed = selectMarkers(srcPed, markers)
      if(verbose)
        message("Transfering marker data from pedigree ", source, " to empty pedigree(s): ", toString(which(empty)))
      x[empty] = lapply(x[empty], transferMarkers, from = srcPed)
    }
  }

  # By default, use all markers
  if(is.null(markers)) {
    nm = vapply(x, nMarkers, FUN.VALUE = 1)
    if(!all(nm == nm[1]))
      stop2("When `markers = NULL`, all pedigrees must have the same number of attached markers: ", nm)
    markers = seq_len(nm[1])
  }

  # Extract names of selected markers
  markernames = name(x[[1]], markers)

  # Fix missing marker names
  if(any(NAnames <- is.na(markernames)))
    markernames[NAnames] = paste0("M", which(NAnames))

  # Fix hypothesis names
  hypnames = names(x)
  if(is.null(hypnames))
    hypnames = paste0("H", seq_along(x))
  else {
    unnamed = hypnames == "" | is.na(hypnames)
    hypnames[unnamed] = paste0("H", which(unnamed))
  }
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
  liks = lapply(x_loopfree,
                function(xx) vapply(markers, function(i) likelihood(xx, marker1 = i),
                                    FUN.VALUE = 1))
  likelihoodsPerMarker = do.call(cbind, liks)

  # LR per marker and total
  LRperMarker = do.call(cbind, lapply(1:length(x), function(j) liks[[j]]/liks[[refIdx]]))

  # Create names for LR comparisons
  colnames(LRperMarker) = paste0(hypnames, ":", hypnames[refIdx])

  # total LR
  LRtotal = apply(LRperMarker, 2, prod)

  # output
  rownames(likelihoodsPerMarker) = rownames(LRperMarker) = markernames

  structure(list(
    LRtotal = LRtotal,
    LRperMarker = LRperMarker,
    likelihoodsPerMarker = likelihoodsPerMarker,
    time = proc.time() - st),
    class = "LRresult")
}

# Alias
LR = kinshipLR


#' @export
print.LRresult = function(x, ...) {
  cat("Total LR:\n")
  print(x$LRtotal)
}
