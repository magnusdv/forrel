#' Likelihood ratios for kinship testing
#'
#' This function computes likelihood ratios (LRs) for a given a list of
#' pedigrees with attached markers. The user must indicate which of the
#' pedigrees is the 'reference', which will be used in the denominator in each
#' LR.
#'
#' @param ... A list of pedigree alternatives. Each alternative should be either a
#'   single `ped` object or a list of such.
#' @param ref A single integer, indicating the index (in the list `x`) of the
#'   reference alternative. This is used in the denominator of each LR. If NULL
#'   (the default), the last pedigree is used as reference.
#' @param markers A vector of integers, indexing which markers should be
#'   included. If NULL (the default) all markers are used.
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
kinshipLR = function(..., ref = NULL, markers = NULL, verbose = FALSE) {
  st = proc.time()

  x = list(...)
  if(length(x) == 1)
    x = x[[1]]

  if(is.ped(x) || length(x) == 1)
    stop2("The input must contain at least two pedigrees")

  if(!all(sapply(x, function(ped) is.ped(ped) || is.pedList(ped))))
    stop2("The input is not a list of pedigrees")

  if(is.null(ref))
    ref = length(x)

  if(verbose)
    message("Using pedigree ", ref, " as reference")

  # Number of attached markers for each pedigree
  nm = vapply(x, nMarkers, FUN.VALUE = 1)

  # Are any of the pedigrees without marker data?
  empty = nm == 0

  if(all(empty))
    stop2("None of the pedigrees has attached marker data")

  # For those that are empty, transfer from the first that has data
  if(any(empty)) {
    srcIdx = which(!empty)[1]
    if(verbose)
      message(sprintf("Transfering marker data from pedigree %d to empty pedigrees: ", srcIdx), which(empty))
    src = x[[srcIdx]]
    for(i in which(empty))
      x[[i]] = transferMarkers(from = src, to = x[[i]])
  }

  # Update marker counts
  nm = vapply(x, nMarkers, FUN.VALUE = 1)

  # Marker names
  if(missing(markers)) {
    if(!all(nm == nm[1]))
      stop2("The pedigrees have different number of markers: ", nm)
    markers = seq_len(nm[1])
  }
  markernames = name(x[[1]], markers)

  # Fix NA names
  if(any(NAnames <- is.na(markernames)))
    markernames[NAnames] = paste0("M", which(NAnames))

  # break all loops (NB: would use rapply, but doesnt work since is.list(ped) = TRUE
  breaklp = function(a) {
    if(is.pedList(a))
      return(lapply(a, breaklp))
    if(is.singleton(a))
      return(a)
    breakLoops(a, verbose = FALSE)
  }

  x_loopfree = lapply(x, breaklp)

  # compute likelihoods
  liks = lapply(x_loopfree, function(xx) vapply(markers, function(i)
      likelihood(xx, marker1 = i), FUN.VALUE = 1))
  likelihoodsPerMarker = do.call(cbind, liks)

  # LR per marker and total
  LRperMarker = do.call(cbind, lapply(1:length(x), function(j) liks[[j]]/liks[[ref]]))

  # ensure output has labels for the hypotheses
  hypnames = colnames(likelihoodsPerMarker)
  if(is.null(hypnames)) {
    hypnames = paste("Hyp", seq_along(x))
    colnames(likelihoodsPerMarker) = hypnames
  }
  colnames(LRperMarker) = hypnames

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
