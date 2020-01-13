#' Likelihood ratios for kinship testing
#'
#' This function computes likelihood ratios (LRs) for a given a list of
#' pedigrees with attached markers. The user must indicate which of the
#' pedigrees is the 'reference', which will be used in the denominator in each
#' LR.
#'
#' @param x A list of pedigree alternatives. Each alternative should be either a
#'   single `ped` object or a list of such.
#' @param ref A single integer, indicating the index (in the list `x`) of the
#'   reference alternative. This is used in the denominator of each LR.
#' @param markers A vector of integers, indexing which markers should be
#'   included. If NULL (the default) all markers are used.
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
#'   * `time` user system and elapsed time
#'
#' @author Magnus Dehli Vigeland and Thore Egeland
#'
#' @examples
#'
#' # Simulate 5 markers for a pair of full sibs
#' set.seed(123)
#' sibs = nuclearPed(children = c("A", "B"))
#' sibs = simpleSim(sibs, N = 5, alleles = 1:4, ids = c("A", "B"))
#'
#' # Create two alternative hypotheses and transfer the simulated genotypes to them
#' halfsibs = relabel(halfSibPed(), old = 4:5, new = c("A", "B"))
#' halfsibs = transferMarkers(sibs, halfsibs)
#'
#' unrel = list(singleton("A"), singleton("B"))
#' unrel = transferMarkers(sibs, unrel)
#'
#' # Compute LR with 'unrelated' as reference
#' res = kinshipLR(list(sibs, halfsibs, unrel), ref = 3)
#' res
#'
#' # Detailed results
#' res$LRperMarker
#' res$likelihoodsPerMarker
#' @export
kinshipLR = function(x, ref, markers) {
  st = proc.time()

  # get marker names
  y = if(is.ped(x[[1]])) x[[1]] else x[[1]][[1]]
  if(missing(markers))
    markers = seq_len(nMarkers(y))
  markernames = sapply(y$MARKERS[markers], attr, "name")
  NAnames = is.na(markernames)
  if(any(NAnames)) markernames[NAnames] = paste0("M", which(NAnames))

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
    time = proc.time()-st),
    class = "LRresult")
}

# Alias
LR = kinshipLR


#' @export
print.LRresult = function(x, ...) {
  cat("Total LR:\n")
  print(x$LRtotal)
}
