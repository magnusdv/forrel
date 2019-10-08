#' Likelihood ratios of pedigree hypotheses
#'
#' This function computes likelihood ratios for a given a list of pedigrees
#' (`ped` and/or `singleton` objects), with genotype data from the same set of
#' markers. Data exported from the 'Familias' software can be analysed by using
#' [Familias2ped()] prior to calling this function.
#'
#' @param x A list of pedigree alternatives. Each alternative should be either a
#'   single `ped` object or a list of such.
#' @param ref A single integer, indicating the index (in the list `x`) of the
#'   reference alternative. This is used in the denominator of each LR.
#' @param markers A vector of integers, indexing which markers should be
#'   included. If NULL (the default) all markers are used.
#'
#' @return A list with entries
#'
#'   * `LR` : Likelihood ratios
#'
#'   * `LRperMarker` : Likelihood ratios for each marker
#'
#'   * `likelihoodsPerSystem` : Likelihoods for each marker
#'
#'   * `time` user system and elapsed time
#'
#' @author Magnus Dehli Vigeland and Thore Egeland
#'
#' @examples
#'
#' # Simulate genotypes for 5 tetra-allelic markers for a pair of full sibs
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
#' LR(list(sibs, halfsibs, unrel), ref=3)
#'
#' @export
LR = function(x, ref, markers) {
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
    breakLoops(a, verbose=F)
  }

  x_loopfree = lapply(x, breaklp)

  # compute likelihoods
  liks = lapply(x_loopfree, function(xx) vapply(markers, function(i)
      likelihood(xx, marker1 = i), FUN.VALUE = 1))
  likelihoodsPerSystem = do.call(cbind, liks)

  # LR per marker and total
  LRperMarker = do.call(cbind, lapply(1:length(x), function(j) liks[[j]]/liks[[ref]]))

  # total LR
  LR = apply(LRperMarker, 2, prod)

  # output
  rownames(likelihoodsPerSystem) = rownames(LRperMarker) = markernames
  list(LR = LR,
       LRperMarker = LRperMarker,
       likelihoodsPerSystem = likelihoodsPerSystem,
       time = proc.time()-st)
}

