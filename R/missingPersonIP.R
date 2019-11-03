#' Inclusion power statistics in missing person cases.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#' @param markers Names or indices of the markers to be included. By default,
#'   all markers.
#' @param nsim A positive integer: the number of simulations
#' @param threshold A numeric vector with one or more positive numbers used
#'   as the likelihood ratio tresholds for inclusion
#' @param disableMutations This parameter determines how mutation models are
#'   treated. If `TRUE`, mutations are disabled for all markers. If `NA` (the
#'   default), mutations are disabled only for those markers with nonzero
#'   baseline likelihood. (In other words: Mutations are NOT disabled if the
#'   reference genotypes are inconsistent with the pedigree.) If `FALSE` no
#'   action is done to disable mutations. Finally, if a vector of marker names
#'   or indices is given, mutations are disabled for these markers exclusively.
#' @param seed A numeric seed for the random number generator (optional)
#' @param verbose A logical, by default TRUE.
#'
#' @return A list with the following entries:
#'
#'   * `LRperSim`: A numeric vector of length `nsim` containing the total LR for
#'   each simulation.
#'
#'   * `ELRperMarker`: The mean LR per marker, over all simulations.
#'
#'   * `ELRtotal`: The mean total LR over all simulations.
#'
#'   * `IP`: A numeric of the same length as `threshold`. For each element of
#'   `threshold`, the fraction of simulations resulting in a LR exceeding the
#'   given number.
#'
#'   * `params`: A list containing the input parameters `missing`, `markers`, `nsim`, `threshold` and
#'   `disableMutations`
#'
#' @examples
#'
#' # Four siblings; the fourth is missing
#' x = nuclearPed(4)
#'
#' # Remaining sibs typed with 5 triallelic markers
#' x = markerSim(x, N = 5, ids = 3:5, alleles = 1:3, seed = 123, verbose = FALSE)
#'
#' # Compute exclusion power statistics
#' missingPersonIP(x, missing = 6, nsim = 5, threshold = c(10, 100))
#'
#' # Compare with genotypes
#' x
#'
#' @importFrom pedprobr likelihood
#' @export
missingPersonIP = function(reference, missing, markers, nsim = 1, threshold = NULL,
                           disableMutations = NA, seed = NULL, verbose = T) {
  st = Sys.time()

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")

  nmark = nMarkers(reference)
  if(nmark == 0)
    stop2("No markers attached to the input reference")

  if(missing(markers)) {
    if(verbose)
      message("Using all ", nmark, " attached markers")
    markers = name(reference, 1:nmark)
    if(anyNA(markers))
      markers = 1:nmark
  }

  # Do any of the markers model mutatinos?
  hasMut = sapply(getMarkers(reference, markers), allowsMutations)

  # For which marker should mutations be disabled?
  disALL = any(hasMut) && isTRUE(disableMutations)
  disGOOD = any(hasMut) && length(disableMutations) == 1 && is.na(disableMutations)
  disSELECT = any(hasMut) && !isFALSE(disableMutations) && !is.null(disableMutations)

  if(disALL)
    disable = markers[hasMut]
  else if(disGOOD) { # disable only if consistent
    refNoMut = reference
    mutmod(refNoMut, markers[hasMut]) = NULL
    liksNoMut = vapply(markers[hasMut], function(i) pedprobr::likelihood(refNoMut, i), 0)
    disable = markers[hasMut][liksNoMut > 0]
  }
  else if(disSELECT)
    disable = whichMarkers(reference, disableMutations)
  else
    disable = NULL

  # Disable mutations in the chosen cases
  if(length(disable) > 0) {
    if(verbose) message("Disabling mutations for marker ", toString(disable))
    mutmod(reference, disable) = NULL
  }

  # Extract markers and set up pedigrees
  midx = whichMarkers(reference, markers)
  ped_related = relabel(reference, old = missing, new = "_POI_")
  ped_unrelated = list(reference, singleton("_POI_"))

  # Raise error if impossible markers
  liks = sapply(midx, function(i) pedprobr::likelihood(reference, i))
  if(any(liks == 0))
    stop2("Marker incompatible with reference pedigree: ", markers[liks == 0],
          "\nThis makes conditional simulations impossible. Exclude the marker from the computation or add a mutation model")

  # Set seed once
  set.seed(seed)

  # Simulate nsim complete profiles of ped_related
  if(verbose)
    message("\nSimulating ", nsim, " profiles...", appendLF = F)

  allsims = profileSim(ped_related, ids = "_POI_", N = nsim, markers = midx)

  if(verbose)
    message("done\nComputing likelihood ratios...", appendLF = F)

  # Compute the exclusion power of each marker
  lrs = vapply(1:nsim, function(i) {
    rel.sim = allsims[[i]]
    unrel.sim = transferMarkers(from = rel.sim, to = ped_unrelated)

    lr = LR(list(rel.sim, unrel.sim), ref = 2)
    lr$LRperMarker[,1]
  }, FUN.VALUE = numeric(length(markers)))

  # Ensure matrix
  if(length(markers) == 1)
    lrs = rbind(lrs, deparse.level = 0)

  rownames(lrs) = markers

  if(verbose)
    message("done")

  # Results
  ELRperMarker = apply(lrs, 1, mean)
  LRperSim = apply(lrs, 2, prod)
  ELRtotal = mean(LRperSim)
  IP = sapply(threshold, function(thr) mean(LRperSim >= thr))
  names(IP) = threshold

  # Timing
  if(verbose)
    message("\nTotal time used: ", format(Sys.time() - st, digits = 3))

  # Lits of input parameters
  params = list(missing = missing, markers = markers,
                nsim = nsim, threshold = threshold,
                disableMutations = disableMutations)

  structure(list(LRperSim = LRperSim, ELRperMarker = ELRperMarker,
       ELRtotal = ELRtotal, IP = IP, params = params), class = "mpIP")
}

print.mpIP = function(x, ...) {
  cat("\n")
  cat("Total ELR:", round(x$ELRtotal, 3), "\n")
  cat("Estimated inclusion powers:\n")
  for(i in seq_along(x$IP))
    cat(sprintf("  P(LR > %s) = %.2g\n", names(x$IP)[i], x$IP[i]))
}


