#' Inclusion power statistics in missing person cases.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#' @param markers Names or indices of the markers to be included. By default,
#'   all markers.
#' @param nsim A positive integer: the number of simulations
#' @param threshold A numeric vector with one or more positive numbers used
#'   as the likelihood ratio tresholds for inclusion
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
#' @importFrom poibin dpoibin
#' @importFrom pedprobr likelihood
#' @export
missingPersonIP = function(reference, missing, markers, nsim = 1, threshold = NULL,
                           seed = NULL, verbose = T) {
  st = proc.time()

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
    message("Simulating ", nsim, " profiles...", appendLF = F)

  allsims = profileSim(ped_related, ids = "_POI_", N = nsim, conditions = midx)

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

  # Return
  list(LRperSim = LRperSim, ELRperMarker = ELRperMarker, ELRtotal = ELRtotal, IP = IP,
       time = proc.time() - st)
}
