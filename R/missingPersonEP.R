#' Exclusion power statistics in missing person cases
#'
#' @param reference A [pedtools::ped()] object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#' @param markers Names or indices of the markers to be included. By default, all markers are included.
#' @param disableMutations A logical (default: TRUE).
#' @param verbose A logical.
#'
#' @examples
#'
#' # 3 siblings; third is missing
#' x = nuclearPed(3)
#'
#' # First two brothers typed with 10 triallelic markers
#' x = markerSim(x, N = 5, ids = 3:4, alleles = 1:3, seed = 3, verbose = FALSE)
#' x
#'
#' # Compute exclusion power statistics
#' missingPersonEP(x, missing = 5)
#'
#' @importFrom poibin dpoibin
#' @importFrom pedprobr likelihood
#' @export
missingPersonEP = function(reference, missing, markers, disableMutations = TRUE, verbose = TRUE) {

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")

  if(missing(markers))
    markers = seq_len(nMarkers(reference))

  # Remove mutation models if indicated
  if(disableMutations)
    mutmod(reference, markers) = NULL

  useMarkers = markers

  # Identify and remove impossible markers
  liks = sapply(markers, function(i) pedprobr::likelihood(reference, i))
  imp = liks == 0
  useMarkers = markers[!imp]
  if(verbose && any(imp)) {
    message("Removing markers with mutations in reference: ", toString(which(imp)))
  }

  # Extract markers and set up pedigrees
  ref = selectMarkers(reference, useMarkers)
  ped_related = relabel(ref, old = missing, new = "_POI_")
  ped_unrelated = list(ref, singleton("_POI_"))

  # Compute the exclusion power of each marker
  ep = vapply(seq_along(useMarkers), function(i) {

    if(verbose) message("Marker ", useMarkers[i], " ... ", appendLF = F)

    this.ep = exclusionPower(ped_claim = ped_related, ped_true = ped_unrelated,
                             ids = "_POI_", markerindex = i, plot = F,
                             verbose = F)
    if(verbose) message("PE = ", this.ep)

    this.ep
  }, FUN.VALUE = 0)

  # Result: EP per marker
  perMarker = numeric(length(markers))
  perMarker[!imp] = ep
  perMarker[imp] = NA
  names(perMarker) = markers

  # Result: Total EP
  tot = 1 - prod(1 - ep)

  # Result: Expected number of exclusions
  meanMis = sum(ep)

  # Result: Distribution of number of mismatches
  # This is a sum of different Bernoulli variables, i.e., Poisson binomial.
  n.nonz = sum(ep > 0)
  distrib = structure(numeric(n.nonz + 1), names = 0:n.nonz)
  if(n.nonz > 0)
    distrib[] = poibin::dpoibin(kk = 0:n.nonz, pp = ep[ep > 0])
  else
    distrib[] = 1

  list(PEperMarker = perMarker, PEtotal = tot,
       meanMismatches = meanMis, mismatchDistrib = distrib)
}
