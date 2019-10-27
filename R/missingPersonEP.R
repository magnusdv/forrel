#' Exclusion power statistics in missing person cases.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#' @param markers Names or indices of the markers to be included. By default,
#'   all markers.
#' @param disableMutations A logical, by default TRUE.
#' @param verbose A logical, by default TRUE.
#'
#' @return A list of four entries:
#'
#'   * `EPperMarker`: A numeric vector containing the exclusion power of each
#'   marker. If the genotypes of a marker are incompatible with the `reference`
#'   pedigree, the corresponding entry is NA.
#'
#'   * `EPtotal`: The total exclusion power, computed as `1 - prod(1 -
#'   EPperMarker, na.rm = T)`
#'
#'   * `expectedMismatch`: The expected number of markers giving exclusion,
#'   computed as `sum(EPperMarker, na.rm = T)`
#'
#'   * `distribMismatch`: The probability distribution of the number of markers
#'   giving exclusion. This is given as a numeric vector of length `n+1`, where
#'   `n` is the number of nonzero element of `EPperMarker`. The vector has
#'   names `0:n`.
#'
#' @examples
#'
#' # Four siblings; the fourth is missing
#' x = nuclearPed(4)
#'
#' # Remaining sibs typed with 4 triallelic markers
#' x = markerSim(x, N = 4, ids = 3:5, alleles = 1:3, seed = 577, verbose = FALSE)
#'
#' # Add marker with inconsistency in reference genotypes
#' # (this should be ignored)
#' badMarker = marker(x, `3` = 1, `4` = 2, `5` = 3)
#' x = addMarkers(x, badMarker)
#'
#' # Compute exclusion power statistics
#' missingPersonEP(x, missing = 6)
#'
#' # Compare with genotypes
#' x
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

  # Extract markers and set up pedigrees
  ref = selectMarkers(reference, markers)
  ped_related = relabel(ref, old = missing, new = "_POI_")
  ped_unrelated = list(ref, singleton("_POI_"))

  # Compute the exclusion power of each marker
  ep = vapply(seq_along(markers), function(i) {

    if(verbose)
      message("Marker ", markers[i], " ... ", appendLF = F)

    # If impossible, return NA
    lik = pedprobr::likelihood(reference, i)
    if(lik == 0) {
      message("Genotypes incompatible with reference pedigree - ignoring marker")
      return(NA_real_)
    }

    # Otherwise, compute EP
    this.ep = exclusionPower(ped_claim = ped_related, ped_true = ped_unrelated,
                             ids = "_POI_", markerindex = i, plot = F,
                             verbose = F)
    if(verbose) message("PE = ", this.ep)
    this.ep
  }, FUN.VALUE = 0)

  names(ep) = markers

  # Total EP
  tot = 1 - prod(1 - ep, na.rm = T)

  # Result: Expected number of exclusions
  expMis = sum(ep, na.rm = T)

  # Result: Distribution of number of mismatches
  # This is a sum of different Bernoulli variables, i.e., Poisson binomial.
  n.nonz = sum(ep > 0, na.rm = T)
  distrib = structure(numeric(n.nonz + 1), names = 0:n.nonz)
  if(n.nonz > 0)
    distrib[] = poibin::dpoibin(kk = 0:n.nonz, pp = ep[!is.na(ep) & ep > 0])
  else
    distrib[] = 1

  list(EPperMarker = ep, EPtotal = tot, expectedMismatch = expMis, distribMismatch = distrib)
}
