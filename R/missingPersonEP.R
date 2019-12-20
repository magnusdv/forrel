#' Exclusion power for missing person cases.
#'
#' This is wrapper of [exclusionPower()] for the special case of a reference
#' family with a single missing member. Some reference members should already be
#' genotyped. The function computes the power to exclude an unrelated
#' individual, i.e. the probability of observing (in a truly unrelated
#' individual) a genotype incompatible with the reference.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#' @param disableMutations This parameter determines how mutation models are
#'   treated. Possible values are as follows:
#'
#'   * `NA` (the default): Mutations are disabled only for those markers whose
#'   known genotypes are consistent with `reference`. This is determined by
#'   temporarily removing all mutation models and checking which markers have
#'   nonzero likelihood.
#'
#'   * `TRUE`: Mutations are disabled for all markers. This will result in an
#'   error if any markers are inconsistent with `reference`.
#'
#'   * `FALSE`: No action is done to disable mutations.
#'
#'   * A vector containing the names or indices of those markers for which
#'   mutations should be disabled.
#' @inheritParams exclusionPower
#'
#' @return The `EPresult` object returned by [exclusionPower()].
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
#' # (this should be ignored by `missingPersonEP()`)
#' badMarker = marker(x, `3` = 1, `4` = 2, `5` = 3)
#' x = addMarkers(x, badMarker)
#'
#' # Compute exclusion power statistics
#' missingPersonEP(x, missing = 6)
#'
#' # With marker names:
#' name(x, 1:5) = paste0("M", 1:5)
#' missingPersonEP(x, missing = 6)
#'
#' @importFrom pedprobr likelihood
#' @export
missingPersonEP = function(reference, missing, markers = NULL, disableMutations = NA, verbose = TRUE) {

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")

  poiLabel = "_POI_"

  relatedPed = relabel(reference, old = missing, new = poiLabel)
  unrelatedPed = list(reference, singleton(poiLabel, sex = getSex(reference, missing)))

  ep = exclusionPower(claimPed = relatedPed, truePed = unrelatedPed, ids = poiLabel,
                      markers = markers, source = "claim", disableMutations = disableMutations,
                      plot = FALSE, verbose = verbose)

  # Change the `ids` entry from "_POI_" to `missing`
  ep$params$ids = missing

  ep
}

