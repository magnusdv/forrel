#' Exclusion power for missing person cases
#'
#' This is a special case of [exclusionPower()] for use in missing person cases.
#' The function computes the probability that a random person is genetically
#' incompatible with the typed relatives of the missing person.
#'
#' This function is identical to [randomPersonEP()], but with different argument
#' names. This makes it consistent with [missingPersonIP()] and the other
#' missing person functions.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#'
#' @inheritParams randomPersonEP
#' @seealso [randomPersonEP()], [exclusionPower()]
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
#' # (by default this is ignored by `missingPersonEP()`)
#' badMarker = marker(x, `3` = 1, `4` = 2, `5` = 3)
#' x = addMarkers(x, badMarker)
#'
#' # Compute exclusion power statistics
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

