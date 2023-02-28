#' Find markers excluding an identification
#'
#' Find markers for which the genotypes of a candidate individual is incompatible with a pedigree
#'
#' @param x A `ped` object or a list of such.
#' @param id A character of length 1; the name of an untyped member of `x`.
#' @param candidate A singleton pedigree, with genotypes for the same markers as `x`.
#' @param removeMut A logical. If TRUE (default), all mutations models are stripped.
#' @return A character vector containing the names of incompatible markers.
#'
#' @examples
#'
#' db = NorwegianFrequencies[1:5]
#'
#' # Pedigree with 3 siblings; simulate data for first two
#' x = nuclearPed(3) |>
#'   profileSim(ids = 3:4, markers = db, seed = 1)
#'
#' # Simulate random person
#' poi = singleton("POI") |>
#'   profileSim(markers = db, seed = 1)
#'
#' # Identify incompatible markers
#' findExclusions(x, id = 5, candidate = poi)   # D21S11
#'
#' # Inspect
#' plotPedList(list(x, poi), marker = "D21S11", hatched = typedMembers)
#'
#' @importFrom pedprobr likelihood
#' @export
findExclusions = function(x, id, candidate, removeMut = TRUE) {

  # Insert candidate at position without erasing existing data
  y = transferMarkers(from = candidate, to = x, erase = FALSE,
                      idsFrom = labels(candidate), idsTo = id)

  if(removeMut)
    y = setMutationModel(y, NULL)

  # Vector of marker names
  mvec = name(y)

  # Likelihood of each marker
  lik = likelihood(y, markers = mvec)

  # Return names of markers with zero likelihood
  mvec[lik == 0]
}
