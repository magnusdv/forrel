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
#' plot(list(x, poi), marker = "D21S11", hatched = typedMembers)
#'
#' @export
findExclusions = function(x, id, candidate, removeMut = TRUE) {

  # Insert candidate at position without erasing existing data
  y = transferMarkers(from = candidate, to = x, erase = FALSE,
                      idsFrom = labels(candidate), idsTo = id)

  inconsistentMarkers(y, names = TRUE, removeMut = removeMut)
}


# Test if genotypes are consistent with ped
# (A better, but slower, alternative to `mendelianCheck()`)
consistentMarkers = function(x, markers = NULL, names = FALSE, removeMut = TRUE) {
  if(!is.null(markers))
    x = selectMarkers(x, markers)

  if(removeMut)
    x = setMutmod(x, model = NULL) # works also with 0 markers

  # Log-likelihood of each marker
  logliks = pedprobr::likelihood(x, logbase = exp(1))
  cons = logliks > -Inf

  # Return logical vector or names of consistent markers
  if(names) name(x)[cons] else cons
}

# Test if genotypes are consistent with ped
# (A better, but slower, alternative to `mendelianCheck()`)
inconsistentMarkers = function(x, markers = NULL, names = FALSE, removeMut = TRUE) {
  if(!is.null(markers))
    x = selectMarkers(x, markers)

  if(removeMut)
    x = setMutmod(x, model = NULL) # works also with 0 markers

  # Log-likelihood of each marker
  logliks = pedprobr::likelihood(x, logbase = exp(1))
  incons = logliks == -Inf

  # Return logical vector or names of consistent markers
  if(names) name(x)[incons] else incons
}
