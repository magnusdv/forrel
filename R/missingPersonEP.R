#' Exclusion power for missing person cases.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#' @inheritParams exclusionPower
#'
#' @return A `mpEP` object, which is essentially a list with the following entries:
#'
#'   * `EPperMarker`: A numeric vector containing the exclusion power of each
#'   marker. If the genotypes of a marker are incompatible with the `reference`
#'   pedigree, the corresponding entry is NA
#'
#'   * `EPtotal`: The total exclusion power, computed as `1 - prod(1 -
#'   EPperMarker, na.rm = TRUE)`
#'
#'   * `expectedMismatch`: The expected number of markers giving exclusion,
#'   computed as `sum(EPperMarker, na.rm = TRUE)`
#'
#'   * `distribMismatch`: The probability distribution of the number of markers
#'   giving exclusion. This is given as a numeric vector of length `n+1`, where
#'   `n` is the number of nonzero element of `EPperMarker`. The vector has names
#'   `0:n`
#'
#'   * `params`: A list containing the input parameters `missing`, `markers` and
#'   `disableMutations`
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
  st = Sys.time()

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")

  relatedPed = relabel(reference, old = missing, new = "_POI_")
  unrelatedPed = list(reference, singleton("_POI_"))

  ep = exclusionPower(claimPed = relatedPed, truePed = unrelatedPed, ids = "_POI_",
                     markers = markers, disableMutations = disableMutations, plot = FALSE,
                     verbose = verbose)

  ep
}

