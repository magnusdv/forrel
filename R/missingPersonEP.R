#' Exclusion power statistics in missing person cases.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing pedigree member.
#' @param markers Names or indices of the markers to be included. By default,
#'   all markers.
#' @param disableMutations This parameter determines how mutation models are
#'   treated. If `TRUE`, mutations are disabled for all markers. If `NA` (the
#'   default), mutations are disabled only for those markers with nonzero
#'   baseline likelihood. (In other words: Mutations are NOT disabled if the
#'   reference genotypes are inconsistent with the pedigree.) If `FALSE` no
#'   action is done to disable mutations. Finally, if a vector of marker names
#'   or indices is given, mutations are disabled for these markers exclusively.
#' @param verbose A logical, by default TRUE.
#'
#' @return A `mpEP` object, which is essentially a list with the following entries:
#'
#'   * `EPperMarker`: A numeric vector containing the exclusion power of each
#'   marker. If the genotypes of a marker are incompatible with the `reference`
#'   pedigree, the corresponding entry is NA
#'
#'   * `EPtotal`: The total exclusion power, computed as `1 - prod(1 -
#'   EPperMarker, na.rm = T)`
#'
#'   * `expectedMismatch`: The expected number of markers giving exclusion,
#'   computed as `sum(EPperMarker, na.rm = T)`
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
missingPersonEP = function(reference, missing, markers, disableMutations = NA, verbose = TRUE) {
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
  if(length(disable) > 0)
    mutmod(reference, disable) = NULL

  # Extract markers and set up pedigrees
  ref = selectMarkers(reference, markers)
  ped_related = relabel(ref, old = missing, new = "_POI_")
  ped_unrelated = list(ref, singleton("_POI_"))
  mseq = seq_along(markers)

  ######################
  ### Setup finished ###
  ######################

  # Baseline likelihoods
  refliks = vapply(mseq, function(i) pedprobr::likelihood(ref, i), 0)

  # Compute the exclusion power of each marker
  ep = vapply(mseq, function(i) {

    if(verbose)
      message("Marker ", markers[i], " ... ", appendLF = F)

    # If impossible, return NA
    if(refliks[i] == 0) {
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
  if(n.nonz > 0) {
    if (requireNamespace("poibin", quietly = TRUE))
      distrib[] = poibin::dpoibin(kk = 0:n.nonz, pp = ep[!is.na(ep) & ep > 0])
    else {
      warning("Package `poibin` not found. Cannot compute the distribution of exclusion counts without this; returning NA's")
      distrib[] = NA_real_
    }
  }
  else
    distrib[] = 1

  # Timing
  if(verbose)
    message("\nTotal time used: ", format(Sys.time() - st, digits = 3))

  # List of input parameters
  params = list(missing = missing, markers = markers, disableMutations = disableMutations)

  structure(list(EPperMarker = ep, EPtotal = tot, expectedMismatch = expMis,
       distribMismatch = distrib, params = params), class = "mpEP")
}

print.mpEP = function(x, ...) {
  cat("\n")
  cat("Total EP:", round(x$EPtotal, 3), "\n")
  cat("Markers with potential mismatch:", sum(x$EPperMarker > 0), "\n")
  cat("Expected number of mismatches:", round(x$expectedMismatch, 3), "\n")
}

