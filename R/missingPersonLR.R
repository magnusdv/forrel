#' Likelihood ratio calculation for missing person identification
#'
#' This is a wrapper function for [kinshipLR()] for the special case of missing
#' person identification. A person of interest (POI) is matched against a
#' reference dataset containing genotypes of relatives of the missing person.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing member of `reference`.
#' @param poi A `singleton` object with genotypes for the person of interest.
#'   The marker names of `poi` must coincide with those of `reference`. The ID
#'   label of `poi` is irrelevant.
#'
#' @return The `LRresult` object returned by [kinshipLR()], but without the
#'   trivial `H2:H2` comparison.
#'
#' @examples
#' #------------------------------------------------
#' # Example: Identification of a missing grandchild
#' #------------------------------------------------
#'
#' set.seed(2509)
#'
#' ### Reference pedigree with missing grandchild (MP)
#' x = relabel(linearPed(2), old = 5, new = "MP")
#'
#' # Simulate reference data for grandmother (5 STR markers)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)[[1]]
#'
#' ### Person of interest 1: Unrelated
#' poi1 = singleton("poi1")
#'
#' # Transfer (empty) markers and simulate genotypes
#' poi1 = transferMarkers(from = x, to = poi1)
#' poi1 = profileSim(poi1, N = 1)[[1]]
#'
#' # Compute LR
#' lr1 = missingPersonLR(x, missing = "MP", poi = poi1)
#' lr1
#' lr1$LRperMarker
#'
#'
#' ### Person of interest 2: The true MP
#'
#' # Simulate MP conditional on reference, and extract as singleton
#' poi2 = profileSim(x, N = 1, ids = "MP")[[1]]
#'
#' # Extract MP as singleton
#' poi2 = subset(poi2, "MP")
#'
#' # Compute LR
#' lr2 = missingPersonLR(x, missing = "MP", poi = poi2)
#' lr2
#' lr2$LRperMarker
#'
#' @export
missingPersonLR = function(reference, missing, poi) {

  if(is.pedList(reference))
    stop2("Argument `reference` must be a connected pedigree, not a list of pedigrees")

  if(!is.ped(reference))
    stop2("Argument `reference` must be a `ped` object")

  if(!is.singleton(poi))
    stop2("Argument `poi` must be a singleton")

  poiName = labels(poi)

  # Hypothesis 1: POI = MP
  H1 = relabel(reference, old = missing, new = poiName)
  H1 = transferMarkers(from = poi, to = H1, erase = FALSE)

  # Hypothesis 2: POI unrelated to MP
  H2 = list(reference, poi)

  # Calculate LR
  lr = kinshipLR(H1, H2, ref = 2, verbose = FALSE)
  lr$LRtotal = lr$LRtotal[1]
  lr$LRperMarker = lr$LRperMarker[,1]
  lr
}

