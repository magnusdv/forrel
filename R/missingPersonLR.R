#' Likelihood ratio calculation for missing person identification
#'
#' This is a wrapper function for [kinshipLR()] for the special case of missing
#' person identification. A person of interest (POI) is matched against a
#' reference dataset containing genotypes of relatives of the missing person.
#'
#' @param reference A `ped` object with attached markers.
#' @param missing The ID label of the missing member of `reference`.
#' @param poi A `singleton` object, or NULL. If NULL, and `missing` is
#'   genotyped, this data is extracted and used as `poi`.
#' @param verbose A logical.
#' @param ... Optional parameters to be passed onto [kinshipLR()].
#'
#' @return The `LRresult` object returned by [kinshipLR()], but without the
#'   trivial `H2:H2` comparison.
#'
#' @examples
#' #------------------------------------------------
#' # Example: Identification of a missing grandchild
#' #------------------------------------------------
#'
#' # Reference pedigree with missing grandchild (MP)
#' x = linearPed(2) |> relabel(old = 5, new = "MP")
#'
#' # Database with 5 STR markers
#' db = NorwegianFrequencies[1:5]
#'
#' # Simulate reference data for grandmother and MP (5 STR markers)
#' x = profileSim(x, markers = db, ids = c(2, "MP"), seed = 2509)
#'
#' # Compute LR
#' lr = missingPersonLR(x, missing = "MP")
#' lr
#' lr$LRperMarker
#'
#' @export
missingPersonLR = function(reference, missing, poi = NULL, verbose = TRUE, ...) {

  if(is.pedList(reference))
    stop2("Argument `reference` must be a connected pedigree, not a list of pedigrees")

  if(!is.ped(reference))
    stop2("Argument `reference` must be a `ped` object")

  # Case I: poi data is attached to MP in the pedigree
  if(is.null(poi)) {

    if(!missing %in% typedMembers(reference))
      stop2(sprintf("Interpreting `%s` as POI, but this individual is not typed", missing))

    if(verbose)
      cat(sprintf("Interpreting `%s` as the person of interest\n", missing))

    # Good to lump early!
    reference = lumpAlleles(reference)

    # Hyp1
    if(verbose)
      cat(sprintf("\nForming H1 from reference:\n  * Renaming `%s` to `POI`\n", missing))
    H1 = relabel(reference, old = missing, new = "POI")

    # Hyp2
    if(verbose)
      cat(sprintf("\nForming H2 from reference:\n  * Extracting `%s` as singleton named `POI`\n", missing))
    poiSingleton = subset(reference, missing) |>
      relabel(old = missing, new = "POI")
    H2 = list(reference, poiSingleton)
  }
  else {
    if(!is.singleton(poi))
      stop2("Argument `poi` must be a singleton or NULL: ", poi)

    # Hypothesis 1: POI = MP
    if(verbose)
      cat(sprintf("\nForming H1 from reference:\n  * Transferring genotypes from POI to '%s'\n", missing))
    poi$ID = missing
    H1 = transferMarkers(from = poi, to = reference, erase = FALSE)

    # Hypothesis 2: POI unrelated to MP
    if(verbose)
      cat("\nForming H2 from reference:\n  * List of reference and poi singleton\n")
    poi$ID = "POI"
    H2 = list(reference, poi)
  }

  if(missing %in% typedMembers(H2)) {
    if(verbose)
      cat(sprintf("  * Removing genotypes from `%s`\n", missing))
    H2 = setAlleles(H2, missing, alleles = 0)
  }

  if(verbose) cat("\n")

  # Calculate LR
  lr = kinshipLR(H1, H2, ref = 2, verbose = verbose, ...)
  lr$LRtotal = lr$LRtotal[1]
  lr$LRperMarker = lr$LRperMarker[,1]
  lr
}

