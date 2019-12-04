#' Missing person power simulations
#'
#' Estimate the exclusion/inclusion power for various selections of available
#' individuals.
#' @param reference A connected `ped` object
#' @param missing The label of the missing person. By default "MP".
#' @param selections A list of pedigree member subsets. In the special case that
#'   all subsets consist of a single individual, `selections` can be given as a
#'   simple vector.
#' @param baseline A logical. If TRUE (default) an *empty* selection, named
#'   "Baseline", is added as the first element of `selection`.
#' @param nProfiles The number of profile simulations for each selection.
#' @param lrSims,thresholdIP Parameters passed onto [missingPersonIP()]
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
#' @param seed A numeric seed for the random number generator (optional)
#' @param verbose A logical
#'
#' @return A list with the same length as `selections`. Each entry has elements
#'   `ep` and `ip`, each of which is a list of length `nProfiles`.
#' @export
#'
#' @examples
#' x = nuclearPed(fa = "Gf", mo = "Gm", children = c("Uncle", "Mother"), sex = 1:2)
#' x = addChildren(x, fa = "Father", mo = "Mother", nch = 3, sex = c(1,2,1),
#'                 id = c("S1", "S2", "MP"))
#' x = addSon(x, "Father", id = "HS")
#'
#' # Brother S1 is already genotyped with a marker with 5 alleles
#' m = marker(x, S1 = 1:2, alleles = 1:4)
#' x = setMarkers(x, m)
#'
#' # Alternatives for additional genotyping
#' sel = list("Father", "S2", "HS", c("Gm", "Uncle"))
#'
#' plot(x, marker = 1, shaded = sel)
#'
#' # Simulate
#' simData = MPPsims(x, selections = sel, nProfiles = 2, lrSims = 2)
#'
#' # Power plot
#' powerPlot(simData, type = 3)
#'
#' ### With  mutations
#' # Create inconsistent marker
#' m2 = m
#' genotype(m2, "Father") = 3
#' x = setMarkers(x, list(m, m2))
#'
#' # Set mutation models for both
#' mutmod(x, 1:2) = list("equal", rate = 0.1)
#'
#' # By default mutations are disabled for consistent markers
#' MPPsims(x, selections = "Father", baseline = FALSE, seed = 123)
#'
#' # Don't disable anything
#' MPPsims(x, selections = "Father", baseline = FALSE, seed = 123,
#'         disableMutations = FALSE)
#'
#' \donttest{
#' # Disable all mutation models. SHOULD GIVE ERROR FOR SECOND MARKER
#' MPPsims(x, selections = "Father", baseline = FALSE, seed = 123,
#'         disableMutations = TRUE)
#' }
MPPsims = function(reference, missing = "MP", selections, baseline = TRUE,
                   nProfiles = 1, lrSims = 1, thresholdIP = NULL,
                   disableMutations = NA, seed = NULL, verbose = TRUE) {

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")

  nMark = nMarkers(reference)
  if(nMark == 0)
    stop2("No markers attached to `reference` pedigree")

  if(!is.list(selections))
    selections = as.list(selections)

  names(selections) = sapply(selections, toString)

  if(baseline)
    selections = c(list(Baseline = NULL), selections)

  ### Mutation disabling
  if(!is.null(disableMutations) && !isFALSE(disableMutations)) {

    # Which of the markers allow mutations?
    hasMut = allowsMutations(reference)

    if(isTRUE(disableMutations))
      disableMutations = which(hasMut)
    else if(identical(disableMutations, NA)) {
      # Disable mutations for consisten markers
      disableTF = logical(nMark)
      disableTF[hasMut] = consistentMarkers(reference, hasMut)
      disableMutations = which(disableTF)
    }
    else {# if numeric or character
      disableMutations = whichMarkers(reference, disableMutations)
    }

    # Disable
    if(length(disableMutations)) {
      if(verbose)
        message("Disabling mutations for markers: ", toString(disableMutations))
      mutmod(reference, disableMutations) = NULL
    }
  }
  # Wrappers (just for simplification)
  epfun = function(x) missingPersonEP(x, missing = missing, disableMutations = FALSE, verbose = FALSE)
  ipfun = function(x) missingPersonIP(x, missing = missing, disableMutations = FALSE, nsim = lrSims,
                                      threshold = thresholdIP, verbose = FALSE)

  set.seed(seed)

  powSims = lapply(names(selections), function(sel) {
    if(verbose)
      message("Selection: ", sel)

    ids = selections[[sel]]

    # Baseline
    if(is.null(ids)) {
      ep0 = epfun(reference)
      ip0 = ipfun(reference)
      return(list(ep = ep0, ip = ip0))
    }

    # Simulate profile for `ids`
    sims = profileSim(reference, ids = ids, N = nProfiles)

    # Compute updated EP and IP for each profile
    ep = lapply(sims, epfun)
    ip = lapply(sims, ipfun)
    list(ep = ep, ip = ip)
  })

  structure(powSims, nProfiles = nProfiles, lrSims = lrSims, thresholdIP = thresholdIP,
            names = names(selections), class = "MPPsim")
}

