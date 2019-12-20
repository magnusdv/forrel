#' Missing person power simulations
#'
#' Estimate the exclusion/inclusion power for various selections of available
#' individuals.
#'
#' @inheritParams missingPersonIP
#' @param reference A connected `ped` object, or a list of pedigrees. In the
#'   latter case, the list must have the same length as `selections`.
#' @param selections A list of pedigree member subsets. In the special case that
#'   all subsets consist of a single individual, `selections` can be given as a
#'   simple vector.
#' @param addBaseline A logical. If TRUE (default) an *empty* selection, named
#'   "Baseline", is added as the first element of `selection`.
#' @param nProfiles The number of profile simulations for each selection.
#' @param lrSims,thresholdIP Parameters passed onto [missingPersonIP()]
#'
#' @return An object of class "MPPsim", which is basically a list with one entry
#'   for each element of `selections`. Each entry has elements `ep` and `ip`,
#'   each of which is a list of length `nProfiles`.
#'
#'   The output object has various attributes reflecting the input. Note that
#'   `reference` and `selection` may differ slightly from the original input, since they
#'   may be modified during the function run. (For instance, a "Baseline" entry
#'   is added to `selection` if `addBaseline` is TRUE.) The crucial point is
#'   that the output attributes correspond exactly to the output data.
#'
#'   * `reference` (always a list, of the same length as the `selections` attribute
#'
#'   * `selections`
#'
#'   * `nProfiles`,`lrSims`,`thresholdIP`,`seed` (as in the input)
#'
#'   * `totalTime` (the total time used)
#'
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
#' MPPsims(x, selections = "Father", addBaseline = FALSE, seed = 123)
#'
#' # Don't disable anything
#' MPPsims(x, selections = "Father", addBaseline = FALSE, seed = 123,
#'         disableMutations = FALSE)
#'
#' \donttest{
#' # Disable all mutation models. SHOULD GIVE ERROR FOR SECOND MARKER
#' MPPsims(x, selections = "Father", addBaseline = FALSE, seed = 123,
#'         disableMutations = TRUE)
#' }
#'
#' # Effect of variable number of alleles
#' y = nuclearPed(father = "fa", child = "MP")
#' mlist = lapply(2:5, function(k) marker(y, alleles = 1:k))
#' y = setMarkers(y, mlist)
#' peds = lapply(1:nMarkers(y), function(i) selectMarkers(y, i))
#' sel = rep("fa", 4)
#' names(sel) = paste(2:5, "alleles")
#' pows = MPPsims(peds, selections = sel, addBaseline = FALSE, lrSims = 10)
#' powerPlot(pows, type = 3)
MPPsims = function(reference, missing = "MP", selections, addBaseline = TRUE,
                   nProfiles = 1, lrSims = 1, thresholdIP = NULL,
                   disableMutations = NA, seed = NULL, verbose = TRUE) {
  st = Sys.time()

  if(!is.list(selections))
    selections = as.list(selections)

  if(is.null(names(selections)))
    names(selections) = sapply(selections, toString)
  if(anyDuplicated(names(selections)))
    stop2("`selections` cannot have duplicated names")

  if(addBaseline)
    selections = c(list(Baseline = NULL), selections)

  if(is.ped(reference)) {
    reference = disableMutationModels(reference, disableMutations, verbose = verbose)
    reference = rep(list(reference), length(selections))
  }
  else {
    if(!is.pedList(reference))
      stop2("`reference` must be a single connected pedigree, or a list of pedigrees")
    else if(length(reference) != length(selections))
      stop2("Incompatible lengths of `reference` and `selections`")
    reference = lapply(reference, disableMutationModels, disable = disableMutations, verbose = FALSE)
  }

  # Wrappers (just for simplification)
  epfun = function(x) missingPersonEP(x, missing = missing, disableMutations = FALSE, verbose = FALSE)
  ipfun = function(x) missingPersonIP(x, missing = missing, disableMutations = FALSE, nsim = lrSims,
                                      threshold = thresholdIP, verbose = FALSE)

  set.seed(seed)

  powSims = lapply(seq_along(selections), function(i) {
    ref = reference[[i]]
    ids = selections[[i]]

    if(verbose)
      message("Selection: ", names(selections)[i])

    # Baseline
    if(is.null(ids)) {
      ep0 = epfun(ref)
      ip0 = ipfun(ref)
      return(list(ep = ep0, ip = ip0))
    }

    # Simulate profile for `ids`
    sims = profileSim(ref, ids = ids, N = nProfiles)

    # Compute updated EP and IP for each profile
    ep = lapply(sims, epfun)
    ip = lapply(sims, ipfun)
    list(ep = ep, ip = ip)
  })

  names(powSims) = names(selections)

  # Timing
  totalTime = format(Sys.time() - st, digits = 3)
  if(verbose)
    message("Total time used: ", totalTime)

  # Return with attributes
  structure(powSims, reference = reference, selections = selections,
            nProfiles = nProfiles, lrSims = lrSims, thresholdIP = thresholdIP,
            seed = seed, totalTime = totalTime, class = "MPPsim")
}

