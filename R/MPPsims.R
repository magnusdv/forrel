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
#' @param ep A logical: Estimate the exclusion power? (Default: TRUE)
#' @param ip A logical: Estimate the inclusion power? (Default: TRUE)
#' @param addBaseline A logical. If TRUE (default) an *empty* selection, named
#'   "Baseline", is added as the first element of `selection`.
#' @param nProfiles The number of profile simulations for each selection.
#' @param lrSims,thresholdIP Parameters passed onto [missingPersonIP()].
#' @param numCores The number of cores used for parallelisation, by default 1.
#'
#' @return An object of class "MPPsim", which is basically a list with one entry
#'   for each element of `selections`. Each entry has elements `ep` and `ip`,
#'   each of which is a list of length `nProfiles`.
#'
#'   The output object has various attributes reflecting the input. Note that
#'   `reference` and `selection` may differ slightly from the original input,
#'   since they may be modified during the function run. (For instance, a
#'   "Baseline" entry is added to `selection` if `addBaseline` is TRUE.) The
#'   crucial point is that the output attributes correspond exactly to the
#'   output data.
#'
#'   * `reference` (always a list, of the same length as the `selections`
#'   attribute
#'
#'   * `selections`
#'
#'   * `nProfiles`,`lrSims`,`thresholdIP`,`seed` (as in the input)
#'
#'   * `totalTime` (the total time used)
#'
#' @examples
#' x = nuclearPed(fa = "Gf", mo = "Gm", children = c("Uncle", "Mother"), sex = 1:2)
#' x = addChildren(x, fa = "Father", mo = "Mother", nch = 3, sex = c(1,2,1),
#'                 id = c("S1", "S2", "MP"))
#' x = addSon(x, "Father", id = "HS")
#'
#' # Brother S1 is already genotyped with a marker with 4 alleles
#' x = addMarker(x, S1 = "1/2", alleles = 1:4)
#'
#' # Alternatives for additional genotyping
#' sel = list("Father", "S2", "HS", c("Gm", "Uncle"))
#'
#' plot(x, marker = 1, hatched = sel)
#'
#' # Simulate
#' simData = MPPsims(x, selections = sel, nProfiles = 2, lrSims = 2)
#'
#' # Power plot
#' powerPlot(simData, type = 3)
#'
#' \donttest{
#' ### With  mutations
#' # Add inconsistent marker
#' x = addMarker(x, S1 = "1/2", Father = "3/3", alleles = 1:4)
#'
#' # Set mutation models for both
#' mutmod(x, 1:2) = list("equal", rate = 0.1)
#'
#' # By default mutations are disabled for consistent markers
#' MPPsims(x, selections = "Father", addBaseline = FALSE)
#'
#' # Don't disable anything
#' MPPsims(x, selections = "Father", addBaseline = FALSE,
#'         disableMutations = FALSE)
#'
#'
#' # Disable all mutation models. SHOULD GIVE ERROR FOR SECOND MARKER
#' # MPPsims(x, selections = "Father", addBaseline = FALSE,
#' #         disableMutations = TRUE)
#'
#' }
#'
#' @importFrom parallel makeCluster stopCluster detectCores parLapply
#'   clusterEvalQ clusterExport clusterSetRNGStream
#' @export
MPPsims = function(reference, missing = "MP", selections, ep = TRUE, ip = TRUE,
                   addBaseline = TRUE, nProfiles = 1, lrSims = 1, thresholdIP = NULL,
                   disableMutations = NA, numCores = 1, seed = NULL, verbose = TRUE) {
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

  # Setup parallelisation
  if(is.na(numCores))
    numCores = max(detectCores() - 1, 1)

  if(numCores > 1) {
    cl = makeCluster(numCores)
    on.exit(stopCluster(cl))
    if(verbose) {
      message("Preparing parallelisation using ", length(cl), " cores")
    }
    clusterEvalQ(cl, library(forrel))
    clusterExport(cl, c("missingPersonEP", "missingPersonIP", "lrSims", "thresholdIP"), envir = environment())
    clusterSetRNGStream(cl, iseed = sample.int(1e6,1))
  }

  powSims = lapply(seq_along(selections), function(i) {
    ref = reference[[i]]
    ids = selections[[i]]

    if(verbose)
      message("Selection: ", names(selections)[i])

    # Baseline
    if(is.null(ids)) {
      ep0 = if(ep) epfun(ref) else NULL
      ip0 = if(ip) ipfun(ref) else NULL
      return(list(ep = ep0, ip = ip0))
    }

    # Simulate profile for `ids`
    sims = profileSim(ref, ids = ids, N = nProfiles, numCores = 1, simplify1 = FALSE, verbose = FALSE)


    # Compute updated EP and IP for each profile
    if(numCores == 1) {
      epRes = if(ep) lapply(sims, epfun) else NULL
      ipRes = if(ip) lapply(sims, ipfun) else NULL
    }
    else {
      epRes = if(ep) parLapply(cl, sims, epfun) else NULL
      ipRes = if(ip) parLapply(cl, sims, ipfun) else NULL
    }
    list(ep = epRes, ip = ipRes)
  })

  names(powSims) = names(selections)

  # Timing
  totalTime = ftime(st)
  if(verbose)
    message("Total time used: ", totalTime)

  # Return with attributes
  structure(powSims, reference = reference, selections = selections,
            nProfiles = nProfiles, lrSims = lrSims, thresholdIP = thresholdIP,
            seed = seed, totalTime = totalTime, class = "MPPsim")
}


#' @export
`[.MPPsim` = function(x, i) {
  structure(unclass(x)[i],
            reference = attr(x, 'reference'),
            selections = attr(x, 'selections')[i],
            nProfiles = attr(x, 'nProfiles'),
            lrSims = attr(x, 'lrSims'),
            thresholdIP = attr(x, 'thresholdIP'),
            seed = "NA (order changed)",
            totalTime = attr(x, 'totalTime'),
            class = "MPPsim")

}
