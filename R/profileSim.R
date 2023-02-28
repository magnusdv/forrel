#' Simulation of complete DNA profiles
#'
#' Simulation of DNA profiles for specified pedigree members. Some pedigree
#' members may already be genotyped; in that case the simulation is conditional
#' on these. The main work of this function is done by [markerSim()].
#'
#' @param x A `ped` object or a list of such.
#' @param N The number of complete simulations to be performed.
#' @param ids A character (or coercible to character) with ID labels indicating
#'   whose genotypes should be simulated.
#' @param markers A list of marker objects, or a vector containing names or
#'   indices referring to markers attached to `x`. By default (`NULL`), all
#'   attached markers are used. The simulations will be conditional on the locus
#'   attributes (allele frequencies, mutation models a.s.o.) and any existing
#'   genotypes in the indicated markers.
#' @param seed An integer seed for the random number generator (optional).
#' @param numCores The number of cores to be used. The default is 1, i.e., no
#'   parallelisation.
#' @param simplify1 A logical, by default TRUE, removing the outer list layer
#'   when `N = 1`. See Value.
#' @param verbose A logical, by default TRUE.
#' @param ... Further arguments passed on to [markerSim()].
#'
#' @return A list of `N` objects similar to `x`, but with simulated genotypes.
#'   Any previously attached markers are replaced by the simulated profiles. If
#'   the indicated markers contained genotypes for some pedigree members, these
#'   are still present in the simulated profiles.
#'
#'   If `N = 1` and `simplify1 = TRUE`, the outer list layer is removed, i.e.,
#'   `profileSim(..., N = 1, simplify1 = T)` is equivalent to `profileSim(..., N
#'   = 1, simplify1 = F)[[1]]`. This is usually the desired object in
#'   interactive use, and works well with piping.
#'
#'   When using `profileSim()` in other functions, it is recommended to add
#'   `simplify1 = FALSE` to safeguard against issues with `N = 1`.
#'
#' @examples
#' # Example with two brothers
#' x = nuclearPed(children = c("B1", "B2"))
#'
#' # Attach two markers; one brother is already genotyped
#' x = addMarker(x, B1 = "1/2", alleles = 1:3)
#' x = addMarker(x, B1 = "1", alleles = 1:4, afreq = (1:4)/10, chrom = "X")
#'
#' # Simulate 3 profiles of B2 conditional on the above
#' profileSim(x, N = 3, ids = "B2")
#'
#' @importFrom parallel makeCluster stopCluster detectCores parLapply
#'   clusterEvalQ clusterExport clusterSetRNGStream
#' @export
profileSim = function(x, N = 1, ids = NULL, markers = NULL, seed = NULL,
                      numCores = 1, simplify1 = TRUE, verbose = TRUE, ...){

  if(!is.ped(x) && !is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

  # Check for linked markers
  if(hasLinkedMarkers(x))
    warning("Linked markers detected. Be aware that this function ignores linkage.\n",
            "(You may want to use `ibdsim2::profileSimIBD()` instead.)", call. = FALSE)

  # Set seed once (instead of passing it to markerSim)
  if(!is.null(seed))
    set.seed(seed)

  # Check that all `ids` are in x
  labs = if(is.ped(x)) labels(x) else unlist(labels(x))
  if(length(err <- setdiff(ids, labs)))
    stop2("Unknown ID label: ", err)

  # If pedlist input: Recurse over components
  if(is.pedList(x)) {
    if(is.marker(markers) || is.markerList(markers))
      stop2("When `x` is a list of pedigrees, `markers` must be a vector of marker names/indices referring to attached markers")

    res_compwise = lapply(x, function(comp)
      profileSim(comp, N = N, ids = if(!is.null(ids)) intersect(ids, labels(comp)),
                 markers = markers, numCores = numCores, simplify1 = FALSE, verbose = verbose, ...))

    # Transpose: Collect j'th sim of each component.
    res = lapply(1:N, function(j) lapply(res_compwise, `[[`, j))

    if(simplify1 && N == 1)
      res = res[[1]]

    return(res)
  }

  # If no markers are indicated, use all attached markers
  if(is.null(markers))
    markers = seq_len(nMarkers(x))

  # If marker object(s), attach
  if(is.marker(markers)) {
    x = setMarkers(x, markers, checkCons = TRUE)
    markers = 1L
  }
  if(is.markerList(markers)) {
    x = setMarkers(x, markers, checkCons = TRUE)
    markers = seq_along(markers)
  }

  if(length(markers) == 0) {
    message("Empty profile; returning `ped` object unchanged")
    return(x)
  }

  # Marker names/positions are lost in the sims - must be re-added
  map = getMap(x, markers = markers, verbose = FALSE)
  useMap = !all(is.na(map))


  ### SIMULATIONS ###

  # Parallelise?
  if(is.na(numCores))
    numCores = max(detectCores() - 1, 1)

  if(numCores > 1) {
    cl = makeCluster(numCores)
    on.exit(stopCluster(cl))
    if(verbose) {
      message("Preparing parallelisation using ", length(cl), " cores")
    }
    clusterEvalQ(cl, library(forrel))
    clusterExport(cl, c("markerSim", "N", "ids"), envir = environment())

    # Random number seed for cluster. NB: Keep this outside of function call
    iseed = sample.int(1e6,1)
    clusterSetRNGStream(cl, iseed = iseed)

    # Iterate over the loci, make N simulations of each.
    sims_markerwise = parLapply(cl, markers, function(pm)
      markerSim(x, N = N, ids = ids, partialmarker = pm, verbose = FALSE))
  }
  else {
    # Iterate over the loci, make N simulations of each.
    sims_markerwise = lapply(markers, function(pm)
        markerSim(x, N = N, ids = ids, partialmarker = pm, verbose = FALSE, ...))
  }

  ### Transpose: Extract i'th marker from each sim above.
  # Output: List of N `ped`s, each with length(markers) attached markers
  sims = lapply(1:N, function(i) {
    mlist = lapply(sims_markerwise, function(y) y$MARKERS[[i]])
    s = setMarkers(x, mlist, checkCons = FALSE)

    # Add names/positions if necessary
    if(useMap)
      s = setMap(s, map, matchNames = FALSE)

    s
  })

  if(simplify1 && N == 1)
    sims = sims[[1]]

  sims
}
