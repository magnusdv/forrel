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
#' @param numCores The number of cores used for parallelisation, by default 1.
#' @param verbose A logical, by default TRUE.
#' @param ... Further arguments passed on to [markerSim()].
#'
#' @return A list of `N` objects similar to `x`, but with simulated genotypes.
#'   Any previously attached markers are replaced by the simulated profiles. If
#'   the indicated markers contained genotypes for some pedigree members, these
#'   are still present in the simulated profiles.
#'
#' @examples
#' # Example with two brothers
#' x = nuclearPed(children = c("B1", "B2"))
#'
#' # Attach two markers; one brother is already genotyped
#' m1 = marker(x, B1 = 1:2, alleles = 1:3)
#' m2 = marker(x, B1 = 1, alleles = 1:4, afreq = (1:4)/10, chrom = "X")
#' x = setMarkers(x, list(m1, m2))
#'
#' # Simulate 3 profiles of B2 conditional on the above
#' profileSim(x, N = 3, ids = "B2")
#'
#' @importFrom parallel makeCluster stopCluster detectCores parLapply
#'   clusterEvalQ clusterExport clusterSetRNGStream
#' @export
profileSim = function(x, N = 1, ids = NULL, markers = NULL, seed = NULL,
                      numCores = 1, verbose = TRUE, ...){

  if(!is.ped(x) && !is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

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
                 markers = markers, numCores = numCores, verbose = verbose, ...))

    # Transpose: Collect j'th sim of each component.
    res = lapply(1:N, function(j) lapply(res_compwise, `[[`, j))
    return(res)
  }

  # If no markers are indicated, use all attached markers
  if(is.null(markers))
    markers = seq_len(nMarkers(x))

  # If single marker object, convert to list
  if(is.marker(markers))
    markers = list(markers)

  if(length(markers) == 0) {
    message("Empty profile; returning `ped` object unchanged")
    return(x)
  }

  # Marker names are lost in the sims - must be re-added
  if (is.markerList(markers))
    mnames = vapply(markers, name, "")
  else if (is.atomic(markers))
    mnames = name(x, markers)
  nonNAs = which(!is.na(mnames))


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
    s = setMarkers(x, mlist)

    # Add names if necessary
    if(length(nonNAs) > 0)
      name(s, nonNAs) = mnames[nonNAs]

    s
  })

  sims
}
