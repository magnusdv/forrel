#' Simulation of DNA profiles
#'
#' Simulation of DNA profiles for specified pedigree members. Some pedigree
#' members may already be genotyped; in that case the simulation is conditional
#' on these. The main work of this function is done by [markerSim()].
#'
#' @param x A `ped` object
#' @param N The number of complete simulations to be performed
#' @param ids A character (or coercible to character) with ID labels indicating
#'   whose genotypes should be simulated
#' @param conditions A list of marker objects, (or an integer vector indicating
#'   markers attached to `x`). If these contain genotypes the simulations will
#'   condition on these. Locus attributes (allele frequencies, mutation models
#'   a.s.o.) are extracted from each marker.
#' @param seed NULL, or a numeric seed for the random number generator
#' @param ... Further arguments passed on to [markerSim()]
#'
#' @return A list of length `N`. Each element is a `ped` object with `K`
#'   attached markers, where `K = length(conditions)`.
#'
#' @examples
#' # Example with two brothers, one of which is already genotyped with 2 markers.
#' x = nuclearPed(children = c("B1", "B2"))
#'
#' m1 = marker(x, B1 = 1:2, alleles = 1:3, afreq = c(.2, .3, .5))
#' m2 = marker(x, B1 = 4, alleles = 1:4, afreq = c(.1,.2,.3,.4))#, chrom="X")
#'
#' # These contain the profile of B1
#' cond = list(m1, m2)
#'
#' # Simulate 3 profiles of B2 conditional on the above
#' profileSim(x, N = 3, ids = "B2", conditions = cond)
#'
#'
#' @export
profileSim = function(x, N = 1, ids = NULL, conditions = NULL, seed = NULL, ...){

  # Set seed once (instead of passing it to markerSim)
  if(!is.null(seed))
    set.seed(seed)

  # If pedlist input: Recurse over components
  if(is.pedList(x)) {
    if(is.marker(conditions) || is.markerList(conditions))
      stop2("When `x` is a list of pedigrees, `conditions` must be a vector of marker names/indices referring to attached markers")

    res_compwise = lapply(x, function(comp)
      profileSim(comp, N = N, ids = if(!is.null(ids)) intersect(ids, labels(comp)),
                 conditions = conditions, ...))

    # Transpose: Collect j'th sim of each component.
    res = lapply(1:N, function(j) lapply(res_compwise, `[[`, j))
    return(res)
  }

  # Ensure conditions is a list
  if(is.marker(conditions))
    conditions = list(conditions)
  nCond = length(conditions)

  # Marker names are lost in the sims - must be re-added
  if (nCond > 0) {
    if (is.markerList(conditions))
      mnames = vapply(conditions, name, "")
    else if (is.atomic(conditions))
      mnames = name(x, conditions)
    nonNAs = which(!is.na(mnames))
  }
  addNames = nCond > 0 && length(nonNAs) > 0


  # Iterate over the loci, make N simulations of each.
  sims_markerwise = lapply(conditions, function(pm)
      markerSim(x, N = N, ids = ids, partialmarker = pm, verbose=F, ...))

  # Transpose: Extract i'th marker from each sim above.
  # Output: List of length N, each with length(conditions) markers
  sims = lapply(1:N, function(i) {
      mlist = lapply(sims_markerwise, function(y) y$markerdata[[i]])
      s = setMarkers(x, mlist)
      if(addNames)
        name(s, nonNAs) = mnames[nonNAs]
      s
    })

  sims
}
