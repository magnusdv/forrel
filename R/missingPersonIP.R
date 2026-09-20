#' LR simulation for missing person cases
#'
#' This function simulates the LR distribution in a missing person case, either for the
#' true missing person or for an unrelated person. The output contains both the total and
#' marker-wise LR of each simulation, as well as various summary statistics. If a specific
#' LR threshold is given, the fraction of simulations exceeding the threshold is computed.
#' When simulating the true missing person (the default), this fraction is referred to
#' as the *inclusion power* (Vigeland et al., 2020).
#'
#' @inheritParams missingPersonEP
#' @param nsim A positive integer: the number of simulations
#' @param threshold A numeric vector with one or more positive numbers used as the
#'   likelihood ratio thresholds for inclusion
#' @param seed An integer seed for the random number generator (optional).
#' @param true Either "missing" (default) or "unrelated", indicating whether profiles are
#'   simulated for the true missing person or for an unrelated person.
#'
#' @return A `mpIP` object, which is essentially a list with the following entries:
#'
#'   * `LRperSim`: A numeric vector of length `nsim` containing the total LR for
#'   each simulation.
#'
#'   * `meanLRperMarker`: The mean LR per marker, over all simulations.
#'
#'   * `meanLR`: The mean total LR over all simulations.
#'
#'   * `meanLogLR`: The mean total `log10(LR)` over all simulations.
#'
#'   * `IP`: A named numeric of the same length as `threshold`. For each element
#'   of `threshold`, the fraction of simulations resulting in an LR exceeding the given
#'   number. With `true = "missing"` this is the inclusion power; with `true =
#'   "unrelated"` it is the corresponding fraction among unrelated persons.
#'
#'   * `params`: A list containing the input parameters `missing`, `markers`,
#'   `nsim`, `threshold`, `disableMutations` and `true`.
#'
#' @references
#' Vigeland MD, Marsico FL, Herrera Piñero M, Egeland T (2020).
#' "Prioritising family members for genotyping in missing person cases:
#' A general approach combining the statistical power of exclusion and inclusion."
#' *FSI: Genetics*, 49, 102376. \doi{10.1016/j.fsigen.2020.102376}
#'
#' @seealso [missingPersonEP()], [missingPersonLR()], [missingPersonPlot()]
#'
#' @examples
#'
#' # Two brothers are looking for their missing sibling.
#' # They are typed with 5 triallelic markers.
#' x = nuclearPed(3) |>
#'   markerSim(N = 5, ids = 3:4, alleles = 1:3, seed = 123, verbose = FALSE)
#'
#' missingPersonPlot(x, missing = 5)
#'
#' nsim = 20 # increase!
#'
#' # Inclusion power statistics
#' ip = missingPersonIP(x, missing = 5, nsim = nsim, threshold = c(10, 100))
#' ip
#' head(ip$LRperSim)
#'
#' # Simulate LRs for a random unrelated person
#' ip2 = missingPersonIP(x, missing = 5, nsim = nsim, threshold = c(10, 100),
#'                       true = "unrelated")
#' ip2
#' head(ip2$LRperSim)
#'
#' # Plot distributions
#' LRpowerPlot(data = list(ip, ip2), threshold = 100)
#'
#' @importFrom pedprobr likelihood
#' @export
missingPersonIP = function(reference, missing, markers, nsim = 1, threshold = NULL,
                           disableMutations = NA, seed = NULL,
                           true = c("missing", "unrelated"), verbose = TRUE) {
  st = Sys.time()

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")

  true = match.arg(true)

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

  # Do any of the markers model mutations?
  hasMut = allowsMutations(reference, markers)

  # For which marker should mutations be disabled?
  disALL = any(hasMut) && isTRUE(disableMutations)
  disGOOD = any(hasMut) && length(disableMutations) == 1 && is.na(disableMutations)
  disSELECT = any(hasMut) && !isFALSE(disableMutations) && !is.null(disableMutations)

  if(disALL)
    disable = markers[hasMut]
  else if(disGOOD) # disable only if consistent
    disable = consistentMarkers(reference, hasMut, names = TRUE)
  else if(disSELECT)
    disable = whichMarkers(reference, disableMutations)
  else
    disable = NULL

  # Disable mutations in the chosen cases
  if(length(disable) > 0) {
    if(verbose) message("Disabling mutations for marker ", toString(disable))
    mutmod(reference, disable) = NULL
  }

  poiLabel = "_POI_"

  # Extract markers and set up pedigrees
  reference = selectMarkers(reference, markers)
  relatedPed = relabel(reference, old = missing, new = poiLabel)
  unrelatedPed = list(reference, singleton(poiLabel, sex = getSex(reference, missing)))

  # Raise error if impossible markers
  imp = inconsistentMarkers(reference, names = TRUE, removeMut = FALSE)
  if(length(imp))
    stop2("Marker incompatible with reference pedigree: ", imp,
          "\nThis makes conditional simulations impossible. Exclude the marker from the computation or add a mutation model")

  # Set seed once
  if(!is.null(seed)) {
    if(.miraiWorkers() > 0L)
      stop2("`seed` is incompatible with mirai workers; set with `mirai::daemons(n, seed = ...)`")
    set.seed(seed)
  }

  # Pedigree used for simulation
  simPed = if(true == "missing") relatedPed else transferMarkers(reference, unrelatedPed[[2]])

  # Simulate nsim complete profiles
  if(verbose) {
    who = if(true == "missing") "the true missing person" else "unrelated person"
    message(sprintf("Simulating %d profile%s for %s...", nsim, pluralise(nsim), who), appendLF = FALSE)
  }

  allsims = profileSim(simPed, ids = "_POI_", N = nsim, simplify1 = FALSE, verbose = FALSE)

  if(verbose)
    message("done\nComputing likelihood ratios...", appendLF = FALSE)

  # Compute log-LR of each marker
  lrs = vapply(allsims, function(s) {
    if(true == "missing") {
      relSim = s
      unrelSim = transferMarkers(from = s, to = unrelatedPed)
    }
    else {
      relSim = transferMarkers(from = s, to = relatedPed, erase = FALSE, matchNames = FALSE)
      unrelSim = list(reference, s)
    }
    lr = kinshipLR(relSim, unrelSim, ref = 2)
    lr$lnLRperMarker[, 1] / log(10)
  }, FUN.VALUE = numeric(length(markers)))

  # Ensure matrix
  if(length(markers) == 1)
    lrs = rbind(lrs, deparse.level = 0)

  rownames(lrs) = markers

  if(verbose)
    message("done")

  # Results
  log10LRperSim = colSums(lrs)
  LRperSim = 10^log10LRperSim
  meanLRperMarker = rowMeans(10^lrs)
  meanLR = mean(LRperSim)
  meanLogLR = mean(log10LRperSim)
  IP = sapply(threshold, function(thr) mean(log10LRperSim >= log10(thr)))
  names(IP) = threshold

  # Timing
  time = Sys.time() - st # included in output
  if(verbose)
    message("Total time used: ", format(time, digits = 3))

  # List of input parameters
  params = list(missing = missing, markers = markers,
                nsim = nsim, threshold = threshold, seed = seed,
                disableMutations = disableMutations, true = true)

  structure(list(LRperSim = LRperSim, meanLRperMarker = meanLRperMarker,
                 meanLR = meanLR, meanLogLR = meanLogLR, IP = IP, params = params,
                 log10LRperSim = log10LRperSim),
            class = c("mpIP", "LRpowerResult"))
}


