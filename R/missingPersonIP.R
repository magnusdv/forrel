#' Inclusion power for missing person cases
#'
#' This function simulates the LR distribution for the true missing person in a
#' reference family. The output contains both the total and marker-wise LR of
#' each simulation, as well as various summary statistics. If a specific LR
#' threshold is given, the _inclusion power_ is computed as the probability that
#' LR exceeds the threshold.
#'
#' @inheritParams missingPersonEP
#' @param nsim A positive integer: the number of simulations
#' @param threshold A numeric vector with one or more positive numbers used as
#'   the likelihood ratio thresholds for inclusion
#' @param seed An integer seed for the random number generator (optional).
#'
#' @return A `mpIP` object, which is essentially a list with the following
#'   entries:
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
#'   of `threshold`, the fraction of simulations resulting in a LR exceeding the
#'   given number.
#'
#'   * `time`: The total computation time.
#'
#'   * `params`: A list containing the input parameters `missing`, `markers`,
#'   `nsim`, `threshold` and `disableMutations`
#'
#' @examples
#'
#' # Four siblings; the fourth is missing
#' x = nuclearPed(4)
#'
#' # Remaining sibs typed with 5 triallelic markers
#' x = markerSim(x, N = 5, ids = 3:5, alleles = 1:3, seed = 123, verbose = FALSE)
#'
#' # Compute inclusion power statistics
#' ip = missingPersonIP(x, missing = 6, nsim = 5, threshold = c(10, 100))
#' ip
#'
#' # LRs from each simulation
#' ip$LRperSim
#'
#' @importFrom pedprobr likelihood
#' @export
missingPersonIP = function(reference, missing, markers, nsim = 1, threshold = NULL,
                           disableMutations = NA, seed = NULL, verbose = TRUE) {
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
  hasMut = allowsMutations(reference, markers)

  # For which marker should mutations be disabled?
  disALL = any(hasMut) && isTRUE(disableMutations)
  disGOOD = any(hasMut) && length(disableMutations) == 1 && is.na(disableMutations)
  disSELECT = any(hasMut) && !isFALSE(disableMutations) && !is.null(disableMutations)

  if(disALL)
    disable = markers[hasMut]
  else if(disGOOD) { # disable only if consistent
    refNoMut = reference
    mutmod(refNoMut, markers[hasMut]) = NULL
    liksNoMut = likelihood(refNoMut, markers = markers[hasMut])
    disable = markers[hasMut][liksNoMut > 0]
  }
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
  midx = whichMarkers(reference, markers)
  relatedPed = relabel(reference, old = missing, new = poiLabel)
  unrelatedPed = list(reference, singleton(poiLabel, sex = getSex(reference, missing)))

  # Raise error if impossible markers
  liks = likelihood(reference, markers = midx)
  if(any(liks == 0))
    stop2("Marker incompatible with reference pedigree: ", markers[liks == 0],
          "\nThis makes conditional simulations impossible. Exclude the marker from the computation or add a mutation model")

  # Set seed once
  set.seed(seed)

  # Simulate nsim complete profiles of relatedPed
  if(verbose)
    message("Simulating ", nsim, " profiles...", appendLF = FALSE)

  allsims = profileSim(relatedPed, ids = "_POI_", N = nsim, markers = midx, simplify1 = FALSE, verbose = FALSE)

  if(verbose)
    message("done\nComputing likelihood ratios...", appendLF = FALSE)

  # Compute LR of each marker
  lrs = vapply(allsims, function(s) {
    unrelSim = transferMarkers(from = s, to = unrelatedPed)

    lr = kinshipLR(list(s, unrelSim), ref = 2)
    lr$LRperMarker[, 1]
  }, FUN.VALUE = numeric(length(markers)))

  # Ensure matrix
  if(length(markers) == 1)
    lrs = rbind(lrs, deparse.level = 0)

  rownames(lrs) = markers

  if(verbose)
    message("done")

  # Results
  LRperSim = apply(lrs, 2, prod)
  meanLRperMarker = apply(lrs, 1, mean)
  meanLR = mean(LRperSim)
  meanLogLR = mean(log10(LRperSim))
  IP = sapply(threshold, function(thr) mean(LRperSim >= thr))
  names(IP) = threshold

  # Timing
  time = Sys.time() - st
  if(verbose)
    message("Total time used: ", format(time, digits = 3))

  # List of input parameters
  params = list(missing = missing, markers = markers,
                nsim = nsim, threshold = threshold, seed = seed,
                disableMutations = disableMutations)

  structure(list(LRperSim = LRperSim, meanLRperMarker = meanLRperMarker,
                 meanLR = meanLR, meanLogLR = meanLogLR, IP = IP,
                 time = time, params = params), class = "LRpowerResult")
}


