#' Simulation of complete DNA profiles
#'
#' Simulation of DNA profiles for specified pedigree members. Some pedigree
#' members may already be genotyped; in that case the simulation is conditional
#' on these. The main work of this function is done by [markerSim()].
#'
#' ## Parallel profile simulation
#'
#' Simulations can be distributed across several R processes using the
#' **mirai** package. To use this, start mirai workers before calling
#' `profileSim()`:
#'
#' ```r
#' mirai::daemons(8, seed = 123)
#'
#' x = nuclearPed()
#' db = NorwegianFrequencies
#' res = profileSim(x, N = 100, markers = db)
#'
#' mirai::daemons(0)
#' ```
#'
#' The call to `mirai::daemons(8, seed = 123)` starts 8 persistent workers and
#' configures reproducible parallel RNG. These workers remain available for
#' later calls until stopped with `mirai::daemons(0)`.
#'
#' If no mirai daemons are set, `profileSim()` runs normally, and the `seed`
#' argument works as before. When mirai daemons are used, set the seed in
#' `mirai::daemons(n, seed = ...)` and omit the `seed` argument to `profileSim()`.
#'
#' @param x A `ped` object or a list of such.
#' @param N The number of complete simulations to be performed.
#' @param ids A character (or coercible to character) of ID labels indicating
#'   whose genotypes should be simulated. Alternatively, a function taking `x`
#'   as input and returning a character vector of ID labels.
#' @param markers Either a vector indicating a subset of markers attached to
#'   `x` (default: all), or a named list of frequency vectors. In the former
#'   case, simulations are conditional on all locus attributes (allele
#'   frequencies, mutation models, etc) and any existing genotypes in the
#'   indicated markers. If `markers` is a list of frequency vectors, these
#'   are attached as new markers to `x` (removing any existing markers) and
#'   used for unconditional simulations.
#' @param loopBreakers A vector containing IDs of individuals to be used as
#'   loop breakers. Can usually be left as NULL, in which case automatic loop
#'   breaking is applied. See [pedtools::breakLoops()].
#' @param seed An integer seed for the random number generator. Used only when
#'   no mirai workers are set. For parallel simulations, set the seed in
#'   [mirai::daemons()].
#' @param numCores Deprecated. Use [mirai::daemons()] to set the number of
#'   workers.
#' @param simplify1 A logical, by default TRUE, removing the outer list layer
#'   when `N = 1`. See Value.
#' @param verbose A logical, by default TRUE.
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
#' # Example pedigree with two brothers
#' x = nuclearPed(children = c("B1", "B2"))
#'
#' ### Simulate profiles using built-in freq database
#' profileSim(x, markers = NorwegianFrequencies[1:3])
#'
#' ### Conditioning on known genotypes for one brother
#'
#' # Attach two SNP markers with genotypes for B1
#' y = x |>
#'   addMarker(B1 = "1/2", alleles = 1:2) |>
#'   addMarker(B1 = "1",   alleles = 1:2, chrom = "X")
#'
#' # Simulate 2 profiles of B2 conditional on the above
#' profileSim(y, N = 2, ids = "B2", seed = 123)
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
profileSim = function(x, N = 1, ids = NULL, markers = NULL, loopBreakers = NULL,
                      seed = NULL, numCores = 1, simplify1 = TRUE, verbose = TRUE) {

  if(!is.ped(x) && !is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

  st = Sys.time()

  if(!missing(numCores))
    warning("`numCores` is deprecated; use `mirai::daemons()` for parallelisation.", call. = FALSE)

  # If `markers` is a list of frequency vectors, attach as new markers
  if(length(markers) && is.list(markers) && !is.marker(markers[[1]]) && is.numeric(markers[[1]])) {
    nms = names(markers)
    if(is.null(nms <- names(markers)))
      stop2("`markers` appears to be a list of frequency vectors, but marker names are missing")

    checkFreqs = vapply(markers, function(m) is.numeric(m) && round(sum(m)) == 1,
                        FUN.VALUE = TRUE, USE.NAMES = FALSE)

    if(!all(checkFreqs))
      stop2("`markers` appears to be a list of frequency vectors, but some entries do not sum to 1: ",
            nms[!checkFreqs])

    x = setMarkers(x, locusAttributes = markers, checkCons = FALSE)
    markers = nms
    if(verbose)
      message(sprintf("Attached %d markers based on frequency database", length(nms)))
  }
  else {
    if(!is.null(markers))
      x = selectMarkers(x, markers)

    if(hasLinkedMarkers(x))
      warning("Linked markers detected. Be aware that this function ignores linkage.\n",
              "(You may want to use `ibdsim2::profileSimIBD()` instead.)", call. = FALSE)
  }

  if(!hasMarkers(x)) {
    if(verbose) message("Empty profile; returning pedigree unchanged")
    return(if(simplify1 && N == 1) x else rep(list(x), N))
  }

  # Check that all `ids` are in x
  if(is.function(ids))
    ids = ids(x)

  labs = labels(x, unlist = TRUE)
  if(length(err <- .mysetdiff(ids, labs)))
    stop2("Unknown ID label: ", err)

  # Parallelisation
  useMirai = mirai::daemons_set()
  nworkers = mirai::info()[["connections"]] %||% 0L

  if(verbose) {
    if(useMirai)
      message(sprintf("Parallelisation: %d mirai %s", nworkers, pluralise("worker", nworkers)))
    else
      message("No parallelisation; use `mirai::daemons(n)` to set up workers")
  }

  if(!is.null(seed)) {
    if(useMirai)
      stop2("`seed` is incompatible with mirai workers; set with `mirai::daemons(n, seed = ...)`")
    else
      set.seed(seed)
  }

  # If pedlist input: Recurse over components
  if(is.pedList(x)) {
    if(is.marker(markers) || is.markerList(markers))
      stop2("When `x` is a list of pedigrees, `markers` must be a vector indicating to attached markers")

    res_compwise = lapply(x, function(comp)
      profileSim(comp, N = N, ids = .myintersect(comp$ID, ids), markers = NULL,
               simplify1 = FALSE, loopBreakers = .myintersect(comp$ID, loopBreakers),
               verbose = FALSE))

    # Transpose: Collect j'th sim of each component.
    res = lapply(1:N, function(j) lapply(res_compwise, `[[`, j))

    if(simplify1 && N == 1)
      res = res[[1]]

    if(verbose)
      message(sprintf("Simulated %d %s in %s", N, pluralise("profile", N),
                      format(Sys.time() - st, digits = 2)))
    return(res)
  }


  # From here: Connected component --------------------------------------------------------------

  Nmark = length(x$MARKERS)

  # Marker names are lost in the sims and must be re-added
  mnames = name(x, markers = seq_len(Nmark))

  progbar = verbose && interactive() && Nmark > 1

  if(useMirai) {
    mm = mirai::mirai_map(1:Nmark, .profileSimMarker,
                          .args = list(x = x, N = N, ids = ids, lb = loopBreakers))
    sims_markerwise = if(progbar) mm[mirai::.stop, mirai::.progress] else mm[mirai::.stop]

  }
  else {
    if(progbar)
      pb = txtProgressBar(min = 0, max = Nmark, style = 3)

    sims_markerwise = lapply(1:Nmark, function(j) {
      a = .profileSimMarker(j, x = x, N = N, ids = ids, lb = loopBreakers)
      if(progbar)
        setTxtProgressBar(pb, j)
      a
    })

    if(progbar) close(pb)
  }

  ### Transpose: Extract i'th marker from each sim above.
  # Output: List of N `ped`s, each with Nmark attached markers
  sims = lapply(1:N, function(i) {
    mlist = lapply(1:Nmark, function(j) {
      m = sims_markerwise[[j]]$MARKERS[[i]]
      attr(m, "name") = mnames[[j]]
      m
    })
    class(mlist) = "markerList"
    setMarkers(x, mlist, checkCons = FALSE)
  })

  if(simplify1 && N == 1)
    sims = sims[[1]]

  if(verbose)
    message(sprintf("Simulated %d %s in %s", N, pluralise("profile", N),
                    format(Sys.time() - st, digits = 2)))
  sims
}

.profileSimMarker = function(j, x, N, ids, lb = NULL)
  markerSim(x, N, ids = ids, partialmarker = j, loopBreakers = lb, verbose = FALSE)
