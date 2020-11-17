#' Pairwise relatedness estimation
#'
#' Estimate the IBD coefficients \eqn{\kappa = (\kappa_0, \kappa_1,
#' \kappa_2)}{(k0, k1, k2)} or the condensed identity coefficients \eqn{\Delta =
#' (\Delta_1, ..., \Delta_9)}{(d1, ..., d9)} between a pair (or several pairs)
#' of pedigree members, using maximum likelihood methods.
#'
#' It should be noted that this procedure estimates the *realised* identity
#' coefficients of each pair, i.e., the actual fractions of the autosomes in
#' each IBD state. These may deviate substantially from the theoretical pedigree
#' coefficients.
#'
#' Maximum likelihood estimation of relatedness coefficients originates with
#' Thompson (1975). Optimisation of \eqn{\kappa} is done in the \eqn{(\kappa_0,
#' \kappa_2)}{(k0, k2)}-plane and restricted to the triangle defined by
#' \deqn{\kappa_0 \ge 0, \kappa_2 \ge 0, \kappa_0 + \kappa_2 \le 1}{k0 >= 0, k2
#' >= 0, k0 + k2 <= 1}. Optimisation of \eqn{\Delta} is done in unit simplex of
#' R^8, using the first 8 coefficients.
#'
#' The implementation uses [constrOptim()], with the "BFGS" method.
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either a vector with ID labels, or a data frame/matrix with two
#'   columns, where each row contains the ID labels of two individuals. The
#'   entries are coerced to characters, and must match uniquely against the ID
#'   labels of `x`. By default, all pairs of members of `x` are included.
#' @param markers A vector with names or indices of markers attached to x,
#'   indicating which markers to include. If NULL (default), all markers are
#'   used.
#' @param param Either "kappa" (default) or "delta"; indicating which set of
#'   coefficients should be estimated.
#' @param start Numeric of length 3 (if `param = "kappa"`) or 9 (if `param =
#'   "delta"`), indicating an initial value of for the optimisation.
#' @param reltol Relative convergence tolerance; passed on to `constrOptim()`.
#'   Defaults to 1e-12,
#' @param returnArgs A logical; for debugging purposes.
#' @param ... Further arguments.
#'
#' @return If `param = "kappa"`: A data frame with 6 columns: `id1`, `id2`, `N`
#'   (the number of markers with no missing alleles), `k0`, `k1` and `k2`.
#'
#'   If `param = "delta"`: A data frame with 12 columns: `id1`, `id2`, `N` (the
#'   number of markers with no missing alleles), `d1`, `d2`, ... `d9`.
#'
#' @author Magnus Dehli Vigeland
#'
#' @seealso [constrOptim()]
#'
#' @references
#'
#' * E. A. Thompson (1975). _The estimation of pairwise relationships._ Annals
#' of Human Genetics 39.
#'
#' * E. A. Thompson (2000). _Statistical Inference from Genetic Data on
#' Pedigrees._ NSF-CBMS Regional Conference Series in Probability and
#' Statistics. Volume 6.
#'
#' @examples
#'
#' ### Example 1: Siblings
#' x = nuclearPed(2)
#'
#' # Simulate 100 markers
#' x = markerSim(x, N = 100, alleles = 1:4, seed = 123, verbose = FALSE)
#'
#' # Estimate Delta (expectation = (0,0,0,0,0,0,1/4,1/2,1/4))
#' ibdEstim(x, ids = 3:4)
#'
#'
#' ### Example 2: Full sib mating
#' y = fullSibMating(1)
#'
#' # Simulate 200 SNP markers
#' y = markerSim(y, N = 1000, alleles = 1:10, seed = 123, verbose = FALSE)
#'
#' # Estimate
#' ibdEstim(y, ids = 5:6, param = "delta")
#'
#'
#' @export
ibdEstim = function(x, ids = typedMembers(x), param = c("kappa", "delta"),
                    markers = NULL, start = NULL, reltol = 1e-12,
                    returnArgs = FALSE, ...) {

  param = match.arg(param)

  if(!is.null(markers))
    x = selectMarkers(x, markers)

  # Argument `ids` should be either a vector or a matrix-like with 2 columns
  if(is.vector(ids))
    ids = .comb2(ids, vec = TRUE)
  else if(is.data.frame(ids))
    ids = as.matrix(ids)
  if(!is.matrix(ids) || nrow(ids) == 0 || ncol(ids) != 2)
      stop2("`ids` must be either a vector of length at least 2, or a matrix-like with 2 columns")
  mode(ids) = "character"

  # Convert to list of pairs
  pairs = lapply(seq_len(nrow(ids)), function(i) ids[i, ])

  # Alleles and frequencies
  alleleData = .getAlleleData2(x, ids = unique.default(unlist(pairs)))

  FUN = switch(param, kappa = .kappaEstim, delta = .deltaEstim)

  # Debugging
  if(returnArgs) {
    if(length(pairs) != 1) stop2("When debugging, please indicate a single pair")
    res = FUN(alleleData[pairs[[1]]], start = start, reltol = reltol, returnArgs = TRUE, ...)
    return(res)
  }

  # Estimate each pair
  res = lapply(pairs, function(pair) {
    FUN(alleleData[pair], start = start, reltol = reltol, ...)
  })

  res = do.call(rbind, res)

  # If all convergence OK, remove column
  if(all(res$convergence == "OK"))
    res$convergence = NULL

  # Return data frame
  res
}

.kappaEstim = function(dat, start, reltol, returnArgs = FALSE, ...) {
  pair = names(dat)

  # Remove missing
  keep = !is.na(dat[[1]]$f1) & !is.na(dat[[2]]$f2)
  if(!all(keep))
    dat = list(lapply(dat[[1]], `[`, keep), lapply(dat[[2]], `[`, keep))

  # Coordinate-wise likelihoods: P(G_j | IBD = i)
  wei = .likelihoodWeights(dat, param = "kappa")

  # Log-likelihood function: Input (k0, k2)
  loglik = function(p) sum(log(.pairwiseLikelihood(p, wei, param = "kappa")))

  # Gradient
  weiRed = wei[c(1,3), ] - rep(wei[2, ], each = 2)
  grad = function(p) weiRed %*% (1/.pairwiseLikelihood(p, wei, param = "kappa"))

  # Optimise
  ML = .maxlik(loglik, grad, param = "kappa", start = start, reltol = reltol, ...)
  par = ML$res$par

  df = data.frame(id1 = pair[1], id2 = pair[2], N = sum(keep),
                   k0 = par[['k0']], k1 = 1 - sum(par), k2 = par[['k2']],
                   convergence = ML$res$message %||% "OK")

  if(returnArgs)
    list(res = df, args = ML$args)
  else
    df
}


.deltaEstim = function(dat, start, reltol, returnArgs = FALSE, ...) {
  pair = names(dat)

  # Remove missing
  keep = !is.na(dat[[1]]$f1) & !is.na(dat[[2]]$f2)
  if(!all(keep))
    dat = list(lapply(dat[[1]], `[`, keep), lapply(dat[[2]], `[`, keep))

  # Coordinate-wise likelihoods: P(G_j | state_i)
  wei = .likelihoodWeights(dat, param = "delta")

  # Log-likelihood function: Input (d1, ..., d8)
  loglik = function(p) sum(log(.pairwiseLikelihood(p, wei, param = "delta")))

  # Gradient
  weiRed = wei[1:8, ] - rep(wei[9, ], each = 8)
  grad = function(p) weiRed %*% (1/.pairwiseLikelihood(p, wei, param = "delta"))

  # Optimise (returns list(res, args))
  ML = .maxlik(loglik, grad, param = "delta", start = start, reltol = reltol, ...)
  par = ML$res$par

  df = data.frame(id1 = pair[1], id2 = pair[2], N = sum(keep),
                   rbind(par), d9 = 1 - sum(par),
                   convergence = ML$res$message %||% "OK")

  if(returnArgs)
    list(res = df, args = ML$args)
  else
    df
}


.pairwiseLikelihood = function(p, likelMat, param = "kappa") {
  # likelMat: matrix (3 or 9) x nMark. Entry (i,j) = P(G_j | state_i)

  L = length(p)

  if(param == "kappa") {
    if(L == 2)
      p = c(p[1], 1 - sum(p), p[2])
    else if(L != 3)
      stop2("`p` must have length 2 or 3 when `param = 'kappa'`")
  }
  else if(param == "delta") {
    if(L == 8)
      p = c(p, 1 - sum(p))
    else if(L != 9)
      stop2("`p` must have length 8 or 9 when `param = 'delta'`")
  }

  as.numeric(p %*% likelMat)
}

# Prepare data for fast computation of IBD likelihood
# Input: ids = a vector of ID labels
# Output: List of lists. For each indiv, the output contains a list of 4 vectors:
#   a1: first allele
#   a2: second allele
#   f1: frequency of first allele
#   f2: frequency of second allele
.getAlleleData2 = function(x, ids) {
  nMark = nMarkers(x)
  nSeq = seq_len(nMark)
  ids = as.character(ids)

  if(is.pedList(x)) {
    pednr = getComponent(x, ids, checkUnique = TRUE)
    if(all(pednr == pednr[1]))
      x = x[[pednr[1]]]
  }

  if(is.pedList(x)) {
    alsMat = do.call(rbind, lapply(x, function(cmp) matrix(unlist(cmp$MARKERS), ncol = 2*nMark)))
    freqlist = lapply(nSeq, function(i) unname(afreq(x, i))) # checks consistency!
  }
  else {
    alsMat = matrix(unlist(x$MARKERS), ncol = 2*nMark)
    freqlist = lapply(x$MARKERS, function(m) attr(m, "afreq")) # faster than the generic afreq()
  }

  rownames(alsMat) = unlist(labels(x))
  alsMat = alsMat[ids, , drop = FALSE]

  maxAlNum = max(lengths(freqlist))
  freqMat = matrix(unlist(lapply(freqlist, `length<-`, maxAlNum),
                          recursive = FALSE, use.names = FALSE),
                   ncol = nMark)

  even = 2 * nSeq
  odd = even - 1

  res = lapply(ids, function(id) {
      a1 = alsMat[id, odd]
      a2 = alsMat[id, even]
      a1nonz = a1 > 0
      a2nonz = a2 > 0
      f1 = f2 = rep(NA_real_, nMark)
      f1[a1nonz] = freqMat[cbind(a1, nSeq)]
      f2[a2nonz] = freqMat[cbind(a2, nSeq)]
      list(a1 = a1, a2 = a2, f1 = f1, f2 = f2)
  })

  names(res) = ids
  res
}


# Input: alleleData = ouput from .getAlleleData() for a pair of indivs
.likelihoodWeights = function(alleleData, param = "kappa") {
  .a = alleleData[[1]]$a1
  .b = alleleData[[1]]$a2
  .c = alleleData[[2]]$a1
  .d = alleleData[[2]]$a2
  pa = alleleData[[1]]$f1
  pb = alleleData[[1]]$f2
  pc = alleleData[[2]]$f1
  pd = alleleData[[2]]$f2

  homoz1 = .a == .b
  homoz2 = .c == .d
  mac = .a == .c
  mbc = .b == .c
  mad = .a == .d
  mbd = .b == .d
  g1.fr = 2^(!homoz1) * pa * pb
  g2.fr = 2^(!homoz2) * pc * pd

  if(param == "kappa") {
    # Matrix with 3 rows: conditional probs given UN, PO, MZ
    w = matrix(0, nrow = 3, ncol = length(.a))
    w[1, ] = g1.fr * g2.fr                  # Prob(g1, g2 | UN)
    w[2, ] = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd)) # Prob(g1, g2 | PO)
    w[3, ] = g1.fr * ((mac & mbd) | (mad & mbc)) # Prob(g1, g2 | MZ)
  }
  else if(param == "delta") {
    # Matrix with 9 rows: Row i are the conditional probs given state Ji
    w = matrix(0, nrow = 9, ncol = length(.a))
    w[1, ] = pa * (homoz1 & homoz2 & mac)   # Prob(g1, g2 | J1)
    w[2, ] = pa * pc * (homoz1 & homoz2)    # Prob(g1, g2 | J2)
    w[3, ] = pa * homoz1 * (mac *pd + mad * pc) * (.5)^homoz2 # Prob(g1, g2 | J3)
    w[4, ] = pa * g2.fr * homoz1            # Prob(g1, g2 | J4)
    w[5, ] = pc * homoz2 * (mac *pb + mbc * pa) * (.5)^homoz1 # Prob(g1, g2 | J5)
    w[6, ] = g1.fr * pc * homoz2            # Prob(g1, g2 | J6)
    w[7, ] = g1.fr * ((mac & mbd) | (mad & mbc)) # Prob(g1, g2 | MZ)
    w[8, ] = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd)) # Prob(g1, g2 | PO)
    w[9, ] = g1.fr * g2.fr                  # Prob(g1, g2 | UN)
  }
  else
    stop2("`param` must be either 'kappa' or 'delta': ", param)

  w
}

#' @importFrom stats constrOptim
.maxlik = function(loglik, grad, start, param = c("kappa", "delta"),
                   reltol, ...) {
  param = match.arg(param)

  if(param == "kappa") {

    # Create 3-dim version first
    if(is.null(start))
      start = c(1/3, 1/3, 1/3)
    else if(length(start) == 2)
      start = c(start[1], 1-start, start[2])
    else if(length(start) != 3)
      stop2("`start` must have length 2 or 3")

    # If necessary, pull inside
    if(any(start == 0))
      start = (1 - 1e-6) * start + 1e-6 * c(1/3, 1/3, 1/3)

    # Project to 2-dim and fix names
    start = start[c(1,3)]
    names(start) = c("k0", "k2")

    # Region: Triangle
    ui = rbind(diag(2), -1)
    ci = c(0,0,-1)
  }
  else {
    # Start point: In R^9 first
    if(is.null(start))
      start = rep(1/9, 9)
    else if(length(start) == 8)
      start = c(start, 1-start)
    else if(length(start) != 9)
      stop2("`start` must have length 8 or 9")

    # If necessary, pull inside
    if(any(start == 0))
      start = (1 - 1e-6) * start + 1e-6 * rep(1/9, 9)

    # Project to R^8 and fix names
    start = start[1:8]
    names(start) = paste0("d", 1:8)

    # Region: unit simplex in R^8
    ui = rbind(diag(8), -1)
    ci = c(rep(0, 8), -1)
  }

  args = list(theta = start, f = loglik, grad = grad, ui = ui, ci = ci,
              control = list(fnscale = -1, maxit = 1000, reltol = reltol), ...)

  res = do.call(constrOptim, args)

  list(res = res, args = args)
}



# NOT YET USED
# This function would be needed to implement the "projected gradient method",
# for optimising over the probability simplex.
# Project any vector in R^D onto the probability simplex.
# Algorithm found in Wang et al. (2013): Projection onto the probability simplex.
simplexProject = function(y) {
  D = length(y)
  u = sort.default(y, decreasing = TRUE, method = "shell")
  v = u + 1/seq_len(D) * (1 - cumsum(u))
  p = match(T, v <= 0, nomatch = D+1) - 1
  lambda = 1/p * (1 - sum(u[1:p]))
  x = y + lambda
  x[x < 0] = 0
  x
}

