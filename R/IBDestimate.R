#' Relatedness estimation
#'
#' Estimate the pairwise IBD coefficients \eqn{(\kappa_0, \kappa_1,
#' \kappa_2)}{(\kappa0, \kappa1, \kappa2)} for specified pairs of pedigree
#' members, using maximum likelihood methods. The optimisation machinery is
#' imported from the `maxLik` package.
#'
#' This function optimises the log-likelihood function first described by
#' Thompson (1975). Optimisation is done in the \eqn{(\kappa_0,
#' \kappa_2)}{(\kappa0, \kappa2)}-plane and restricted to the probability
#' triangle defined by \deqn{\kappa_0 \ge 0, \kappa_2 \ge 0, \kappa_0 + \kappa_2
#' \le 1.}{\kappa0 \ge 0, \kappa2 \ge 0, \kappa0 + \kappa2 \le 1.}
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either a vector with ID labels, or a data frame/matrix with two
#'   columns, where each row contains the ID labels of two individuals. The
#'   entries are coerced to characters, and must match uniquely against the ID
#'   labels of `x`. If `ids` is a vector, it is converted to a matrix containing
#'   all pairs. By default, all individuals of `x` are included.
#' @param markers A vector with names or indices of markers attached to x,
#'   indicating which markers to include. If NULL (default), all markers are
#'   used.
#' @param start Numeric of length 2, indicating the initial value of
#'   \eqn{(\kappa_0, \kappa_2)}{(\kappa0, \kappa2)} in the optimisation (passed
#'   on to `maxLik`).
#' @param tol A single numeric: the optimising tolerance value; passed on to
#'   `maxLik()`).
#' @param contourPlot A logical. If TRUE, contours of the log-likelihood
#'   function are plotted overlaying the IBD triangle.
#' @param levels A numeric vector of levels at which to draw contour lines. If
#'   NULL (default), levels are chosen automatically. (This option is ignored
#'   unless `contourPlot = TRUE`.)
#'
#' @return A data frame with 6 columns: `id1`, `id2`, `N` (the number of markers
#'   with no missing alleles), `kappa0`, `kappa1` and `kappa2`.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [maxLik::maxLik()], [showInTriangle()], [checkPairwise()]
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
#' ids = c("sib1", "sib2")
#'
#' # Pedigree
#' x = nuclearPed(children = ids)
#'
#' # Simulate 100 markers
#' x = markerSim(x, N = 100, alleles = 1:4, seed = 123, verbose = FALSE)
#'
#' # Estimate IBD coefficients (exact = (0.25, 0.5, 0.25))
#' est = IBDestimate(x, ids = ids)
#'
#' # Show the result in the IBD triangle
#' showInTriangle(est, labels = TRUE)
#'
#' # Contour plot (just a few markers to save time)
#' IBDestimate(x, ids = ids, markers = 1:10,
#'             contourPlot = TRUE, levels = -(43:50))
#'
#' ### Example 2: Unrelated singletons
#' y = list(singleton(1), singleton(2))
#'
#' # Simulate 200 SNP markers
#' y = markerSim(y, N = 200, alleles = 1:2, verbose = FALSE)
#'
#' # Estimate
#' IBDestimate(y, ids = 1:2)
#'
#'
#' @importFrom maxLik maxLik
#' @importFrom graphics contour
#' @export
IBDestimate = function(x, ids = NULL, markers = NULL,
                       start = c(0.99,0.001), tol = 1e-7,
                       contourPlot = FALSE, levels = NULL) {

  if(!is.null(markers))
    x = selectMarkers(x, markers)

  # ARgument `ids` should be either a vector or a matrix-like with 2 columns
  if(is.null(ids))
    ids = unlist(labels(x))
  if(length(ids) < 2)
    stop2("`ids` must be either a vector of length at least 2, or a matrix/data.frame with 2 columns")

  if(is.vector(ids))
    ids = .comb2(ids)
  else if(is.data.frame(ids))
    ids = as.matrix(ids)

  if(!is.matrix(ids) && ncol(ids) == 2)
    stop2("`ids` must be either a vector of length at least 2, or a matrix/data.frame with 2 columns")

  # Convert to list of pairs
  pairs = lapply(1:nrow(ids), function(i) ids[i, ])

  # Optimisation constaints: IBD triangle
  constraints = list(ineqA = matrix(c(1,0,-1,0,1,-1), nrow = 3, ncol = 2),
                     ineqB = c(0,0,1))

  res = lapply(pairs, function(pair) {
    dat = .getAlleleData(x, pair)
    amat = dat$alleleMat
    fmat = dat$freqMat

    # Log-likelihood function
    loglik_FUN = function(k) sum(log(.IBDlikelihoodFAST(k, amat, fmat)))

    # Optimise
    ML = maxLik(loglik_FUN, start = start, constraints = constraints, tol = tol)
    est = ML$estimate
    data.frame(id1 = pair[1],
               id2 = pair[2],
               N = ncol(amat),
               kappa0 = est[1],
               kappa1 = 1 - sum(est),
               kappa2 = est[2],
               stringsAsFactors = FALSE)
  })

  res.df = do.call(rbind, res)

  if(contourPlot) { # Not optimised for speed
    if(nrow(ids) > 1)
      stop2("Contour plots require `ids` to be a single pair of individuals")
    pair = ids[1, ]
    n = 51
    k0 = seq(0, 1, length.out = n)
    k2 = seq(0, 1, length.out = n)

    logliks = matrix(NA_real_, ncol = n, nrow = n)
    for(i in 1:n) for(j in 1:n)
      logliks[i,j] <- .IBDlikelihood(x, pair, kappa = c(k0[i], k2[j]))

    if(is.null(levels)) {
      mx = max(logliks)
      levels = pretty(c(mx - 20, mx), n = 8, min.n = 4)
      levels = unique.default(c(floor(mx), floor(mx) + .5, levels)) # include closest integer
    }

    IBDtriangle()
    showInTriangle(res.df, new = FALSE)
    contour(k0, k2, z = logliks, add = TRUE, levels = levels)
  }

  res.df
}


.IBDlikelihood = function(x, ids, kappa, log = TRUE, total = TRUE) {
  if(length(ids) != 2)
    stop2("`ids` must have length 2")

  if(!is.numeric(kappa))
    stop2("`kappa` must be numeric")

  if (!length(kappa) %in% 2:3)
    stop2("`kappa` must have length 2 or 3")

  kappa02 = if(length(kappa) == 3) kappa[c(1, 3)] else kappa

  dat = .getAlleleData(x, ids)
  liks = .IBDlikelihoodFAST(kappa02, dat$alleleMat, dat$freqMat)

  # Correction of negative values (caused by rounding errors)
  liks[liks < 0] = 0

    if(log)
    liks = log(liks)

  if(total)
    if(log) sum(liks) else prod(liks)
  else
    liks
}

# Prepare data for fast computation of IBD likelihood
# Output: 2 matrices w/ 1 col per marker and 4 rows (id1-a1, id1-a2, id2-a1, id2-a2)
#   alleleMat: internal allele indices
#   probMat: allele frequencies corresponding to entries in alleleMat
# Input: ids = a pair of ID labels
.getAlleleData = function(x, ids) {
  pednr = getComponent(x, ids, checkUnique = TRUE)
  pednr1 = pednr[1]
  pednr2 = pednr[2]

  if(pednr1 == pednr2) {
    ped = if(is.ped(x)) x else x[[pednr1]]
    idsInt = internalID(ped, ids)
    A = vapply(ped$MARKERS, function(m) {
      als = c(m[idsInt[1], ], m[idsInt[2], ])
      frq = if(all(als > 0)) attr(m, 'afreq')[als] else rep_len(NA_real_, 4)
      c(als, frq)
    }, FUN.VALUE = numeric(8))
  }
  else {
    ped1 = x[[pednr1]]
    ped2 = x[[pednr2]]
    idsInt1 = internalID(ped1, ids[1])
    idsInt2 = internalID(ped2, ids[2])

    A = vapply(seq_len(nMarkers(x)), function(i) {
      m1 = ped1$MARKERS[[i]]
      m2 = ped2$MARKERS[[i]]
      als = c(m1[idsInt1, ], m2[idsInt2, ])
      frq = if(all(als > 0)) attr(m1, 'afreq')[als] else rep_len(NA_real_, 4)
      c(als, frq)
    }, FUN.VALUE = numeric(8))
  }

  # Missing data: Check any of the freq rows
  miss = is.na(A[5, ])

  # Split alleles and frequencies
  alleleMat = A[1:4, !miss, drop = FALSE]
  mode(alleleMat) = "integer"
  freqMat = A[5:8, !miss, drop = FALSE]

  list(alleleMat = alleleMat, freqMat = freqMat)
}


.IBDlikelihoodFAST = function(kappa02, alleleMat, freqMat) {
  ### Fast computation of kappa likelihoods, given alleles/freqs for two individuals
  # k: numeric of length 2 = (kappa0, kappa2)
  .a = alleleMat[1,]
  .b = alleleMat[2,]
  .c = alleleMat[3,]
  .d = alleleMat[4,]
  pa = freqMat[1,]
  pb = freqMat[1,]
  pc = freqMat[1,]
  pd = freqMat[1,]

  homoz1 = .a == .b
  homoz2 = .c == .d
  mac = .a == .c
  mbc = .b == .c
  mad = .a == .d
  mbd = .b == .d
  g1.fr = 2^(!homoz1) * pa * pb
  g2.fr = 2^(!homoz2) * pc * pd

  # Prob(g1, g2 | unrelated)
  UN = g1.fr * g2.fr

  # Prob(g1, g2 | parent-offspring)
  PO = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd))

  # Prob(g1, g2 | monozygotic twins)
  MZ = g1.fr * ((mac & mbd) | (mad & mbc))

  # return likelihoods (Thompson)
  k0 = kappa02[1]
  k2 = kappa02[2]
  k1 = 1 - k0 - k2
  k0 * UN + k1 * PO + k2 * MZ
}
