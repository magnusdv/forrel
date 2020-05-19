#' Relatedness estimation
#'
#' Estimate the pairwise IBD coefficients \eqn{(\kappa_0, \kappa_1,
#' \kappa_2)}{(\kappa0, \kappa1, \kappa2)} for specified pairs of pedigree
#' members, using maximum likelihood methods. The optimisation machinery is
#' imported from the `maxLik` package.
#'
#' This function optimises the log-likelihood function first described in
#' (Thompson, 1975). Optimisation is done in the \eqn{(\kappa_0,
#' \kappa_2)}{(\kappa0, \kappa2)}-plane and restricted to the probability
#' triangle defined by \deqn{\kappa_0 \ge 0, \kappa_2 \ge 0, \kappa_0 + \kappa_2
#' \le 1.}{\kappa0 \ge 0, \kappa2 \ge 0, \kappa0 + \kappa2 \le 1.}
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either a vector with ID labels, or a data frame/matrix with two
#'   columns, where each row contains the ID labels of two individuals. The entries
#'   are coerced to characters, and must match uniquely against the ID labels of
#'   `x`. If `ids` is a vector, it is converted to a matrix containing all pairs.
#'   By default, all individuals of `x` are included.
#' @param markers A numeric indicating which marker(s) to include. If NULL
#'   (default), all markers are used.
#' @param start Numeric of length 2, indicating the initial value of
#'   \eqn{(\kappa_0, \kappa_2)}{(\kappa0, \kappa2)} in the optimisation (passed
#'   on to `maxLik`).
#' @param tol A single numeric: the optimising tolerance value; passed on to
#'   `maxLik()`).
#'
#' @return A data frame with 6 columns: `ID1`, `ID2`, `N` (the number of markers
#'   with no missing alleles), `k0`, `k1` and `k2`.
#' @author Magnus Dehli Vigeland
#' @seealso [maxLik::maxLik()], [showInTriangle()]
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
#' x = nuclearPed(children = c("sib1", "sib2"))
#'
#' # Simluate 200 equifrequent SNPs
#' x = markerSim(x, N = 200, alleles = 1:2, verbose = FALSE)
#'
#' # Estimate IBD coefficients (exact = (0.25, 0.5, 0.25))
#' est = IBDestimate(x, ids = c("sib1", "sib2"))
#' showInTriangle(est, labels = TRUE)
#'
#' ### Example 2: Unrelated singletons
#' y = list(singleton(1), singleton(2))
#' y = markerSim(y, N = 200, alleles = 1:2, verbose = FALSE)
#'
#' IBDestimate(y, ids = 1:2)
#'
#'
#' @importFrom maxLik maxLik
#' @export
IBDestimate = function(x, ids = NULL, markers = NULL, start = c(0.99,0.001), tol = 1e-7) {
  if(is.ped(x))
    x = list(x)
  else if(!is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

  if(is.null(markers))
    markers = seq_len(nMarkers(x))

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

  ids_df_list = lapply(seq_len(nrow(ids)), function(i) pedlistMembership(x, ids[i,]))

  # Optimize the above function in the triangle
  constraints = list(ineqA = matrix(c(1,0,-1,0,1,-1), nrow = 3, ncol = 2),
                     ineqB = c(0,0,1))

  res = lapply(ids_df_list, function(ids_df) {
    A = IBDest_getAlleleData(x, ids_df, markers)
    # Remove markers with missing alleles
    if(any(miss <- apply(A, 2, anyNA)))
        A = A[, !miss, drop = FALSE]

    # Likelihood function
    a=A[1,]; b=A[2,]; cc=A[3,]; d=A[4,]; pa=A[5,]; pb=A[6,]; pc=A[7,]; pd=A[8,]
    loglik_FUN = function(k) sum(log(.IBDlikelihood(k,a,b,cc,d,pa,pb,pc,pd)))

    # Optimise
    ML = maxLik(loglik_FUN, start = start, constraints = constraints, tol = tol)
    est = ML$estimate
    data.frame(ID1 = ids_df$id[1],
               ID2 = ids_df$id[2],
               N = ncol(A),
               k0 = est[1],
               k1 = 1 - sum(est),
               k2 = est[2],
               stringsAsFactors = FALSE)
  })

  do.call(rbind, res)
}

# Match ID labels against a list of pedigrees
pedlistMembership = function(pedl, ids) {
    ped_match = vapply(pedl,
                       function(p) ids %in% labels(p),
                       FUN.VALUE = logical(length(ids)))

    # Convert to matrix if length(ids == 1)
    ped_match = rbind(ped_match, deparse.level = 0)

    if(any(nomatch <- rowSums(ped_match) == 0))
      stop2("IDs not found in pedigrees: ", ids[nomatch])
    if(any(multimatch <- rowSums(ped_match) > 1))
      stop2("IDs matching multiple pedigrees: ", ids[multimatch])

    pednr = apply(ped_match, 1, which)

    data.frame(pednr = pednr, id = ids, stringsAsFactors = FALSE)
}


IBDest_getAlleleData = function(x, ped_id_df, markers = NULL) {
  stopifnot(is.pedList(x), is.data.frame(ped_id_df))
  pednr1 = ped_id_df$pednr[1]
  pednr2 = ped_id_df$pednr[2]

  # Match IDs with orig.ids.
  id1_int = internalID(x[[pednr1]], ped_id_df$id[1])
  id2_int = internalID(x[[pednr2]], ped_id_df$id[2])

  # Collect alleles and frequencies in a matrix with 8 rows
  # (first 4 alleles a,b,c,d followed by their freqs), one column per marker
  # Note that alleles are internal integers

  if(pednr1 == pednr2) {
    ped = x[[pednr1]]
    A = vapply(ped$MARKERS[markers], function(m) {
      als = c(m[id1_int,], m[id2_int,])
      frq = rep_len(NA_real_, length(als))
      frq[als > 0] = attr(m, 'afreq')[als] # works, since 0's in als are ignored when indexing
      c(als, frq)
    }, FUN.VALUE = numeric(8))
  }
  else {
    ped1 = x[[pednr1]]
    ped2 = x[[pednr2]]
    A = vapply(markers, function(i) {
      m1 = ped1$MARKERS[[i]]
      m2 = ped2$MARKERS[[i]]
      als = c(m1[id1_int,], m2[id2_int,])
      frq = rep_len(NA_real_, length(als))
      frq[als > 0] = attr(m1, 'afreq')[als] # works, since 0's in als are ignored when indexing
      c(als, frq)
    }, FUN.VALUE = numeric(8))
  }

  return(A)
}


.IBDlikelihood = function(k, a, b, cc, d, pa, pb, pc, pd) {
  ### Vectorized function for computing kappa likelihoods, given genotypes for two related individuals
  # k: numeric of length 2 = (kappa0, kappa2)
  # a: vector of positive integers (allele 1 of individual 1)
  # b: vector of positive integers (allele 2 of individual 1)
  # cc: vector of positive integers (allele 1 of individual 2) (avoid overloading of 'c')
  # d: vector of positive integers (allele 2 of individual 2)
  # pa, pb, pc, pd: numeric vectors with frequencies of the above alleles.
  #
  # Note that all input vectors (except k) have the same length = #markers.
  homoz1 = a == b
  homoz2 = cc == d
  mac = a == cc
  mbc = b == cc
  mad = a == d
  mbd = b == d
  g1.fr = 2^(!homoz1) * pa * pb
  g2.fr = 2^(!homoz2) * pc * pd

  # Prob(g1, g2 | unrelated)
  UN = g1.fr * g2.fr

  # Prob(g1, g2 | parent-offspring)
  PO = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd))

  # Prob(g1, g2 | monozygotic twins)
  MZ = g1.fr * ((mac & mbd) | (mad & mbc))

  # return likelihoods (Thompson)
  k[1] * UN + (1-k[1]-k[2]) * PO + k[2] * MZ
}
