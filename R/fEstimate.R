#' Estimation of inbreeding coefficients
#'
#' Estimate the inbreeding coefficient \eqn{f} of pedigree members from their
#' genotypes. The default method uses maximum likelihood, but moment estimators
#' are also available.
#'
#' The ML estimator assumes independent markers and known allele frequencies.
#' For a single marker with allele frequencies \eqn{p_a}, the genotype
#' likelihood \eqn{P(G \mid f)} is:
#'
#' \deqn{P(G = a/a \mid f) = f p_a + (1-f)p_a^2,} \deqn{P(G = a/b \mid f) =
#' 2(1-f)p_a p_b,\quad a \neq b.}
#'
#' For multiple markers, the full likelihood is the product over markers. The
#' function uses `optimise()` to maximise the log-likelihood over \eqn{[0,1]}.
#'
#' The `"simple"` method is a method-of-moments estimator based on observed and
#' expected homozygosity. Let \eqn{L} be the number of markers, \eqn{N_H} the
#' number of observed homozygous genotypes, and \eqn{H^* = \sum_1^L \sum_a
#' p_a^2} the expected number under Hardy--Weinberg equilibrium. The estimator
#' is then:
#'
#' \deqn{\hat{f} = \dfrac{N_H - H^*}{L - H^*}.}
#'
#' The `"ritland"` method implements the estimator of Ritland (1996). For marker
#' \eqn{j}, let \eqn{n_j} be the number of alleles, and define \eqn{T_j} to be
#' \eqn{1/p_a} if the genotype is \eqn{a/a}, and 0 otherwise. The estimator is
#' then:
#'
#' \deqn{\hat{f} = \dfrac{\sum_{j=1}^L (T_j - 1)}{\sum_{j=1}^L (n_j - 1)}
#'       = \dfrac{\sum_{\text{hom}} 1/p_a - L}{A - L},}
#'
#' where \eqn{A = \sum n_j} is the total number of alleles across markers. Note
#' that both this and the "simple" estimator may produce estimates outside the
#' interval \eqn{[0,1]}.
#'
#' The MLE approach is described explicitly in Section 3.2 of Thompson (1986).
#' Ritland (1996) introduces the Ritland estimator, and also mentions both the
#' simple estimator and the MLE. Note that Ritlands formulation of the MLE is
#' slightly different (but equivalent) to Thompson's.
#'
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector with ID labels, by default, all genotyped members of `x`.
#' @param method A character, either "mle" (maximum-likelihood estimate;
#'   default), "simple" (simple method-of-moments estimator), or "ritland" (see
#'   Details).
#' @param markers A vector with names or indices of markers attached to x,
#'   indicating which markers to include. By default, all markers are included.
#'   Missing genotypes are removed before estimation.
#' @param verbose A logical.
#'
#' @return A named numeric vector.
#'
#' @seealso [ibdEstimate()]
#'
#' @references
#'
#' * Thompson, E. A. (1986). Pedigree Analysis in Human Genetics.
#' _Johns Hopkins University Press_.
#' * Ritland, K. (1996). Estimators for pairwise relatedness and individual
#' inbreeding coefficients. _Genetical Research_.
#'
#' @examples
#'
#' # Simulate genotype for 100 SNPs, for a child of full sibs
#' x = nuclearPed(2, sex = 1:2) |> addSon(3:4) |>
#'   markerSim(N = 100, ids = 5, alleles = 1:2, afreq = c(0.3, 0.7), seed = 123)
#' x
#'
#' # Estimate inbreeding coefficient (pedigree: 0.25)
#' fEstimate(x) # MLE
#' fEstimate(x, method = "simple")
#' fEstimate(x, method = "ritland")
#'
#'
#' #-----------------------------
#' # Compare different estimators
#' #-----------------------------
#'
#' # Simulate 200 DNA profiles (35 STR markers) for a child of half-sibs
#' x = halfSibPed(sex2 = 2) |> addSon(4:5) |>
#'   profileSim(N = 200, ids = 6, markers = NorwegianFrequencies, seed = 123)
#'
#' # Estimates using all three methods
#' fhat = data.frame(
#'   mle = sapply(x, fEstimate, method = "mle"),
#'   sim = sapply(x, fEstimate, method = "simple"),
#'   rit = sapply(x, fEstimate, method = "ritland")
#' )
#'
#' # Mean estimates (expected: 0.125)
#' apply(fhat, 2, mean)
#'
#' # RMSE
#' apply(fhat, 2, function(v) sqrt(mean((v - 0.125)^2)))
#'
#' @export
fEstimate = function(x, ids = typedMembers(x), method = c("mle", "simple", "ritland"),
                     markers = NULL, verbose = FALSE) {
  st = Sys.time()
  method = match.arg(method)

  if(!is.null(markers))
    x = selectMarkers(x, markers)

  ids = as.character(ids)
  if(!all(ids %in% typedMembers(x)))
    stop2("Untyped pedigree member: ", setdiff(ids, typedMembers(x)))

  # Alleles and frequencies
  alleleData = .prepAlleleData(x, ids = ids)

  if(verbose) {
    s = "Estimating inbreeding coefficients for %d individuals using %d markers (method: %s)"
    message(sprintf(s, length(ids), nMarkers(x), method))
  }

  if(method == "mle") {
    res = lapply(alleleData, .fest_mle)
  }
  else if(method == "simple") {
    db = getFreqDatabase(x)
    res = lapply(alleleData, .fest_simple, db = db)
  }
  else if(method == "ritland") {
    nals = nAlleles(x)
    res = lapply(alleleData, .fest_ritland, nals = nals)
  }

  res = unlist(res)

  if(verbose)
    message("Total time: ", ftime(st))

  res
}


#' @importFrom stats optimise
.fest_mle = function(dat) {

  dat = .removeMiss1(dat)
  n = length(dat$a1)
  het = dat$a1 != dat$a2
  p1 = dat$f1
  p2 = dat$f2
  p1sq = p1^2

  loglik = function(f) {
    lik = f * p1 + (1-f) * p1sq  #hom
    lik[het] = (1-f) * 2 * p1[het] * p2[het]
    sum(log(lik), na.rm = TRUE)
  }
  opt = optimise(loglik, interval = c(0,1), maximum = TRUE)
  opt$maximum
}

.fest_simple = function(dat, db) {
  if(length(dat$miss)) {
    db = db[-dat$miss]
    dat = .removeMiss1(dat)
  }
  L = length(db)
  if(L != length(dat$a1))
    stop2("Length mismatch in .fest_simple")

  # Expected homozygosity under HWE
  Hstar = sum(unlist(db)^2)

  # Observed number of homozygotes
  nHom = sum(dat$a1 == dat$a2)

  # Inbreeding estimate
  (nHom - Hstar) / (L - Hstar)
}

.fest_ritland = function(dat, nals) {
  if(length(dat$miss)) {
    nals = nals[-dat$miss]
    dat = .removeMiss1(dat)
  }
  L = length(nals)
  if(L != length(dat$a1))
    stop2("Length mismatch in .fest_ritland")

  ishom = dat$a1 == dat$a2
  numer = sum(1/dat$f1[ishom]) - L
  denom = sum(nals - 1)
  numer / denom
}

.removeMiss1 = function(dat) {
  miss = dat$miss
  if(!length(miss))
    return(dat)
  dat$a1 = dat$a1[-miss]
  dat$a2 = dat$a2[-miss]
  dat$f1 = dat$f1[-miss]
  dat$f2 = dat$f2[-miss]
  dat$miss = integer(0)
  dat
}


# Special likelihood functions ------------------------------------------


#' Parent-child pedigree likelihoods
#'
#' Fast likelihood computation of pedigrees in which the only genotyped
#' individuals are a (noninbred) parent and their (possibly inbred) child. In
#' this case the pedigree information is fully captured by the inbreeding
#' coefficient \eqn{f} of the child, enabling efficient calculations by
#' vectorising over the markers.
#'
#' @param x A `ped` object with exactly two typed members.
#' @param parent ID of parent; identified automatically if missing.
#' @param child ID of child; identified automatically if missing.
#' @param f Inbreeding coefficient of child, calculated with
#'   `ribd::inbreeding()` if missing.
#' @param check A logical indicating whether to check parent and child in `x`.
#'
#' @returns A numeric vector with per-marker likelihoods P(parent & child
#'   genotypes | f).
#'
#' @examples
#'
#' db = NorwegianFrequencies[1:10]
#'
#' # Simulate genotypes for a child of first cousins, and his mother
#' x = cousinPed(1, child = TRUE) |> profileSim(ids = 8:9, markers = db)
#'
#' plot(x, hatched = typedMembers)
#'
#' liks = parentChildLikelihood(x)
#'
#' # Compare with standard likelihood calculation, requiring loop breaking etc.
#' stopifnot(all.equal(liks, pedprobr::likelihood(x)))
#'
#' # bench::mark(parentChildLikelihood(x), pedprobr::likelihood(x))
#'
#' @export
#' @importFrom ribd inbreeding
parentChildLikelihood = function(x, parent = NULL, child = NULL, f = NULL, check = TRUE) {
  if(check || is.null(parent) || is.null(child)) {
    dat = .checkParentChild(x, parent, child)
    parent = dat$parent
    child = dat$child
  }

  if(is.null(f))
    f = inbreeding(x, child)

  d = .prepAlleleData(x, c(parent, child))
  m = d[[parent]]
  c = d[[child]]

  moHom = m$a1 == m$a2
  chHom = c$a1 == c$a2

  m1c1 = c$a1 == m$a1
  m1c2 = c$a2 == m$a1
  m2c1 = c$a1 == m$a2
  m2c2 = c$a2 == m$a2

  # Mother genotype under HWE (per marker)
  pm = 2*m$f1*m$f2
  pm[moHom] = m$f1[moHom]^2

  # P(father -> child allele | mother genotype); "IBD with one of mother's alleles" has prob 2f
  pt1 = (1 - 2*f)*c$f1 + f * (m1c1 + m2c1)
  pt2 = (1 - 2*f)*c$f2 + f * (m1c2 + m2c2)

  # Maternal transmission weights: (1,0) if homozygous; (0.5,0.5) otherwise
  w1 = 0.5 + moHom/2
  w2 = 1 - w1

  # Child genotype given maternal transmitted allele u = m1 (pc1) or u = m2 (pc2)
  pc1 = m1c1 * pt2 + m1c2 * pt1
  pc2 = m2c1 * pt2 + m2c2 * pt1
  pc1[chHom] = m1c1[chHom] * pt1[chHom]
  pc2[chHom] = m2c1[chHom] * pt1[chHom]

  pm * (w1*pc1 + w2*pc2)
}


.checkParentChild = function(x, parent = NULL, child = NULL) {
    typed = typedMembers(x)
    hasPar = !is.null(parent)
    hasCh = !is.null(child)

    if(length(typed) != 2)
      stop2("The pedigree must contain exactly two typed members")

    if(hasPar && !(parent %in% typed))
      stop2("The indicated parent is not typed")
    if(hasCh && !(child %in% typed))
      stop2("The indicated child is not typed")

    if(!hasPar && !hasCh) {
      if(typed[1] %in% parents(x, typed[2]))
        return(list(parent = typed[1], child = typed[2]))
      if(typed[2] %in% parents(x, typed[1]))
        return(list(parent = typed[2], child = typed[1]))
      stop2("Typed members are not in a parent-child relationship")
    }

    if(!hasPar) {
      parent = intersect(typed, parents(x, child))
      if(!length(parent))
        stop2("The indicated child has no typed parents")
      return(list(parent = parent, child = child))
    }

    if(!hasCh) {
      child = intersect(typed, children(x, parent))
      if(!length(child))
        stop2("Could not uniquely identify child of parent '", parent, "'")
      return(list(parent = parent, child = child))
    }

    if(!(parent %in% parents(x, child)))
      stop2("The indicated parent is not a parent of the indicated child")

    list(parent = parent, child = child)
}
