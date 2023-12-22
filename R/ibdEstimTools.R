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
    list(a1 = a1, a2 = a2, f1 = f1, f2 = f2, miss = which(!(a1nonz & a2nonz)))
  })

  names(res) = ids
  res
}

.removeMissing = function(dat) {
  anymiss = lapply(dat, function(d) d$miss) |>
    unlist(recursive = FALSE) |>
    unique.default()

  if(length(anymiss)) {
    dat = lapply(dat, function(d) {
      d[[1]] = d[[1]][-anymiss]
      d[[2]] = d[[2]][-anymiss]
      d[[3]] = d[[3]][-anymiss]
      d[[4]] = d[[4]][-anymiss]
      d
    })
  }
  dat
}

# Input: alleleData = output from .getAlleleData() for a pair of indivs
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


stationary = function(p, grad, tol) {
  # print(grad,digits = 12)
  # print(sapply(seq_along(p), function(i) {pi = p; pi[i] = p[i] - 1; grad %*% pi}))
  for(i in seq_along(p)) {
    pi = p
    pi[i] = p[i] - 1
    if(grad %*% pi < -tol)
      return(FALSE)
  }
  TRUE
}

# Project any vector in R^D onto the probability simplex.
# Algorithm found in Wang et al. (2013): Projection onto the probability simplex.
simplexProject = function(y) {
  D = length(y)
  if(D == 3)
    u = .sort3(y)
  else
    u = sort.int(y, decreasing = TRUE, method = "shell")
  v = u + 1/seq_len(D) * (1 - cumsum(u))
  p = match(TRUE, v <= 0, nomatch = D + 1) - 1
  lambda = 1/p * (1 - sum(u[1:p]))
  x = y + lambda
  x[x < 0] = 0
  x
}

.sort3 = function(y) {
  if(y[1] < y[2]) {tmp = y[1]; y[1] = y[2]; y[2] = tmp}
  if(y[2] < y[3]) {tmp = y[2]; y[2] = y[3]; y[3] = tmp}
  if(y[1] < y[2]) {tmp = y[1]; y[1] = y[2]; y[2] = tmp}
  y
}





##### Moved from the old IBDestimate

# TODO: These should be merged into the newer functions and removed.

# Used in checkPairwise()
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
  # TODO: use afreq() to extract frqs. This checks consistency across components!

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
  pb = freqMat[2,]
  pc = freqMat[3,]
  pd = freqMat[4,]

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

