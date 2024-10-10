# Prepare data for fast computation of IBD likelihood
# Input: ids = a vector of ID labels
# Output: List of lists. For each indiv, the output contains a list of 5 vectors:
#   a1: first allele
#   a2: second allele
#   f1: frequency of first allele
#   f2: frequency of second allele
#   miss: Indices of markers with missing data
.prepAlleleData = function(x, ids) {
  nMark = nMarkers(x)
  nSeq = seq_len(nMark)
  ids = as.character(ids)
  isList = is.pedList(x)

  if(isList) {
    pednr = getComponent(x, ids, checkUnique = TRUE)
    if(all(pednr == pednr[1])) {
      x = x[[pednr[1]]]
      isList = FALSE
    }
  }

  if(isList) {
    alsMat = do.call(rbind, lapply(x, function(cmp) matrix(unlist(cmp$MARKERS), ncol = 2*nMark)))

    # freqlist = lapply(nSeq, function(i) unname(afreq(x, i))) # previous, slow but secure
    x1 = x[[1]] # faster to use only 1st comp
    freqlist = lapply(x1$MARKERS, function(m) attr(m, "afreq"))
  }
  else {
    alsMat = matrix(unlist(x$MARKERS), ncol = 2*nMark)
    freqlist = lapply(x$MARKERS, function(m) attr(m, "afreq")) # faster than the generic afreq()
  }

  rownames(alsMat) = labels(x, unlist = TRUE)
  alsMat = alsMat[ids, , drop = FALSE]

  maxAlNum = max(lengths(freqlist))
  freqMat = unlist(lapply(freqlist, `length<-`, maxAlNum), recursive = FALSE, use.names = FALSE)
  dim(freqMat) = c(maxAlNum, nMark)

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

# Version of .prepAlleleData for interaction with e.g. profileSimParametric
# als: list of vectors a,b,c,d; where a/b and c/d are the genotypes
# freqList: e.g. NorwegianFreqs
.prepAlleleData2 = function(als, freqList) {
  alsMat = do.call(rbind, als)
  freqMat = vapply(seq_along(freqList),
                   function(i) freqList[[i]][alsMat[, i]],
                   FUN.VALUE = numeric(4))

  list(list(a1 = als[[1]], a2 = als[[2]], f1 = freqMat[1, ], f2 = freqMat[2, ]),
       list(a1 = als[[3]], a2 = als[[4]], f1 = freqMat[3, ], f2 = freqMat[4, ]))
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
      d$miss = integer(0)
      d
    })
  }
  dat
}

# Input: alleleData = output from .prepAlleleData() for a pair of indivs
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

