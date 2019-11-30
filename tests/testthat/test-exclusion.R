context("Exclusion power")

quickEP = function(claim, true, ids, afreq = NULL, ...) {
  als = if(is.null(afreq)) NULL else seq_along(afreq)
  exclusionPower(claim, true, ids, alleles = als, afreq = afreq,
                 plot = F, verbose = F, ...)$EPtotal
}

test_that("EP works in empty paternity case", {
  claim = nuclearPed(1)
  true = list(singleton(1), singleton(3))
  ids = c(1, 3)
  afr = c(.1, .9)

  # Autosomal
  ep_aut = quickEP(claim, true, ids, afr)
  expect_equal(ep_aut, 2 * afr[1]^2 * afr[2]^2)

  # X: father/son (no exclusion)
  expect_equal(quickEP(claim, true, ids, afr, Xchrom = T), 0)

  # X: father/daughter
  claim2 = nuclearPed(1, sex = 2)
  true2 = list(singleton(1), singleton(3, sex = 2))

  ep_X = quickEP(claim2, true2, ids, afr, Xchrom = T)
  expect_equal(ep_X, afr[1] * afr[2]) # p^2*q + p*q^2
})

test_that("EP works in empty pat-case with added singletons", {
  claim = list(nuclearPed(1), singleton(4))
  true = list(singleton(1), singleton(3), singleton(4))
  ids = c(1, 3, 4)
  afr = c(.1, .9)

  # Autosomal
  expect_equal(quickEP(claim, true, ids, afr), 2 * afr[1]^2 * afr[2]^2)

  # X
  expect_equal(quickEP(claim, true, ids, afr, Xchrom = T), 0)

  # X: father/daughter
  claim2 = list(nuclearPed(1, sex = 2), singleton(4))
  true2 = list(singleton(1), singleton(3, sex = 2), singleton(4))
  expect_equal(quickEP(claim2, true2, ids, afr, Xchrom = T), afr[1] * afr[2])
})

test_that("EP works in paternity case with child typed", {
  claim = nuclearPed(1, sex = 2)
  true = list(singleton(1), singleton(3))
  afr = c(.5, .3, .2)

  m = mX = marker(claim, `3` = 1, alleles = 1:3, afreq = afr)
  chrom(mX) = 23L
  claim = setMarkers(claim, list(m, mX))

  # Autosomal
  ep_aut = quickEP(claim, true, ids = 1, markers = 1)
  expect_equal(ep_aut, sum(afr[-1])^2)

  # X
  ep_X = quickEP(claim, true, ids = 1, markers = 2)
  expect_equal(ep_X, sum(afr[-1]))
})

test_that("EP works in paternity case with parents typed", {
  claim = nuclearPed(1, sex = 2)
  true = list(singleton(1), singleton(3, sex = 2))
  afr = c(.5, .3, .2)

  m = mX = marker(claim, `1` = 1, `2` = 2, alleles = 1:3, afreq = afr)
  chrom(mX) = 23L
  claim = setMarkers(claim, list(m, mX))

  # Autosomal
  ep_aut = quickEP(claim, true, ids = 3, markers = 1)
  expect_equal(ep_aut, 1 - 2*afr[1]*afr[2])

  # X
  ep_X = quickEP(claim, true, ids = 3, markers = 2)
  expect_equal(ep_X, 1 - 2*afr[1]*afr[2])
})
