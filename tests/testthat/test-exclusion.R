context("Exclusion power")

quickPE = function(claim, true, ids, afreq)
  exclusionPower(claim, true, ids, alleles=1:length(afreq), afreq = afreq, plot = F)

test_that("PE works in paternity case", {
  claim = nuclearPed(1)
  true = list(singleton(1), singleton(3))
  ids = c(1,3)
  afr = c(.1, .9)
  pe = quickPE(claim, true, ids, afr)
  expect_equal(pe, sum(2 * afr[1]^2 * afr[2]^2))
})

test_that("PE in pat. case with extra singleton in claim AND true ", {
  claim = list(nuclearPed(1), singleton(4))
  true = list(singleton(1), singleton(3), singleton(4))
  ids = c(1, 3, 4)
  afr = c(.1, .9)
  pe = quickPE(claim, true, ids, afr)
  expect_equal(pe, sum(2 * afr[1]^2 * afr[2]^2))
})

test_that("PE in pat. case with extra indiv: singleton (claim) / sib (true)", {
  claim = list(nuclearPed(1), singleton(4))
  true = list(singleton(1), singleton(3), singleton(4))
  ids = c(1, 3, 4)
  afr = c(.1, .9)
  pe = quickPE(claim, true, ids, afr)
  expect_equal(pe, sum(2 * afr[1]^2 * afr[2]^2))
})

