context("Exclusion power")

test_that("exclusion power works in paternity case", {
  claim = nuclearPed(1)
  true = list(singleton(1), singleton(3))
  ids = c(1,3)
  als = c(.1, .9)
  pe = exclusionPower(claim, true, ids, alleles=1:length(als), afreq = als)
  expect_equal(pe, sum(2 * als[1]^2 * als[2]^2))
})

test_that("exclusion power works in paternity case with extra singleton", {
  claim = list(nuclearPed(1), singleton(4))
  true = list(singleton(1), singleton(3), singleton(4))
  ids = c(1,3,4)
  als = c(.1, .9)
  pe = exclusionPower(claim, true, ids, alleles=1:length(als), afreq = als)
  expect_equal(pe, sum(2 * als[1]^2 * als[2]^2))
})

