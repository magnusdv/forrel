context("simple simulation")

library(pedtools)

test_that("simpleSim() runs in trivial example", {
  x = nuclearPed(1)
  y = simpleSim(x, N=1, alleles=1:2, verbose=F)
  expect_is(y, "ped")
  expect_equal(nMarkers(y), 1)
  expect_true(all(y$markerdata[[1]] > 0))
})

test_that("markerSim() runs in simple example", {
  x = nuclearPed(1)
  m = marker(x, '3'=1, alleles=1:2)
  y = markerSim(x, N=1, partialmarker=m, verbose=F)
  expect_is(y, "ped")
  expect_equal(nMarkers(y), 1)
  expect_equal(genotype(y, 1, id=3), c("1", "1"))
  expect_true(all(y$markerdata[[1]] > 0))
})

test_that("markerSim() works with partial with mutations", {
  x = nuclearPed(1)
  m = marker(x, '1'=1, '2'=1, alleles=1:2, mutmod = "eq", rate=1)
  y = markerSim(x, N=1, partialmarker=m, verbose=F)
  expect_equal(genotype(y, 1, id=3), c("2", "2"))
})

test_that("markerSim() works with explicit mutation args", {
  x = nuclearPed(1)
  m = pedmut::mutationModel("custom", matrix = matrix(c(0,0,1,1), ncol=2, dimnames=list(1:2,1:2)))
  y = markerSim(x, N=1, alleles = 1:2, afreq = c(.99, .01), mutmod = m, verbose=F)
  expect_equal(genotype(y, 1, id=3), c("2", "2"))
})
