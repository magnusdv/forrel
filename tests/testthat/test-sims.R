context("Marker simulation")

library(pedtools)

test_that("simpleSim() runs in trivial example", {
  x = nuclearPed(1)
  y = simpleSim(x, N=1, alleles=1:2, verbose=F)
  expect_is(y, "ped")
  expect_equal(nMarkers(y), 1)
  expect_true(all(y$MARKERS[[1]] > 0))
})

test_that("markerSim() runs in simple example", {
  x = nuclearPed(1)
  m = marker(x, '3'=1, alleles=1:2)
  y = markerSim(x, N=1, partialmarker=m, verbose=F)
  expect_is(y, "ped")
  expect_equal(nMarkers(y), 1)
  expect_equal(genotype(y, 1, id=3), c("1", "1"))
  expect_true(all(y$MARKERS[[1]] > 0))
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

test_that("profileSim() keeps marker names", {
  x = nuclearPed(1)
  m = marker(x)
  x = setMarkers(x, list(m, m))
  name(x, 1:2) = c("m1", "m2")

  s = profileSim(x, N = 1)[[1]]
  expect_identical(name(s, 1:2), c("m1", "m2"))

  s2 = profileSim(x, N = 1, markers = list(m, m))[[1]]
  expect_identical(name(s2, 1:2), rep(NA_character_, 2))
})

test_that("profileSim() treats pedlists as expected", {
  x = singleton(1)
  x = setMarkers(x, marker(x, alleles = 1:5, name = "M"))
  y = relabel(x, 2)
  SEED = 777

  sim_pedlist = profileSim(list(x, y), N = 3, markers = "M", seed = SEED)

  set.seed(SEED)
  sim_compwise = list(profileSim(x, N = 3, markers = "M"),
                      profileSim(y, N = 3, markers = "M"))

  # Check third sim
  expect_identical(sim_pedlist[[3]], lapply(sim_compwise, `[[`, 3))
})

test_that("markerSim() works with peds in non-standard ordering", {
x = reorderPed(nuclearPed(2), 4:1)
m = marker(x, "3" = 1, "4" = 2) # parents must be 1:2

y = markerSim(x, partialmarker = m, verbose = F)
expect_identical(genotype(y, id = 1, marker = 1), c("1", "2"))
})

test_that("markerSim() works in looped pedigree 1", {
  x = addChildren(linearPed(2), 5,2,1)
  x = setMarkers(x, marker(x, "5" = 1, alleles = 1:2, afreq = c(0.001, 0.999)))
  set.seed(123)
  y = markerSim(x, partialmarker = 1, verbose = FALSE)
  expect_identical(genotype(y, id = 3, marker = 1), c("1", "2"))
  expect_identical(genotype(y, id = 4, marker = 1), c("1", "2"))
})

test_that("markerSim() works in looped pedigree 2", {
  x = addChildren(linearPed(2), 5,2,1)
  x = setMarkers(x, marker(x, "6" = 1, "2" = c(0,2)))
  # plot(x,1)
  y1 = markerSim(x,partialmarker = 1, verbose = FALSE)
  expect_identical(genotype(y1, id = 2, marker = 1), c("1", "2")) # Forced

  x = addParents(x, 2, 10, 11, verbose = FALSE)
  x = setMarkers(x, marker(x, "6" = 1, "10" = 2))
  # plot(x,1)
  y2 = markerSim(x, partialmarker = 1, verbose = FALSE)
  expect_identical(genotype(y2, id = 2, marker = 1), c("1", "2")) # Forced!

  y3 = markerSim(x, partialmarker = 1, loopBreaker = "5", verbose = FALSE)
  expect_identical(genotype(y3, id = 2, marker = 1), c("1", "2")) # Forced!
})
