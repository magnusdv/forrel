
mSim = function(...) markerSim(..., verbose = FALSE)
pSim = function(...) profileSim(..., numCores = 1, verbose = FALSE)
sSim = function(...) simpleSim(..., verbose = FALSE)

test_that("simpleSim() runs in trivial example", {
  x = nuclearPed(1)
  y = sSim(x, N=1, alleles=1:2)
  expect_is(y, "ped")
  expect_equal(nMarkers(y), 1)
  expect_true(all(y$MARKERS[[1]] > 0))
})

test_that("markerSim() runs in simple example", {
  x = nuclearPed(1)
  m = marker(x, '3'=1, alleles=1:2)
  y = mSim(x, N=1, partialmarker=m)
  expect_is(y, "ped")
  expect_equal(nMarkers(y), 1)
  expect_equal(genotype(y, 1, id=3), c("1", "1"))
  expect_true(all(y$MARKERS[[1]] > 0))
})

test_that("markerSim() works with partial with mutations", {
  x = nuclearPed(1)
  m = marker(x, '1'=1, '2'=1, alleles=1:2, mutmod = "eq", rate=1)
  y = mSim(x, N=1, partialmarker=m)
  expect_equal(genotype(y, 1, id=3), c("2", "2"))
})

test_that("markerSim() works with explicit mutation args", {
  x = nuclearPed(1)
  m = pedmut::mutationModel("custom", matrix = matrix(c(0,0,1,1), ncol=2, dimnames=list(1:2,1:2)))
  y = mSim(x, N=1, alleles = 1:2, afreq = c(.99, .01), mutmod = m)
  expect_equal(genotype(y, 1, id=3), c("2", "2"))
})

test_that("profileSim() keeps marker names", {
  x = nuclearPed(1)
  m = marker(x)
  x = setMarkers(x, list(m, m))
  name(x, 1:2) = c("m1", "m2")

  s = pSim(x, N = 1)[[1]]
  expect_identical(name(s, 1:2), c("m1", "m2"))

  s2 = pSim(x, N = 1, markers = list(m, m))[[1]]
  expect_identical(name(s2, 1:2), rep(NA_character_, 2))
})

test_that("profileSim() treats pedlists as expected", {
  x = singleton(1)
  x = setMarkers(x, marker(x, alleles = 1:5, name = "M"))
  y = relabel(x, 2)
  SEED = 777

  sim_pedlist = pSim(list(x, y), N = 3, markers = "M", seed = SEED)

  set.seed(SEED)
  sim_compwise = list(pSim(x, N = 3, markers = "M"),
                      pSim(y, N = 3, markers = "M"))

  # Check third sim
  expect_identical(sim_pedlist[[3]], lapply(sim_compwise, `[[`, 3))
})

test_that("markerSim() works with peds in non-standard ordering", {
x = reorderPed(nuclearPed(2), 4:1)
m = marker(x, "3" = 1, "4" = 2) # parents must be 1:2

y = mSim(x, partialmarker = m)
expect_identical(genotype(y, id = 1, marker = 1), c("1", "2"))
})

test_that("markerSim() works in looped pedigree 1", {
  x = addChildren(linearPed(2), 5,2,1)
  x = setMarkers(x, marker(x, "5" = 1, alleles = 1:2, afreq = c(0.001, 0.999)))

  y = mSim(x, partialmarker = 1, seed = 123)
  expect_identical(genotype(y, id = 3, marker = 1), c("1", "2"))
  expect_identical(genotype(y, id = 4, marker = 1), c("1", "2"))
})

test_that("markerSim() works in looped pedigree 2", {
  x = addChildren(linearPed(2), 5,2,1)
  x = setMarkers(x, marker(x, "6" = 1, "2" = c(0,2)))
  # plot(x,1)
  y1 = mSim(x, partialmarker = 1)
  expect_identical(genotype(y1, id = 2, marker = 1), c("1", "2")) # Forced

  x = addParents(x, 2, 10, 11, verbose = FALSE)
  x = setMarkers(x, marker(x, "6" = 1, "10" = 2))
  # plot(x,1)
  y2 = mSim(x, partialmarker = 1)
  expect_identical(genotype(y2, id = 2, marker = 1), c("1", "2")) # Forced!

  y3 = mSim(x, partialmarker = 1, loopBreaker = "5")
  expect_identical(genotype(y3, id = 2, marker = 1), c("1", "2")) # Forced!
})

test_that("markerSim() works in looped pedigree 3", {
  x = cousinPed(0, child = T)
  x = relabel(x, letters[1:5])
  x = setMarkers(x, marker(x, c = 1:2))

  y1 = mSim(x, partial = 1, loopBreaker = "c", seed = 1234)
  expect_equal(as.numeric(getAlleles(y1)), c(1,2,1,2,1,2,2,2,2,2))

  y2 = mSim(x, partial = 1, loopBreaker = "d", seed = 1234)
  expect_equal(as.numeric(getAlleles(y2)), c(2,1,1,1,1,2,1,2,2,2))
})

