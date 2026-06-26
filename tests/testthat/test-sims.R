
mSim = function(...) markerSim(..., verbose = FALSE)
pSim = function(...) profileSim(..., verbose = FALSE)
sSim = function(...) simpleSim(..., verbose = FALSE)

test_that("simpleSim() runs in trivial example", {
  x = nuclearPed(1)
  y = sSim(x, N=1, alleles=1:2)
  expect_is(y, "ped")
  expect_equal(nMarkers(y), 1)
  expect_true(all(y$MARKERS[[1]] > 0))
})

test_that("simpleSim() works with explicit mutation args", {
  x = nuclearPed(1)
  mut = structure(c(0, 0, 1, 1), dim = c(2L, 2L), dimnames = list(c("1", "2"), c("1", "2")),
                  model = "custom", class = "mutationMatrix")
  mod = structure(list(female = mut, male = mut), sexEqual = TRUE,
                  alwaysLumpable = TRUE, class = "mutationModel")
  y = sSim(x, N=1, alleles = 1:2, afreq = c(.99, .01), mutmod = mod)
  expect_equal(genotype(y, 1, id=3), c("2", "2"))
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

test_that("profileSim() keeps marker names", {
  x = nuclearPed(1)
  m = marker(x)
  x = setMarkers(x, list(m, m))
  name(x, 1:2) = c("m1", "m2")

  s = pSim(x, N = 1)
  expect_identical(name(s, 1:2), c("m1", "m2"))

  s2 = pSim(x, N = 1, markers = list(m, m))
  expect_identical(name(s2, 1:2), rep(NA_character_, 2))
})


test_that("profileSim() treats pedlists as expected", {
  x = singleton(1) |> addMarker(alleles = 1:5, name = "M")
  y = relabel(x, 2)

  set.seed(777)
  sim_pedlist = pSim(list(x, y), N = 3, markers = "M")

  set.seed(777)
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
  x = linearPed(2) |> addSon(c(2,5)) |>
    addMarker("5" = "1/1", alleles = 1:2, afreq = c(0.001, 0.999))

  y = mSim(x, partialmarker = 1, seed = 123)
  expect_identical(genotype(y, id = 3, marker = 1), c("1", "2"))
  expect_identical(genotype(y, id = 4, marker = 1), c("1", "2"))
})

test_that("markerSim() works in looped pedigree 2", {
  x = linearPed(2) |> addSon(c(2,5)) |> addMarker("6" = "1/1", "2" = "0/2")
  # plot(x, mark =1)

  y1 = mSim(x, partialmarker = 1)
  expect_identical(genotype(y1, id = 2, marker = 1), c("1", "2")) # Forced

  x2 = x |> addParents(2, 10, 11, verbose = FALSE) |> addMarker("6" = 1, "10" = 2)
  # plot(x2, mark = 2)

  y2 = mSim(x2, partialmarker = 2)
  expect_identical(genotype(y2, id = 2, marker = 1), c("1", "2")) # Forced!

  y3 = mSim(x2, partialmarker = 2, loopBreaker = "5")
  expect_identical(genotype(y3, id = 2, marker = 1), c("1", "2")) # Forced!
})

test_that("markerSim() works in looped pedigree 3", {
  x = cousinPed(0, child = T) |> relabel(new = letters[1:5]) |> addMarker(c = 1:2)

  y1 = mSim(x, partial = 1, loopBreaker = "c", seed = 1234)
  expect_equal(as.numeric(getAlleles(y1)), c(1,2,1,2,1,2,2,2,2,2))

  y2 = mSim(x, partial = 1, loopBreaker = "d", seed = 1234)
  expect_equal(as.numeric(getAlleles(y2)), c(2,1,1,1,1,2,1,2,2,2))
})


# Added regression tests 2026 -----------------------------------------------------------------


test_that("markerSim() indexes founder inbreeding after subsetting", {
  x = halfSibPed()
  founderInbreeding(x, 1) = 1
  m = marker(x, `4` = "1/1", alleles = 1:2)

  y = mSim(x, N = 40, ids = 3, partialmarker = m, seed = 1)
  geno = as.matrix(as.data.frame(y)[internalID(y, 3), -(1:4)])

  expect_true(any(geno == "1/2"))
})

test_that("markerSim() supports founder loop breakers", {

  # Pure marriage loop, can only be broken by founder
  x = quadHalfFirstCousins() |> removeIndividuals(9:10, verbose = F) |>
    addMarker(`5` = "1/1", alleles = 1:2, afreq = c(0.0001, 0.9999))

  # plot(x, marker = 1, arr = T)

  y = mSim(x, N = 3, partialmarker = 1, seed = 2)

  # Both parents of 5 should carry allele "1"
  g = getGenotypes(y, parents(y, 5)) |> as.character()
  expect_all_true(g == "1/2")
})

test_that("markerSim() on pedlists uses all ids by default", {
  x = singletons(1:2)
  y = mSim(x, N = 3, alleles = 1:2, seed = 3)
  z = mSim(x, N = 3, ids = 2, alleles = 1:2, seed = 3)

  expect_equal(nMarkers(y), 3)
  expect_equal(typedMembers(y), c("1", "2"))
  expect_equal(typedMembers(z), "2")
})

test_that("markerSim() on pedlists forwards mutation arguments", {
  x = list(nuclearPed(1), relabel(nuclearPed(1), 4:6))

  y = mSim(x, N = 2, alleles = 1:2, afreq = c(1, 0), mutmod = "equal", rate = 1, seed = 4)
  expect_false(is.null(mutmod(y[[1]],1)))
  expect_false(is.null(mutmod(y[[2]],1)))

  childGenos = as.character(getGenotypes(y, c(3,6)))
  expect_all_true(childGenos == "2/2")
})

test_that("profileSim() returns N outputs when markers are empty", {
  x = nuclearPed(1) |> addMarker(alleles = 1:2)

  y = suppressMessages(pSim(x, N = 3, markers = integer(0)))

  expect_equal(length(y), 3)
  expect_identical(y, list(x, x, x))
})

test_that("profileSim() forwards loopBreakers", {
  x = halfCousinPed(1, child = TRUE) |>
    addMarker(`10` = "1/1", alleles = 1:2)

  a = pSim(x, loopBreakers = 2, seed = 1)
  b = pSim(x, loopBreakers = 4, seed = 1)
  expect_true(!identical(a, b))
})

test_that(".mutate() respects source-specific mutation rows", {
  M = rbind(c(0, 1, 0),
            c(0, 0, 1),
            c(1, 0, 0))

  expect_identical(.mutate(c(1L, 2L, 3L, 1L, 2L), M),
                   c(2L, 3L, 1L, 2L, 3L))
})

test_that("markerSim() applies mutation model in gene drops", {
  x = nuclearPed(1)

  y = mSim(x, N = 3, alleles = 1:2, afreq = c(1, 0), mutmod = "equal", rate = 1, seed = 1)

  expect_equal(as.character(getGenotypes(y, ids = 3)), rep("2/2", 3))
})
