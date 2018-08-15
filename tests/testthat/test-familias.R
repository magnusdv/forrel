context("Familias conversion")

test_that("connectedComponents() works", {
  w = data.frame(id = 1:7, fid = c(0,0,1,0,0,0,0), mid = c(0,0,2,0,0,7,0))
  w = w[sample(1:7), ]
  comps = connectedComponents(w$id, w$fid, w$mid)
  comps_sorted = comps[order(sapply(comps, '[[', 1))]
  expect_equal(comps_sorted, list(1:3, 4, 5, 6:7))

})

test_that("Familias2ped() converts fullsib pedigree", {

  famped = list(id = c(1,3,5,2,4,6,7),
                findex = c(0, 1, 2, 0, 1, 2, 3),
                mindex = c(0, 4, 5, 0, 4, 5, 6),
                sex = c("male", "male", "male", "female", "female", "female", "female"))
  class(famped) = "FamiliasPedigree"

  ped = Familias2ped(famped, NULL, NULL)
  expect_identical(ped, reorderPed(addChildren(fullSibMating(2),5,6,1,sex=2), famped$id))

})

test_that("Familias2ped() converts singletons", {

  famped = list(id = 3:1, findex = c(0,0,0),  mindex = c(0,0,0),  sex = c("male", "male", "female"))
  class(famped) = "FamiliasPedigree"

  ped = Familias2ped(famped, NULL, NULL)
  expect_identical(ped, list('1' = singleton(3, famid=1),
                             '2' = singleton(2, famid=2),
                             '3' = singleton(1, sex=2, famid=3)))

})

test_that("Familias2ped() reverses pedlikCompare:::ped2Familias()", {

  p2f = pedlikCompare:::ped2Familias
  f2p = function(x) Familias2ped(x$pedigree, datamatrix = x$datamatrix, loci = x$loci)

  # Random ped with 3 founders and 12 matings; simulate 2 markers
  x = randomPed(12,3)
  x = markerSim(x, N=2, alleles=1:3, afreq=c(0.1, 0.4, 0.5), verbose=F)

  # convert back and forth
  y = f2p(p2f(x))

  if(is.ped(y)) # don't test if disconnected (it will fail)
    expect_identical(x, y)
})
