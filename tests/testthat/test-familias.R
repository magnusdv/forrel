context("Familias conversion")

test_that("connectedComponents() works", {
  w = data.frame(id = 1:7, fid = c(0,0,1,0,0,0,0), mid = c(0,0,2,0,0,7,0))
  w = w[sample(1:7), ]
  comps = connectedComponents(w$id, w$fid, w$mid)
  comps_sorted = comps[order(sapply(comps, '[[', 1))]
  comps_sorted = lapply(comps_sorted, sort)
  expect_equal(comps_sorted, list(1:3, 4, 5, 6:7))
})

test_that("Familias2ped() converts fullsib pedigree", {

  famped = list(id = c(1,3,5,2,4,6,7),
                findex = c(0, 1, 2, 0, 1, 2, 3),
                mindex = c(0, 4, 5, 0, 4, 5, 6),
                sex = c("male", "male", "male", "female", "female", "female", "female"))
  class(famped) = "FamiliasPedigree"

  ped = Familias2ped(famped, NULL, NULL)
  expect_identical(ped, reorderPed(addChildren(fullSibMating(1),5,6,1,sex=2), famped$id))

})

test_that("Familias2ped() converts a list of singletons", {

  famped = list(id = 3:1, findex = c(0,0,0),  mindex = c(0,0,0),  sex = c("male", "male", "female"))
  class(famped) = "FamiliasPedigree"

  ped = Familias2ped(famped, NULL, NULL)
  expect_equivalent(ped, list('_comp1' = singleton(3, famid='_comp1'),
                              '_comp2' = singleton(2, famid='_comp2'),
                              '_comp3' = singleton(1, sex=2, famid='_comp3')))

  datamatrix = data.frame(m1.1 = c(1,"A",NA), m1.2 = c(1,"B",1), row.names = 1:3)
  expect_error(Familias2ped(famped, datamatrix, NULL), NA)
})

test_that("Familias2ped() converts a single singleton", {

  famped = list(id = 1, findex = 0,  mindex = 0,  sex = "male")
  class(famped) = "FamiliasPedigree"

  datamatrix = data.frame(m1.1 = "A", m1.2 = "B", row.names = 1)
  ped = Familias2ped(famped, datamatrix, NULL)

  true = singleton(1)
  true = addMarkers(true, marker(true, name = "m1", `1` = c("A", "B")))

  expect_identical(ped, true)
})

test_that("Familias2ped() reverses pedlikCompare:::ped2Familias()", {
  skip_on_cran()
  skip_if_not_installed("Familias")
  p2f = pedlikCompare:::ped2Familias
  f2p = function(x) Familias2ped(x$pedigree, datamatrix = x$datamatrix, loci = x$loci)

  x = nuclearPed(1)
  x = markerSim(x, N=1, alleles=1:2, verbose=F)
  name(x, 1) = "m1"

  y = f2p(p2f(x)) # convert back and forth
  expect_identical(x, y)


  # Random ped with 3 founders and 12 matings; simulate 2 markers
  x = randomPed(12, 3, seed = 123)
  x = markerSim(x, N=2, alleles=1:3, afreq=c(0.1, 0.4, 0.5), verbose=F)
  name(x, 1:2) = paste0("m", 1:2)

  # convert back and forth
  y = f2p(p2f(x))

  if(is.ped(y)) # don't test if disconnected (it will fail)
    expect_identical(x, y)
})

test_that("readFamiliasLoci() works", {
  skip_on_cran()
  skip_if_not_installed("Familias")
  loc = Familias::FamiliasLocus(frequencies = c(0.5, 0.5),
                                MutationModel = "equal",
                                MutationRate = 0.1)
  m = readFamiliasLoci(loc)[[1]]
  expect_identical(as.vector(m$mutmod$male),
                             as.vector(loc$maleMutationMatrix))
  expect_identical(dimnames(m$mutmod$male),
                             dimnames(loc$maleMutationMatrix))


})

test_that("Complete Familias2ped example", {
  skip_on_cran()
  skip_if_not_installed("Familias")
  data(NorwegianFrequencies, package = "Familias")
  loci = lapply(NorwegianFrequencies[1:2], Familias::FamiliasLocus)

  ped = Familias::FamiliasPedigree(
    id = c('mother', 'daughter', 'AF'),
    dadid = c(NA, 'AF', NA),
    momid = c(NA, 'mother', NA),
    sex   = c('female', 'female', 'male'))

  datamatrix = data.frame(
    TH01.1 = c(NA, 8, NA),
    TH01.2 = c(NA, 9.3, NA),
    row.names = ped$id)

  expect_identical(Familias2ped(ped, datamatrix, loci[2]),
                   Familias2ped(ped, datamatrix, loci, matchLoci = T))

  datamatrix2 = cbind(datamatrix,
                     D3S1358.1 = c(10, 11, 12),
                     D3S1358.2 = c(12, 11, 10))

  expect_identical(Familias2ped(ped, datamatrix2, loci, matchLoci = T),
                   Familias2ped(ped, datamatrix2, loci[2:1], matchLoci = T))

})

test_that("Another complete Familias example", {
  skip_on_cran()
  skip_if_not_installed("Familias")
  # Familias pedigree which gave bug in Familias2ped
  ped1 = Familias::FamiliasPedigree(id = c("AF", "CH"), dadid = c(NA, "AF"),
                                    momid = c(NA,NA), sex = c("male", "male"))
  pedigrees = list(ped1)

  locus1 = Familias::FamiliasLocus(frequencies=c(0.1, 0.9), name="S", allelenames=c("1","2"))
  loci = list(locus1)

  datamatrix = data.frame(locus1.1 = c("1","2"), locus1.2 = c("2","2"))
  rownames(datamatrix) = ped1$id

  # Convert
  x = Familias2ped(pedigrees, datamatrix, loci)
  x1 = x[[1]]

  # Should equal this:
  y = nuclearPed(fa = "AF", mo = "added_1", child = "CH")
  y = setMarkers(y, marker(y, AF = 1:2, CH = 2, name = "S", afreq = c(`1`=.1, `2`=.9)))
  y = reorderPed(y, c(1,3,2))

  expect_identical(x1, y)
})
