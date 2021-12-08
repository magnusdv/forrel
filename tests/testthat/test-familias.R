
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

  true = singleton(1) |> addMarker(name = "m1", "1" = "A/B")

  expect_identical(ped, true)
})

