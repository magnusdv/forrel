context("Familias conversion")

test_that("connectedComponents() works", {
  w = data.frame(ID=1:7, FID=c(0,0,1,0,0,0,0), MID=c(0,0,2,0,0,7,0))
  w = w[sample(1:7), ]
  comps = connectedComponents(w$ID, w$FID, w$MID)
  comps_sorted = comps[order(sapply(comps, '[[', 1))]
  expect_equal(comps_sorted, list(1:3,4,5,6:7))

})
