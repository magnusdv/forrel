library(pedprobr)

if(!checkMerlin())
  skip("Merlin not installed")

test_that("kinshipLR() works with two linked markers", {
  x = linearPed(2)
  m = marker(x, geno = c("1/1", NA, "1/2", NA, "1/1"), alleles = 1:10)
  x = setMarkers(x, list(m, m))

  unrel = extractSingletons(x, c(1,3,5))

  rho = 0.25
  map = data.frame(CHROM = 1, MARKER = NA, CM = c(0, pedprobr:::haldane(rho = rho)))
  lr = kinshipLR(x, unrel, linkageMap = map)$LRtotal[[1]]
  lr2 = likelihood2(x, 1, 2, rho = rho)/likelihood2(unrel, 1, 2, rho = rho)
  expect_equal(signif(lr, 3), signif(lr2, 3))
})


test_that("kinshipLR() works with two linked markers on X", {
  x = linearPed(2, sex = 2:1)
  m = marker(x, geno = c("1/1", NA, NA, "1/2", "1/1"), alleles = 1:10, chrom = "X")
  x = setMarkers(x, list(m, m))

  unrel = extractSingletons(x, c(1,4,5))

  rho = 0.25
  map = data.frame(CHROM = "X", MARKER = NA, CM = c(0, pedprobr:::haldane(rho = rho)))
  lr = kinshipLR(x, unrel, linkageMap = map)$LRtotal[[1]]
  lr2 = likelihood2(x, 1, 2, rho = rho)/likelihood2(unrel, 1, 2, rho = rho)
  expect_equal(signif(lr, 2), signif(lr2, 2))
})
