context("Founder inbreeding")

library(pedtools)

test_that("simpleSim() works on inbred singleton", {
  x = singleton(1)
  founderInbreeding(x, 1) = 1
  y = markerSim(x, N = 5, alleles = 1:4, verbose = F)
  geno = as.matrix(as.data.frame(y)[1, -(1:4)])

  # all should be homozygous
  expect_true(all(geno %in% c("1/1", "2/2", "3/3", "4/4")))
})

test_that("markerSim() does conditional sims with inbred founders", {
  x = nuclearPed(1)
  founderInbreeding(x, 1:2) = 1
  m = marker(x, '3' = 1:2, alleles = 1:4) # heterozygous child

  y = markerSim(x, N=5, partial = m, verbose = F)

  # parental genotypes
  geno = as.matrix(as.data.frame(y)[1:2, -(1:4)])

  # all should be homozygous
  expect_true(all(geno %in% c("1/1", "2/2")))
})
