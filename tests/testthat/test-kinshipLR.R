
test_that("kinshipLR() catches input errors", {
  s = singleton(1)
  expect_error(kinshipLR(s), "The input must contain at least two pedigrees")
  expect_error(kinshipLR(list(s)), "The input must contain at least two pedigrees")
  expect_error(kinshipLR(list(s, 1)), "The input is not a list of pedigrees")
  expect_error(kinshipLR(list(s, s)), "None of the pedigrees")
  expect_error(kinshipLR(s, s), "None of the pedigrees")
  expect_error(kinshipLR(s, s, source = 1), "The source pedigree has no attached markers")
  expect_error(kinshipLR(list(s, s), source = 1), "The source pedigree has no attached markers")

  s1 = addMarker(s)
  expect_error(kinshipLR(s1, s1, ref = 3), "Invalid value for `ref`")
  expect_error(kinshipLR(s1, s1, ref = 0), "Invalid value for `ref`")
  expect_error(kinshipLR(s1, s1, ref = "unrel"), "Invalid value for `ref`")

  expect_error(kinshipLR(a=s, a=s1), "Duplicated hypothesis name")
  expect_error(kinshipLR(s, H1=s1), "Duplicated hypothesis name")

  s2 = addMarker(s1)
  expect_error(kinshipLR(s1, s2), "The pedigrees have different")
})

test_that("kinshipLR() computes correctly in paternity case", {
  H1 = nuclearPed(fa = "fa", child = "ch") |>
    addMarker(name = "M1", fa = "A/A", ch = "A/A", afreq = c(A=0.05, B=0.95)) |>
    addMarker(name = "M2", fa = "A/A", ch = "A/A", afreq = c(A=0.1, B=0.9))
  H2 = singletons(c("fa", "ch"))

  expect_equal(kinshipLR(H1, H2, source = 1, verbose = F)$LRtotal[[1]], 20*10)
  expect_equal(kinshipLR(H1, H2, markers = "M1", verbose = F)$LRtotal[[1]], 20)
  expect_equal(kinshipLR(H1, H2, markers = -1, source = 1, verbose = F)$LRtotal[[1]], 10)
})


test_that("kinshipLR() returns complete unlinked output", {
  Hpat = nuclearPed(fa = "fa", child = "ch") |>
    addMarker(name = "M1", fa = "A/A", ch = "A/A", afreq = c(A = .05, B = .95)) |>
    addMarker(name = "M2", fa = "A/A", ch = "A/A", afreq = c(A = .10, B = .90))

  Hunrel = singletons(c("fa", "ch"))
  res = kinshipLR(Hpat = Hpat, Hunrel = Hunrel, source = "Hpat", verbose = FALSE)

  expect_is(res, "LRresult")
  expect_identical(names(res), c("LRtotal", "LRperMarker", "likelihoodsPerMarker", "time"))
  expect_identical(names(res$LRtotal), c("Hpat:Hunrel", "Hunrel:Hunrel"))
  expect_identical(rownames(res$LRperMarker), c("M1", "M2"))
  expect_identical(colnames(res$LRperMarker), names(res$LRtotal))
  expect_identical(rownames(res$likelihoodsPerMarker), c("M1", "M2"))
  expect_identical(colnames(res$likelihoodsPerMarker), c("Hpat", "Hunrel"))

  expect_equal(res$LRperMarker[, "Hpat:Hunrel"], c(M1 = 20, M2 = 10))
  expect_equal(res$LRperMarker[, "Hunrel:Hunrel"], c(M1 = 1, M2 = 1))
  expect_equal(res$LRtotal, c(`Hpat:Hunrel` = 200, `Hunrel:Hunrel` = 1))
  expect_equal(apply(res$LRperMarker, 2, prod), res$LRtotal)
  expect_equal(res$likelihoodsPerMarker[, "Hpat"], pedprobr::likelihood(Hpat), check.names = F)
})


test_that("kinshipLR() keeps undefined impossible-data ratios", {
  ok = singletons(1:3) |>
    addMarker("1" = "1/1", "2" = "1/1", "3" = "2/2", alleles = 1:2)
  bad = nuclearPed()

  res = kinshipLR(bad1 = bad, ok = ok, bad2 = bad,
                  ref = "bad2", source = "ok", verbose = FALSE)

  expect_equal(res$likelihoodsPerMarker[1, c("bad1", "bad2")], c(bad1 = 0, bad2 = 0))
  expect_gt(res$likelihoodsPerMarker[1, "ok"], 0)

  expect_true(is.nan(res$LRperMarker[1, "bad1:bad2"]))
  expect_true(is.infinite(res$LRperMarker[1, "ok:bad2"]))
  expect_true(is.nan(res$LRperMarker[1, "bad2:bad2"]))

  expect_true(is.nan(res$LRtotal["bad1:bad2"]))
  expect_true(is.infinite(res$LRtotal["ok:bad2"]))
  expect_true(is.nan(res$LRtotal["bad2:bad2"]))
})


test_that("quickLR() keeps LRresult layout", {
  x = nuclearPed(fa = "fa", child = "ch", sex = 1) |>
    addMarker(fa = "1/1", ch = "1/1", afreq = c("1" = .1, "2" = .9))

  r = quickLR(x, ids = c("fa", "ch"), test = c("pat", "sib"))

  expect_is(r, "LRresult")
  expect_identical(names(r$LRtotal), c("pat:unrel", "sib:unrel"))
  expect_identical(colnames(r$LRperMarker), names(r$LRtotal))
  expect_equal(r$LRperMarker[1, "pat:unrel"], 10)
})
