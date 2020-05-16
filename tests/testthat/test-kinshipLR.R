context("LR calculations")

test_that("kinshipLR() catches input errors", {
  s = singleton(1)
  expect_error(kinshipLR(list(s), ref = 1), "The input must contain at least two pedigrees")
  expect_error(kinshipLR(list(s, s), ref = 1), "None of the pedigrees")
})

test_that("kinshipLR() computes correctly in paternity case", {
  H1 = nuclearPed(fa = "fa", child = "ch")
  H2 = list(singleton("fa"), singleton("ch"))

  m = marker(H1, fa = "A/A", ch = "A/A", afreq = c(A=0.05, B=0.95))
  H1 = setMarkers(H1, m)

  lr = kinshipLR(list(H1, H2))$LRtotal[[1]]
  expect_equal(lr, 20)
})
