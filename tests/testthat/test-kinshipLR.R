
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
  expect_error(kinshipLR(s1, s2), "When `markers = NULL`, all pedigrees must have the same number of attached markers")
})

test_that("kinshipLR() computes correctly in paternity case", {
  H1 = nuclearPed(fa = "fa", child = "ch") |>
    addMarker(fa = "A/A", ch = "A/A", afreq = c(A=0.05, B=0.95))

  H2 = list(singleton("fa"), singleton("ch"))

  lr = kinshipLR(H1, H2)$LRtotal[[1]]
  expect_equal(lr, 20)
})

