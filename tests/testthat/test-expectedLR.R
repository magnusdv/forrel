test_that("expectedLR() agrees with exact paternity formulas", {
  num = nuclearPed(father = "fa", mother = "mo", child = "ch")
  den = singletons(c("fa", "mo", "ch"), sex = c(1, 2, 1))

  p = 0.9
  m = marker(num, mo = "1/2", afreq = c("1" = p, "2" = 1-p))

  expect_equal(expectedLR(num, den, ids = "ch", marker = m),
               1/(8*p*(1-p)) + 1/2)

  genotype(m, id = "fa") = "1/1"
  expect_equal(expectedLR(num, den, ids = "ch", marker = m),
               1/(8*p*(1-p)) + 1/(4*p^2))
})


test_that("expectedLR() handles zero-probability genotype combinations", {
  H = nuclearPed(fa = "fa", child = "ch")
  m = marker(H, afreq = c("1" = 0.5, "2" = 0.5))

  expect_equal(expectedLR(H, H, ids = c("fa", "ch"), marker = m), 1)
})


test_that("expectedLR() checks sex for X-linked markers", {
  H1 = singleton("A", sex = 1)
  H2 = singleton("A", sex = 1)
  truth = singleton("A", sex = 2)
  m = marker(H1, afreq = c("1" = 0.5, "2" = 0.5), chrom = 23)

  expect_error(expectedLR(H1, H2, truePed = truth, ids = "A", marker = m),
               "sex of `ids` must agree")
})

