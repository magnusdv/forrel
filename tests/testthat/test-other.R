
test_that("findExclusions() returns marker names only", {
  x = nuclearPed(fa = "fa", mo = "mo", child = "ch") |>
    addMarker(name = "M1", fa = "1/1", mo = "1/1", alleles = 1:2) |>
    addMarker(name = "M2", fa = "1/1", mo = "1/2", alleles = 1:2)

  cand = singleton("cand") |>
    addMarker(name = "M1", cand = "2/2", alleles = 1:2) |>
    addMarker(name = "M2", cand = "1/2", alleles = 1:2)

  expect_identical(findExclusions(x, id = "ch", candidate = cand), "M1")
})


test_that("rankProfiles() returns ranked profile output", {
  x = nuclearPed(1) |>
    addMarker(name = "M", `1` = "1/1", `2` = "2/2", alleles = 1:2)

  r = rankProfiles(x, id = 3)

  expect_identical(names(r), c("profiles", "best", "marginal1", "marginal2"))
  expect_identical(unname(r$best), "1/2")
  expect_equal(r$marginal1, c(`1/2` = 1))
  expect_true(is.na(r$marginal2))
  expect_equal(r$profiles$likelihood, 1)
})

