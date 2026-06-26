test_that("ibdEstimate() returns complete kappa output", {
  x = nuclearPed(2) |>
    markerSim(N = 6, ids = 3:4, alleles = 1:3, seed = 10, verbose = FALSE)

  k = ibdEstimate(x, ids = 3:4, maxval = TRUE, verbose = FALSE)
  co = unlist(as.data.frame(k)[paste0("k", 0:2)], use.names = FALSE)

  expect_is(k, "ibdEst")
  expect_identical(names(k), c("id1", "id2", "N", "k0", "k1", "k2", "maxloglik"))
  expect_equal(nrow(k), 1)
  expect_identical(as.character(k$id1), "3")
  expect_identical(as.character(k$id2), "4")
  expect_equal(k$N, 6)
  expect_true(all(is.finite(co)))
  expect_true(all(co >= -1e-8))
  expect_equal(sum(co), 1, tolerance = 1e-8)
  expect_true(is.finite(k$maxloglik))
})

test_that("ibdEstimate() returns one row per requested pair", {
  x = nuclearPed(3) |>
    markerSim(N = 4, ids = 3:5, alleles = 1:3, seed = 11, verbose = FALSE)

  ids = matrix(c(3, 4, 3, 5), ncol = 2, byrow = TRUE)
  k = ibdEstimate(x, ids = ids, verbose = FALSE)
  df = as.data.frame(k)
  co = as.matrix(df[paste0("k", 0:2)])

  expect_equal(nrow(k), 2)
  expect_identical(as.character(df$id1), c("3", "3"))
  expect_identical(as.character(df$id2), c("4", "5"))
  expect_equal(df$N, c(4, 4))
  expect_equal(rowSums(co), c(1, 1), tolerance = 1e-8)
})

test_that("ibdEstimate() returns complete delta output", {
  x = fullSibMating(1) |>
    markerSim(N = 4, ids = 5:6, alleles = 1:3, seed = 1, verbose = FALSE)

  d = ibdEstimate(x, ids = 5:6, param = "delta", verbose = FALSE)
  co = unlist(as.data.frame(d)[paste0("d", 1:9)], use.names = FALSE)

  expect_is(d, "ibdEst")
  expect_identical(names(d), c("id1", "id2", "N", paste0("d", 1:9)))
  expect_equal(nrow(d), 1)
  expect_equal(d$N, 4)
  expect_true(all(is.finite(co)))
  expect_true(all(co >= -1e-8))
  expect_equal(sum(co), 1, tolerance = 1e-8)
})
