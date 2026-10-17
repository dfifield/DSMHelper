# The surfaces issue #87 is about, as 02.15 measured them on 2026-10-15
# (Atl IMRP) and 2026-09-05 (NL_EXPL_DRL_RA).
test_that("each arm of the rule catches what the other misses", {
  expect_identical(
    is.broken.surface(
      max_over_obs = c(1.31, 21.6, 75.4, 127.8, 500,  NA),
      n_over_limit = c(0,    34,   799,  118,   0,    5),
      ref.mult = 100),
    c(FALSE,              # worst sound surface, Atl IMRP NOGA Winter nb
      TRUE,               # Murres Spring tw: missed by the 100x rule alone
      TRUE,               # NL Petrels Winter nb: missed by the 100x rule alone
      TRUE,               # HERG Winter tw: both arms fire
      TRUE,               # past the ceiling multiple, under the flat limit
      TRUE))              # no observed reference, still over the limit
})

test_that("NA in both arms is not broken, and the result is never NA", {
  out <- is.broken.surface(c(NA, 0.5, NA), c(NA, NA, 0), ref.mult = 100)
  expect_identical(out, c(FALSE, FALSE, FALSE))
})

test_that("zero-length input gives zero-length output", {
  expect_identical(is.broken.surface(numeric(0), numeric(0), 100), logical(0))
})

test_that("mismatched lengths stop", {
  expect_error(is.broken.surface(c(1, 2), 0, 100))
})
