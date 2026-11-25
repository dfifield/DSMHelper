# Cross-validation helpers for the sdmTMB track (Atlantic-DSM #130, step 05.00).
# Pure functions: no sdmTMB needed.

test_that("the fishing term is appended once and keeps the model", {
  f <- add.sdmtmb.fishing.term(abundance.est ~ 0 + Season + platform + s(depth, by = Season, bs = "ts"))
  expect_true("fishing" %in% all.vars(f))
  expect_true("Season:fishing" %in% attr(stats::terms(f), "term.labels"))
  # No global intercept before or after.
  expect_identical(attr(stats::terms(f), "intercept"), 0L)
  expect_true('s(depth, by = Season, bs = "ts")' %in% attr(stats::terms(f), "term.labels"))
  expect_error(add.sdmtmb.fishing.term(f), "already uses")
})

test_that("folds merge evenly and every row keeps one fold", {
  expect_identical(combine.cv.folds(1:10, 5), c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L))
  f <- rep(1:10, times = 3)
  m <- combine.cv.folds(f, 5)
  expect_length(m, length(f))
  # Rows in the same original fold stay together.
  expect_true(all(tapply(m, f, function(z) length(unique(z))) == 1))
  expect_identical(combine.cv.folds(1:10, 10), 1:10)
  expect_error(combine.cv.folds(1:10, 3), "evenly")
  expect_error(combine.cv.folds(c(1, 2, 4), 2), "every value")
})

test_that("ranking follows 02.05's two-margin rule", {
  r <- rank.cv.candidates(data.frame(candidate = c("tw", "dl", "dg"), score = c(0.50, 0.40, 0.60)),
                          tol = 0.05, abs.tol = 0.05)
  expect_identical(r$best, "dl")
  expect_identical(r$runner_up, "tw")
  expect_true(r$decisive)                     # 20% and 0.10 apart
  # Clears the relative margin but not the absolute one.
  r2 <- rank.cv.candidates(data.frame(candidate = c("a", "b"), score = c(0.02, 0.04)), 0.05, 0.05)
  expect_false(r2$decisive)
  # A non-finite score cannot win.
  r3 <- rank.cv.candidates(data.frame(candidate = c("a", "b"), score = c(NaN, 0.9)), 0.05, 0.05)
  expect_identical(r3$best, "b")
  # One candidate: nothing compared.
  expect_false(rank.cv.candidates(data.frame(candidate = "a", score = 1), 0.05, 0.05)$decisive)
})

test_that("the range is the longest within tolerance of the best", {
  s <- data.frame(range_km = c(75, 150, 300), score = c(0.60, 0.40, 0.41))
  r <- choose.sdmtmb.range(s, tol = 0.05, abs.tol = 0.05, extend = c(50, 450))
  expect_identical(r$best_range_km, 150)
  expect_identical(r$range_km, 300)           # tied with 150, and longer
  expect_identical(r$tied, "150, 300")
  expect_false(r$decisive)
  # 300 is the top of the ladder and 450 is untried, so extend.
  expect_true(r$at_edge)
  expect_identical(r$extend_to, 450)
  # Once 450 has been tried and is clearly worse, 300 is no longer at an edge
  # that can be extended.
  s2 <- rbind(s, data.frame(range_km = 450, score = 0.70))
  r2 <- choose.sdmtmb.range(s2, 0.05, 0.05, extend = c(50, 450))
  expect_identical(r2$range_km, 300)
  expect_false(r2$at_edge)
  # A clear winner in the middle: decisive, not at an edge.
  r3 <- choose.sdmtmb.range(data.frame(range_km = c(75, 150, 300), score = c(0.9, 0.4, 0.9)),
                            0.05, 0.05, extend = c(50, 450))
  expect_true(r3$decisive)
  expect_identical(r3$range_km, 150)
  expect_false(r3$at_edge)
  expect_true(is.na(r3$extend_to))
  # The bottom edge extends downwards, to the largest untried value below it.
  r4 <- choose.sdmtmb.range(data.frame(range_km = c(75, 150), score = c(0.3, 0.9)),
                            0.05, 0.05, extend = c(25, 50, 450))
  expect_identical(r4$extend_to, 50)
})
