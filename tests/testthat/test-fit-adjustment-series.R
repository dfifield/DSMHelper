test_that("fit.adjustment.series fits every order, keeps nested AICs consistent and returns the AIC-best", {
  skip_if_not_installed("Distance")
  set.seed(113)
  # Binned half-normal distances, like ECSAS ship data
  x <- abs(rnorm(600, sd = 0.1))
  x <- x[x <= 0.3]
  cuts <- c(0, 0.05, 0.1, 0.2, 0.3)
  bin <- cut(x, cuts, include.lowest = TRUE, labels = FALSE)
  dat <- data.frame(distbegin = cuts[bin], distend = cuts[bin + 1],
                    object = seq_along(x), size = 1)

  m <- fit.adjustment.series(dat, "hn", "cos", max_adjustments = 2,
                             ds.args = list(truncation = 0.3))
  s <- attr(m, "adj.search")
  expect_equal(s$nadj, 0:2)
  # Each order nests the one before, so at its optimum AIC rises by at most 2
  expect_true(all(diff(s$AIC) <= 2.1))
  expect_true(all(s$flag == ""))
  expect_equal(m$ddf$criterion, min(s$AIC), tolerance = 1e-3)

  # The per-key cap and a fixed nadj
  m1 <- fit.adjustment.series(dat, "hn", "cos", max_adjustments = 2,
                              max.adj.by.key = c(hn = 1),
                              ds.args = list(truncation = 0.3))
  expect_equal(attr(m1, "adj.search")$nadj, 0:1)
  m2 <- fit.adjustment.series(dat, "hn", "cos", max_adjustments = 2, nadj = 1,
                              ds.args = list(truncation = 0.3))
  expect_equal(length(m2$ddf$ds$aux$ddfobj$adjustment$order), 1)
})
