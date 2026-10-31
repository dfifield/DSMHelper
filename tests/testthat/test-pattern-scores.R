test_that("a correctly specified surface scores near 0 whatever its level", {
  set.seed(1)
  mu <- rexp(20000, 0.5)
  y <- rpois(20000, mu)
  right <- score.pattern.calibration(y, mu)
  expect_lt(right$score, 0.05)
  expect_equal(right$total_OE, sum(y) / sum(mu))

  # Wrong level, right shape: still near 0, and the level shows in total_OE.
  scaled <- score.pattern.calibration(y, 1.5 * mu)
  expect_equal(scaled$score, right$score, tolerance = 1e-12)
  expect_equal(scaled$total_OE, right$total_OE / 1.5, tolerance = 1e-12)
})

test_that("an inflated hotspot is penalised, and in the top bin", {
  set.seed(2)
  mu <- rexp(20000, 0.5)
  y <- rpois(20000, mu)
  # Contrast exaggerated: the highest-predicted rows are inflated.
  hot <- mu * ifelse(mu > stats::quantile(mu, 0.9), 2, 1)
  right <- score.pattern.calibration(y, mu)
  bad <- score.pattern.calibration(y, hot)
  expect_gt(bad$score, 5 * right$score)
  expect_lt(bad$bins$OE[bad$bins$bin == 1], 0.8)
  expect_gt(bad$bins$OE[bad$bins$bin == max(bad$bins$bin)], 1.1)
})

test_that("bins carry equal predicted mass and number from the top", {
  pred <- c(50, rep(1, 90))   # the big row carries 36% of the total
  b <- score.pattern.calibration(rep(1, 91), pred, n.bins = 5)$bins
  # The single big row opens bin 1 and holds it alone.
  expect_equal(b$bin[1], 1L)
  expect_equal(b$rows[1], 1L)
  expect_equal(sum(b$rows), 91L)
})

test_that("an empty bin counts as max.ratio, not infinity", {
  # 10 rows predicted high with no birds, 90 low rows carrying all of them.
  pred <- c(rep(9, 10), rep(1, 90))
  obs  <- c(rep(0, 10), rep(2, 90))
  r <- score.pattern.calibration(obs, pred, n.bins = 2, max.ratio = 10)
  expect_equal(r$bins$O[1], 0)
  expect_true(is.finite(r$score))
  # Bin 1 is capped at log(10); bin 2 is off by 2x after rescaling.
  expect_equal(r$score, mean(c(log(10), abs(log(r$bins$OE[2])))))
  expect_lte(r$score, log(10))
})

test_that("degenerate input returns NA rather than failing", {
  expect_true(is.na(score.pattern.calibration(c(0, 0), c(1, 2))$score))
  expect_true(is.na(score.pattern.calibration(numeric(), numeric())$score))
})

test_that("cell ranks follow the surface and ignore sparse cells", {
  set.seed(3)
  n <- 20000
  x <- runif(n, 0, 2e5)
  y <- runif(n, 0, 2e5)
  mu <- exp(x / 5e4)
  obs <- rpois(n, mu)
  good <- score.cell.rank(obs, mu, x, y, rep(1, n), seq_len(n), cell.m = 25000)
  expect_gt(good$rho, 0.9)
  expect_equal(good$n_cells, 64L)
  flat <- score.cell.rank(obs, rev(sort(mu))[rank(x)], x, y, rep(1, n), seq_len(n))
  expect_lt(flat$rho, 0)
  # A season splits the cells.
  s <- score.cell.rank(obs, mu, x, y, rep(1, n), seq_len(n),
                       season = rep(c("a", "b"), n / 2))
  expect_equal(s$n_cells, 128L)
  # min.segments counts distinct segments, not rows.
  few <- score.cell.rank(obs, mu, x, y, rep(1, n), rep(1, n), min.segments = 2)
  expect_equal(few$n_cells, 0L)
  expect_true(is.na(few$rho))
})

test_that("pooling uses the kept folds and returns one row per finalist", {
  set.seed(4)
  n <- 4000
  mk <- function(fam, f) {
    mu <- rexp(n)
    tibble::tibble(species = "SP", family_key = fam, modname = paste0("m_", fam),
                   fold = rep(1:4, each = n / 4), row = seq_len(n),
                   obs = rpois(n, mu), pred = mu * f,
                   x = runif(n, 0, 1e5), y = runif(n, 0, 1e5),
                   area = 1, season = "Fall", seg = as.character(seq_len(n)))
  }
  h <- dplyr::bind_rows(mk("nb", 1), mk("tw", 1))
  p <- pool.family.cv.pattern(h)
  expect_equal(nrow(p), 2L)
  expect_true(all(p$n_folds == 4L))
  expect_true(all(c("OE_bin1", "OE_bin5") %in% names(p)))
  k <- pool.family.cv.pattern(h, keep = tibble::tibble(species = "SP", fold = 1:3))
  expect_true(all(k$n_folds == 3L))
  expect_equal(nrow(pool.family.cv.pattern(h[0, ])), 0L)
})

test_that("provenance records the scoring version, so old CV parts go stale", {
  p <- get.family.cv.provenance(10, 1e5, 200, 1)
  expect_equal(p$scores, "crps+pattern+heldout")
  expect_equal(p$pattern_bins, 5)
  old <- p[setdiff(names(p), c("scores", "pattern_bins", "pattern_cell_m",
                               "pattern_min_segments"))]
  f <- tempfile(); writeLines("x", f)
  expect_match(saved.result.stale.reason(f, f, old, p), "settings differ")
})
