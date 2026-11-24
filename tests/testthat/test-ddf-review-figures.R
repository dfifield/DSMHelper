# Per-season fit, review figures and term-drop table (Atlantic-DSM #168).

# A binned ship-like data set with two seasons whose detection differs:
# Spring birds are seen less far out, so a model without Season misfits both.
sim.binned <- function(seed = 168, n = 1500) {
  set.seed(seed)
  season <- sample(c("Fall", "Spring"), n, TRUE)
  sd <- ifelse(season == "Fall", 0.12, 0.06)
  x <- abs(stats::rnorm(n, sd = sd))
  keep <- x <= 0.3
  cuts <- c(0, 0.05, 0.1, 0.2, 0.3)
  bin <- cut(x[keep], cuts, include.lowest = TRUE, labels = FALSE)
  data.frame(distbegin = cuts[bin], distend = cuts[bin + 1], object = seq_len(sum(keep)),
             size = 1, Season = factor(season[keep]))
}
fit.ds <- function(dat, formula = ~1, key = "hn") {
  suppressWarnings(suppressMessages(
    Distance::ds(dat, key = key, formula = formula, adjustment = NULL,
                 truncation = 0.3, quiet = TRUE)))
}

test_that("ddf.season.fit's expected counts add up to ddf.gof()'s, split by season", {
  skip_if_not_installed("Distance")
  dat <- sim.binned()
  for (f in list(~1, ~Season)) {
    m <- fit.ds(dat, f)
    sf <- ddf.season.fit(m$ddf)
    g <- mrds::ddf.gof(m$ddf, qq = FALSE)$chisquare$chi1
    expect_equal(unname(colSums(sf$expected)), as.numeric(g$expected), tolerance = 1e-6)
    expect_equal(unname(colSums(sf$obs)), as.numeric(g$observed))
    expect_equal(sf$seasons, c("Fall", "Spring"))
    expect_equal(sum(sf$n), nrow(dat))
    # each season's expected counts sum to its n (probabilities given detection)
    expect_equal(unname(rowSums(sf$expected)), unname(sf$n), tolerance = 1e-6)
    # 2 seasons x 3 free cells, minus the parameters
    expect_equal(sf$df, 2 * 3 - length(m$ddf$par))
    # the curve is g, so it starts near 1
    expect_true(all(abs(sf$curve[, 1] - 1) < 0.02))
  }
  # Without Season the model misfits both seasons, with it it does not:
  # the per-season chi-square says so where the pooled one cannot.
  chi0 <- ddf.season.fit(fit.ds(dat, ~1)$ddf)$chi_total
  chi1 <- ddf.season.fit(fit.ds(dat, ~Season)$ddf)$chi_total
  expect_gt(chi0, 5 * chi1)
})

test_that("ddf.season.fit's season curve averages to the expected bars", {
  skip_if_not_installed("Distance")
  sf <- ddf.season.fit(fit.ds(sim.binned(), ~Season, key = "hr")$ddf, curve.n = 301)
  br <- sf$breaks; W <- max(br)
  # mean of the curve over a bin (trapezoid) equals that bin's expected bar
  # on the g scale, expected / width * p * W / n
  for (s in 1:2) {
    for (j in 1:4) {
      in.bin <- sf$curve_x >= br[j] & sf$curve_x <= br[j + 1]
      xx <- sf$curve_x[in.bin]; yy <- sf$curve[s, in.bin]
      mean.g <- sum(diff(xx) * (head(yy, -1) + tail(yy, -1)) / 2) / diff(range(xx))
      bar <- sf$expected[s, j] / (br[j + 1] - br[j]) * sf$p[s] * W / sf$n[s]
      expect_equal(mean.g, unname(bar), tolerance = 0.01)
    }
  }
})

test_that("ddf.review.figures writes the PNGs, caches, and returns the statistics", {
  skip_if_not_installed("Distance")
  res <- withr::local_tempdir()
  folder <- file.path(res, "BLKI", "DF Summaries_ECSAS_Ship_W_D")
  dir.create(folder, recursive = TRUE)
  dat <- sim.binned(n = 600)
  model <- fit.ds(dat, ~Season)
  rd <- file.path(folder, "AIC_100.0_hn.Season.RData")
  save(model, file = rd)
  cands <- tibble::tibble(species = "BLKI", ddftype = "ECSAS_Ship_W_D",
                          model = c("hn.Season", "hr"), failed = c(FALSE, TRUE),
                          rdata = c(rd, NA))
  out.dir <- file.path(res, "ddf_review_files")
  cache <- file.path(res, "figs.rds")
  f1 <- ddf.review.figures(cands, out.dir, cache, workers = 1)
  expect_equal(nrow(f1), 1)
  expect_equal(f1$fig_season, "ddf_review_files/BLKI/ECSAS_Ship_W_D/hn.Season_season.png")
  expect_true(file.exists(file.path(res, f1$fig_season)))
  expect_true(file.exists(file.path(res, f1$fig_df)))
  expect_true(file.exists(file.path(res, f1$fig_p)))
  expect_equal(names(f1$season_p_num[[1]]), c("Fall", "Spring"))
  expect_equal(f1$season_df, 2 * 3 - 2)
  expect_equal(f1$fig_note, "")
  # second call: everything cached, nothing redrawn
  expect_no_message(f2 <- ddf.review.figures(cands, out.dir, cache, workers = 1))
  expect_equal(f2$season_chi, f1$season_chi)
})

test_that("ddf.term.drop compares each covariate candidate with its nested ones", {
  cands <- tibble::tibble(
    species = "BLKI", ddftype = "ECSAS_Ship_W_D",
    model = c("hr", "hr.Season", "hr.size", "hr.Season.size", "hn.Season", "hr.cos"),
    key = c("hr", "hr", "hr", "hr", "hn", "hr"),
    formula = c("~1", "~ Season", "~ size", "~ Season + size", "~ Season", "~1"),
    adj = c(NA, NA, NA, NA, NA, "cos"), failed = FALSE,
    AIC = c(120, 105, 118, 100, 110, 115), BIC = c(124, 117, 126, 116, 122, 123),
    p = c(0.40, 0.42, 0.41, 0.421, 0.39, 0.38), n = 1000,
    season_p_num = list(c(Fall = 0.4, Spring = 0.4), c(Fall = 0.5, Spring = 0.34),
                    c(Fall = 0.41, Spring = 0.41), c(Fall = 0.501, Spring = 0.341),
                    c(Fall = 0.45, Spring = 0.33), NULL))
  out <- ddf.term.drop(cands)
  d <- setNames(out$drop_text, out$model)
  m <- setNames(out$drop_min_change, out$model)
  # hr.Season.size: dropping size barely moves p (0.421 -> 0.42), dropping
  # Season moves the seasons a lot (0.501 -> 0.41)
  expect_match(d[["hr.Season.size"]], "&minus;size: &Delta;AIC \\+5.0, &Delta;BIC \\+1.0, p -0.2%")
  expect_match(d[["hr.Season.size"]], "&minus;Season: &Delta;AIC \\+18.0")
  expect_lt(m[["hr.Season.size"]], 0.01)
  # hr.Season: its only nested candidate is hr (not hn.Season or hr.cos)
  expect_match(d[["hr.Season"]], "&minus;Season: &Delta;AIC \\+15.0")
  expect_gt(m[["hr.Season"]], 0.1)
  expect_equal(d[["hr"]], "")
  expect_true(is.na(m[["hr.cos"]]))
  # flag.ddf.candidates turns the small change into term-no-p
  ch <- data.frame(species = "BLKI", ddftype = "ECSAS_Ship_W_D", model = "hr.Season", notes = "")
  out$delta_aic <- out$AIC - min(out$AIC)
  f <- setNames(flag.ddf.candidates(out, ch)$candidates$flags, out$model)
  expect_match(f[["hr.Season.size"]], "term-no-p")
  expect_false(grepl("term-no-p", f[["hr.Season"]]))
})

test_that("ddf.term.drop reports birds and moved from per-detection p", {
  # Three detections: flocks of 1, 1 and 10. Without size every p is 0.5; with
  # size the big flock's p rises to 0.8 and the singles' fall to 0.4. With
  # size the birds are N = 2.5 + 2.5 + 12.5 = 17.5, without N = 2 + 2 + 20 = 24.
  obs.with <- data.frame(object = 1:3, size = c(1, 1, 10), p = c(0.4, 0.4, 0.8))
  obs.wo   <- data.frame(object = 3:1, size = c(10, 1, 1), p = c(0.5, 0.5, 0.5))  # shuffled order
  cands <- tibble::tibble(
    species = "NOFU", ddftype = "ECSAS_Ship_W_D",
    model = c("hr.size", "hr"), key = "hr", formula = c("~ size", "~1"),
    adj = NA, failed = FALSE, AIC = c(100, 110), BIC = c(104, 112), p = c(0.48, 0.5),
    obs_p = list(obs.with, obs.wo))
  out <- ddf.term.drop(cands)
  d <- out$drop_text[out$model == "hr.size"]
  # dropping size: birds 24 / 17.5 - 1 = +37.1%
  expect_match(d, "birds \\+37.1%")
  # shares with size: 1/7, 1/7, 5/7; without: 1/12, 1/12, 10/12.
  # moved = 0.5 * (2 * (1/7 - 1/12) + (10/12 - 5/7)) = 0.5 * (10/84 + 10/84) = 10/84
  expect_match(d, sprintf("moved %.1f%%", 100 * 10 / 84))
  expect_equal(out$drop_min_change[out$model == "hr.size"], 24 / 17.5 - 1)
  # a mismatched pair (different detections) gives no birds measures, not a guess
  cands$obs_p[[2]] <- data.frame(object = 4:6, size = 1, p = 0.5)
  expect_false(grepl("birds", ddf.term.drop(cands)$drop_text[1]))
  # without obs_p the text is as before
  expect_false(grepl("birds", ddf.term.drop(cands[, names(cands) != "obs_p"])$drop_text[1]))
})

test_that("ddf.obs.p reads each detection's object, size and p", {
  skip_if_not_installed("Distance")
  set.seed(1)
  n <- 300
  dat <- data.frame(object = seq_len(n), size = sample(1:5, n, replace = TRUE),
                    distance = abs(rnorm(n, 0, 0.12)))
  dat <- dat[dat$distance < 0.3, ]
  model <- suppressWarnings(suppressMessages(
    Distance::ds(dat, key = "hn", formula = ~ size, truncation = 0.3, quiet = TRUE)))
  td <- tempfile("obsp"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  f <- file.path(td, "AIC_1_hn.size.RData")
  save(model, file = f)
  cands <- tibble::tibble(species = "X", ddftype = "Ship_W_D", model = "hn.size",
                          failed = FALSE, adj = NA, rdata = f)
  op <- ddf.obs.p(cands, file.path(td, "cache.rds"), workers = 1)
  o <- op$obs_p[[1]]
  expect_equal(o$object, model$ddf$data$object)
  expect_equal(o$p, as.numeric(stats::fitted(model$ddf)))
  # second call comes from the cache
  expect_true(file.exists(file.path(td, "cache.rds")))
  expect_equal(ddf.obs.p(cands, file.path(td, "cache.rds"), workers = 1)$obs_p[[1]], o)
})

test_that("shape>3 counts shape parameters only, gof-df0 and gof-season fire", {
  cands <- tibble::tibble(
    species = "BLKI", ddftype = "ECSAS_Ship_W_D",
    model = c("hr.Season.size.Windspeed", "hn.cos", "hr"),
    failed = FALSE, AIC = c(100, 110, 120), delta_aic = c(0, 10, 20), p = c(0.4, 0.4, 0.4),
    # 7 parameters, 5 of them covariate coefficients; hn + 2 adjustments = 3
    n_par = c(7L, 3L, 2L), n_cov = c(5L, 0L, 0L),
    obs = list(c(10, 8, 6, 2), c(10, 8, 6, 2), c(10, 8, 6, 2)),
    season_gof_p = c(0.001, NA, 0.5), n = 300)
  ch <- data.frame(species = "BLKI", ddftype = "ECSAS_Ship_W_D", model = "hr", notes = "")
  f <- setNames(flag.ddf.candidates(cands, ch)$candidates$flags, cands$model)
  expect_false(grepl("shape>3", f[["hr.Season.size.Windspeed"]]))
  expect_match(f[["hr.Season.size.Windspeed"]], "gof-df0")     # 4 - 1 - 7 < 0
  expect_match(f[["hr.Season.size.Windspeed"]], "gof-season")
  expect_match(f[["hn.cos"]], "gof-df0")                        # saturated: 4 - 1 - 3 = 0
  expect_false(grepl("gof-df0", f[["hr"]]))
  expect_false(grepl("gof-season", f[["hr"]]))
})
