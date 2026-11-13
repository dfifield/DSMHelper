# sdmTMB DSM functions (Atlantic-DSM #130), on a small synthetic data set: a
# 200 km square with a land block in the middle, four seasons, two platform
# levels. Each fit takes seconds.

seasons <- c("Fall", "Spring", "Summer", "Winter")

make.sim.data <- function(n = 1200, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 0, 2e5), y = stats::runif(n, 0, 2e5),
                  Season = sample(seasons, n, replace = TRUE),
                  platform = sample(c("WD", "FD"), n, replace = TRUE),
                  segment.area = stats::runif(n, 0.5, 1.5))
  d$depth <- d$x / 1000 + stats::rnorm(n, 0, 5)
  mu <- exp(-1 + 0.01 * d$depth - 0.5 * (d$platform == "FD")) * d$segment.area
  d$abundance.est <- stats::rbinom(n, 1, 0.3) * stats::rgamma(n, shape = 2, rate = 2 / mu)
  d
}
domain <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(
  cbind(c(0, 2e5, 2e5, 0, 0), c(0, 0, 2e5, 2e5, 0))))))
land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(
  cbind(c(80, 120, 120, 80, 80), c(20, 20, 180, 180, 20))))))

test_that("the formula has the season intercepts, platform, ts smooths and the fishing term", {
  f <- sdmtmb.dsm.formula(c("depth", "sst"), fishing = TRUE)
  expect_identical(all.vars(f), c("abundance.est", "Season", "platform", "depth", "sst", "fishing"))
  expect_match(deparse(f, width.cutoff = 500), 's(depth, by = Season, bs = "ts")', fixed = TRUE)
  expect_match(deparse(f, width.cutoff = 500), "Season:fishing", fixed = TRUE)
  expect_false(grepl("fishing", deparse(sdmtmb.dsm.formula(), width.cutoff = 500)))
})

test_that("prepare.sdmtmb.data adds km coordinates, ordered seasons and log area", {
  d <- prepare.sdmtmb.data(make.sim.data(10), seasons)
  expect_equal(d$X, d$x / 1000)
  expect_identical(levels(d$Season), seasons)
  expect_identical(d$season_num, as.integer(d$Season))
  expect_equal(d$log_area, log(d$segment.area))
  bad <- make.sim.data(10)
  bad$Season[1] <- "Autumn"
  expect_error(prepare.sdmtmb.data(bad, seasons), "Autumn")
})

skip_if_not_installed("sdmTMB")
skip_if_not_installed("fmesher")
skip_if_not_installed("sdmTMBextra")

dat <- prepare.sdmtmb.data(make.sim.data(), seasons)
fm <- build.sdmtmb.mesh(dat, domain, max.edge.km = 25)
mesh <- bind.sdmtmb.mesh(dat, fm, land)

test_that("the mesh is in km and the barrier marks the land block", {
  expect_true(max(fm$loc[, 1]) < 1000)
  expect_gt(length(mesh$barrier_triangles), 0)
  expect_null(bind.sdmtmb.mesh(dat, fm)$barrier_triangles)
})

test_that("fix.smooths.off fixes exactly the chosen smooths' SDs and coefficients", {
  params <- list(ln_smooth_sigma = matrix(0.5, 3, 2), b_smooth = matrix(1, 9, 2))
  off <- matrix(c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE), 3, 2)
  sm <- fix.smooths.off(params, c(0L, 3L, 6L), off, list(), list())
  expect_identical(is.na(matrix(as.integer(sm$map$ln_smooth_sigma), 3, 2)), off)
  b_off <- matrix(as.integer(sm$map$b_smooth), 9, 2)
  expect_identical(which(is.na(b_off[, 1])), 4:6)
  expect_true(all(is.na(b_off[, 2])))
  expect_true(all(sm$start$b_smooth[4:6, 1] == 0) && all(sm$start$b_smooth[-(4:6), 1] == 1))
  expect_equal(sm$start$ln_smooth_sigma[2, 1], log(1e-6))
})

fit_tw <- fit.sdmtmb.dsm(sdmtmb.dsm.formula("depth"), dat, mesh, "tw", range.km = 60)

test_that("tweedie: the range stays fixed and the fit carries its attributes", {
  rp <- sdmTMB::tidy(fit_tw, "ran_pars")
  expect_equal(rp$estimate[rp$term == "range"], rep(60, sum(rp$term == "range")), tolerance = 1e-6)
  expect_gt(attr(fit_tw, "fit_min"), 0)
  expect_type(attr(fit_tw, "sanity_failed"), "character")
})

test_that("a warm start from another range lands on the new range and the cold fit's optimum", {
  form <- sdmtmb.dsm.formula("depth")
  cold <- fit.sdmtmb.dsm(form, dat, mesh, "tw", range.km = 90)
  warm <- fit.sdmtmb.dsm(form, dat, mesh, "tw", range.km = 90, warm.start = fit_tw)
  rp <- sdmTMB::tidy(warm, "ran_pars")
  expect_equal(unique(rp$estimate[rp$term == "range"]), 90, tolerance = 1e-6)
  # Same model, same data: the optimum must not depend on where it started.
  expect_equal(as.numeric(stats::logLik(warm)), as.numeric(stats::logLik(cold)), tolerance = 1e-4)
  expect_error(fit.sdmtmb.dsm(form, dat, mesh, "tw", range.km = 90, warm.start = fit_tw,
                              previous.fit = fit_tw), "not both")
})

test_that("spatial.varying is passed through", {
  # One field shared by all seasons instead of one per season: a single
  # spatial_varying column, so a single field SD in the random-effect pars.
  f1 <- fit.sdmtmb.dsm(sdmtmb.dsm.formula("depth"), dat, mesh, "tw", range.km = 60,
                       spatial.varying = ~ 1)
  expect_identical(ncol(f1$tmb_data$z_i), 1L)
  expect_identical(ncol(fit_tw$tmb_data$z_i), 4L)
})

test_that("the surface sums platform levels the way get_index does", {
  grid <- expand.grid(cell = 1:50, Season = seasons, stringsAsFactors = FALSE)
  set.seed(2)
  grid$x <- rep(stats::runif(50, 0, 2e5), 4)
  grid$y <- rep(stats::runif(50, 0, 2e5), 4)
  grid$depth <- grid$x / 1000
  grid$area <- 4
  s <- get.sdmtmb.surface(fit_tw, grid, seasons, cell.col = "cell", nsim = 20)
  expect_identical(nrow(s$cells), nrow(grid))
  expect_identical(s$cells$cell, grid$cell)
  # The plug-in sum over cells and platform levels is get_index's plain total.
  expect_equal(s$totals$N_median_sum, s$totals$N_plain, tolerance = 1e-8)
  # Rescaled cells sum to the bias-corrected total.
  resc <- tapply(s$cells$N_rescaled, s$cells$Season, sum)
  expect_equal(as.vector(resc[s$totals$Season]), s$totals$N, tolerance = 1e-8)
  expect_true(all(s$cells$lo <= s$cells$hi))
  expect_true(all(s$totals$lo_draws < s$totals$hi_draws))
  expect_length(s$notes, 0)
  expect_error(get.sdmtmb.surface(fit_tw, transform(grid, platform = "WD"), seasons), "platform")
  expect_error(get.sdmtmb.surface(fit_tw, rbind(grid, grid[1, ]), seasons, cell.col = "cell"), "more than one row")
})

test_that("delta-lognormal: the positive part's smooths are fixed off", {
  # The first fit warns of a non-positive-definite Hessian; the refit is the point.
  fit_dl <- suppressWarnings(fit.sdmtmb.dsm(sdmtmb.dsm.formula("depth"), dat, mesh, "dl", range.km = 60))
  lss_map <- matrix(fit_dl$tmb_map$ln_smooth_sigma, ncol = 2)
  expect_true(all(is.na(lss_map[, 2])))
  # The occurrence part's smooths stay free unless a refit fixed collapsed ones;
  # on this data it does (a non-positive-definite Hessian), so both paths run.
  if (is.null(attr(fit_dl, "refit_note"))) expect_false(any(is.na(lss_map[, 1]))) else
    expect_match(attr(fit_dl, "refit_note"), "collapsed smooth")
  expect_false(all(is.na(lss_map[, 1])))
  expect_true(all(matrix(fit_dl$parlist$b_smooth, ncol = 2)[, 2] == 0))
  for (k in 1:2) {
    rp <- sdmTMB::tidy(fit_dl, "ran_pars", model = k)
    expect_equal(unique(rp$estimate[rp$term == "range"]), 60, tolerance = 1e-6)
  }
})
