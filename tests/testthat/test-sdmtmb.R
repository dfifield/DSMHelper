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

test_that("a fixed range gives the same fit and surface as 0.41.0", {
  # 0.42.0 added range.km = NULL and random-intercept prediction; with a fixed
  # range and no random intercept nothing may change. Reference values from
  # DSMHelper 0.41.0 on this same data, mesh and call (sdmTMB 1.2.0, TMB
  # 1.9.x, R 4.6.1, 2026-11-15): log-likelihood, the first six fixed
  # parameters, and the seasonal totals with and without bias correction.
  expect_equal(as.numeric(stats::logLik(fit_tw)), -1056.8840778245, tolerance = 1e-8)
  expect_equal(unname(fit_tw$model$par[1:6]),
               c(-1.6669772525, -1.7378645435, -1.5665760242, -1.6320518621, 0.5064366246, 12.7924570265),
               tolerance = 1e-6)
  expect_equal(attr(fit_tw, "range_km"), 60, tolerance = 1e-6)
  grid <- expand.grid(cell = 1:50, Season = seasons, stringsAsFactors = FALSE)
  set.seed(2)
  grid$x <- rep(stats::runif(50, 0, 2e5), 4)
  grid$y <- rep(stats::runif(50, 0, 2e5), 4)
  grid$depth <- grid$x / 1000
  grid$area <- 4
  s <- get.sdmtmb.surface(fit_tw, grid, seasons, cell.col = "cell")
  expect_equal(s$totals$N, c(132.26928727, 112.47607336, 129.69491764, 137.21921789), tolerance = 1e-6)
  expect_equal(s$totals$N_plain, c(131.21015367, 111.50019030, 128.34065368, 136.06093161), tolerance = 1e-6)
})

test_that("range.km = NULL estimates the range and records it", {
  free <- fit.sdmtmb.dsm(sdmtmb.dsm.formula("depth"), dat, mesh, "tw", range.km = NULL)
  rp <- sdmTMB::tidy(free, "ran_pars")
  r <- unique(rp$estimate[rp$term == "range"])
  expect_length(r, 1)
  expect_equal(attr(free, "range_km"), r)
  # Estimated, so not left at the fixed value of the other fits.
  expect_false(isTRUE(all.equal(r, 60, tolerance = 1e-3)))
  # Estimating a parameter can only raise the maximised likelihood.
  expect_gte(as.numeric(stats::logLik(free)), as.numeric(stats::logLik(fit_tw)) - 1e-6)
})

test_that("a random intercept is predicted at its mean level", {
  # Eight groups (as cruise-days in #142) with a real group effect.
  d2 <- dat
  d2$grp <- factor(rep_len(1:8, nrow(d2)))
  set.seed(3)
  d2$abundance.est <- d2$abundance.est * exp(stats::rnorm(8, 0, 0.8))[d2$grp]
  f.re <- stats::as.formula(paste(deparse(sdmtmb.dsm.formula("depth"), width.cutoff = 500L), "+ (1 | grp)"))
  fit_re <- fit.sdmtmb.dsm(f.re, d2, mesh, "tw", range.km = 60)
  rp <- sdmTMB::tidy(fit_re, "ran_pars")
  sig <- rp$estimate[rp$term == "sd__(Intercept)" & rp$group_name == "grp"]
  expect_length(sig, 1)
  grid <- expand.grid(cell = 1:20, Season = seasons, stringsAsFactors = FALSE)
  grid$x <- rep(seq(1e4, 1.9e5, length.out = 20), 4)
  grid$y <- 5e4
  grid$depth <- grid$x / 1000
  grid$area <- 4
  # The grid has no grp column: the surface must not need one.
  s <- get.sdmtmb.surface(fit_re, grid, seasons, cell.col = "cell", nsim = 10, bias.correct = FALSE)
  # By hand: each platform level predicted at the median level (re_form_iid =
  # NA), times exp(sigma^2 / 2), summed over the levels, times area.
  g <- prepare.sdmtmb.data(grid[rep(seq_len(nrow(grid)), 2), ], seasons)
  g$platform <- factor(rep(levels(fit_re$data$platform), each = nrow(grid)), levels = levels(fit_re$data$platform))
  p <- stats::predict(fit_re, newdata = g, offset = rep(0, nrow(g)), type = "response", re_form_iid = NA)
  by.hand <- rowsum(p$est * exp(sig^2 / 2) * g$area, rep(seq_len(nrow(grid)), 2))[, 1]
  expect_equal(unname(s$cells$N), unname(by.hand), tolerance = 1e-8)
  # get_index()'s plain total is scaled by the same factor, so it still equals
  # the cells' sum; the draws were taken (scaled per draw, no level noise).
  expect_equal(s$totals$N_plain, s$totals$N_median_sum, tolerance = 1e-8)
  expect_false(anyNA(s$cells$lo))
  expect_match(paste(s$notes, collapse = " "), "random intercept")
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
  # Second-order cells: (1 + v/2) >= 1, so never below the plug-in, and they
  # sum to N_so_sum.
  expect_true(all(s$cells$N_so >= s$cells$N))
  so <- tapply(s$cells$N_so, s$cells$Season, sum)
  expect_equal(as.vector(so[s$totals$Season]), s$totals$N_so_sum, tolerance = 1e-8)
  expect_length(s$notes, 0)
  # Without draws there is no second-order surface.
  s0 <- get.sdmtmb.surface(fit_tw, grid, seasons, cell.col = "cell", bias.correct = FALSE)
  expect_true(all(is.na(s0$cells$N_so)))
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

test_that("a collapsed season field is fixed off in the refit (0.48.0)", {
  # The first delta-lognormal fit on this data has a non-positive-definite
  # Hessian (test above), so the refit runs. The threshold is set just above
  # the second-smallest fitted field SD (over both parts; several are tiny on
  # this synthetic data), so exactly two fields count as collapsed: they must
  # come back fixed at SD 1e-6 with all their values 0, the others left free.
  f0 <- suppressWarnings(fit.sdmtmb.dsm(sdmtmb.dsm.formula("depth"), dat, mesh, "dl", range.km = 60,
                                        refit.collapsed = FALSE))
  expect_false(isTRUE(f0$sd_report$pdHess))
  sig0 <- matrix(f0$tmb_obj$report(f0$tmb_obj$env$last.par.best)$sigma_Z, ncol = 2)
  thr <- sort(sig0)[2] * (1 + 1e-6)          # catches the two smallest
  fit <- suppressWarnings(fit.sdmtmb.dsm(sdmtmb.dsm.formula("depth"), dat, mesh, "dl", range.km = 60,
                                         collapsed.field.sd = thr))
  want <- sig0 < thr
  expect_identical(sum(want), 2L)
  expect_match(attr(fit, "refit_note"), "2 collapsed field\\(s\\) \\[")
  expect_match(attr(fit, "refit_note"), "Season[A-Za-z]+ \\(component [12]\\)")
  # Fixed: map NA, SD 1e-6, every vertex value 0. Free: map not NA.
  ltz_map <- matrix(fit$tmb_map$ln_tau_Z, ncol = 2)
  expect_identical(is.na(ltz_map), want)
  sig1 <- matrix(fit$tmb_obj$report(fit$tmb_obj$env$last.par.best)$sigma_Z, ncol = 2)
  expect_equal(sig1[want], rep(1e-6, 2), tolerance = 1e-6)
  zeta <- array(fit$parlist$zeta_s, dim = c(length(fit$parlist$zeta_s) / 8, 4, 2))
  for (w in which(want)) {
    j <- (w - 1) %% 4 + 1
    k <- (w - 1) %/% 4 + 1
    expect_true(all(zeta[, j, k] == 0))
  }
  # collapsed.field.sd = 0 switches the field check off: smooths only, as before.
  fit_nf <- suppressWarnings(fit.sdmtmb.dsm(sdmtmb.dsm.formula("depth"), dat, mesh, "dl", range.km = 60,
                                            collapsed.field.sd = 0))
  expect_null(fit_nf$tmb_map$ln_tau_Z)
})

test_that("a fishing term with no detection in a season is fixed, as if its rows were dropped (0.49.0)", {
  # About 10% of rows fishing-affected, and every Spring one of them a zero:
  # SeasonSpring:fishing then has no finite estimate (Atl IMRP Gulls had 160
  # such rows and no detection on them).
  d <- make.sim.data()
  set.seed(2)
  d$fishing <- stats::rbinom(nrow(d), 1, 0.1)
  d$abundance.est[d$fishing == 1 & d$Season == "Spring"] <- 0
  dd <- prepare.sdmtmb.data(d, seasons)
  mesh_f <- bind.sdmtmb.mesh(dd, fm, land)
  form <- sdmtmb.dsm.formula("depth", fishing = TRUE)
  fit <- suppressWarnings(fit.sdmtmb.dsm(form, dd, mesh_f, "tw", range.km = 60))
  X <- fit$tmb_data$X_ij[[1]]
  j <- which(colnames(X) == "SeasonSpring:fishing")
  expect_length(j, 1)
  # Fixed at -20 (map NA), the other seasons' fishing terms estimated.
  expect_identical(attr(fit, "fixed_coefs"), "SeasonSpring:fishing (component 1) = -20")
  expect_true(is.na(fit$tmb_map$b_j[j]))
  expect_equal(fit$parlist$b_j[j], -20)
  others <- grep(":fishing$", colnames(X))
  expect_false(any(is.na(fit$tmb_map$b_j[setdiff(others, j)])))

  # Same estimates as dropping those rows: refit without them (new mesh
  # binding, since the projection matrix is per row). The Spring fishing
  # column is then all zero, so it is fixed at 0 there. The platform and the
  # other fishing coefficients must agree to well within the optimiser's
  # tolerance. The Season intercepts may not: dropping rows also changes the
  # depth smooth's basis, which s() builds from the data (its centring
  # constraint is a mean over the rows), and that moves the intercepts by
  # about 0.007 on this data, so they get a looser 0.02.
  keep <- !(dd$fishing == 1 & dd$Season == "Spring")
  dk <- dd[keep, ]
  fit_drop <- suppressWarnings(fit.sdmtmb.dsm(form, dk, bind.sdmtmb.mesh(dk, fm, land), "tw", range.km = 60))
  expect_identical(attr(fit_drop, "fixed_coefs"), "SeasonSpring:fishing (component 1) = 0")
  ints <- grep("^Season[A-Za-z]+$", colnames(X))
  rest <- setdiff(seq_len(ncol(X)), c(ints, j))
  expect_equal(fit$parlist$b_j[rest], fit_drop$parlist$b_j[rest], tolerance = 1e-4)
  expect_true(all(abs(fit$parlist$b_j[ints] - fit_drop$parlist$b_j[ints]) < 0.02))

  # Delta-lognormal: also fixed at 0 in the positive part, which has no data
  # for it; no flat direction is left from it in either part.
  fit_dl <- suppressWarnings(fit.sdmtmb.dsm(form, dd, mesh_f, "dl", range.km = 60))
  expect_true(all(c("SeasonSpring:fishing (component 1) = -20", "SeasonSpring:fishing (component 2) = 0") %in%
                    attr(fit_dl, "fixed_coefs")))
  expect_true(is.na(fit_dl$tmb_map$b_j2[j]))

  # fix.separated = FALSE estimates it, as before 0.49.0.
  fit_off <- suppressWarnings(fit.sdmtmb.dsm(form, dd, mesh_f, "tw", range.km = 60, fix.separated = FALSE))
  expect_length(attr(fit_off, "fixed_coefs"), 0)
  # sdmTMB always fills a map for b_j itself, so check nothing in it is fixed.
  expect_false(any(is.na(fit_off$tmb_map$b_j)))
})
