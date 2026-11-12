# sdmTMB density surface models (Atlantic-DSM issue #130)
#
# The spec settled in #130's Phase 1 trials (Atlantic-DSM
# Docs/issue_130_sdmtmb.md): the response is dsm's Horvitz-Thompson
# abundance.est with log segment area as the offset; one Matern field per
# season (spatial_varying = ~ 0 + Season) in place of s(x, y, by = Season); the
# range fixed rather than estimated; a barrier mesh so the field does not reach
# across land; and the covariate smooths kept as bs = "ts". Coordinates are in
# km throughout, so the range is in km.
#
# sdmTMB, fmesher and sdmTMBextra are Suggests, not Imports: sdmTMBextra is on
# GitHub only, and the mgcv pipeline does not need any of them.

#' Stop unless the sdmTMB packages are installed
#'
#' @param pkgs Package names.
#' @return `TRUE`, invisibly.
#' @examples
#' \dontrun{
#' check.sdmtmb.installed("sdmTMB")
#' }
#' @noRd
check.sdmtmb.installed <- function(pkgs = c("sdmTMB", "fmesher")) {
  checkmate::expect_character(pkgs, min.len = 1)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("the sdmTMB track needs ", paste(missing, collapse = ", "),
         "; see renv.lock (sdmTMBextra is pbs-assess/sdmTMBextra on GitHub)", call. = FALSE)
  invisible(TRUE)
}

#' Land polygons for an sdmTMB barrier mesh
#'
#' ETOPO altitude above 0 m around the domain, polygonised, projected to the
#' domain's CRS, simplified, and converted to km with the CRS removed, which is
#' what `sdmTMBextra::add_barrier_mesh()` needs when the mesh is in km.
#'
#' @param domain.sf The domain (sf), normally the buffered training study
#'   area, in `segProj` metres.
#' @param etopo.file Path to the ETOPO netCDF (`etopo180.nc`).
#' @param margin.deg Degrees added around the domain's lon/lat box before
#'   cropping, so land just outside it is included.
#' @param simplify.m Simplification tolerance in metres.
#' @return sf polygons of land, coordinates in km, no CRS.
#' @examples
#' \dontrun{
#' land <- build.land.polygons(training.study.area.buffered,
#'                             file.path(predLayerDir, "NetCDF", "etopo180.nc"))
#' }
#' @export
build.land.polygons <- function(domain.sf, etopo.file, margin.deg = 2, simplify.m = 1000) {
  checkmate::expect_class(domain.sf, "sf")
  checkmate::expect_file_exists(etopo.file)
  checkmate::expect_number(margin.deg, lower = 0)
  checkmate::expect_number(simplify.m, lower = 0)
  bb <- sf::st_bbox(sf::st_transform(domain.sf, 4326))
  alt <- terra::rast(etopo.file)
  alt <- terra::crop(alt, terra::ext(bb[["xmin"]] - margin.deg, bb[["xmax"]] + margin.deg,
                                     bb[["ymin"]] - margin.deg, bb[["ymax"]] + margin.deg))
  land <- terra::classify(alt > 0, cbind(0, NA)) %>%
    terra::as.polygons(dissolve = TRUE) %>%
    sf::st_as_sf() %>%
    sf::st_transform(sf::st_crs(domain.sf)) %>%
    sf::st_simplify(dTolerance = simplify.m) %>%
    sf::st_make_valid()
  sf::st_geometry(land) <- sf::st_geometry(land) / 1000
  sf::st_crs(land) <- NA
  land
}

#' The fmesher mesh for an sdmTMB DSM
#'
#' Vertices at the unique segment locations, refined to `max.edge.km` inside
#' the domain and to 4x that in an outer extension, with a cutoff of a quarter
#' of the edge. The mesh depends only on the locations and the domain, so it
#' can be built once and cached; [bind.sdmtmb.mesh()] attaches it to a data
#' set. On #130's trials: 8,625 vertices for NL_EXPL_DRL_RA and 6,515 for
#' Atl IMRP at 50 km, following the domain's area rather than the row count.
#'
#' @param locs Data frame with `x` and `y` in `segProj` metres (segments or
#'   segment x ddftype rows; duplicates are dropped).
#' @param domain.sf The domain (sf), normally the buffered training study area.
#' @param max.edge.km Largest triangle edge inside the domain, km.
#' @return An `fm_mesh_2d` in km.
#' @examples
#' \dontrun{
#' fm <- build.sdmtmb.mesh(segdata, training.study.area.buffered, 50)
#' }
#' @export
build.sdmtmb.mesh <- function(locs, domain.sf, max.edge.km) {
  check.sdmtmb.installed("fmesher")
  checkmate::expect_data_frame(locs)
  checkmate::expect_names(names(locs), must.include = c("x", "y"))
  checkmate::expect_class(domain.sf, "sf")
  checkmate::expect_number(max.edge.km, lower = 0, finite = TRUE)
  domain <- sf::st_coordinates(sf::st_geometry(domain.sf))[, c("X", "Y")] / 1000
  xy <- unique(cbind(locs$x, locs$y) / 1000)
  fmesher::fm_mesh_2d_inla(loc = xy, loc.domain = domain,
                           max.edge = c(1, 4) * max.edge.km, offset = c(1, 4) * max.edge.km,
                           cutoff = max.edge.km / 4)
}

#' Attach a cached mesh to a data set, with the barrier if land is given
#'
#' `sdmTMB::make_mesh()` binds the mesh to the data's rows, so this is done
#' once per data set (seconds), while [build.sdmtmb.mesh()] is done once.
#' With `land`, `sdmTMBextra::add_barrier_mesh()` marks triangles whose centre
#' is on land, which then get `range.fraction` x the range.
#'
#' @param dat Data with `X` and `Y` in km ([prepare.sdmtmb.data()]).
#' @param fm The mesh from [build.sdmtmb.mesh()].
#' @param land Land polygons from [build.land.polygons()], or `NULL` for no
#'   barrier.
#' @param range.fraction Share of the range a land triangle keeps (0.2, the
#'   sdmTMBextra default, used in every #130 trial).
#' @return An `sdmTMBmesh`.
#' @examples
#' \dontrun{
#' mesh <- bind.sdmtmb.mesh(dat, fm, land)
#' }
#' @export
bind.sdmtmb.mesh <- function(dat, fm, land = NULL, range.fraction = 0.2) {
  check.sdmtmb.installed(c("sdmTMB", if (!is.null(land)) "sdmTMBextra"))
  checkmate::expect_data_frame(dat)
  checkmate::expect_names(names(dat), must.include = c("X", "Y"))
  checkmate::expect_class(fm, "fm_mesh_2d")
  if (!is.null(land)) checkmate::expect_class(land, "sf")
  checkmate::expect_number(range.fraction, lower = 0, upper = 1)
  mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), mesh = fm)
  if (!is.null(land))
    mesh <- sdmTMBextra::add_barrier_mesh(mesh, land, range_fraction = range.fraction,
                                          proj_scaling = 1, plot = FALSE)
  mesh
}

#' The columns an sdmTMB DSM needs
#'
#' Adds `X` and `Y` (km), `Season` as a factor in `season.levels` order,
#' `season_num` (the `time` variable), drops unused `platform` levels, and
#' adds `log_area` where `segment.area` is present.
#'
#' @param d Segment rows or a prediction grid, with `x`, `y` (metres) and
#'   `Season`.
#' @param season.levels Season names in order (`season.names`).
#' @return `d` with the columns added.
#' @examples
#' \dontrun{
#' dat <- prepare.sdmtmb.data(dat, season.names)
#' }
#' @export
prepare.sdmtmb.data <- function(d, season.levels) {
  checkmate::expect_data_frame(d)
  checkmate::expect_names(names(d), must.include = c("x", "y", "Season"))
  checkmate::expect_character(season.levels, min.len = 1, unique = TRUE)
  if (!all(as.character(d$Season) %in% season.levels))
    stop("Season values not in season.levels: ",
         paste(setdiff(unique(as.character(d$Season)), season.levels), collapse = ", "))
  d$X <- d$x / 1000
  d$Y <- d$y / 1000
  d$Season <- factor(as.character(d$Season), levels = season.levels)
  d$season_num <- as.integer(d$Season)
  if ("platform" %in% names(d)) d$platform <- droplevels(factor(d$platform))
  if ("segment.area" %in% names(d)) d$log_area <- log(d$segment.area)
  d
}

#' The formula of an sdmTMB DSM
#'
#' The mgcv candidates with `s(x, y, by = Season)` removed, since the
#' per-season fields replace it: a Season intercept, `platform`, one smooth
#' per covariate by Season, and `Season:fishing` for fishing-affected species
#' (issue #107).
#'
#' @param covars Covariates with a smooth each.
#' @param fishing Add `Season:fishing`.
#' @param response Response column.
#' @param bs Smooth basis; `"ts"` shrinks a smooth out entirely (#130 trial 2).
#' @return A formula.
#' @examples
#' sdmtmb.dsm.formula()
#' @export
sdmtmb.dsm.formula <- function(covars = c("depth", "depth.g", "sst"), fishing = FALSE,
                               response = "abundance.est", bs = "ts") {
  checkmate::expect_character(covars, any.missing = FALSE, unique = TRUE)
  checkmate::expect_flag(fishing)
  checkmate::expect_string(response)
  checkmate::expect_string(bs)
  rhs <- c("0 + Season + platform",
           if (length(covars)) sprintf('s(%s, by = Season, bs = "%s")', covars, bs))
  if (fishing) rhs <- c(rhs, "Season:fishing")
  stats::as.formula(paste(response, "~", paste(rhs, collapse = " + ")), env = globalenv())
}

#' An sdmTMB family from its short key
#'
#' @param key `"tw"` tweedie, `"nb"` nbinom2, `"dg"` delta-gamma, `"dl"`
#'   delta-lognormal (both deltas with `type = "poisson-link"`).
#' @return The family object.
#' @examples
#' \dontrun{
#' sdmtmb.family("tw")
#' }
#' @export
sdmtmb.family <- function(key) {
  check.sdmtmb.installed("sdmTMB")
  checkmate::expect_choice(key, c("tw", "nb", "dg", "dl"))
  switch(key,
         tw = sdmTMB::tweedie(link = "log"),
         nb = sdmTMB::nbinom2(link = "log"),
         dg = sdmTMB::delta_gamma(type = "poisson-link"),
         dl = sdmTMB::delta_lognormal(type = "poisson-link"))
}

#' Fix chosen smooths at "off" in a start/map pair
#'
#' A smooth is off when its SD is fixed at `log(1e-6)` and its coefficients at
#' 0. Used for a delta model's positive part (sdmTMB 1.2.0 refuses different
#' formulas in the two parts when they contain smooths) and for smooths whose
#' SD collapsed towards 0, which leave the Hessian not positive definite.
#'
#' @param params `tmb_params` from a build or a fit's `parlist` (supplies
#'   the shapes and the starting values kept for the other smooths).
#' @param b.smooth.start `tmb_data$b_smooth_start`: 0-based first row of each
#'   smooth's coefficients in `b_smooth`.
#' @param off Logical matrix, smooths x components: which to fix.
#' @param start,map The start and map lists to add to.
#' @return List with `start` and `map`.
#' @noRd
fix.smooths.off <- function(params, b.smooth.start, off, start, map) {
  lss <- matrix(params$ln_smooth_sigma, ncol = ncol(off))
  bsm <- matrix(params$b_smooth, ncol = ncol(off))
  lss_map <- if (is.null(map$ln_smooth_sigma)) matrix(seq_along(lss), ncol = ncol(off)) else
    matrix(as.integer(map$ln_smooth_sigma), ncol = ncol(off))
  bsm_map <- if (is.null(map$b_smooth)) matrix(seq_along(bsm), ncol = ncol(off)) else
    matrix(as.integer(map$b_smooth), ncol = ncol(off))
  size <- diff(c(b.smooth.start, nrow(bsm)))
  for (k in seq_len(ncol(off))) for (j in which(off[, k])) {
    rows <- b.smooth.start[j] + seq_len(size[j])
    lss[j, k] <- log(1e-6)
    lss_map[j, k] <- NA
    bsm[rows, k] <- 0
    bsm_map[rows, k] <- NA
  }
  start$ln_smooth_sigma <- lss
  start$b_smooth <- bsm
  map$ln_smooth_sigma <- factor(lss_map)
  map$b_smooth <- factor(bsm_map)
  list(start = start, map = map)
}

#' Fit an sdmTMB DSM to #130's spec
#'
#' One Matern field per season (`spatial_varying = ~ 0 + Season`, `time =
#' "season_num"`, no shared spatial or spatiotemporal field), the range fixed
#' at `range.km` by mapping `ln_kappa`, and the offset as a column. For a delta
#' family the positive part's smooths are fixed off unless
#' `positive.smooths = TRUE`, since sdmTMB will not take a simpler formula
#' there. If the Hessian is not positive definite and some smooth SDs are below
#' `collapsed.sd`, the model is refitted once with those smooths fixed off,
#' which restored it in the #134 simulation; draws and the bias-corrected index
#' need a positive definite Hessian.
#'
#' @param formula From [sdmtmb.dsm.formula()].
#' @param dat Data from [prepare.sdmtmb.data()].
#' @param mesh From [bind.sdmtmb.mesh()], bound to `dat`.
#' @param family Family key for [sdmtmb.family()].
#' @param range.km The fixed range, km (150 for NL_EXPL_DRL_RA, #130 trial 1b).
#' @param offset Offset column (log segment area).
#' @param positive.smooths Keep a delta model's positive-part smooths.
#' @param refit.collapsed Refit once with collapsed smooths fixed off when the
#'   Hessian is not positive definite.
#' @param collapsed.sd Smooth SD below which a smooth counts as collapsed.
#' @param previous.fit A fit of the same model to the same rows, to warm-start
#'   from (`sdmTMB(previous_fit = )`); only the response may differ.
#' @param nthreads TMB threads (`sdmTMBcontrol(parallel = )`).
#' @return The fitted `sdmTMB`, with attributes `fit_min` (minutes, all fits
#'   included), `refit_note` (`NULL` if no refit) and `sanity_failed`.
#' @examples
#' \dontrun{
#' fit <- fit.sdmtmb.dsm(sdmtmb.dsm.formula(), dat, mesh, "tw", range.km = 150)
#' }
#' @export
fit.sdmtmb.dsm <- function(formula, dat, mesh, family = "tw", range.km,
                           offset = "log_area", positive.smooths = FALSE,
                           refit.collapsed = TRUE, collapsed.sd = 1e-3,
                           previous.fit = NULL, nthreads = 1L) {
  check.sdmtmb.installed("sdmTMB")
  checkmate::expect_formula(formula)
  checkmate::expect_data_frame(dat)
  checkmate::expect_names(names(dat), must.include = c("Season", "season_num", offset))
  checkmate::expect_class(mesh, "sdmTMBmesh")
  checkmate::expect_choice(family, c("tw", "nb", "dg", "dl"))
  checkmate::expect_number(range.km, lower = 0, finite = TRUE)
  checkmate::expect_flag(positive.smooths)
  checkmate::expect_flag(refit.collapsed)
  checkmate::expect_number(collapsed.sd, lower = 0)
  if (!is.null(previous.fit)) checkmate::expect_class(previous.fit, "sdmTMB")
  checkmate::expect_count(nthreads, positive = TRUE)
  t0 <- Sys.time()
  n_comp <- if (family %in% c("dg", "dl")) 2L else 1L
  # ln_kappa is (spatial, spatiotemporal) x component.
  start <- list(ln_kappa = matrix(log(sqrt(8) / range.km), nrow = 2, ncol = n_comp))
  map <- list(ln_kappa = factor(rep(NA, 2 * n_comp)))
  fit.it <- function(start, map, do.fit = TRUE) {
    args <- list(formula = formula, data = dat, mesh = mesh, family = sdmtmb.family(family),
                 offset = offset, time = "season_num", spatial = "off", spatiotemporal = "off",
                 spatial_varying = ~ 0 + Season, do_fit = do.fit,
                 control = sdmTMB::sdmTMBcontrol(start = start, map = map, parallel = nthreads))
    if (!is.null(previous.fit) && do.fit) args$previous_fit <- previous.fit
    do.call(sdmTMB::sdmTMB, args)
  }
  if (n_comp == 2L && !positive.smooths) {
    built <- fit.it(start, map, do.fit = FALSE)
    n_sm <- length(built$tmb_data$b_smooth_start)
    if (n_sm > 0) {
      off <- matrix(FALSE, nrow = n_sm, ncol = 2)
      off[, 2] <- TRUE
      sm <- fix.smooths.off(built$tmb_params, built$tmb_data$b_smooth_start, off, start, map)
      start <- sm$start
      map <- sm$map
    }
  }
  fit <- fit.it(start, map)
  refit_note <- NULL
  if (refit.collapsed && !isTRUE(fit$sd_report$pdHess)) {
    lss_hat <- matrix(fit$parlist$ln_smooth_sigma, ncol = n_comp)
    lss_map <- if (is.null(map$ln_smooth_sigma)) matrix(TRUE, nrow(lss_hat), n_comp) else
      matrix(!is.na(map$ln_smooth_sigma), ncol = n_comp)
    off <- exp(lss_hat) < collapsed.sd & lss_map
    if (any(off)) {
      sm <- fix.smooths.off(fit$parlist, fit$tmb_data$b_smooth_start, off, start, map)
      refit_note <- sprintf("refitted with %d collapsed smooth(s) fixed off; first fit's Hessian not positive definite",
                            sum(off))
      fit <- fit.it(sm$start, sm$map)
    }
  }
  for (k in seq_len(n_comp)) {
    rp <- sdmTMB::tidy(fit, "ran_pars", model = k)
    r <- rp$estimate[rp$term == "range"]
    if (!length(r) || any(abs(r / range.km - 1) > 1e-6))
      stop(sprintf("range meant to be fixed at %g km, but component %d reports %s", range.km, k,
                   paste(signif(r, 4), collapse = ", ")))
  }
  san <- sdmTMB::sanity(fit, silent = TRUE)
  attr(fit, "sanity_failed") <- setdiff(names(san)[!unlist(san)], "all_ok")
  attr(fit, "refit_note") <- refit_note
  attr(fit, "fit_min") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  fit
}

#' Predicted surface, seasonal totals and their uncertainty from an sdmTMB DSM
#'
#' The grid is replicated over the model's platform levels and predicted at
#' offset 0, so each row is a density, and a cell's abundance sums density x
#' area over the levels, as `dsm.pred()` does. Returned per cell:
#'
#' * `N`, the plug-in surface. Each row is exp(linear predictor), the median of
#'   a log-normal, not the mean.
#' * `N_rescaled`, `N` scaled per season so it sums to the bias-corrected
#'   total: the median's shape at the mean's level. #130 trial 6 found per-cell
#'   means from draws unusable, so this is the proposed map product.
#' * with `nsim > 0`, `lo`/`hi`, the 2.5% and 97.5% quantiles of the per-cell
#'   abundance over draws, each draw summed over platform levels first.
#'
#' Per season: `get_index()` with and without bias correction, and with draws
#' the quantiles of each draw's total. #130 trial 6b found the draws' spread of
#' the log total 1.4-4.2x `get_index()`'s standard error on NL ATPU; which is
#' right is open, so both are returned. Each uncertainty step can fail where
#' the surface does not (it needs the joint precision), so each runs on its
#' own: a failure leaves its columns `NA` and is recorded in `notes`.
#'
#' @param fit Fitted `sdmTMB` ([fit.sdmtmb.dsm()]).
#' @param grid One row per cell x season, with no `platform` column: `x`, `y`
#'   (metres), `Season`, `area` (km²), the covariates, and `fishing = 0` for a
#'   fishing-affected species.
#' @param season.levels Season names in order (`season.names`).
#' @param cell.col Column identifying a cell, or `NULL` to identify cells by
#'   `x` and `y`.
#' @param nsim Number of draws; 0 for none.
#' @param bias.correct Which `get_index()` totals to compute: any of `TRUE`
#'   and `FALSE`.
#' @return List: `cells` (the grid's cell id, `Season`, `N`, `N_rescaled`,
#'   `lo`, `hi`), `totals` (per season: `N` and `se_log` bias-corrected,
#'   `N_plain` and `se_log_plain`, `N_median_sum`, `lo_draws`, `hi_draws`),
#'   `notes`, `pred_min`.
#' @examples
#' \dontrun{
#' s <- get.sdmtmb.surface(fit, grid, season.names, nsim = 100)
#' }
#' @export
get.sdmtmb.surface <- function(fit, grid, season.levels, cell.col = NULL, nsim = 0L,
                               bias.correct = c(TRUE, FALSE)) {
  check.sdmtmb.installed("sdmTMB")
  checkmate::expect_class(fit, "sdmTMB")
  checkmate::expect_data_frame(grid)
  checkmate::expect_names(names(grid), must.include = c("x", "y", "Season", "area", cell.col),
                          disjunct.from = "platform")
  checkmate::expect_character(season.levels, min.len = 1, unique = TRUE)
  checkmate::expect_count(nsim)
  checkmate::expect_logical(bias.correct, any.missing = FALSE, min.len = 1, unique = TRUE)
  t0 <- Sys.time()
  lv <- levels(fit$data$platform)
  if (is.null(lv)) stop("the fitted model has no platform factor")
  key <- if (is.null(cell.col)) paste(grid$x, grid$y, grid$Season) else paste(grid[[cell.col]], grid$Season)
  if (anyDuplicated(key)) stop("grid has more than one row per cell x season")
  g <- grid[rep(seq_len(nrow(grid)), times = length(lv)), , drop = FALSE]
  g$platform <- factor(rep(lv, each = nrow(grid)), levels = lv)
  g <- prepare.sdmtmb.data(g, season.levels)
  g$platform <- factor(as.character(g$platform), levels = lv)
  row_of <- rep(seq_len(nrow(grid)), times = length(lv))
  off0 <- rep(0, nrow(g))
  season <- as.character(grid$Season)
  notes <- character()
  step <- function(name, expr) tryCatch(expr, error = function(e) {
    notes <<- c(notes, sprintf("%s: %s", name, conditionMessage(e)))
    NULL
  })

  p <- stats::predict(fit, newdata = g, offset = off0, type = "response")
  n_cell <- rowsum(p$est * g$area, row_of, reorder = TRUE)[, 1]
  rm(p)
  index <- function(bc) {
    ix <- sdmTMB::get_index(fit, newdata = g, area = g$area, offset = off0, bias_correct = bc)
    ix$Season <- season.levels[ix$season_num]
    ix
  }
  ix_bc <- if (TRUE %in% bias.correct) step("get_index bias-corrected", index(TRUE))
  ix_plain <- if (FALSE %in% bias.correct) step("get_index", index(FALSE))
  lo <- hi <- rep(NA_real_, nrow(grid))
  tot_d <- NULL
  if (nsim > 0) {
    draws <- step("draws", stats::predict(fit, newdata = g, offset = off0, nsim = nsim))
    if (!is.null(draws)) {
      Nd <- rowsum(exp(draws) * g$area, row_of, reorder = TRUE)
      rm(draws)
      q <- t(apply(Nd, 1, stats::quantile, probs = c(0.025, 0.975)))
      lo <- q[, 1]
      hi <- q[, 2]
      tot_d <- rowsum(Nd, season)
    }
  }
  pick <- function(ix, s, col) if (is.null(ix)) NA_real_ else ix[[col]][ix$Season == s]
  totals <- purrr::map_dfr(intersect(season.levels, unique(season)), function(s) {
    tibble::tibble(Season = s,
                   N = pick(ix_bc, s, "est"), se_log = pick(ix_bc, s, "se"),
                   N_plain = pick(ix_plain, s, "est"), se_log_plain = pick(ix_plain, s, "se"),
                   N_median_sum = sum(n_cell[season == s]),
                   lo_draws = if (is.null(tot_d)) NA_real_ else unname(stats::quantile(tot_d[s, ], 0.025)),
                   hi_draws = if (is.null(tot_d)) NA_real_ else unname(stats::quantile(tot_d[s, ], 0.975)))
  })
  scale <- totals$N / totals$N_median_sum
  cells <- tibble::tibble(Season = factor(season, levels = season.levels), N = n_cell,
                          N_rescaled = n_cell * scale[match(season, totals$Season)], lo = lo, hi = hi)
  id <- if (is.null(cell.col)) tibble::tibble(x = grid$x, y = grid$y) else
    stats::setNames(tibble::tibble(grid[[cell.col]]), cell.col)
  cells <- dplyr::bind_cols(id, cells)
  list(cells = cells, totals = totals, notes = notes,
       pred_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
