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
  # requireNamespace() loads a package's namespace without attaching it and
  # returns FALSE (rather than erroring) when it is not installed, so this
  # collects every missing package at once instead of stopping at the first.
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  # Packages in Suggests are not installed automatically with DSMHelper, so a
  # missing one has to be reported here, with where to get it; otherwise the
  # first sdmTMB:: call would fail with a less helpful message.
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

  # The ETOPO grid is in lon/lat (EPSG:4326), so find the domain's extent in
  # lon/lat to know which part of the (global) grid to read.
  bb <- sf::st_bbox(sf::st_transform(domain.sf, 4326))

  # Crop the global grid to that box plus a margin. The margin matters because
  # the mesh extends beyond the domain (build.sdmtmb.mesh()'s outer ring is
  # 4 x max.edge.km wide), and land out there must be a barrier too, or the
  # field could route around a coastline through the mesh's outer edge.
  alt <- terra::rast(etopo.file)
  alt <- terra::crop(alt, terra::ext(bb[["xmin"]] - margin.deg, bb[["xmax"]] + margin.deg,
                                     bb[["ymin"]] - margin.deg, bb[["ymax"]] + margin.deg))

  # Turn the elevation grid into land polygons:
  #  - `alt > 0` is a logical raster: 1 (TRUE) on land, 0 (FALSE) at sea;
  #  - classify() maps 0 to NA, so only land cells keep a value;
  #  - as.polygons(dissolve = TRUE) merges adjacent land cells into one
  #    polygon per connected landmass (NA cells produce no polygon);
  #  - project to the domain's CRS (segProj, metres), where distances mean km;
  #  - simplify, so the coastline has a manageable number of vertices
  #    (simplify.m metres of tolerance, much finer than any mesh triangle);
  #  - st_make_valid() repairs self-intersections that simplification can
  #    create, which would otherwise break the point-in-polygon tests below.
  land <- terra::classify(alt > 0, cbind(0, NA)) %>%
    terra::as.polygons(dissolve = TRUE) %>%
    sf::st_as_sf() %>%
    sf::st_transform(sf::st_crs(domain.sf)) %>%
    sf::st_simplify(dTolerance = simplify.m) %>%
    sf::st_make_valid()

  # The mesh is in km (see build.sdmtmb.mesh()), so divide the coordinates by
  # 1000. Dividing an sfc geometry by a number scales every coordinate. The
  # CRS is then dropped, because a km version of segProj is not a defined CRS,
  # and add_barrier_mesh() compares the polygons with the mesh's plain (CRS-
  # free) coordinates.
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

  # The domain boundary as a plain matrix of vertices, in km. st_coordinates()
  # returns one row per vertex, with X and Y plus ring/part index columns
  # (L1, L2, ...); only X and Y are wanted. fmesher uses these points to
  # decide where the inner (fine) part of the mesh must reach.
  domain <- sf::st_coordinates(sf::st_geometry(domain.sf))[, c("X", "Y")] / 1000

  # The data locations, in km. segdata has one row per segment x ddftype, so
  # every location appears several times; unique() keeps one copy of each, as
  # duplicates would only slow the mesh build.
  xy <- unique(cbind(locs$x, locs$y) / 1000)

  # Build the triangulation (the same arguments every #130 trial used):
  #  - loc: the data locations are used as vertices (thinned by `cutoff`), so
  #    the mesh is finest where there is data;
  #  - max.edge = c(inner, outer): the longest triangle edge allowed inside the
  #    domain, and in the outer extension. The outer extension (4x coarser)
  #    moves the mesh boundary away from the data, because the SPDE field's
  #    variance is distorted near a mesh boundary;
  #  - offset = c(inner, outer): how far the inner and outer regions reach
  #    beyond the domain;
  #  - cutoff: vertices closer than this are merged, so clustered segments
  #    (along a ship track, every few km) do not produce thousands of tiny
  #    triangles. A quarter of the edge length (12.5 km at 50 km).
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

  # make_mesh() wraps the prebuilt triangulation (`mesh = fm`, so it does not
  # build a new one) and computes the projection matrix A from it: one row per
  # data row, giving the weights that interpolate the field at that row's
  # location from the three vertices of the triangle it falls in. That matrix
  # is why this step is per data set: A is tied to `dat`'s rows and order.
  mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), mesh = fm)

  # The barrier (Bakka et al. 2019). add_barrier_mesh() tests each triangle's
  # centre against the land polygons, records the land triangles in
  # mesh$barrier_triangles, and builds the matrices that give those triangles
  # range.fraction x the field's range. Correlation then dies out quickly
  # across land instead of passing through it, so density on one side of a
  # peninsula no longer leaks to the other.
  #  - proj_scaling = 1: the mesh and the land polygons are both in km, so no
  #    unit conversion is needed between them;
  #  - plot = FALSE: the function can draw the mesh; not wanted in a package.
  # Land triangles are kept in the mesh, not cut out as holes: a hole gives the
  # field a hard boundary at the coast, which inflates its variance there.
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

  # A Season value missing from season.levels would silently become NA in the
  # factor() call below, and those rows would then be dropped or mis-handled
  # by sdmTMB. Stop instead, naming the offending values.
  if (!all(as.character(d$Season) %in% season.levels))
    stop("Season values not in season.levels: ",
         paste(setdiff(unique(as.character(d$Season)), season.levels), collapse = ", "))

  # Coordinates in km. The mesh is in km, and so the fitted range is in km: in
  # metres the range parameter would be ~1e5 and its log-scale optimisation
  # badly scaled. Upper case so the metre columns x and y are kept as they are.
  d$X <- d$x / 1000
  d$Y <- d$y / 1000

  # Season as a factor with a fixed level order. as.character() first, so an
  # existing factor with a different level order is re-levelled by name rather
  # than by its integer codes.
  d$Season <- factor(as.character(d$Season), levels = season.levels)

  # sdmTMB's `time` argument needs an integer-like column. Each season's field
  # is a spatially varying coefficient on its Season dummy, but sdmTMB still
  # wants a time column to index the seasons in its output (get_index() returns
  # one row per time value), so season_num is 1, 2, 3, 4 in season.levels order.
  d$season_num <- as.integer(d$Season)

  # A platform level with no rows (e.g. a ddftype this species never had)
  # would be an all-zero dummy column, whose coefficient is not identified.
  # droplevels() removes unused levels.
  if ("platform" %in% names(d)) d$platform <- droplevels(factor(d$platform))

  # The offset: log of each segment's surveyed area (km^2), so the linear
  # predictor models density and exp(prediction) at offset 0 is birds per km^2.
  # Prediction grids have no segment.area, so this is skipped for them.
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

  # The right-hand side, as text pieces joined with " + " below:
  #  - "0 + Season": no global intercept, one intercept per season instead.
  #    Each season gets its own level, which a by = Season smooth needs because
  #    mgcv centres each factor-by smooth on zero;
  #  - "platform": the ddftype level (WD, FD, WN, FN), additive on the log
  #    scale, so each platform's copy of a segment gets its own multiplier;
  #  - one smooth per covariate per season, e.g.
  #      s(depth, by = Season, bs = "ts")
  #    sdmTMB turns each into penalised random effects through
  #    mgcv::smooth2random(); "ts" puts a penalty on the null space too, so a
  #    smooth the data do not support shrinks to nothing (#130 trial 2).
  #    There is no s(x, y): the per-season Matern fields are the spatial term.
  rhs <- c("0 + Season + platform",
           if (length(covars)) sprintf('s(%s, by = Season, bs = "%s")', covars, bs))

  # Fishing-affected species (issue #107): one log-multiplier per season for
  # segments surveyed while the ship was fishing. Predictions set fishing = 0.
  if (fishing) rhs <- c(rhs, "Season:fishing")

  # env = globalenv(): a formula carries the environment it was made in, and
  # sdmTMB evaluates variables in it. This function's environment holds
  # nothing useful, and keeping it would also keep it alive inside every
  # saved fit.
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
  # All four use a log link on the combined expected density:
  #  - tweedie and nbinom2 are single-part models with a log link;
  #  - the deltas are hurdle models (occurrence x positive density). With
  #    type = "poisson-link" (Thorson 2018) both parts share one log-linear
  #    predictor for the combined density, so the log-area offset enters both
  #    parts and get_index() sums densities exactly as for tweedie. The default
  #    delta link (logit occurrence) would not take the offset that way.
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
  # How TMB fixes parameters. `start` gives each parameter's starting value,
  # and `map` says which are estimated: map$<name> is a factor with one entry
  # per element of that parameter. Elements sharing a level are estimated as
  # one value, and an NA element is not estimated at all; it stays at its
  # start value. So fixing a smooth "off" means: start its SD at ~0 and its
  # coefficients at 0, and set their map entries to NA.
  #
  # The two parameters involved, both stored by sdmTMB as matrices with one
  # column per model component (1 for tweedie, 2 for a delta model):
  #  - ln_smooth_sigma: log SD of each smooth's random coefficients, one row
  #    per smooth (12 rows = 3 covariates x 4 seasons);
  #  - b_smooth: the smooths' random coefficients, all smooths stacked in one
  #    column (9 rows per ts smooth with k = 10, so 108 rows).
  # params holds them as vectors or matrices; matrix(ncol = ) gives them the
  # (rows x components) shape either way, filling column by column.
  lss <- matrix(params$ln_smooth_sigma, ncol = ncol(off))
  bsm <- matrix(params$b_smooth, ncol = ncol(off))

  # The existing map for each, as integer codes in the same shape. If nothing
  # is mapped yet, every element is estimated separately, which is the code
  # sequence 1, 2, 3, ... (all distinct). If a map exists (e.g. the delta
  # positive part's smooths already fixed off), keep it: its NA entries stay
  # NA and its codes are reused, so earlier fixes are not undone.
  lss_map <- if (is.null(map$ln_smooth_sigma)) matrix(seq_along(lss), ncol = ncol(off)) else
    matrix(as.integer(map$ln_smooth_sigma), ncol = ncol(off))
  bsm_map <- if (is.null(map$b_smooth)) matrix(seq_along(bsm), ncol = ncol(off)) else
    matrix(as.integer(map$b_smooth), ncol = ncol(off))

  # Number of coefficient rows of each smooth. b.smooth.start holds each
  # smooth's first row in b_smooth, 0-based (it comes from the TMB C++ side),
  # so the differences between successive starts, closing with the total row
  # count, are the sizes. E.g. starts 0, 9, 18 with 27 rows give sizes 9, 9, 9.
  size <- diff(c(b.smooth.start, nrow(bsm)))

  # For every (smooth j, component k) marked in `off`:
  for (k in seq_len(ncol(off))) for (j in which(off[, k])) {
    # 1-based rows of smooth j in b_smooth: a 0-based start plus 1..size is
    # exactly the 1-based rows start+1 .. start+size.
    rows <- b.smooth.start[j] + seq_len(size[j])
    # SD fixed at 1e-6, effectively zero. Not exactly zero, because log(0) is
    # -Inf, and TMB still evaluates the random-effect density with this SD.
    lss[j, k] <- log(1e-6)
    lss_map[j, k] <- NA
    # Coefficients fixed at 0, so the smooth contributes nothing.
    bsm[rows, k] <- 0
    bsm_map[rows, k] <- NA
  }

  # Write the results into the start and map lists. factor() of an integer
  # matrix gives one level per distinct non-NA code, in column-major order,
  # which is the element order TMB expects. NAs stay NA.
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
  # Timed from here, so fit_min includes the build and any refit.
  t0 <- Sys.time()

  # A delta model has two components (occurrence and positive density), each
  # with its own linear predictor, fields and parameters; the others have one.
  # Several parameters below are matrices with one column per component.
  n_comp <- if (family %in% c("dg", "dl")) 2L else 1L

  # ---- Fix the Matern range ----
  # The range is fixed rather than estimated: estimated, it follows the mesh
  # down (70 -> 44 -> 28 km as the mesh was refined, #130 trial 5) and the
  # totals never settle. sdmTMB parameterises the range through
  # kappa, with range = sqrt(8) / kappa (the distance at which correlation
  # falls to about 0.13), so kappa = sqrt(8) / range and the parameter
  # ln_kappa = log(sqrt(8) / range.km).
  # ln_kappa is a 2 x n_comp matrix: row 1 for the spatial field (which here
  # is also what the spatial_varying fields use), row 2 for the spatiotemporal
  # field (switched off here, but the parameter still exists). Both rows are
  # given the fixed value and mapped to NA, so neither is estimated.
  start <- list(ln_kappa = matrix(log(sqrt(8) / range.km), nrow = 2, ncol = n_comp))
  map <- list(ln_kappa = factor(rep(NA, 2 * n_comp)))

  # ---- One sdmTMB call, used for the build, the fit and any refit ----
  # Defined once so the three calls cannot drift apart. The model:
  #  - offset: the log segment area column, so the model is for density;
  #  - spatial = "off", spatiotemporal = "off": no field shared by all seasons
  #    and no field per time step;
  #  - spatial_varying = ~ 0 + Season: one Matern field per season, each
  #    multiplying its Season dummy (1 in that season's rows, 0 elsewhere), so
  #    each season has its own independent spatial surface. This replaces
  #    mgcv's s(x, y, by = Season);
  #  - time = "season_num": indexes the seasons for get_index();
  #  - do_fit = FALSE builds the TMB inputs without optimising (seconds),
  #    which is used below to read the smooths' layout;
  #  - parallel: TMB's OpenMP threads.
  fit.it <- function(start, map, do.fit = TRUE) {
    args <- list(formula = formula, data = dat, mesh = mesh, family = sdmtmb.family(family),
                 offset = offset, time = "season_num", spatial = "off", spatiotemporal = "off",
                 spatial_varying = ~ 0 + Season, do_fit = do.fit,
                 control = sdmTMB::sdmTMBcontrol(start = start, map = map, parallel = nthreads))
    # Warm start: previous_fit makes sdmTMB start every parameter, fixed and
    # random, at that fit's estimates. Used by the detection-function
    # uncertainty refits, where only the response changes slightly. Not for
    # the do_fit = FALSE build, which does not optimise.
    if (!is.null(previous.fit) && do.fit) args$previous_fit <- previous.fit
    do.call(sdmTMB::sdmTMB, args)
  }

  # ---- Delta models: switch off the positive part's smooths ----
  # The positive part was meant to have no smooths (trial 3b: they shrank to
  # nothing anyway, and their collapsing SDs made the Hessian non-positive
  # definite). sdmTMB 1.2.0 refuses a different formula for the two parts of a
  # delta model when the formula has smooths, so both parts get the full
  # formula, and the positive part's smooths are then fixed at zero through
  # start and map.
  if (n_comp == 2L && !positive.smooths) {
    # Build without fitting, to get the parameter layout: tmb_params has the
    # shapes, and tmb_data$b_smooth_start gives where each smooth's
    # coefficients begin.
    built <- fit.it(start, map, do.fit = FALSE)
    n_sm <- length(built$tmb_data$b_smooth_start)
    if (n_sm > 0) {
      # off: smooths x components. Column 1 (occurrence) all FALSE, column 2
      # (positive part) all TRUE, so every positive-part smooth is fixed off.
      off <- matrix(FALSE, nrow = n_sm, ncol = 2)
      off[, 2] <- TRUE
      sm <- fix.smooths.off(built$tmb_params, built$tmb_data$b_smooth_start, off, start, map)
      start <- sm$start
      map <- sm$map
    }
  }

  # ---- The fit ----
  fit <- fit.it(start, map)

  # ---- Refit if the Hessian is not positive definite ----
  # A smooth whose SD heads to 0 (shrunk out by its ts penalty) has a log SD
  # heading to -Inf, where the likelihood is flat. That makes the Hessian of
  # the fixed effects singular (not positive definite: sd_report$pdHess is
  # FALSE), and then standard errors, draws and the bias-corrected index are
  # unreliable or fail. The smooth contributes nothing, so fixing it at zero
  # changes no estimate but removes the flat direction. This is done once only;
  # a Hessian still not positive definite after it is left for the caller to
  # see in sanity_failed.
  refit_note <- NULL
  if (refit.collapsed && !isTRUE(fit$sd_report$pdHess)) {
    # Fitted log SDs, as smooths x components.
    lss_hat <- matrix(fit$parlist$ln_smooth_sigma, ncol = n_comp)
    # Which SDs were actually estimated (TRUE) rather than already fixed by
    # the map (e.g. the delta positive part above). Only estimated ones are
    # candidates; the fixed ones sit at 1e-6 and would all count as collapsed.
    lss_map <- if (is.null(map$ln_smooth_sigma)) matrix(TRUE, nrow(lss_hat), n_comp) else
      matrix(!is.na(map$ln_smooth_sigma), ncol = n_comp)
    # Collapsed: estimated, and SD below collapsed.sd (1e-3 by default). A
    # shrunk-out smooth's SD is typically 1e-4 to 1e-5 (trial 2), and a smooth
    # in use has an SD of order 0.1-10, so the threshold separates them widely.
    off <- exp(lss_hat) < collapsed.sd & lss_map
    if (any(off)) {
      # Start from this fit's estimates (parlist), so the refit is warm, and
      # add the newly fixed smooths to the existing start and map.
      sm <- fix.smooths.off(fit$parlist, fit$tmb_data$b_smooth_start, off, start, map)
      refit_note <- sprintf("refitted with %d collapsed smooth(s) fixed off; first fit's Hessian not positive definite",
                            sum(off))
      fit <- fit.it(sm$start, sm$map)
    }
  }

  # ---- Check that the range really was fixed ----
  # A mistake in the start/map above would let the range be estimated with no
  # error from sdmTMB, and every later result would quietly depend on it. So
  # read the range back from each component's random-effect parameters and
  # require it to equal range.km to 1 part in a million.
  for (k in seq_len(n_comp)) {
    rp <- sdmTMB::tidy(fit, "ran_pars", model = k)
    r <- rp$estimate[rp$term == "range"]
    if (!length(r) || any(abs(r / range.km - 1) > 1e-6))
      stop(sprintf("range meant to be fixed at %g km, but component %d reports %s", range.km, k,
                   paste(signif(r, 4), collapse = ", ")))
  }

  # ---- Record diagnostics on the fit ----
  # sanity() returns a list of named TRUE/FALSE checks plus an overall
  # "all_ok". Keep the names of the failed checks only. Expect
  # "se_magnitude_ok" to fail routinely: it flags the large standard errors of
  # shrunk smooths' log SDs, which is shrinkage working (trial 2).
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

  # ---- Replicate the grid over platform levels ----
  # The model has a platform term, so a prediction is for one platform level
  # (WD, FD, WN, FN). Each level is a disjoint part of the birds (water vs
  # flying, distance vs no-distance; see CLAUDE.md's three-axes section), so a
  # cell's density is the SUM over the levels the model was fitted with, as
  # dsm.pred() does for mgcv. Take the levels from the fit itself: they vary
  # by species, and predicting a subset would silently under-report.
  lv <- levels(fit$data$platform)
  if (is.null(lv)) stop("the fitted model has no platform factor")

  # Each grid row must be one cell in one season, or the per-cell sums below
  # would merge rows. key identifies a row by cell (coordinates or cell.col)
  # and season.
  key <- if (is.null(cell.col)) paste(grid$x, grid$y, grid$Season) else paste(grid[[cell.col]], grid$Season)
  if (anyDuplicated(key)) stop("grid has more than one row per cell x season")

  # g is the grid stacked length(lv) times: rows 1..n are level 1, rows
  # n+1..2n level 2, and so on (rep(times = ) repeats the whole row sequence;
  # rep(each = ) gives each level for a whole block of n rows, matching it).
  g <- grid[rep(seq_len(nrow(grid)), times = length(lv)), , drop = FALSE]
  g$platform <- factor(rep(lv, each = nrow(grid)), levels = lv)
  g <- prepare.sdmtmb.data(g, season.levels)
  # prepare.sdmtmb.data() drops unused platform levels. All levels are used
  # here, but re-impose the fit's exact level set and order anyway, because
  # predict() builds the design matrix from the factor's levels and a
  # mismatch would shift columns.
  g$platform <- factor(as.character(g$platform), levels = lv)

  # row_of[i] is the original grid row of stacked row i: 1..n, 1..n, ...
  # rowsum(..., row_of) then adds the platform copies of each cell back into
  # one value per original grid row.
  row_of <- rep(seq_len(nrow(grid)), times = length(lv))

  # Offset 0 for every prediction row. sdmTMB would otherwise look for the
  # offset column (log segment area), which the grid does not have. With
  # offset 0, exp(prediction) is a density (birds per km^2).
  off0 <- rep(0, nrow(g))
  season <- as.character(grid$Season)

  # ---- Error handling for the uncertainty steps ----
  # The plug-in surface below needs only the point estimates, but the index
  # and the draws need the joint precision matrix, and can fail (e.g. a
  # Hessian that is not positive definite). step() runs one of them, and on an
  # error records "<step name>: <message>" in `notes` and returns NULL, so a
  # failed step leaves its columns NA instead of losing the whole surface.
  # `<<-` assigns to `notes` in this function's environment, not inside step().
  notes <- character()
  step <- function(name, expr) tryCatch(expr, error = function(e) {
    notes <<- c(notes, sprintf("%s: %s", name, conditionMessage(e)))
    NULL
  })

  # ---- Plug-in surface ----
  # type = "response" gives exp(linear predictor) in p$est: the density at the
  # parameter estimates. For a log-normal quantity this is the median, not the
  # mean, so these cells sum to less than the bias-corrected total.
  # Abundance per stacked row = density x cell area (km^2; coastal cells are
  # smaller than the nominal cell). rowsum() sums the platform copies of each
  # cell: it groups by row_of and returns a 1-column matrix in row_of order
  # (reorder = TRUE sorts the groups, which here is 1..n, the grid order), so
  # [, 1] gives a vector aligned with grid.
  p <- stats::predict(fit, newdata = g, offset = off0, type = "response")
  n_cell <- rowsum(p$est * g$area, row_of, reorder = TRUE)[, 1]
  # Free the prediction (one row per stacked grid row, several columns):
  # hundreds of thousands of rows on Atl IMRP.
  rm(p)

  # ---- Seasonal totals with get_index() ----
  # get_index() sums density x area over the rows of each time value (season)
  # and gives a standard error for the log total by the generalised delta
  # method. With bias_correct = TRUE, TMB's epsilon method corrects the total
  # from the median towards the mean (needed: uncorrected totals ran
  # 1.3-1.65x low on NL ATPU, trial 5). It is slower: about 6 min on NL, more
  # for delta families. The plain version is kept for the mgcv comparison.
  # newdata/offset/area are passed directly; predict(return_tmb_object = TRUE)
  # is deprecated in sdmTMB 1.2.0.
  index <- function(bc) {
    ix <- sdmTMB::get_index(fit, newdata = g, area = g$area, offset = off0, bias_correct = bc)
    # get_index() labels its rows by season_num (1, 2, ...); translate back.
    ix$Season <- season.levels[ix$season_num]
    ix
  }
  # Each version only if asked for (bias.correct can hold TRUE, FALSE or both);
  # an index not asked for, or one that failed, is NULL.
  ix_bc <- if (TRUE %in% bias.correct) step("get_index bias-corrected", index(TRUE))
  ix_plain <- if (FALSE %in% bias.correct) step("get_index", index(FALSE))

  # ---- Optional draws: per-cell and per-total intervals ----
  # Default: no draws, so the interval columns are NA.
  lo <- hi <- rep(NA_real_, nrow(grid))
  tot_d <- NULL
  if (nsim > 0) {
    # predict(nsim = ) draws parameters from their joint (Gaussian)
    # approximation and returns a matrix of linear predictors: one row per
    # stacked grid row, one column per draw (log density, since offset is 0).
    draws <- step("draws", stats::predict(fit, newdata = g, offset = off0, nsim = nsim))
    if (!is.null(draws)) {
      # exp() to density, x area to abundance, then sum the platform copies of
      # each cell (rowsum works on every column at once): Nd is cells x draws.
      # Platform levels must be summed within each draw, before any quantile,
      # because a quantile of a sum is not the sum of quantiles.
      Nd <- rowsum(exp(draws) * g$area, row_of, reorder = TRUE)
      rm(draws)
      # Per cell, the 2.5% and 97.5% quantiles across draws. apply() over rows
      # returns a 2 x cells matrix, so t() makes it cells x 2.
      q <- t(apply(Nd, 1, stats::quantile, probs = c(0.025, 0.975)))
      lo <- q[, 1]
      hi <- q[, 2]
      # Seasonal total in each draw: sum the cells of each season.
      # rowsum(Nd, season) is seasons x draws, with seasons as row names.
      tot_d <- rowsum(Nd, season)
    }
  }

  # ---- Assemble the per-season table ----
  # pick() reads one column of an index for season s, or NA if that index was
  # not computed or failed (NULL).
  pick <- function(ix, s, col) if (is.null(ix)) NA_real_ else ix[[col]][ix$Season == s]
  # One row per season present in the grid, in season.levels order.
  # N/se_log: bias-corrected total and the standard error of its log;
  # N_plain/se_log_plain: the same without bias correction (equal to the sum
  # of the plug-in cells); N_median_sum: the plug-in cells summed directly;
  # lo_draws/hi_draws: 95% interval of the total from the draws.
  totals <- purrr::map_dfr(intersect(season.levels, unique(season)), function(s) {
    tibble::tibble(Season = s,
                   N = pick(ix_bc, s, "est"), se_log = pick(ix_bc, s, "se"),
                   N_plain = pick(ix_plain, s, "est"), se_log_plain = pick(ix_plain, s, "se"),
                   N_median_sum = sum(n_cell[season == s]),
                   lo_draws = if (is.null(tot_d)) NA_real_ else unname(stats::quantile(tot_d[s, ], 0.025)),
                   hi_draws = if (is.null(tot_d)) NA_real_ else unname(stats::quantile(tot_d[s, ], 0.975)))
  })
  # ---- Per-cell table ----
  # The rescaled surface: each season's plug-in cells multiplied by one factor,
  # bias-corrected total / sum of plug-in cells, so it keeps the plug-in map's
  # shape but sums to the bias-corrected total. The factor is per season, so
  # match() looks up each cell's season in `totals`. If the bias-corrected
  # index failed, N is NA, and so is N_rescaled.
  scale <- totals$N / totals$N_median_sum
  cells <- tibble::tibble(Season = factor(season, levels = season.levels), N = n_cell,
                          N_rescaled = n_cell * scale[match(season, totals$Season)], lo = lo, hi = hi)
  # Identify each cell as the caller did: by its coordinates, or by cell.col.
  # setNames() names the column after cell.col, avoiding rlang's `:=`.
  # Every vector above is in grid row order, so the id columns line up.
  id <- if (is.null(cell.col)) tibble::tibble(x = grid$x, y = grid$y) else
    stats::setNames(tibble::tibble(grid[[cell.col]]), cell.col)
  cells <- dplyr::bind_cols(id, cells)
  list(cells = cells, totals = totals, notes = notes,
       pred_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
