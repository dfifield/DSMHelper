# Over-water ("geodesic") coordinates for a spatial smooth that respects land
# barriers - Atlantic-DSM issue #127.
#
# s(x, y) smooths on straight-line distance, so density surveyed on one side of
# a peninsula leaks onto the other: Northumberland Strait into the Bay of Fundy,
# Placentia Bay into Trinity Bay across the Avalon isthmus. These functions
# follow Miller & Wood (2014, "Finite area smoothing with generalized distance
# splines", Environ. Ecol. Stat. 21:715-731):
#
#   1. shortest over-water path lengths on a fine water raster
#      (build.water.raster(), landmark.distances());
#   2. landmark MDS (de Silva & Tenenbaum 2004) to coordinates geo1, geo2, ...
#      whose Euclidean distances approximate those paths (lmds.fit(),
#      lmds.embed());
#   3. a Duchon spline on those coordinates in place of s(x, y)
#      (geodesic.smooth(), geo.formula()).
#
# The coordinates are ordinary covariates, so the model stays an ordinary dsm:
# detection-function variance, prediction, cross-validation and the delta
# method are unchanged. create.geodesic.embedding() and add.geodesic.coords()
# are the two calls a pipeline makes.


#' Fine water raster: the physical sea around a study area
#'
#' Water is where `depth` is positive (below sea level), OR inside `area`.
#' Neither alone is right, measured on Atl IMRP:
#'
#' - A study-area polygon alone fails twice. A buffered one spreads over land
#'   and fills in straits, bays and peninsulas (a 25 km buffer erased
#'   Northumberland Strait and PEI). An unbuffered one has political holes (the
#'   French EEZ corridor south of Newfoundland) that would become false
#'   barriers.
#' - ETOPO alone misses inshore water: it has Lake Melville a few metres above
#'   sea level, and at 2 km it closes inner bays, leaving 329 prediction cells
#'   more than 4 km from water, 26 of them up to 186 km away.
#'
#' So `area` should be the EXACT (unbuffered) training polygon, drawn on the
#' coastline; it supplies the inshore water and the depth raster supplies the
#' sea outside it. Paths may then run through sea outside the study area -
#' through Cabot Strait or the Strait of Belle Isle - because that is how the
#' water connects.
#'
#' The lattice nests in `template`'s: `template` extended by `margin.km`, then
#' disaggregated by an integer factor.
#'
#' Water cut off from the main body has no finite path to anything. A detached
#' patch holding any of `need.xy` and covering at least `join.min.km2` is
#' water the study area includes but whose connection is narrower than the fine
#' cell - Lake Melville is a separate piece of the Atl IMRP polygon, and ETOPO
#' closes the Narrows at Rigolet - so it is joined to the main body along the
#' shortest straight gap, where the narrows should be. Every other detached
#' patch is dropped, NOT joined: on Atl IMRP joining by size alone cut
#' straight-line canals of 190-270 km across land, which is the bleed this
#' exists to remove; and a 1-cell sliver holding one prediction-cell centroid
#' would be joined across 21 km of land, where snapping the point
#' (snap.to.water()) moves it 21 km.
#'
#' @param depth SpatRaster of depth, positive below sea level, any CRS.
#' @param area sf polygon, the exact training study area.
#' @param template SpatRaster defining the lattice and CRS (genpredrast).
#' @param cell.km Target fine cell size, km.
#' @param margin.km How far past `template`'s extent water is kept, km.
#' @param need.xy Two-column matrix of points, in `template`'s CRS, that must be
#'   reachable - the prediction cells.
#' @param join.min.km2 Smallest detached patch that is joined, km^2.
#' @return SpatRaster: 1 for water, NA for land, with attribute `joins`: one row
#'   per join, giving the gap's two ends and its length, or NULL.
#' @examples
#' \dontrun{
#' water <- build.water.raster(depth, training.study.area, genpredrast, 2, 300,
#'                             as.matrix(predgrid[, c("x", "y")]))
#' }
#' @export
build.water.raster <- function(depth, area, template, cell.km, margin.km, need.xy,
                               join.min.km2 = 100) {
  checkmate::expect_class(depth, "SpatRaster")
  checkmate::expect_class(area, "sf")
  checkmate::expect_class(template, "SpatRaster")
  checkmate::expect_number(cell.km, lower = 0.1)
  checkmate::expect_number(margin.km, lower = 0)
  checkmate::expect_matrix(need.xy, mode = "numeric", ncols = 2, any.missing = FALSE)
  checkmate::expect_number(join.min.km2, lower = 0)

  res.km <- mean(terra::res(template)) / 1000
  fact <- max(1L, as.integer(round(res.km / cell.km)))
  fine <- terra::extend(terra::rast(template), ceiling(margin.km / res.km)) %>%
    terra::disagg(fact = fact)
  d <- terra::project(depth, fine, method = "bilinear")
  if (anyNA(terra::values(d)))
    warning("build.water.raster: ", sum(is.na(terra::values(d))), " fine cells have no ",
            "depth - the depth raster does not cover the margin. They are treated as land.",
            immediate. = TRUE)
  inside <- terra::rasterize(terra::vect(sf::st_transform(area, terra::crs(fine))),
                             fine, field = 1, background = 0)
  water <- terra::ifel((!is.na(d) & d > 0) | inside == 1, 1, NA)

  p <- terra::patches(water, directions = 8)
  sizes <- terra::freq(p)
  main <- sizes$value[which.max(sizes$count)]
  pv <- terra::values(p)[, 1]
  main.xy <- terra::xyFromCell(p, which(pv == main))
  need.patch <- unique(stats::na.omit(pv[terra::cellFromXY(p, need.xy)]))
  min.cells <- join.min.km2 / prod(terra::res(p) / 1000)
  joined <- sizes[sizes$value != main & sizes$value %in% need.patch &
                    sizes$count >= min.cells, ]
  joins <- NULL
  for (v in joined$value) {
    pc <- which(pv == v)
    if (!length(pc)) next   # already absorbed by an earlier canal
    pxy <- terra::xyFromCell(p, pc)
    nn <- FNN::get.knnx(main.xy, pxy, k = 1)
    i <- which.min(nn$nn.dist[, 1])
    from <- pxy[i, ]
    to <- main.xy[nn$nn.index[i, 1], ]
    n.step <- ceiling(sqrt(sum((to - from)^2)) / (mean(terra::res(p)) / 2)) + 1
    canal <- terra::cellFromXY(p, cbind(seq(from[1], to[1], length.out = n.step),
                                        seq(from[2], to[2], length.out = n.step)))
    pv[c(pc, canal)] <- main
    joins <- rbind(joins, data.frame(patch_cells = length(pc),
                                     gap_km = nn$nn.dist[i, 1] / 1000,
                                     x_from = from[1], y_from = from[2],
                                     x_to = to[1], y_to = to[2]))
    message(sprintf("build.water.raster: joined a %d-cell patch to the main body across a %.1f km gap",
                    length(pc), nn$nn.dist[i, 1] / 1000))
  }
  dropped <- sizes[sizes$value != main & !sizes$value %in% joined$value, ]
  message(sprintf(paste0("build.water.raster: %d x %d at %.2f km, %d water cells; ",
                         "%d patch(es) joined; %d cells in %d small patch(es) dropped"),
                  nrow(water), ncol(water), mean(terra::res(water)) / 1000,
                  sum(pv == main, na.rm = TRUE), nrow(joined), sum(dropped$count),
                  nrow(dropped)))
  out <- terra::rast(p)
  terra::values(out) <- ifelse(!is.na(pv) & pv == main, 1, NA)
  attr(out, "joins") <- joins
  out
}

#' Nearest water cell for each point
#'
#' A segment midpoint can sit on land (a straight segment cutting a headland),
#' and a clipped coastal prediction cell's centroid can fall outside it, so
#' points are snapped to the nearest water cell centre when not already on one.
#'
#' @param water SpatRaster from build.water.raster().
#' @param xy Two-column numeric matrix in the raster's CRS.
#' @return List: `cell` (water cell index per point) and `snap_km`.
#' @examples
#' \dontrun{
#' s <- snap.to.water(water, as.matrix(segdata[, c("x", "y")]))
#' }
#' @export
snap.to.water <- function(water, xy) {
  checkmate::expect_class(water, "SpatRaster")
  checkmate::expect_matrix(xy, mode = "numeric", ncols = 2, any.missing = FALSE)

  wv <- terra::values(water)[, 1]
  cell <- terra::cellFromXY(water, xy)
  on.water <- !is.na(cell) & !is.na(wv[cell])
  snap <- numeric(nrow(xy))
  if (any(!on.water)) {
    wc <- which(!is.na(wv))
    nn <- FNN::get.knnx(terra::xyFromCell(water, wc), xy[!on.water, , drop = FALSE], k = 1)
    cell[!on.water] <- wc[nn$nn.index[, 1]]
    snap[!on.water] <- nn$nn.dist[, 1] / 1000
  }
  list(cell = cell, snap_km = snap)
}

#' Spread landmarks evenly over the water inside a polygon
#'
#' k-means on the water cell centres, each centre then snapped to its nearest
#' water cell, so landmark density follows water area and every basin gets its
#' share. Only water inside `within` is used: the landmarks are where points
#' will be embedded, and sea outside the training area is there for paths to
#' cross, not to be described.
#'
#' @param water SpatRaster from build.water.raster().
#' @param within sf polygon to place landmarks in.
#' @param n Number of landmarks.
#' @param seed Random seed.
#' @return Integer vector of water cell indices.
#' @examples
#' \dontrun{
#' lm.cells <- choose.landmarks(water, training.study.area, 400, 127)
#' }
#' @export
choose.landmarks <- function(water, within, n, seed) {
  checkmate::expect_class(water, "SpatRaster")
  checkmate::expect_class(within, "sf")
  checkmate::expect_count(n, positive = TRUE)
  checkmate::expect_count(seed)

  inside <- terra::rasterize(terra::vect(sf::st_transform(within, terra::crs(water))),
                             water, field = 1)
  wc <- which(!is.na(terra::values(water)[, 1]) & !is.na(terra::values(inside)[, 1]))
  if (length(wc) < n)
    stop("choose.landmarks: only ", length(wc), " water cells inside `within`, fewer than n = ", n)
  wxy <- terra::xyFromCell(water, wc)
  # Convergence does not matter here - the centres only need to be spread out -
  # so the iteration-limit warning is not worth reporting.
  km <- withr::with_seed(seed, suppressWarnings(
    stats::kmeans(wxy, centers = n, iter.max = 50, algorithm = "Lloyd")))
  unique(wc[FNN::get.knnx(wxy, km$centers, k = 1)$nn.index[, 1]])
}

#' Over-water path length from each landmark to every water cell
#'
#' `terra::gridDist()` from each landmark in turn, travelling through water
#' cells only; land (NA) cells are barriers.
#'
#' @param water SpatRaster from build.water.raster().
#' @param lm.cells Landmark cell indices.
#' @return SpatRaster, one layer per landmark, km; NA on land.
#' @examples
#' \dontrun{
#' lmdist <- landmark.distances(water, lm.cells)
#' }
#' @export
landmark.distances <- function(water, lm.cells) {
  checkmate::expect_class(water, "SpatRaster")
  checkmate::expect_integerish(lm.cells, lower = 1, any.missing = FALSE, min.len = 1)

  t0 <- Sys.time()
  layers <- lapply(seq_along(lm.cells), function(i) {
    # `* 1` makes a new raster; a plain copy would share terra's pointer
    tgt <- water * 1
    tgt[lm.cells[i]] <- 0
    if (i %% 50 == 0)
      message(sprintf("landmark.distances: %d / %d, %.1f min", i, length(lm.cells),
                      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    terra::gridDist(tgt, target = 0) / 1000
  })
  out <- terra::rast(layers)
  names(out) <- paste0("lm", seq_along(lm.cells))
  out
}

#' Landmark MDS: fit on the landmarks
#'
#' Classical MDS of the landmark-to-landmark over-water distances (de Silva &
#' Tenenbaum 2004). Only positive eigenvalues can be used. Over-water distance
#' is not Euclidean, so some eigenvalues are negative, and the share of the
#' spectrum they hold is reported as a measure of how non-Euclidean the water
#' is (9.5% on Atl IMRP).
#'
#' @param D L x L matrix of over-water distances between landmarks, km.
#' @param max.dims Largest embedding dimension wanted.
#' @return List: `mu` (column means of squared distances), `vectors`, `values`.
#' @examples
#' \dontrun{
#' fit <- lmds.fit(as.matrix(terra::extract(lmdist, lm.cells)), 6)
#' }
#' @export
lmds.fit <- function(D, max.dims) {
  checkmate::expect_matrix(D, mode = "numeric", any.missing = FALSE)
  checkmate::expect_count(max.dims, positive = TRUE)
  if (nrow(D) != ncol(D)) stop("lmds.fit: D must be square")

  D2 <- ((D + t(D)) / 2)^2
  n <- nrow(D2)
  H <- diag(n) - 1 / n
  e <- eigen(-0.5 * H %*% D2 %*% H, symmetric = TRUE)
  if (any(e$values[seq_len(max.dims)] <= 0))
    stop("lmds.fit: fewer than ", max.dims, " positive eigenvalues")
  neg <- sum(abs(e$values[e$values < 0])) / sum(abs(e$values))
  message(sprintf(paste0("lmds.fit: %d landmarks; negative eigenvalues hold %.1f%% of the ",
                         "spectrum; first %d positive hold %.1f%% of the positive part"),
                  n, 100 * neg, max.dims,
                  100 * sum(e$values[seq_len(max.dims)]) / sum(e$values[e$values > 0])))
  list(mu = colMeans(D2), vectors = e$vectors[, seq_len(max.dims), drop = FALSE],
       values = e$values[seq_len(max.dims)])
}

#' Landmark MDS: place points from their distances to the landmarks
#'
#' The triangulation of de Silva & Tenenbaum (2004): a point with squared
#' distances delta to the landmarks sits at 1/2 L# (mu - delta).
#'
#' @param fit Output of lmds.fit().
#' @param Dp N x L matrix of over-water distances from each point to each
#'   landmark, km.
#' @param d Number of dimensions.
#' @return N x d matrix `geo1`, ..., in metres to match x and y.
#' @examples
#' \dontrun{
#' Y <- lmds.embed(fit, Dp, 4)
#' }
#' @export
lmds.embed <- function(fit, Dp, d) {
  checkmate::expect_list(fit)
  checkmate::expect_matrix(Dp, mode = "numeric", any.missing = FALSE,
                           ncols = length(fit$mu))
  checkmate::expect_count(d, positive = TRUE)
  if (d > length(fit$values))
    stop("lmds.embed: d = ", d, " but the fit holds ", length(fit$values), " dimensions")

  Lsharp <- sweep(fit$vectors[, seq_len(d), drop = FALSE], 2,
                  sqrt(fit$values[seq_len(d)]), "/")
  Y <- 0.5 * sweep(-(Dp^2), 2, fit$mu, "+") %*% Lsharp
  colnames(Y) <- paste0("geo", seq_len(d))
  Y * 1000
}

#' Path lengths from points to every landmark
#'
#' @param lmdist SpatRaster from landmark.distances().
#' @param water SpatRaster from build.water.raster().
#' @param xy Two-column coordinate matrix.
#' @return List: `D` (N x L, km) and `snap_km`.
#' @examples
#' \dontrun{
#' pd <- point.landmark.dist(lmdist, water, as.matrix(predgrid[, c("x", "y")]))
#' }
#' @export
point.landmark.dist <- function(lmdist, water, xy) {
  checkmate::expect_class(lmdist, "SpatRaster")
  checkmate::expect_class(water, "SpatRaster")
  checkmate::expect_matrix(xy, mode = "numeric", ncols = 2, any.missing = FALSE)

  s <- snap.to.water(water, xy)
  D <- as.matrix(terra::extract(lmdist, s$cell))
  if (anyNA(D))
    stop("point.landmark.dist: NA path length - a point snapped to a cell with no water path")
  list(D = D, snap_km = s$snap_km)
}

#' Straight-line distance from points to landmarks, km
#'
#' @param xy N x 2 coordinates, metres.
#' @param lm.xy L x 2 landmark coordinates, metres.
#' @return N x L matrix, km.
#' @examples
#' euclid.to.landmarks(cbind(0, 0), cbind(3000, 4000))
#' @export
euclid.to.landmarks <- function(xy, lm.xy) {
  checkmate::expect_matrix(xy, mode = "numeric", ncols = 2)
  checkmate::expect_matrix(lm.xy, mode = "numeric", ncols = 2)
  sqrt(outer(xy[, 1], lm.xy[, 1], "-")^2 + outer(xy[, 2], lm.xy[, 2], "-")^2) / 1000
}

#' How well geo-space distance reproduces the over-water path
#'
#' Compares, for every point and every landmark, the Euclidean distance between
#' their embedded coordinates with the true over-water path length.
#'
#' @param fit Output of lmds.fit().
#' @param Dp N x L path lengths, km.
#' @param E N x L straight-line distances, km.
#' @param dims Dimensions to assess.
#' @param lm.geo Landmark coordinates, as lmds.embed() of the landmarks, metres.
#' @param radius.km Pairs closer than this in a straight line are "local".
#' @param detour A local pair whose path is more than `detour` times its
#'   straight line is a "barrier" pair - the ones the method exists for.
#' @return Tibble, one row per dimension: `r2_all`, `med_rel_err_local`,
#'   `med_rel_err_barrier`, `n_barrier_pairs`.
#' @examples
#' \dontrun{
#' embedding.fidelity(fit, pd$D, E, 2:6, lm.geo, 100, 2)
#' }
#' @export
embedding.fidelity <- function(fit, Dp, E, dims, lm.geo, radius.km, detour) {
  checkmate::expect_matrix(Dp, mode = "numeric")
  checkmate::expect_matrix(E, mode = "numeric", nrows = nrow(Dp), ncols = ncol(Dp))
  checkmate::expect_integerish(dims, lower = 1, min.len = 1)
  checkmate::expect_number(radius.km, lower = 0)
  checkmate::expect_number(detour, lower = 1)

  local <- E < radius.km
  barrier <- local & Dp > detour * E
  purrr::map_dfr(dims, function(d) {
    Y <- lmds.embed(fit, Dp, d) / 1000
    L <- lm.geo[, seq_len(d), drop = FALSE] / 1000
    G <- sqrt(pmax(0, outer(rowSums(Y^2), rep(1, nrow(L))) +
                     outer(rep(1, nrow(Y)), rowSums(L^2)) - 2 * Y %*% t(L)))
    rel <- abs(G / pmax(Dp, 1e-6) - 1)
    tibble::tibble(
      dims = d,
      r2_all = stats::cor(as.vector(G), as.vector(Dp))^2,
      med_rel_err_local = stats::median(rel[local & Dp > 0]),
      med_rel_err_barrier = if (any(barrier)) stats::median(rel[barrier]) else NA_real_,
      n_barrier_pairs = sum(barrier))
  })
}

#' Barrier exposure of each point
#'
#' Where a straight-line smooth CAN bleed: the share of landmarks within
#' `radius.km` in a straight line whose over-water path is more than `detour`
#' times as long.
#'
#' @param Dp N x L path lengths, km.
#' @param E N x L straight-line distances, km.
#' @param radius.km Straight-line radius, km.
#' @param detour Path-to-straight-line ratio that counts as separated by land.
#' @return Numeric vector, NA where no landmark is in range.
#' @examples
#' barrier.exposure(matrix(c(10, 300), 1), matrix(c(10, 40), 1), 100, 2)
#' @export
barrier.exposure <- function(Dp, E, radius.km, detour) {
  checkmate::expect_matrix(Dp, mode = "numeric")
  checkmate::expect_matrix(E, mode = "numeric", nrows = nrow(Dp), ncols = ncol(Dp))
  checkmate::expect_number(radius.km, lower = 0)
  checkmate::expect_number(detour, lower = 1)
  near <- E < radius.km & E > 0
  n <- rowSums(near)
  ifelse(n > 0, rowSums(near & Dp > detour * E) / n, NA_real_)
}

#' The over-water spatial smooth term, as formula text
#'
#' A Duchon spline in `dims` dimensions with m = c(1, (dims - 1) / 2). A Duchon
#' spline needs m + s > dims / 2 and s < dims / 2; this choice has
#' m + s = (dims + 1) / 2 and a null space of constants only, so after a
#' factor-`by` centring nothing is left unpenalised - the property
#' `bs = "ts"` gives `s(x, y)`. For dims = 2 it is the familiar m = c(1, 0.5).
#'
#' @param dims Embedding dimension.
#' @param k Basis size.
#' @param by Name of the `by` variable, or NULL.
#' @return Character string, e.g.
#'   `s(geo1, geo2, geo3, geo4, bs = "ds", m = c(1, 1.5), k = 100, by = Season)`.
#' @examples
#' geodesic.smooth(4, 100, "Season")
#' @export
geodesic.smooth <- function(dims, k, by = NULL) {
  checkmate::expect_count(dims, positive = TRUE)
  checkmate::expect_count(k, positive = TRUE)
  checkmate::expect_string(by, null.ok = TRUE)
  sprintf('s(%s, bs = "ds", m = c(1, %s), k = %d%s)',
          paste0("geo", seq_len(dims), collapse = ", "),
          format((dims - 1) / 2), as.integer(k),
          if (is.null(by)) "" else paste0(", by = ", by))
}

#' Replace a formula's s(x, y) term with the over-water smooth
#'
#' Keeps the term's `k` and `by`, so the over-water smooth has the same basis
#' size and per-season structure as the term it replaces.
#'
#' @param f Formula with exactly one `s(x, y, ...)` term.
#' @param dims Embedding dimension.
#' @return The new formula.
#' @examples
#' geo.formula(y ~ Season + s(x, y, bs = "ts", k = 100, by = Season), 4)
#' @export
geo.formula <- function(f, dims) {
  checkmate::expect_formula(f)
  checkmate::expect_count(dims, positive = TRUE)
  labs <- attr(stats::terms(f), "term.labels")
  is.xy <- vapply(labs, function(l) {
    cl <- str2lang(l)
    is.call(cl) && identical(cl[[1]], as.name("s")) && length(cl) >= 3 &&
      identical(cl[[2]], as.name("x")) && identical(cl[[3]], as.name("y"))
  }, logical(1))
  if (sum(is.xy) != 1)
    stop("geo.formula: expected exactly one s(x, y, ...) term, found ", sum(is.xy))
  old <- str2lang(labs[is.xy])
  k <- if (is.null(old$k)) 100 else eval(old$k)
  by <- if (is.null(old$by)) NULL else deparse(old$by)
  stats::update(f, stats::as.formula(paste(". ~ . -", labs[is.xy], "+",
                                           geodesic.smooth(dims, k, by))))
}

#' Build and save the over-water embedding for a SubProject
#'
#' The once-per-SubProject part: water raster, landmarks, path lengths from
#' every landmark, and the landmark MDS. Writes `water.tif`,
#' `landmark_dist.tif` and `embedding.rds` to `dir`; load them back with
#' load.geodesic.embedding().
#'
#' @param depth SpatRaster of depth, positive below sea level.
#' @param area sf polygon, the exact training study area.
#' @param template SpatRaster, genpredrast.
#' @param need.xy Prediction-cell coordinates; see build.water.raster().
#' @param dir Directory to write to (created if absent).
#' @param cell.km Fine cell size, km.
#' @param margin.km Water kept this far past `template`, km.
#' @param n.landmarks Number of landmarks.
#' @param seed Seed for landmark placement.
#' @param max.dims Dimensions kept in the MDS fit.
#' @return Invisibly, the list load.geodesic.embedding() returns.
#' @examples
#' \dontrun{
#' create.geodesic.embedding(depth, training.study.area, genpredrast,
#'                           as.matrix(predgrid[, c("x", "y")]), geodesicDir)
#' }
#' @export
create.geodesic.embedding <- function(depth, area, template, need.xy, dir,
                                      cell.km = 2, margin.km = 300,
                                      n.landmarks = 400, seed = 127, max.dims = 6) {
  checkmate::expect_string(dir)
  checkmate::expect_count(n.landmarks, positive = TRUE)
  checkmate::expect_count(max.dims, positive = TRUE)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  water <- build.water.raster(depth, area, template, cell.km, margin.km, need.xy)
  lm.cells <- choose.landmarks(water, area, n.landmarks, seed)
  lmdist <- landmark.distances(water, lm.cells)
  fit <- lmds.fit(as.matrix(terra::extract(lmdist, lm.cells)), max.dims)
  terra::writeRaster(water, file.path(dir, "water.tif"), overwrite = TRUE)
  terra::writeRaster(lmdist, file.path(dir, "landmark_dist.tif"), overwrite = TRUE,
                     datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  saveRDS(list(cells = lm.cells, xy = terra::xyFromCell(water, lm.cells), fit = fit,
               joins = attr(water, "joins"),
               settings = list(cell.km = cell.km, margin.km = margin.km,
                               n.landmarks = n.landmarks, seed = seed, max.dims = max.dims)),
          file.path(dir, "embedding.rds"))
  message(sprintf("create.geodesic.embedding: %d landmarks, %.1f min, written to %s",
                  length(lm.cells), as.numeric(difftime(Sys.time(), t0, units = "mins")), dir))
  invisible(load.geodesic.embedding(dir))
}

#' Load an embedding written by create.geodesic.embedding()
#'
#' @param dir The directory it was written to.
#' @return List: `water`, `lmdist` (SpatRasters) and everything in
#'   `embedding.rds` (`cells`, `xy`, `fit`, `joins`, `settings`).
#' @examples
#' \dontrun{
#' emb <- load.geodesic.embedding(geodesicDir)
#' }
#' @export
load.geodesic.embedding <- function(dir) {
  checkmate::expect_directory_exists(dir)
  files <- file.path(dir, c("water.tif", "landmark_dist.tif", "embedding.rds"))
  checkmate::expect_file_exists(files)
  c(list(water = terra::rast(files[1]), lmdist = terra::rast(files[2])), readRDS(files[3]))
}

#' Add over-water coordinates geo1, ..., geo<dims> to a data frame
#'
#' Each distinct (x, y) is embedded once, so segment data carrying one copy per
#' ddftype costs no more than one copy. Points off the water raster are snapped
#' to the nearest water cell first; the distances moved are summarised in a
#' message and returned in attribute `snap_km`, one per distinct position.
#'
#' @param df Data frame with columns `x` and `y` in the embedding's CRS.
#' @param emb Output of load.geodesic.embedding().
#' @param dims Number of coordinates to add.
#' @return `df` with `geo1`, ..., `geo<dims>` (metres) added or replaced.
#' @examples
#' \dontrun{
#' segdata <- add.geodesic.coords(segdata, emb, GEODESIC_DIMS)
#' }
#' @export
add.geodesic.coords <- function(df, emb, dims) {
  checkmate::expect_data_frame(df, min.rows = 1)
  checkmate::expect_names(names(df), must.include = c("x", "y"))
  checkmate::expect_list(emb)
  checkmate::expect_count(dims, positive = TRUE)

  xy <- cbind(as.numeric(df$x), as.numeric(df$y))
  key <- paste(xy[, 1], xy[, 2])
  u <- !duplicated(key)
  pd <- point.landmark.dist(emb$lmdist, emb$water, xy[u, , drop = FALSE])
  Y <- lmds.embed(emb$fit, pd$D, dims)[match(key, key[u]), , drop = FALSE]
  message(sprintf("add.geodesic.coords: %d positions; %d snapped onto water (median %.1f km, max %.1f km)",
                  sum(u), sum(pd$snap_km > 0),
                  if (any(pd$snap_km > 0)) stats::median(pd$snap_km[pd$snap_km > 0]) else 0,
                  max(pd$snap_km)))
  old <- grep("^geo[0-9]+$", names(df), value = TRUE)
  if (length(old)) df <- df[, setdiff(names(df), old), drop = FALSE]
  out <- cbind(df, as.data.frame(Y))
  attr(out, "snap_km") <- pd$snap_km
  out
}
