# Over-water coordinates (Atlantic-DSM #127), on a synthetic coast: a 300 km
# square of sea with a 20 km wide peninsula reaching down from the north edge to
# y = 60 km, and a 20 x 20 km lake cut off from the sea by 10 km of land.

toy.coast <- function() {
  crs <- "EPSG:3857"
  sq <- function(x0, x1, y0, y1)
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  sea <- sq(0, 3e5, 0, 3e5)
  peninsula <- sq(1.4e5, 1.6e5, 6e4, 3.1e5)
  lake.land <- sq(2.2e5, 2.6e5, 2.2e5, 2.6e5)     # land block ...
  lake <- sq(2.3e5, 2.5e5, 2.3e5, 2.5e5)          # ... with the lake inside it
  area <- sf::st_difference(sf::st_difference(sea, peninsula), lake.land)
  area <- sf::st_sf(geometry = sf::st_sfc(sf::st_union(area, lake), crs = crs))
  template <- terra::rast(xmin = 0, xmax = 3e5, ymin = 0, ymax = 3e5, res = 1e4, crs = crs)
  # depth from the polygon itself, everything outside it land, so the margin
  # adds nothing and the lake is detached
  depth <- terra::rasterize(terra::vect(area), terra::disagg(template, 5),
                            field = 50, background = -10)
  list(area = area, template = template, depth = depth)
}

test_that("water raster: the peninsula is land and a needed lake is joined", {
  tc <- toy.coast()
  lake.pt <- cbind(2.4e5, 2.4e5)
  w <- suppressMessages(build.water.raster(tc$depth, tc$area, tc$template,
                                           cell.km = 2, margin.km = 0, need.xy = lake.pt))
  expect_true(is.na(terra::extract(w, cbind(1.5e5, 2e5))[1, 1]))      # peninsula
  expect_identical(terra::extract(w, lake.pt)[1, 1], 1)               # lake kept
  expect_equal(nrow(attr(w, "joins")), 1)
  expect_lt(attr(w, "joins")$gap_km, 14)   # 10 km of land, ~12 between cell centres

  # not needed: the lake is dropped, not joined
  w2 <- suppressMessages(build.water.raster(tc$depth, tc$area, tc$template, 2, 0,
                                            need.xy = cbind(5e4, 5e4)))
  expect_true(is.na(terra::extract(w2, lake.pt)[1, 1]))
  expect_null(attr(w2, "joins"))
})

test_that("embedding separates the two sides of the peninsula and keeps one side intact", {
  tc <- toy.coast()
  w <- suppressMessages(build.water.raster(tc$depth, tc$area, tc$template, 2, 0,
                                           need.xy = cbind(5e4, 5e4)))
  lm.cells <- choose.landmarks(w, tc$area, 60, 1)
  lmdist <- suppressMessages(landmark.distances(w, lm.cells))
  fit <- suppressMessages(lmds.fit(as.matrix(terra::extract(lmdist, lm.cells)), 4))
  pts <- rbind(c(1.3e5, 2.5e5), c(1.7e5, 2.5e5), c(1.3e5, 2.0e5))
  pd <- point.landmark.dist(lmdist, w, pts)
  Y <- lmds.embed(fit, pd$D, 3) / 1000
  across <- sqrt(sum((Y[1, ] - Y[2, ])^2))    # 40 km apart, ~420 km by water
  same <- sqrt(sum((Y[1, ] - Y[3, ])^2))      # 50 km apart on the same side
  expect_gt(across, 300)
  expect_lt(abs(same - 50), 15)
})

test_that("landmark MDS reproduces a Euclidean configuration exactly", {
  set.seed(3)
  X <- matrix(runif(40), 20) * 100           # km
  D <- as.matrix(stats::dist(X))
  fit <- suppressMessages(lmds.fit(D, 2))
  Y <- lmds.embed(fit, D, 2) / 1000
  expect_equal(as.matrix(stats::dist(Y)), D, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("points on land are snapped to the nearest water", {
  tc <- toy.coast()
  w <- suppressMessages(build.water.raster(tc$depth, tc$area, tc$template, 2, 0,
                                           need.xy = cbind(5e4, 5e4)))
  s <- snap.to.water(w, rbind(c(5e4, 5e4), c(1.5e5, 2e5)))
  expect_identical(s$snap_km[1], 0)
  expect_gt(s$snap_km[2], 8)
  expect_lt(s$snap_km[2], 12)
})

test_that("add.geodesic.coords embeds each position once and replaces old columns", {
  tc <- toy.coast()
  dir <- withr::local_tempdir()
  emb <- suppressMessages(create.geodesic.embedding(
    tc$depth, tc$area, tc$template, need.xy = cbind(5e4, 5e4), dir = dir,
    margin.km = 0, n.landmarks = 40, max.dims = 4))
  df <- data.frame(x = c(5e4, 5e4, 2e5), y = c(5e4, 5e4, 1e5), geo1 = 0)
  out <- suppressMessages(add.geodesic.coords(df, emb, 3))
  expect_named(out, c("x", "y", "geo1", "geo2", "geo3"))
  expect_identical(unname(unlist(out[1, 3:5])), unname(unlist(out[2, 3:5])))
  expect_false(identical(out$geo1[1], 0))
  expect_identical(load.geodesic.embedding(dir)$settings$n.landmarks, 40)
})

test_that("the smooth term and the formula rewrite", {
  expect_identical(geodesic.smooth(4, 100, "Season"),
                   's(geo1, geo2, geo3, geo4, bs = "ds", m = c(1, 1.5), k = 100, by = Season)')
  expect_identical(geodesic.smooth(2, 50), 's(geo1, geo2, bs = "ds", m = c(1, 0.5), k = 50)')
  f <- geo.formula(y ~ Season + s(x, y, bs = "ts", k = 100, by = Season) + s(depth) +
                     offset(off.set), 3)
  labs <- attr(stats::terms(f), "term.labels")
  expect_false(any(grepl("s(x, y", labs, fixed = TRUE)))
  expect_true(any(grepl("geo3", labs) & grepl("by = Season", labs) & grepl("k = 100", labs)))
  expect_true(grepl("offset(off.set)", paste(deparse(f), collapse = ""), fixed = TRUE))
  expect_error(geo.formula(y ~ s(depth), 3), "exactly one")
})

test_that("barrier exposure counts only land-separated neighbours within range", {
  Dp <- matrix(c(10, 300, 500), 1)
  E <- matrix(c(10, 40, 150), 1)
  expect_equal(barrier.exposure(Dp, E, radius.km = 100, detour = 2), 0.5)
  expect_true(is.na(barrier.exposure(matrix(500, 1), matrix(150, 1), 100, 2)))
})
