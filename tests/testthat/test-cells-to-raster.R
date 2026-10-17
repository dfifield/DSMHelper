# Issue #90. Prediction rasters were written by rasterising the clipped cell
# polygons onto a template built from their extent: a different lattice from
# genpredrast, and polygon rasterisation drops any coastal sliver that misses its
# cell's centre. On Atl IMRP that lost 1,122 of 16,004 cells. These tests build a
# small grid the way 00.03 does - a deliberately non-round lattice, cells clipped
# by a "coast", centroids of the clipped pieces as x/y - and check that every
# cell lands in its own raster cell and that density x area gives back the
# abundance.

# A 6 x 5 lattice whose cells are not a round size, as genpredrast's are not
# (#91), clipped by a diagonal coast so some cells become slivers.
clipped.grid <- function() {
  template <- terra::rast(xmin = 0, xmax = 59950, ymin = 0, ymax = 49960,
                          ncols = 6, nrows = 5, crs = "EPSG:3857")
  terra::values(template) <- 1
  cells <- sf::st_as_sf(terra::as.polygons(template, dissolve = FALSE))
  coast <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(-1, -1), c(60000, -1), c(60000, 30000), c(-1, 51000), c(-1, -1)))),
    crs = sf::st_crs(cells))
  grid <- suppressWarnings(sf::st_intersection(cells, coast))
  grid <- grid[as.numeric(sf::st_area(grid)) > 1, ]
  ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(grid)))
  grid$x <- ctr[, 1]
  grid$y <- ctr[, 2]
  grid$area <- as.numeric(sf::st_area(grid)) / 1e6
  list(template = template, grid = grid)
}

test_that("every clipped cell lands in its own raster cell", {
  g <- clipped.grid()
  r <- cells.to.raster(g$grid, "area", g$template)
  expect_true(terra::compareGeom(r, g$template))
  expect_equal(sum(is.finite(terra::values(r))), nrow(g$grid))
  # The slivers are the point: polygon rasterisation loses them.
  lost <- terra::rasterize(terra::vect(g$grid), g$template, field = "area")
  expect_lt(sum(is.finite(terra::values(lost))), nrow(g$grid))
})

test_that("density x area summed over the raster is the abundance", {
  g <- clipped.grid()
  g$grid$NHat <- seq_len(nrow(g$grid)) * 3.7
  g$grid$Dens <- g$grid$NHat / g$grid$area
  dens <- cells.to.raster(g$grid, "Dens", g$template)
  area <- cells.to.raster(g$grid, "area", g$template)
  expect_equal(sum(terra::values(dens * area), na.rm = TRUE), sum(g$grid$NHat))
})

test_that("two rows in one cell, or a row off the template, stop", {
  g <- clipped.grid()
  expect_error(cells.to.raster(rbind(g$grid, g$grid[1, ]), "area", g$template),
               "share a template cell")
  off <- g$grid[1, ]
  off$x <- 1e7
  expect_error(cells.to.raster(off, "area", g$template), "outside the template")
})

test_that("the abundance summaries weight density by cell area", {
  g <- clipped.grid()
  g$grid$NHat <- 10
  g$grid$Dens <- g$grid$NHat / g$grid$area
  dens <- cells.to.raster(g$grid, "Dens", g$template)
  area <- cells.to.raster(g$grid, "area", g$template)
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(dens, f)

  conc <- summarise.prediction.concentration(f, "SP", "Summer", area.raster = area,
                                             top_n = 1)
  expect_equal(conc$total, sum(g$grid$NHat), tolerance = 1e-5)
  # Every cell holds 10 birds, so one cell is 1/n of the total whatever its area.
  expect_equal(conc$pct_top_n, 100 / nrow(g$grid), tolerance = 1e-5)

  cells <- tibble::tibble(cell = which(is.finite(terra::values(dens)[, 1])),
                          type = "analogue", mic_name = "depth")
  ab <- summarise.extrapolation.abundance(cells, dens, area)
  expect_equal(ab$total, sum(g$grid$NHat))
})
