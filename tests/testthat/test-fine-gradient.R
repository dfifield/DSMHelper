# fine.gradient(): clamp.lower clamps the returned covariate as well as the
# gradient's input (0.50.0). A synthetic "depth" field in a projected CRS, with
# a block of "land" (negative depth = height above sea level) on one side.

test_that("clamp.lower clamps the covariate and leaves the gradient unchanged (0.50.0)", {
  skip_if_not_installed("grec")
  # fine.gradient() reads three project globals from the calling environment
  # chain; set them in the global environment for the test and remove them after.
  old <- mget(c("GRADIENT_NA_FILL_PASSES", "SOBEL_KERNEL_GAIN", "GRADIENT_AGG_FUN"),
              envir = globalenv(), ifnotfound = list(NULL))
  assign("GRADIENT_NA_FILL_PASSES", 0, envir = globalenv())
  assign("SOBEL_KERNEL_GAIN", 8, envir = globalenv())
  assign("GRADIENT_AGG_FUN", "median", envir = globalenv())
  on.exit(for (v in names(old)) {
    if (is.null(old[[v]])) rm(list = v, envir = globalenv()) else assign(v, old[[v]], envir = globalenv())
  })

  # Source: 1 km cells over 40 x 40 km, depth rising from -60 (land, west) to
  # +140 m (east). Coarse grid: 10 km cells over the same extent, so N = 10.
  crs <- "EPSG:3347"
  src <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40000, ymin = 0, ymax = 40000, crs = crs)
  terra::values(src) <- rep(seq(-60, 140, length.out = 40), times = 40)   # row-major: x varies fastest
  names(src) <- "depth"
  coarse <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 40000, ymin = 0, ymax = 40000, crs = crs)

  un <- suppressMessages(fine.gradient(src, coarse, coarse.km = 10, mid.lat = 47))
  cl <- suppressMessages(fine.gradient(src, coarse, coarse.km = 10, mid.lat = 47, clamp.lower = 0))

  # Without the clamp the westmost coarse cells are negative; with it none is.
  expect_lt(min(terra::values(un$covar), na.rm = TRUE), 0)
  expect_gte(min(terra::values(cl$covar), na.rm = TRUE), 0)
  # Cells that were already >= 0 are unchanged (clamping commutes with the median).
  pos <- terra::values(un$covar) >= 0
  expect_equal(terra::values(cl$covar)[pos], terra::values(un$covar)[pos])
  # The clamped run's gradient is what the gradient always was with clamp.lower:
  # computed on the clamped field. Recompute that by hand-clamping the source.
  src.cl <- terra::clamp(src, lower = 0, values = TRUE)
  ref <- suppressMessages(fine.gradient(src.cl, coarse, coarse.km = 10, mid.lat = 47))
  expect_equal(terra::values(cl$grad), terra::values(ref$grad))
})
