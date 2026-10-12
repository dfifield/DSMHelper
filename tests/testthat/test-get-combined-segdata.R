# get.combined.segdata() reads ddftype_levels from the global environment,
# where analysis_settings.R puts it.
local_ddftype_levels <- function(env = parent.frame()) {
  assign("ddftype_levels",
         c("SWD", "SFD", "AWD", "AFD", "SWN", "SFN", "AWN", "AFN"),
         envir = globalenv())
  withr::defer(rm("ddftype_levels", envir = globalenv()), envir = env)
}

make_segdata <- function(base, survey, types, est_dens) {
  n <- length(types)
  sf::st_sf(
    Sample.Label = paste(base, types, sep = "_"),
    estDens = est_dens,
    estAbund = est_dens * 2,
    rawCount = seq_len(n),
    SurveyType = survey,
    x = 1, y = 1, LatStart = 47, LongStart = -52,
    Season = "Fall",
    Date = as.Date("2020-10-01"),
    platform = types,
    geometry = sf::st_sfc(rep(list(sf::st_point(c(-52, 47))), n), crs = 4326)
  )
}

test_that("segments with different copy counts are each collapsed to one row", {
  local_ddftype_levels()

  # create.dsm.data() drops copies for ddftypes with no observations, so a
  # mixed Ship/Aerial species can have 4 ship copies beside 1 aerial copy
  # (Atl IMRP ATPU: aerial keeps only AWD).
  segdata <- rbind(
    make_segdata("ECSAS_shp_1", "Ship", c("SWD", "SFD", "SWN", "SFN"), c(1, 2, 3, 4)),
    make_segdata("ECSAS_air_2", "Aerial", "AWD", 5),
    make_segdata("ECSAS_shp_3", "Ship", c("SWD", "SFD", "SWN", "SFN"), c(10, 20, 30, 40))
  )

  res <- get.combined.segdata(segdata)

  expect_equal(nrow(res), 3)
  expect_equal(res$Sample.Label, c("ECSAS_air_2", "ECSAS_shp_1", "ECSAS_shp_3"))
  expect_equal(res$estDens, c(5, 10, 100))
  expect_equal(res$estAbund, c(10, 20, 200))
  expect_equal(res$rawCount, c(1, 10, 10))
  expect_equal(res$SurveyType, c("Aerial", "Ship", "Ship"))
  expect_false("platform" %in% names(res))
  expect_s3_class(res, "sf")
  expect_equal(sum(res$estAbund), sum(segdata$estAbund))
})

test_that("uniform copy counts give the same result as the pre-0.24.1 stride", {
  local_ddftype_levels()

  # The structure before create.dsm.data() dropped empty copies: every segment
  # carries all four ddftypes of its survey type, the empty ones as zeros.
  segdata <- rbind(
    make_segdata("ECSAS_shp_1", "Ship", c("SWD", "SFD", "SWN", "SFN"), c(1, 2, 3, 4)),
    make_segdata("ECSAS_air_2", "Aerial", c("AWD", "AFD", "AWN", "AFN"), c(5, 0, 0, 0))
  )

  res <- get.combined.segdata(segdata)

  expect_equal(res$Sample.Label, c("ECSAS_air_2", "ECSAS_shp_1"))
  expect_equal(res$estDens, c(5, 10))
  expect_equal(res$rawCount, c(10, 10))
})
