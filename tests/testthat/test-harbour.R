# Tests for the harbour functions (Atlantic-DSM #173, DSMHelper 0.51.0), all on
# synthetic data.

test_that("harbour.mixture.rows splits partly-harbour rows and keeps row order", {
  # Three cells: no harbour, 30% harbour, all harbour; areas 100 km2.
  df <- data.frame(id = 1:3, harbour.share = c(0, 0.3, 1), area = 100)
  m <- harbour.mixture.rows(df, df$area, group = c("a", "b", "c"))
  # Rows 1-3 are the originals in order; one copy is appended for cell 2.
  expect_equal(nrow(m$df), 4)
  expect_equal(m$source.row, c(1, 2, 3, 2))
  expect_equal(m$df$id, c(1, 2, 3, 2))
  # harbour is 0 for the no-harbour and partly-harbour originals, 1 for the
  # whole-harbour cell and for the copy.
  expect_equal(m$df$harbour, c(0, 0, 1, 1))
  # Offsets: the partly-harbour cell's area is split 70 / 30, and every
  # cell's offsets still add up to its area.
  expect_equal(m$off.set, c(100, 70, 100, 30))
  expect_equal(as.numeric(rowsum(m$off.set, m$source.row)), rep(100, 3))
  # The copy inherits its row's group, and the group levels are kept.
  expect_equal(as.character(m$group), c("a", "b", "c", "b"))
  expect_equal(levels(m$group), c("a", "b", "c"))
})

test_that("summing mixture-row predictions gives (1 - s) N0 + s N1", {
  # A log-linear model with a harbour multiplier of exp(2): N = area x
  # exp(b0 + b1 harbour). The mixture rows, predicted and summed per cell, must
  # equal (1 - s) N(harbour = 0) + s N(harbour = 1) with the full area.
  b0 <- -1; b1 <- 2
  pred <- function(d, off) off * exp(b0 + b1 * d$harbour)
  df <- data.frame(harbour.share = c(0, 0.25, 0.6, 1), area = c(100, 80, 50, 100))
  m <- harbour.mixture.rows(df, df$area)
  got <- rowsum(pred(m$df, m$off.set), m$source.row)[, 1]
  s <- df$harbour.share
  want <- (1 - s) * df$area * exp(b0) + s * df$area * exp(b0 + b1)
  expect_equal(unname(got), want)
})

test_that("harbour.mixture.rows refuses shares outside [0, 1] or missing", {
  expect_error(harbour.mixture.rows(data.frame(harbour.share = c(0.2, 1.5)), c(1, 1)), "harbour.share")
  expect_error(harbour.mixture.rows(data.frame(harbour.share = c(0.2, NA)), c(1, 1)), "harbour.share")
  expect_error(harbour.mixture.rows(data.frame(x = 1), 1))
})

test_that("harbour.share gives each cell's share of area within the radius", {
  # A 6 x 6 grid of 10 km square cells in a metric CRS (EPSG:3347), and one
  # port at the grid's centre (a cell corner). With a 10 km radius the four
  # cells around the centre are each a quarter-circle of radius 10 km inside a
  # 10 x 10 km square: share pi / 4. A far corner cell has share 0. The total
  # harbour area must be the whole circle, pi x 10^2 km2.
  crs <- sf::st_crs(3347)
  origin <- c(7e6, 1.5e6)
  sq <- lapply(0:35, function(i) {
    x0 <- origin[1] + (i %% 6) * 1e4; y0 <- origin[2] + (i %/% 6) * 1e4
    sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1e4, y0), c(x0 + 1e4, y0 + 1e4), c(x0, y0 + 1e4), c(x0, y0))))
  })
  cells <- sf::st_sf(id = 0:35, geometry = sf::st_sfc(sq, crs = crs))
  centre <- sf::st_transform(sf::st_sfc(sf::st_point(origin + 3e4), crs = crs), 4326)
  ports <- data.frame(lon = sf::st_coordinates(centre)[1], lat = sf::st_coordinates(centre)[2])
  s <- harbour.share(cells, ports, 10)
  # Cells 14, 15, 20, 21 touch the centre (row 2-3, column 2-3).
  expect_equal(s[c(15, 16, 21, 22)], rep(pi / 4, 4), tolerance = 0.01)
  expect_equal(s[1], 0)
  expect_equal(sum(s * 100), pi * 100, tolerance = 0.01)
  expect_true(all(s >= 0 & s <= 1))
  # A geographic CRS is refused: areas and the radius need metres.
  expect_error(harbour.share(sf::st_transform(cells, 4326), ports, 10), "projected")
})

test_that("distance.to.ports.km measures to the nearest port", {
  # Two ports 1 degree of latitude apart on the same meridian; a point at one
  # of them is 0 km away, and a point halfway is about 55.6 km from either.
  ports <- data.frame(lon = c(-60, -60), lat = c(45, 46))
  d <- distance.to.ports.km(c(-60, -60), c(45, 45.5), ports)
  expect_equal(d[1], 0, tolerance = 1e-6)
  expect_equal(d[2], 55.6, tolerance = 0.01)
})

test_that("find.survey.ports keeps start/stop clusters near land used by enough cruises", {
  # Synthetic land: altitude 10 m west of -60, -50 m east of it (water), on a
  # 1-arc-minute lon/lat raster as ETOPO's.
  r <- terra::rast(xmin = -62, xmax = -58, ymin = 44, ymax = 46, resolution = 1 / 60, crs = "EPSG:4326")
  terra::values(r) <- ifelse(terra::xFromCell(r, seq_len(terra::ncell(r))) < -60, 10, -50)
  # Four cruises. Each starts and ends at a harbour just east of the coast
  # (-59.98, 45, about 1.6 km from land) and steams offshore between; cruise 4
  # also stops overnight (a 14 h gap) far offshore at -58.5, which is not a
  # port because it is over 3 km from land.
  mk <- function(cid, lons, hours) data.frame(
    CruiseID = cid, Date = as.Date("2020-06-01") + hours %/% 24,
    StartTime = as.POSIXct("1899-12-30 00:00", tz = "UTC") + (hours %% 24) * 3600,
    LatStart = 45, LongStart = lons, PlatformActivity = 1)
  w <- rbind(mk(1, c(-59.98, -59.5, -59.0, -59.98), c(0, 2, 4, 6)),
             mk(2, c(-59.98, -59.4, -59.97), c(0, 3, 6)),
             mk(3, c(-59.97, -59.3, -59.98), c(0, 3, 6)),
             mk(4, c(-59.98, -58.5, -58.5, -59.98), c(0, 4, 18, 22)))
  p <- find.survey.ports(w, r)
  # One port, at the harbour, from all four cruises; the offshore stop is gone.
  expect_equal(nrow(p), 1)
  expect_equal(p$cruises, 4)
  expect_equal(p$lon, -59.98, tolerance = 1e-3)
  # With min.cruises = 5 nothing qualifies.
  expect_equal(nrow(find.survey.ports(w, r, min.cruises = 5)), 0)
})

test_that("add.harbour.term and model.uses.harbour", {
  specs <- tibble::tibble(modname = "m", formula = list(y ~ Season + s(x)))
  out <- add.harbour.term(specs)
  expect_true("Season:harbour" %in% attr(stats::terms(out$formula[[1]]), "term.labels"))
  expect_true("Season:harbour" %in% attr(stats::terms(add.sdmtmb.harbour.term(y ~ Season)), "term.labels"))
  # A fitted model with and without the term.
  set.seed(1)
  d <- data.frame(y = stats::rpois(60, 3), Season = factor(rep(c("a", "b"), 30)), harbour = rep(0:1, each = 30))
  expect_true(model.uses.harbour(stats::glm(y ~ Season + Season:harbour, poisson, d)))
  expect_false(model.uses.harbour(stats::glm(y ~ Season, poisson, d)))
})

test_that("harbour.mixture.rows works when no cell is partly harbour", {
  # All shares 0 or 1: no copies, the same rows, harbour set from the share.
  df <- data.frame(harbour.share = c(0, 1, 0))
  m <- harbour.mixture.rows(df, c(10, 20, 30))
  expect_equal(nrow(m$df), 3)
  expect_equal(m$df$harbour, c(0, 1, 0))
  expect_equal(m$off.set, c(10, 20, 30))
  expect_equal(m$source.row, 1:3)
})
