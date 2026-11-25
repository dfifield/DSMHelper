# Harbour effect for gulls (Atlantic-DSM #173, DSMHelper 0.51.0).
#
# Gulls, Herring Gulls and Great Black-backed Gulls are many times denser near
# ports than offshore, by an amount that differs by season (Atl IMRP: Herring
# Gull 4.6x in Spring to 69x in Summer, within 10 km). Without a term for it the
# model explains those counts with depth, SST or the spatial smooth, and carries
# them to shallow or coastal cells far from any port (the mgcv Spring Gulls
# surface, about 3x too high on Atl IMRP). These functions:
#
#   find.survey.ports()     the ports, found from where survey ships start,
#                           stop or dock (Dave chose survey-derived ports only,
#                           2026-11-24: a list of named places flagged water
#                           near communities that ships only pass)
#   distance.to.ports.km()  each point's distance to the nearest port
#   harbour.share()         each prediction cell's share of water within the
#                           harbour radius
#   add.harbour.term(), add.sdmtmb.harbour.term()
#                           + Season:harbour, as add.fishing.term() adds fishing
#   model.uses.harbour()    whether a fitted model has the term
#   harbour.mixture.rows()  the prediction rows for the exact mixture
#
# How a cell is predicted (Dave's option 1, 2026-11-24): a cell whose share s
# of its water lies within the radius holds s of its area at the harbour
# density and 1 - s at the ordinary one, so its abundance is
#   (1 - s) x N(harbour = 0) + s x N(harbour = 1).
# harbour.mixture.rows() writes that as two prediction rows with the cell's
# area split between them, grouped as one cell. dsm.pred() sums each group's
# predictions; get.per.cell.var(group = ) and get.per.cell.var.posterior() sum
# the group inside the variance, keeping the covariance of the two rows (they
# come from the same coefficients); get.sdmtmb.surface() sums them per cell
# within each draw. Prediction and variance therefore see the same rows.
#
# Unlike fishing (#107), the harbour term is predicted at its true value, not 0:
# the gulls at a port are there whether or not a survey ship is.

#' Ports, found from where survey ships start, stop or dock
#'
#' A "port event" is a ship watch where a cruise's effort starts or stops: each
#' cruise's first and last watch, the watches either side of a break of
#' `gap.hours` or more within a cruise, and any watch coded PlatformActivity 8
#' ("At dock"). Events more than `land.km` from land are dropped (nights at
#' sea, rigs and stations are not ports; land is ETOPO altitude > 0). The rest
#' are grouped by single linkage at `cluster.km`, and a group whose events come
#' from at least `min.cruises` cruises is a port. Its position is the median of
#' its events' coordinates. On Atl IMRP this gives 22 ports, and the named
#' ports in ECSAS's lkpPorts that survey ships actually use lie 0.3-7.4 km from
#' them (Atlantic-DSM #173).
#'
#' @param watches Ship watches, one row per watch, before any speed or
#'   activity filter (the raw pull): `CruiseID`, `Date`, `StartTime`,
#'   `LatStart`, `LongStart`, `PlatformActivity`. `StartTime` may be a
#'   date-time whose date part is a placeholder (Access stores times as
#'   1899-12-30 hh:mm:ss); only its time of day is used.
#' @param etopo.file The ETOPO netCDF (`altitude` variable, lon/lat), or a
#'   `SpatRaster` of altitude in lon/lat (used by the tests).
#' @param gap.hours A break this long within a cruise counts as a port call.
#' @param land.km Events further than this from land are dropped.
#' @param cluster.km Single-linkage distance that joins events into one port.
#' @param min.cruises A group needs events from this many cruises.
#' @param crs.m A projected CRS in metres for the clustering (default Statistics
#'   Canada Lambert, EPSG:3347, which covers every SubProject).
#' @return Tibble, one row per port: `port` (an integer id, most-used first),
#'   `lat`, `lon`, `events`, `cruises`, `at_dock` (events coded 8).
#' @examples
#' \dontrun{
#' ports <- find.survey.ports(ship.raw, file.path(predLayerDir, "NetCDF", "etopo180.nc"))
#' }
#' @export
find.survey.ports <- function(watches, etopo.file, gap.hours = 12, land.km = 3,
                              cluster.km = 5, min.cruises = 3, crs.m = 3347) {
  checkmate::expect_data_frame(watches, min.rows = 1)
  checkmate::expect_names(names(watches), must.include = c("CruiseID", "Date", "StartTime", "LatStart",
                                                           "LongStart", "PlatformActivity"))
  if (!inherits(etopo.file, "SpatRaster")) checkmate::expect_file_exists(etopo.file)
  checkmate::expect_number(gap.hours, lower = 0)
  checkmate::expect_number(land.km, lower = 0)
  checkmate::expect_number(cluster.km, lower = 0)
  checkmate::expect_count(min.cruises, positive = TRUE)

  # ---- Each watch's start as one date-time ----
  # Date carries the day and StartTime the time of day (its own date part is
  # Access's 1899-12-30 placeholder), so the two are pasted together. UTC is
  # used only so that no daylight-saving shift can create or hide a 12 h gap.
  w <- watches %>%
    dplyr::distinct(.data$CruiseID, .data$Date, .data$StartTime, .data$LatStart, .data$LongStart,
                    .keep_all = TRUE) %>%
    dplyr::mutate(t = as.POSIXct(paste(as.Date(.data$Date), format(as.POSIXct(.data$StartTime), "%H:%M:%S")),
                                 tz = "UTC")) %>%
    dplyr::filter(!is.na(.data$LatStart), !is.na(.data$LongStart), !is.na(.data$t)) %>%
    dplyr::arrange(.data$CruiseID, .data$t)

  # ---- Port events ----
  # Within each cruise (rows in time order): the first and last watch, any
  # watch with a gap of gap.hours or more before or after it, and any watch
  # coded 8. coalesce(, 0) treats the missing gap before the first watch (and
  # after the last) as no gap; those two are kept by row_number() anyway.
  ev <- w %>%
    dplyr::group_by(.data$CruiseID) %>%
    dplyr::mutate(gap_before = as.numeric(difftime(.data$t, dplyr::lag(.data$t), units = "hours")),
                  gap_after = as.numeric(difftime(dplyr::lead(.data$t), .data$t, units = "hours"))) %>%
    dplyr::filter(dplyr::row_number() == 1 | dplyr::row_number() == dplyr::n() |
                    dplyr::coalesce(.data$gap_before, 0) >= gap.hours |
                    dplyr::coalesce(.data$gap_after, 0) >= gap.hours |
                    .data$PlatformActivity %in% 8) %>%
    dplyr::ungroup()
  if (!nrow(ev)) stop("find.survey.ports: no port events")

  # ---- Distance from each event to land ----
  # ETOPO cropped to the events' box plus 1 degree, so the nearest land of an
  # event near the box edge is still inside it. classify() turns water (0) into
  # NA and keeps land (1); terra::distance() then gives every cell its distance
  # in metres to the nearest non-NA (land) cell, computed on the ellipsoid for
  # a lon/lat raster. extract() reads it at each event.
  et <- if (inherits(etopo.file, "SpatRaster")) etopo.file else terra::rast(etopo.file, subds = "altitude")
  bb <- terra::ext(min(ev$LongStart) - 1, max(ev$LongStart) + 1, min(ev$LatStart) - 1, max(ev$LatStart) + 1)
  et <- terra::crop(et, bb)
  land <- terra::classify(et > 0, cbind(0, NA))
  d_land <- terra::distance(land)
  ev$d_land_km <- terra::extract(d_land, cbind(ev$LongStart, ev$LatStart))[, 1] / 1000
  ev <- dplyr::filter(ev, !is.na(.data$d_land_km), .data$d_land_km <= land.km)
  if (nrow(ev) < 2) stop("find.survey.ports: fewer than 2 port events within ", land.km, " km of land")

  # ---- Single-linkage clusters ----
  # Projected to metres so dist() is a distance in metres; hclust(single) then
  # cutree(h = ) joins any two events closer than cluster.km, directly or
  # through a chain of events.
  xy <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(ev, coords = c("LongStart", "LatStart"), crs = 4326), crs.m))
  ev$cl <- stats::cutree(stats::hclust(stats::dist(xy), method = "single"), h = cluster.km * 1000)

  # ---- Ports: clusters used by at least min.cruises cruises ----
  ev %>%
    dplyr::group_by(.data$cl) %>%
    dplyr::summarise(events = dplyr::n(), cruises = dplyr::n_distinct(.data$CruiseID),
                     at_dock = sum(.data$PlatformActivity %in% 8),
                     lat = stats::median(.data$LatStart), lon = stats::median(.data$LongStart),
                     .groups = "drop") %>%
    dplyr::filter(.data$cruises >= min.cruises) %>%
    dplyr::arrange(dplyr::desc(.data$cruises), dplyr::desc(.data$events)) %>%
    dplyr::mutate(port = dplyr::row_number()) %>%
    dplyr::select("port", "lat", "lon", "events", "cruises", "at_dock")
}

#' Distance from each point to the nearest port, km
#'
#' @param lon,lat Point coordinates, decimal degrees (WGS84).
#' Measured on the ellipsoid (geodesic), not in a projection: a conformal one
#' such as EPSG:3347 stretches distances by about 1.3% at 45°N (a 55.6 km
#' distance reads 56.3 km), which would move a 10 km radius by about 130 m.
#'
#' @param lon,lat Point coordinates, decimal degrees (WGS84).
#' @param ports From [find.survey.ports()]: `lat`, `lon`.
#' @return Numeric vector of km, one per point.
#' @examples
#' \dontrun{
#' segdata$d_port_km <- distance.to.ports.km(segdata$LongStart, segdata$LatStart, ports)
#' }
#' @export
distance.to.ports.km <- function(lon, lat, ports) {
  checkmate::expect_numeric(lon, any.missing = FALSE)
  checkmate::expect_numeric(lat, any.missing = FALSE, len = length(lon))
  checkmate::expect_data_frame(ports, min.rows = 1)
  checkmate::expect_names(names(ports), must.include = c("lat", "lon"))
  # Both sets of points in WGS84 lon/lat. For geographic coordinates sf
  # measures on the sphere or ellipsoid (s2, or lwgeom's ellipsoidal distance
  # when s2 is off), so st_nearest_feature() finds each point's nearest port
  # and st_distance(by_element = TRUE) the distance to it, in metres.
  p <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  q <- sf::st_as_sf(as.data.frame(ports), coords = c("lon", "lat"), crs = 4326)
  nn <- sf::st_nearest_feature(p, q)
  as.numeric(sf::st_distance(p, q[nn, ], by_element = TRUE)) / 1000
}

#' Each prediction cell's share of water within the harbour radius
#'
#' The share of each cell polygon's area that lies within `radius.km` of any
#' port. The prediction grid's polygons are the analysis-grid cells clipped to
#' the prediction study area (water only, coastal cells smaller), so this is a
#' share of the cell's water.
#'
#' @param cells `sf` polygons (the prediction grid, one row per cell), in a
#'   projected CRS in metres (segProj).
#' @param ports From [find.survey.ports()].
#' @param radius.km The harbour radius.
#' @return Numeric vector in \[0, 1\], one per cell.
#' @examples
#' \dontrun{
#' predgrid$harbour.share <- harbour.share(predgrid, ports, HARBOUR_RADIUS_KM)
#' }
#' @export
harbour.share <- function(cells, ports, radius.km) {
  checkmate::expect_class(cells, "sf")
  checkmate::expect_data_frame(ports, min.rows = 1)
  checkmate::expect_number(radius.km, lower = 0)
  if (sf::st_is_longlat(cells)) stop("harbour.share: cells must be in a projected CRS in metres")
  # The harbour zone: a circle of radius.km around every port, in the cells'
  # CRS, merged into one (multi)polygon so overlapping circles count once.
  zone <- sf::st_as_sf(as.data.frame(ports), coords = c("lon", "lat"), crs = 4326) %>%
    sf::st_transform(sf::st_crs(cells)) %>%
    sf::st_buffer(radius.km * 1000) %>%
    sf::st_union()
  geom <- sf::st_geometry(cells)
  share <- numeric(length(geom))
  # Only cells touching the zone need the intersection, which is the slow
  # part; every other cell has share 0. st_intersects() returns, for each
  # cell, the indices of the zone parts it touches (here 0 or 1).
  hit <- lengths(sf::st_intersects(geom, zone)) > 0
  # For each touching cell, the area of its intersection with the zone over
  # its own area. Cell by cell, because st_intersection() of a whole sfc
  # drops elements whose intersection is empty and can return a line for a
  # cell that only touches the zone's edge; one cell at a time, an empty or
  # zero-area intersection simply gives 0. Only a few hundred cells touch the
  # zone (22 ports x a few cells on Atl IMRP), so the loop is cheap.
  for (i in which(hit)) {
    inter <- sf::st_intersection(geom[i], zone)
    a_in <- if (length(inter)) sum(as.numeric(sf::st_area(inter))) else 0
    share[i] <- min(1, a_in / as.numeric(sf::st_area(geom[i])))
  }
  share
}

#' Add Season:harbour to every candidate model formula
#'
#' One log-multiplier per season for segments near a port, as
#' [add.fishing.term()] adds Season:fishing. Per season, because the harbour
#' effect differs by season (breeding against non-breeding for Herring and
#' Great Black-backed Gulls, and a different species mix in the Gulls group;
#' Atlantic-DSM #173).
#'
#' @param mod.specs Data frame of candidate models with a list column
#'   `formula` (dsm.mod.specs).
#' @return `mod.specs` with `+ Season:harbour` added to every formula.
#' @examples
#' \dontrun{
#' if (species %in% harbourAffectedSpecies) dsm.mod.specs <- add.harbour.term(dsm.mod.specs)
#' }
#' @export
add.harbour.term <- function(mod.specs) {
  checkmate::expect_data_frame(mod.specs)
  checkmate::expect_list(mod.specs$formula, types = "formula")
  mod.specs$formula <- lapply(mod.specs$formula, function(f) stats::update(f, . ~ . + Season:harbour))
  mod.specs
}

#' Add Season:harbour to one sdmTMB formula
#'
#' @param formula An sdmTMB model formula.
#' @return The formula with `+ Season:harbour`.
#' @examples
#' \dontrun{
#' f <- add.sdmtmb.harbour.term(f)
#' }
#' @export
add.sdmtmb.harbour.term <- function(formula) {
  checkmate::expect_formula(formula)
  stats::update(formula, . ~ . + Season:harbour)
}

#' Does a fitted model have the harbour term?
#'
#' @param model A fitted mgcv/dsm model or an sdmTMB fit.
#' @return TRUE when `harbour` is one of the formula's variables.
#' @examples
#' \dontrun{
#' model.uses.harbour(dsm_final)
#' }
#' @export
model.uses.harbour <- function(model) {
  # An sdmTMB fit keeps its formula as a list (one per part of a delta model);
  # mgcv's formula() returns one formula. all.vars() lists every variable name
  # on either side.
  f <- if (inherits(model, "sdmTMB")) model$formula[[1]] else stats::formula(model)
  "harbour" %in% all.vars(f)
}

#' The prediction rows for the exact harbour mixture
#'
#' Expands a prediction grid so that a cell with harbour share s is predicted
#' as (1 - s) of its area at harbour = 0 plus s of its area at harbour = 1:
#'   s = 0      one row, harbour = 0, the full offset;
#'   s = 1      one row, harbour = 1, the full offset;
#'   0 < s < 1  two rows: the original with harbour = 0 and offset x (1 - s),
#'              and a copy with harbour = 1 and offset x s.
#' The copies are appended after the original rows, so rows 1..nrow(df) keep
#' their order, and `source.row` says which original row each output row comes
#' from. Summing a model's predictions over `source.row` gives each cell's
#' mixture abundance; giving get.per.cell.var() the expanded `group` keeps the
#' two rows' covariance in the cell's variance.
#'
#' @param df Prediction rows with a `harbour.share` column (and, for the
#'   variance, every model predictor). An `sf` object loses its geometry.
#' @param off.set The offset (area, km2) of each row: a vector as long as
#'   `nrow(df)`.
#' @param group Optional grouping of the rows (get.per.cell.var's `group`, e.g.
#'   the platform copies of a cell); the copies inherit their row's group.
#' @return List: `df` (the expanded rows, with `harbour` set), `off.set`
#'   (scaled), `group` (NULL if none was given) and `source.row`.
#' @examples
#' \dontrun{
#' m <- harbour.mixture.rows(predgrid, predgrid$area)
#' n <- rowsum(predict(model, m$df, off.set = m$off.set), m$source.row)
#' }
#' @export
harbour.mixture.rows <- function(df, off.set, group = NULL) {
  checkmate::expect_data_frame(df, min.rows = 1)
  checkmate::expect_names(names(df), must.include = "harbour.share")
  checkmate::expect_numeric(off.set, len = nrow(df), any.missing = FALSE)
  if (!is.null(group) && length(group) != nrow(df))
    stop("harbour.mixture.rows: group has ", length(group), " elements, df has ", nrow(df), " rows")
  # An sf grid loses its geometry: the rows are only predicted from, and a
  # plain data frame binds and subsets without sf's sticky geometry column.
  if (inherits(df, "sf")) df <- sf::st_drop_geometry(df)
  s <- df$harbour.share
  if (anyNA(s) || any(s < 0 | s > 1)) stop("harbour.mixture.rows: harbour.share must be in [0, 1] with no NA")
  # The original rows: harbour = 1 only where the whole cell is harbour,
  # otherwise 0, and the offset scaled to the non-harbour part.
  base <- df
  base$harbour <- as.numeric(s >= 1)
  off_base <- ifelse(s >= 1, off.set, off.set * (1 - s))
  # The copies, for cells partly in the zone: harbour = 1 and the harbour part
  # of the area. which() gives their original row numbers.
  part <- which(s > 0 & s < 1)
  extra <- df[part, , drop = FALSE]
  # rep(1, length(part)), not a bare 1: when no cell is partly in the zone,
  # extra has 0 rows, and assigning a length-1 value to a 0-row data frame
  # is an error ("replacement has 1 row, data has 0").
  extra$harbour <- rep(1, length(part))
  # Both are plain data frames now, with the same columns.
  out <- rbind(base, extra)
  list(df = out, off.set = c(off_base, off.set[part] * s[part]),
       group = if (is.null(group)) NULL else {
         # factor() levels are kept, so the results still come back in the
         # caller's group order.
         g <- if (is.factor(group)) group else factor(group)
         factor(c(as.character(g), as.character(g[part])), levels = levels(g))
       },
       source.row = c(seq_len(nrow(df)), part))
}
