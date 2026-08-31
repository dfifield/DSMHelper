
#' Augment a transect row with WhatCount values from its constituent watches
#'
#' For each row in a transects data frame, determines how many distinct
#' \code{WhatCount} values appear across the watches it aggregates and stores
#' both the values and their count.  Note: performance is poor when called via
#' \code{map()} on large data frames.
#'
#' @param row Single-row data frame for one transect, with a \code{Watches}
#'   column of comma-separated watch IDs.
#' @param watches Data frame of watch records with \code{WatchID} and
#'   \code{WhatCount} columns.
#' @return \code{row} with \code{WhatCounts} (list column of unique values)
#'   and \code{nWhatCount} (integer count) appended.
#' @export
get.what.counts <- function(row, watches){
  wtch.ids <- stringr::str_split_fixed(row$Watches, ",", n = Inf) %>%
    stringr::str_trim()
  watches %<>% dplyr::filter(WatchID %in% wtch.ids)
  stopifnot(length(wtch.ids) == nrow(watches))
  row$WhatCounts <- watches$WhatCount %>%
    unique %>%
    na.omit %>%
    as.vector %>%
    list
  row$nWhatCount <- length(row$WhatCounts[[1]])
  row
}

#' Convert a bearing in degrees to an ECSAS cardinal direction code
#'
#' Maps a numeric bearing to the integer codes used in \code{lkpDirections}
#' in the ECSAS database.  Vectorised via \code{\link[base]{Vectorize}}.
#'
#' @param deg Numeric bearing in degrees (0--360), or \code{NA}.
#' @return Integer direction code (1--9), or \code{NA}.
#' @export
convert.to.cardinal.code <- Vectorize(function(deg) {
  if (is.na(deg))
    NA
  else if ((deg > 337.5 && deg <= 360) || (deg >= 0 && deg <= 22.5))
    2
  else if ((deg > 22.5) && (deg <= 67.5))
    3
  else if ((deg > 67.5) && (deg <= 112.5))
    4
  else if ((deg > 112.5) && (deg <= 157.5))
    5
  else if ((deg > 157.5) && (deg <= 202.5))
    6
  else if ((deg > 202.5) && (deg <= 247.5))
    7
  else if ((deg > 247.5) && (deg <= 292.5))
    8
  else if ((deg > 292.5) && (deg <= 337.5))
    9
  else
    1
})

#' Find which species group a species alpha code belongs to
#'
#' Searches the project global \code{spec.grps} list and returns the name of
#' the group containing \code{species}.  Vectorised over \code{species}.
#'
#' @param species Character string (or vector) of species alpha codes.
#' @details
#' Raises an error if \code{species} belongs to multiple groups.
#'
#' @return Character string group name, or \code{NA} if the species is not
#'   found in any group.
#' @export
find.spec.grp <- Vectorize(function(species) {
  if (is.na(species))
    return(NA)

  res <- names(spec.grps[grepl(species, spec.grps)])

  # No match
  if (length(res) == 0)
    return(NA)

  # Too many matches
  if (length(res) > 1)
    stop(sprintf(
      "%s belongs to more than one species group: %s",
      species,
      paste(res, collapse = ", ")
    ))

  res
})

#' Plot annual survey effort within a study area
#'
#' Filters transects to those whose geometry falls within a given study area
#' polygon and produces two plots: a faceted map of transect lines by year and
#' a bar chart of total annual effort. Also prints a summary table. Transects
#' crossing the study area boundary are included but not clipped (see issue #5).
#'
#' @param transects An sf object of LINESTRING transects with at least columns
#'   \code{Date} (Date) and \code{Effort} (numeric, km). Any CRS is accepted;
#'   coordinates are transformed to WGS84 internally.
#' @param sa An sf polygon defining the study area. Any CRS is accepted.
#' @param sa_label Character string used in plot titles, captions, and the
#'   summary table caption. Defaults to \code{"study area"}.
#' @param buf Numeric. Buffer in decimal degrees added around the study area
#'   bounding box when setting map extents. Must be finite and >= 0. Defaults
#'   to \code{0.5}.
#' @param species Character string identifying the species or group, used in
#'   plot titles. Defaults to \code{""}.
#' @return A data frame of annual effort (columns: \code{Year},
#'   \code{total_effort_km}, \code{n_transects}), returned invisibly. The
#'   summary table and both plots are printed as side effects.
#' @examples
#' \dontrun{
#' plot_annual_effort(the.data$transects, mcp_95, sa_label = "95% MCP")
#' plot_annual_effort(the.data$transects, ccman, sa_label = "min. concave polygon")
#' }
#' @export
plot_annual_effort <- function(transects, sa, sa_label = "study area", buf = 0.5,
                               species = "") {
  checkmate::expect_class(transects, "sf")
  checkmate::expect_names(names(transects), must.include = c("Date", "Effort"))
  checkmate::expect_class(sa, "sf")
  checkmate::expect_string(sa_label, min.chars = 1)
  checkmate::expect_number(buf, lower = 0, finite = TRUE)
  checkmate::expect_string(species)

  sa_ll <- sf::st_transform(sa, 4326) %>% sf::st_make_valid()

  countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  transects_ll <- sf::st_transform(transects, 4326) %>%
    dplyr::mutate(Year = lubridate::year(Date))

  # concaveman polygons can fail s2 spherical validity; use planar GEOS instead
  old_s2 <- sf::sf_use_s2(FALSE)
  transects_in_sa <- sf::st_filter(transects_ll, sa_ll)
  sf::sf_use_s2(old_s2)

  effort_by_year <- transects_in_sa %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(total_effort_km = sum(Effort, na.rm = TRUE),
              n_transects = dplyr::n())

  print(knitr::kable(effort_by_year,
                     caption = paste("Annual effort within", sa_label)))

  bbox <- sf::st_bbox(sa_ll)
  p_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = countries, fill = "grey85", color = "grey60", linewidth = 0.3) +
    ggplot2::geom_sf(data = sa_ll, fill = NA, color = "blue", linewidth = 0.6) +
    ggplot2::geom_sf(data = transects_in_sa, color = "red", linewidth = 0.7,
            alpha = 1.0) +
    ggplot2::facet_wrap(~ Year, ncol = 4) +
    ggplot2::coord_sf(
      xlim = c(bbox["xmin"] - buf, bbox["xmax"] + buf),
      ylim = c(bbox["ymin"] - buf, bbox["ymax"] + buf)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
    ggplot2::labs(
      title = paste(species, "- Annual effort distribution within", sa_label),
      caption = paste("Lines = transects within species study area (red outline).",
                      "\nNote: transects crossing the boundary are not yet clipped",
                      "(see issue #5).")
    )
  print(p_map)

  p_bar <- ggplot2::ggplot(effort_by_year, ggplot2::aes(x = factor(Year), y = total_effort_km)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(x = "Year", y = "Total effort (km)",
         title = paste(species, "- Total annual survey effort within", sa_label)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  print(p_bar)

  invisible(effort_by_year)
}


#' Validate and reconcile distance and interval distance columns
#'
#' Ensures an observation data frame does not simultaneously carry a
#' \code{distance} column and non-NA \code{distbegin}/\code{distend} columns,
#' which would confuse \code{\link[Distance]{ds}}.  All-NA interval columns
#' are dropped; otherwise \code{distance} is dropped after confirming it
#' equals \code{(distbegin + distend) / 2}.
#'
#' @param data Data frame of observation data.
#' @return \code{data} with the redundant column(s) removed.
#' @export
check.distdata.cols <- function(data) {
  if (!is.null(data$distance) && !is.null(data$distbegin) && !is.null(data$distend)) {
    if (all(is.na(data$distbegin)) && all(is.na(data$distend))) {
      message("check.distdata.cols: removing distbegin and distend columns that are all NA from data")
      data <- dplyr::select(data, -distbegin, -distend)
    } else if (isTRUE(all.equal(data$distance, (data$distbegin + data$distend) / 2))) {
      message(
        "check.distdata.cols: removing distance column because data has distbegin and distend columns"
      )
      data <- dplyr::select(data, -distance)
    } else {
      stop(
        paste0(
          "check.distdata.cols: data has both distance and non-NA distbegin/distend ",
          "columns but data$distance != (data$distbegin + data$distend)/2"
        )
      )
    }
  }
  data
}



#' Flag watches with suspicious GPS positions
#'
#' Wraps \code{ECSAS.find.suspicious.posn()} and annotates results with
#' whether each flagged watch was previously recorded as fixed in the
#' database.  Results are written to a CSV in \code{rel.folder}.
#'
#' @param alldat Data frame of all watch data.
#' @param rel.folder Path relative to \code{here()} where the output CSV is
#'   written.
#' @param leave Optional vector of watch IDs to exclude from flagging; passed
#'   to \code{ECSAS.find.suspicious.posn()}.
#' @param fixed_db Optional character vector of watch IDs already recorded as
#'   fixed in the database.
#' @param quiet If \code{TRUE} (default), suppresses printing of problem IDs
#'   and diagnostic histograms.
#' @param filename Name of the output CSV file.
#' @param ... Additional arguments passed to
#'   \code{ECSAS.find.suspicious.posn()}.
#' @return Data frame of suspicious watches; a 0-row data frame if none were
#'   found (previously \code{NULL}, which made callers such as
#'   \code{table(probs$Program)} fail with "nothing to tabulate").
#' @export
check.problem.data = function(alldat,
                              rel.folder,
                              leave = NULL,
                              fixed_db = NULL,
                              quiet = TRUE,
                              filename = "suspicious_posns.csv",
                              ...) {

  # Look for watches with problematic watches
  probs <- ECSASconnect::ECSAS.find.suspicious.posn(alldat, leave = leave, ...) %>%
    dplyr::mutate(prev_fixed = WatchID %in% fixed_db)

  # NB: > 0, not > 1. With exactly one problem watch the old test took the
  # "none found" branch, silently discarding that watch and reporting none.
  if (nrow(probs) > 0) {
    message (sprintf("%d problem watches found. %d of these were flagged as previously fixed in db",
                     nrow(probs), sum(probs$prev_fixed)))

    if (!quiet) {
      message("Problem watchIDs:")
      print(probs$WatchID)

      message("Problem watchIDs that were flagged as previously fixed in db:")
      print(dplyr::filter(probs, prev_fixed == TRUE)$WatchID)

      hist(probs$dist_diff_km)
      hist(probs$pct_diff)
    }

    probs %>%
      dplyr::select(CruiseID, WatchID, ObserverName, PlatformName, Date, StartTime, EndTime,
                    LatStart, LongStart, LatEnd, LongEnd, WatchLenKm, PlatformSpeed, CalcDurMin,
                    prev_fixed, dist_dr_km, dist_geo_km, dist_diff_km, pct_diff) %T>%
      readr::write_csv(file = here::here(rel.folder, filename))
  } else {
    # Leave probs as the 0-row data frame rather than setting it to NULL, so
    # callers can use nrow()/table() on the result without special-casing.

    # None found - clean up old files
    message("No problem watches found.")

    # create empty file or truncate if it exists
    file.create(file.path(here::here(rel.folder, filename)))

    #remove shapefile
    # file.remove(list.files(ShapeDir, pattern = paste0(filename, "\\..*"), full.names = T))
  }

  probs
}



#' Measure how much survey effort a study-area clip discards, by buffer width
#'
#' \code{\link{create.survey.data}} clips effort to a polygon. Where that
#' polygon's boundary is a straight administrative line drawn through surveyed
#' water, the clip cuts continuous coverage dead and leaves the DSM's spatial
#' smooth unconstrained across the cut; on a log link that has produced
#' predictions of 1.7e10 birds/km^2 (issue #33). The remedy is to clip to a
#' \emph{buffered} study area and predict on the exact one, which raises the
#' question this function exists to answer: how wide should the buffer be?
#'
#' It answers with two measurements rather than a rule of thumb.
#'
#' \strong{Recovery.} For each candidate width, how many watches and how many
#' kilometres of effort sit in the ring - outside the study area, inside the
#' buffer - reported \emph{by season}, because the failure this addresses is
#' seasonal and an annual total hides it. The marginal column is the one to
#' read: the width to choose is where recovery stops growing, not where it
#' stops.
#'
#' \strong{Boundary classification.} A buffer is only needed where the boundary
#' cuts through water that was surveyed on both sides. Points are sampled along
#' the boundary every \code{boundary.spacing.km} and each is classified by
#' whether there is effort within the width on the inside, the outside, or
#' both. Both is a \emph{cut}, and cuts are what run away. Inside-only is a
#' boundary the survey simply stopped at; neither is unsurveyed boundary. A
#' study area whose edge is all coastline should come back with no cuts at all,
#' and can honestly take a buffer of zero.
#'
#' Give it the watch table \strong{as it stands immediately before
#' \code{create.survey.data()}} - filtered, but not yet clipped. Measuring the
#' unfiltered database pull would count effort the analysis would never have
#' used; measuring anything after the clip cannot see past the cut at all.
#'
#' The widths worth exploring are bounded by what the database query retrieved.
#' The project derives its extraction bounding box from the study area expanded
#' by \code{BUFFER_MEASURE_MAX_KM} for exactly this reason - asking here for a
#' width beyond that measures the edge of the query, not the edge of the data.
#'
#' @param watches Data frame of watches, filtered but \strong{not} clipped,
#'   carrying longitude, latitude, effort and date columns (see the
#'   \code{*.col} arguments).
#' @param study.area \code{sf} polygon of the exact study area, any CRS.
#' @param widths.km Ascending numeric vector of candidate buffer widths, km.
#' @param season.def Named list of season definitions in the project's
#'   \code{month * 100 + day} form, i.e. one element of the \code{seasons} list
#'   (e.g. \code{seasons[["Petrels"]]}). \code{NULL} skips the seasonal split.
#' @param proj.crs A CRS in metres to do the geometry in. \code{NULL} (default)
#'   uses \code{study.area}'s own CRS, which must then be projected.
#' @param boundary.spacing.km Spacing of the points sampled along the boundary
#'   for the cut classification.
#' @param min.width.km Floor on the recommended width. Defaults to 25 - the
#'   distance within which residual variograms on this project reach 95\% of
#'   their plateau for four of five species (see the note beside
#'   \code{family.cv.block.size.m}). A buffer shorter than the range over which
#'   residuals stay correlated cannot anchor the smooth.
#' @param saturation.frac Recovery is called saturated at the first width from
#'   which onward every step recovers less than this fraction of the most
#'   productive step's per-km rate. Judged on \code{rel_per_km}, and it requires
#'   the condition to hold for every wider step as well - a single quiet step in
#'   the middle of a rising trend is not saturation.
#' @param lon.col,lat.col,effort.col,date.col Column names in \code{watches}.
#' @param watch.id.col Column identifying a watch. \code{watches} is reduced to
#'   one row per distinct value before anything is counted, because the natural
#'   thing to hand this function - \code{ecsas.ship.dat} - is the
#'   \emph{observation}-level table, where a watch appears once per observation
#'   recorded on it. Summing \code{effort.col} over that would multiply each
#'   watch's effort by its bird count, which is both wrong and biased towards
#'   exactly the busy watches a buffer question cares about. \code{NULL} skips
#'   the reduction, for a table that is already one row per watch.
#' @param lonlat.crs CRS of \code{lon.col}/\code{lat.col}. Default 4326.
#' @return Invisibly, a list with \code{widths} (one row per width per season),
#'   \code{boundary} (one row per sampled boundary point per width),
#'   \code{cuts} (contiguous runs of cut boundary, longest first),
#'   \code{cut.fraction}, \code{boundary.length.km}, \code{saturating.km} and
#'   \code{recommended.km}. Prints a summary.
#' @examples
#' \dontrun{
#' # In 00.01_Extract_data.Rmd, immediately before create.survey.data():
#' measure.study.area.buffer(ecsas.ship.dat, study.area,
#'                           season.def = seasons[["Petrels"]])
#' }
#' @export
measure.study.area.buffer <- function(watches,
                                      study.area,
                                      widths.km = seq(0, 150, by = 10),
                                      season.def = NULL,
                                      proj.crs = NULL,
                                      boundary.spacing.km = 10,
                                      min.width.km = 25,
                                      saturation.frac = 0.05,
                                      lon.col = "LongStart",
                                      lat.col = "LatStart",
                                      effort.col = "WatchLenKm",
                                      date.col = "Date",
                                      watch.id.col = "WatchID",
                                      lonlat.crs = 4326) {

  coll <- checkmate::makeAssertCollection()
  checkmate::assert_data_frame(watches, add = coll)
  checkmate::assert_class(study.area, "sf", add = coll)
  checkmate::assert_numeric(widths.km, lower = 0, min.len = 1, any.missing = FALSE,
                            sorted = TRUE, unique = TRUE, add = coll)
  checkmate::assert_list(season.def, null.ok = TRUE, add = coll)
  checkmate::assert_number(boundary.spacing.km, lower = 0.1, add = coll)
  checkmate::assert_number(min.width.km, lower = 0, add = coll)
  checkmate::assert_number(saturation.frac, lower = 0, upper = 1, add = coll)
  checkmate::assert_string(lon.col, add = coll)
  checkmate::assert_string(lat.col, add = coll)
  checkmate::assert_string(effort.col, add = coll)
  checkmate::assert_string(date.col, add = coll)
  checkmate::assert_string(watch.id.col, null.ok = TRUE, add = coll)
  checkmate::assert_subset(c(lon.col, lat.col, effort.col), names(watches), add = coll)
  checkmate::reportAssertions(coll)

  ###--------------------------------------------------------------------------
  ### Reduce to one row per watch
  #
  # See ?watch.id.col: the table this is normally handed is observation-level,
  # so without this every watch's effort is counted once per bird seen on it.
  if (!is.null(watch.id.col) && watch.id.col %in% names(watches)) {
    n.before <- nrow(watches)
    watches <- watches[!duplicated(watches[[watch.id.col]]), , drop = FALSE]
    if (nrow(watches) < n.before)
      message(sprintf(
        "measure.study.area.buffer: %d rows reduced to %d distinct %s.",
        n.before, nrow(watches), watch.id.col))
  }

  ###--------------------------------------------------------------------------
  ### Geometry set-up
  #
  # Everything below is metres, so the working CRS must be projected:
  # st_buffer() and st_line_sample() do not behave on a geographic one.
  if (is.null(proj.crs)) proj.crs <- sf::st_crs(study.area)
  sa <- sf::st_transform(study.area, proj.crs)
  if (isTRUE(sf::st_is_longlat(sa)))
    stop("measure.study.area.buffer: proj.crs must be a projected CRS in metres; ",
         "study.area is geographic and no proj.crs was given.")
  # st_zm() because a study area read from a shapefile can carry Z or M
  # coordinates - NL_EXPL_DRL_RA's does - and every GEOS predicate below then
  # fails with "GEOS does not support XYM or XYZM geometries".
  sa <- sf::st_union(sf::st_zm(sf::st_geometry(sa)))

  keep <- !is.na(watches[[lon.col]]) & !is.na(watches[[lat.col]])
  if (sum(keep) < nrow(watches))
    message(sprintf(
      "measure.study.area.buffer: dropping %d of %d watches with no position.",
      sum(!keep), nrow(watches)))
  w <- watches[keep, , drop = FALSE]

  pts <- w %>%
    sf::st_as_sf(coords = c(lon.col, lat.col), crs = sf::st_crs(lonlat.crs),
                 remove = FALSE) %>%
    sf::st_transform(proj.crs)

  effort <- w[[effort.col]]
  effort[is.na(effort)] <- 0

  if (!is.null(season.def) && date.col %in% names(w)) {
    # season.levels explicitly, so this does not depend on a season.names
    # global being in scope - it is a diagnostic and gets run from odd places.
    season <- assign.season(sf::st_drop_geometry(pts), season.def,
                            datefield = date.col,
                            season.levels = names(season.def))$Season
    season <- as.character(season)
  } else {
    season <- rep("All", nrow(pts))
  }
  season[is.na(season)] <- "unassigned"

  ###--------------------------------------------------------------------------
  ### Distance of every watch to the study area
  #
  # st_distance() to the polygon is 0 for anything inside it, so a single
  # distance vector answers every width at once. That is why the widths are not
  # looped over st_buffer(): each buffer would cost a full geometry operation
  # to learn something this already knows.
  inside <- lengths(sf::st_intersects(pts, sa)) > 0
  dist.m <- rep(0, nrow(pts))
  if (any(!inside))
    dist.m[!inside] <- as.numeric(sf::st_distance(pts[!inside, ], sa))

  message(sprintf(
    "measure.study.area.buffer: %d watches, %d inside the study area, %d outside.",
    nrow(pts), sum(inside), sum(!inside)))

  ###--------------------------------------------------------------------------
  ### Recovery by width and season
  seasons.seen <- sort(unique(season))
  widths <- purrr::map_dfr(widths.km, function(km) {
    ring <- !inside & dist.m <= km * 1000
    purrr::map_dfr(c("All", seasons.seen), function(s) {
      sel <- ring & (s == "All" | season == s)
      tibble::tibble(width_km  = km,
                     season    = s,
                     n_watches = sum(sel),
                     effort_km = sum(effort[sel]))
    })
  }) %>%
    dplyr::arrange(season, width_km) %>%
    dplyr::group_by(season) %>%
    dplyr::mutate(step_km          = width_km - dplyr::lag(width_km),
                  marginal_watches = n_watches - dplyr::lag(n_watches, default = 0),
                  marginal_effort  = effort_km - dplyr::lag(effort_km, default = 0),
                  # Per km of extra buffer, NOT per step. A step's raw gain is
                  # proportional to its width, so comparing raw gains across an
                  # uneven width grid says only which steps were widest - it
                  # showed a false saturation at 30 km on NL_EXPL_DRL_RA purely
                  # because the steps either side of it were 5 km and 10 km.
                  effort_per_km    = marginal_effort / step_km,
                  # As a fraction of the strongest per-km recovery seen. This is
                  # what saturation is judged on: recovery has saturated when
                  # widening further buys little compared with what widening
                  # bought at its most productive.
                  rel_per_km       = effort_per_km / max(effort_per_km, na.rm = TRUE)) %>%
    dplyr::ungroup()

  ###--------------------------------------------------------------------------
  ### Boundary classification
  #
  # Cast to LINESTRING first: a multi-ring study area is one MULTILINESTRING,
  # and sampling it as a single feature would space the points by total length
  # rather than per ring. density = points per metre handles rings of different
  # lengths without any arithmetic here.
  bnd <- sf::st_cast(sf::st_boundary(sa), "LINESTRING")
  bnd.len.m <- sum(as.numeric(sf::st_length(bnd)))
  bpts <- bnd %>%
    sf::st_line_sample(density = 1 / (boundary.spacing.km * 1000)) %>%
    sf::st_cast("POINT") %>%
    sf::st_sf(geometry = .)
  bpts$bid <- seq_len(nrow(bpts))
  bll <- sf::st_coordinates(sf::st_transform(bpts, 4326))
  bpts$lon <- bll[, "X"]
  bpts$lat <- bll[, "Y"]

  pts.in  <- pts[inside, ]
  pts.out <- pts[!inside, ]

  boundary <- purrr::map_dfr(widths.km[widths.km > 0], function(km) {
    d <- km * 1000
    tibble::tibble(
      width_km  = km,
      bid       = bpts$bid,
      lon       = bpts$lon,
      lat       = bpts$lat,
      n_inside  = lengths(sf::st_is_within_distance(bpts, pts.in,  dist = d)),
      n_outside = lengths(sf::st_is_within_distance(bpts, pts.out, dist = d))) %>%
      dplyr::mutate(is_cut = n_inside > 0 & n_outside > 0)
  })

  # Contiguous runs, so one long cut is reported as one cut rather than as a
  # scatter of points. NB runs are contiguous in sampling order, which follows
  # each ring in turn - a run spanning two rings would be split, which is the
  # conservative direction.
  #
  # watches_outside is a sum over the run's sample points, so a watch within
  # reach of several of them is counted several times. It ranks runs against
  # each other and nothing else; the recovery table is where the real counts
  # are.
  runs <- boundary %>%
    dplyr::group_by(width_km) %>%
    dplyr::arrange(bid, .by_group = TRUE) %>%
    dplyr::mutate(run = cumsum(is_cut != dplyr::lag(is_cut, default = FALSE))) %>%
    dplyr::filter(is_cut) %>%
    dplyr::group_by(width_km, run) %>%
    dplyr::summarise(n_points        = dplyr::n(),
                     length_km       = dplyr::n() * boundary.spacing.km,
                     lon_min         = min(lon),
                     lon_max         = max(lon),
                     lat_min         = min(lat),
                     lat_max         = max(lat),
                     watches_outside = sum(n_outside),
                     .groups = "drop") %>%
    dplyr::arrange(width_km, dplyr::desc(length_km))

  cut.frac <- boundary %>%
    dplyr::group_by(width_km) %>%
    dplyr::summarise(pct_boundary_cut = 100 * mean(is_cut),
                     km_boundary_cut  = sum(is_cut) * boundary.spacing.km,
                     .groups = "drop")

  ###--------------------------------------------------------------------------
  ### Recommendation
  #
  # The smallest width past which recovery stops growing materially, floored at
  # min.width.km. This is a reading of the table, not a decision - the constant
  # is set by hand in analysis_settings.R with the numbers written beside it.
  # Saturated means every width from here out recovers little per km, not just
  # this one - a single quiet step in the middle of a rising trend is not
  # saturation, and on NL_EXPL_DRL_RA there is one.
  all.w <- dplyr::filter(widths, season == "All", width_km > 0)
  quiet <- !is.na(all.w$rel_per_km) & all.w$rel_per_km < saturation.frac
  stays.quiet <- rev(cumprod(rev(as.numeric(quiet)))) > 0
  saturating.km  <- if (any(stays.quiet)) min(all.w$width_km[stays.quiet]) else NA_real_
  recommended.km <- if (is.na(saturating.km)) NA_real_ else
    max(saturating.km, min.width.km)

  ###--------------------------------------------------------------------------
  ### Report
  cat(sprintf("\nStudy area boundary: %.0f km, sampled every %.0f km (%d points).\n",
              bnd.len.m / 1000, boundary.spacing.km, nrow(bpts)))
  cat("\nEffort recovered in the ring (outside the study area, inside the buffer):\n\n")
  print(as.data.frame(dplyr::filter(widths, season == "All")), row.names = FALSE)
  if (length(seasons.seen) > 1) {
    cat("\nBy season:\n\n")
    print(as.data.frame(dplyr::filter(widths, season != "All")), row.names = FALSE)
  }
  cat("\nBoundary classified as a cut (effort within the width on BOTH sides):\n\n")
  print(as.data.frame(cut.frac), row.names = FALSE)
  if (nrow(runs)) {
    cat("\nLongest contiguous cuts, per width:\n\n")
    print(as.data.frame(dplyr::slice_head(dplyr::group_by(runs, width_km), n = 3)),
          row.names = FALSE)
  } else {
    cat("\nNo stretch of boundary has effort on both sides at any width tried.\n",
        "This study area does not have the failure mode issue #33 is about,\n",
        "and a buffer of 0 is defensible for it.\n", sep = "")
  }
  cat(sprintf("\nRecovery saturates at %s km (marginal gain < %.0f%%); recommended %s km, floor %.0f km.\n",
              ifelse(is.na(saturating.km), "no width tried", format(saturating.km)),
              100 * saturation.frac,
              ifelse(is.na(recommended.km), "undetermined", format(recommended.km)),
              min.width.km))
  cat("Set STUDY_AREA_BUFFER_KM by hand from the tables above, not from that line.\n\n")

  invisible(list(widths             = widths,
                 boundary           = boundary,
                 cuts               = runs,
                 cut.fraction       = cut.frac,
                 boundary.length.km = bnd.len.m / 1000,
                 saturating.km      = saturating.km,
                 recommended.km     = recommended.km))
}


#' Build observation and watch tables from raw ECSAS data
#'
#' Takes raw ECSAS data, filters observations, clips both watches and
#' observations to \code{clip.area}, and optionally assembles watches into
#' transects.  Returns a list with elements \code{distdata}, \code{watches},
#' and (if \code{create_transects = TRUE}) \code{transects}.
#'
#' \strong{The area effort is clipped to need not be the area results are
#' reported over, and deliberately is not on some SubProjects.} Where a study
#' area boundary is a straight administrative line through surveyed water,
#' clipping at it cuts continuous coverage dead and leaves the spatial smooth
#' unconstrained across the cut - on a log link that has produced predictions
#' of 1.7e10 birds/km^2 (issue #33). Passing a buffered polygon here keeps the
#' effort just outside the boundary, which anchors the smooth, while the
#' prediction grid stays on the exact study area so totals still describe the
#' assessment area. See the study-area buffer section of the project
#' \code{CLAUDE.md}, and \code{\link{measure.study.area.buffer}} for choosing
#' the width.
#'
#' @param raw.dat Data frame of raw ECSAS records.
#' @param dataset Character string identifying the dataset (e.g.
#'   \code{"ECSAS"}, \code{"SOMEC"}).
#' @param file.prefix Filename prefix used to name output shapefiles
#'   (e.g. \code{"ECSAS.ship"}).
#' @param inproj EPSG code or CRS object for the input coordinates.
#' @param outproj EPSG code or CRS object for output data; defaults to the
#'   project global \code{segProj}.
#' @param saveshp If \code{TRUE} (default), write watches and observations
#'   shapefiles to \code{ShapeDir}.
#' @param create_transects If \code{TRUE}, aggregate watches on the same day,
#'   ship, observer, and direction into transects.
#' @param intransect.only If \code{TRUE} (default), retain only observations
#'   where \code{InTransect == TRUE}.
#' @param clip.area \code{sf} polygon that watches and observations are clipped
#'   to.  \code{NULL} (the default) falls back to \code{study.area} from the
#'   calling environment, which is the historical behaviour.
#' @return Named list with elements \code{distdata}, \code{watches}, and
#'   optionally \code{transects}.
#' @export
create.survey.data <- function(raw.dat = NULL,
                               dataset = NULL,
                               file.prefix = NULL,
                               inproj = 4326,
                               outproj = segProj,
                               saveshp = TRUE,
                               create_transects = FALSE,
                               intransect.only = TRUE,
                               clip.area = NULL) {

  # Fall back to the caller's study.area so existing callers are unaffected.
  # Resolved before validation so the assertion reports on what is actually
  # used, whichever it came from.
  if (is.null(clip.area)) {
    if (!exists("study.area", envir = parent.frame()))
      stop("create.survey.data: clip.area is NULL and there is no study.area ",
           "in the calling environment to fall back to.")
    clip.area <- get("study.area", envir = parent.frame())
  }

  coll = checkmate::makeAssertCollection()
  checkmate::assert_data_frame(raw.dat, add = coll)
  checkmate::assert(
    checkmate::check_class(clip.area, "sf"),
    add = coll
  )
  checkmate::assert(
    checkmate::check_string(file.prefix),
    add = coll
  )
  checkmate::reportAssertions(coll)


  # create watches
  message("Creating watches...")
  # NB: WhatCount must be kept. The aerial watches built in 00.01 include it, so
  # dropping it here leaves the ship/SOMEC watches one column short and the
  # rbind() in combine_all_data fails with "numbers of columns of arguments do
  # not match".
  keep_cols <- c("SurveyType", "TransectID", "Program", "CruiseID", "WatchID",
                 "ObserverName", "TransFarEdge", "DistMeth", "Date", "StartTime",
                 "EndTime", "LatStart", "LongStart", "LatEnd", "LongEnd",
                 "WatchLenKm", "ObsHeight", "CalcDurMin", "TotalWidthKm",
                 "WhatCount")

  # These three only needed for making transects
  if (create_transects)
    keep_cols <- c(keep_cols, c("PlatformDir", "PlatformName", "PlatformSpeed"))

  watches <- raw.dat %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    dplyr::mutate(
      Sample.Label = WatchID,
      # For ECSAS, TransectSides is 1 for ship, and assigned elsewhere in
      # Extract_data.Rmd for air. For SOMEC, both aerial and ship data pass
      # through here and all aerial surveys are 2 sided whereas ship are 1.
      TransectSides = dplyr::case_when(
        dataset == "SOMEC" & SurveyType == "Aerial" ~ 2,
        .default = 1)
    ) %>%
    dplyr::distinct()

  # Collapse watches recorded once per observer into a single row, keeping the
  # observer names.
  #
  # NB: do NOT simply drop ObserverName for SOMEC. That was an earlier stopgap
  # to make the watch data consistent, but it discards who observed and leaves
  # the column absent entirely, so the ship/SOMEC watches no longer share a
  # column set with the aerial ones and rbind() fails downstream.
  watches <- watches %>%
    dplyr::group_by(dplyr::across(-ObserverName)) %>%
    dplyr::summarise(
      ObserverName = paste(unique(ObserverName), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(CruiseID, Sample.Label, ObserverName, Date, StartTime)

  # Make sure we didn't lose any watches while collapsing observers
  missing_ids <- setdiff(raw.dat$WatchID, watches$WatchID)
  if (length(missing_ids) != 0)
    stop("create.survey.data: lost the following WatchIDs after collapsing ",
         "watches with multiple observers: ",
         paste(missing_ids, collapse = ", "))

  # make sure all data for a watch is consistent. Find rows with duplicate watchIDs
  # and remove these watches
  dups <- watches %>%
    dplyr::group_by(WatchID) %>%
    dplyr::summarise(nrows = dplyr::n()) %>%
    dplyr::filter(nrows > 1) %>%
    dplyr::pull("WatchID")

  # Remove watches with inconsistent watch info
  if (length(dups) > 0) {
    warning(paste0(sprintf("Removing %d watches due to inconsistent watch info between rows: ", length(dups)),
                   paste(dups, collapse = ", ")), immediate. = TRUE)
    watches <- dplyr::filter(watches, !(WatchID %in% dups))
  }

  # clip to clip.area, which is the study area buffered by
  # STUDY_AREA_BUFFER_KM where the SubProject sets one - see the roxygen above.
  #
  # NB: rmapshaper::ms_clip() throws "Not compatible with STRSXP: [type=list]"
  # rather than returning an empty result when nothing overlaps, so check for
  # overlap first. This happens for a whole dataset when the clip area has no
  # coverage by this survey type at all.
  watch_pts <- watches %>%
    sf::st_as_sf(coords = c("LongStart", "LatStart"), crs = sf::st_crs(inproj)) %>%
    sf::st_transform(sf::st_crs(4326)) %>% # for ms_clip below
    dplyr::select(WatchID) # just keep WatchID
  clip_area_4326 <- clip.area %>% sf::st_transform(sf::st_crs(4326))

  if (nrow(sf::st_filter(watch_pts, clip_area_4326)) == 0) {
    warning(sprintf(
      "create.survey.data: no %s watches overlap clip.area - returning no watches",
      dataset), immediate. = TRUE)
    watches <- watch_pts %>%
      dplyr::slice(0) %>%
      dplyr::left_join(watches, by = "WatchID") %>%
      sf::st_transform(outproj)
  } else {
    watches <- watch_pts %>%
      rmapshaper::ms_clip(clip_area_4326) %>%   # do the clipping -
      dplyr::left_join(watches, by = "WatchID") %>%  # add other cols back in
      sf::st_transform(outproj)
  }


  # Create transects: combine watches on on same day, same ship, same observer
  # and same direction into transects.
  if (create_transects) {
    message("Creating transects...")
    transects <- ECSASconnect::ECSAS.create.transects(sf::st_drop_geometry(watches)) %>%
      sf::st_as_sf()

    # Modify watches and set the Sample.Label for each watch
    watches %<>% ECSASconnect::ECSAS.add.sample.label(sf::st_drop_geometry(transects))

    # Cconvert Watches list column in transects object into vector of watchID's
    # contained in each transect since st_write (and other downstream code?)
    # can't deal with list cols
    transects %<>%
      dplyr::mutate(wtchs = unname(unlist(
        split(., .$Sample.Label) %>% purrr::map( ~ unlist(.x$Watches) %>%
                                                   paste(collapse = ", "))
      )),
      length = sf::st_length(.)) %>%
      dplyr::select(-Watches) %>%
      dplyr::rename(Watches = wtchs)


    if (saveshp) {
      st.write.if.any(transects,
                      layer = paste0(file.prefix, "_transects.shp"),
                      what = "transects")
    }
  }


  # Save watches as shapefile
  if (saveshp) {
    st.write.if.any(watches,
                    layer = paste0(file.prefix, "_watches.shp"),
                    what = "watches")
  }

  # no longer remember why this was desirable
  watches <- dplyr::mutate(watches,
                           StartTime = as.character(StartTime),
                           EndTime = as.character(EndTime))

  # Create Obs: remove non-birds, convert distances to km,  windforce is
  # converted to windspeed if necessary, convert InTransect to T/F, remove ship
  # followers, remove zeros (ie. watches where Count == NA),
  # only keep Flyswim == W or F (not L or S). Add in sample.Label,
  # rename columns and select ones of interest
  message("Creating observations...")

  obs <- raw.dat %>%
    dplyr::mutate(
      Region.Label = 1,
      InTransect = dplyr::case_when(InTransect == -1 ~ TRUE,
                                    InTransect == 0 ~ FALSE,
                                    TRUE ~ NA),
      Distance = Distance / 1000
    ) %>%
    dplyr::filter(
      if (isTRUE(intransect.only))
        InTransect == TRUE
      else
        InTransect %in% c(TRUE, FALSE),
      WatchID %in% watches$WatchID,
      !is.na(Count),
      !is.na(DistMeth),
      # !is.na(Distance), # Now using these with dummy_ddf so don't filter out
      FlySwim %in% c("W", "F"),
      Class == "Bird",
      is.na(Association) | Association != 18
    ) %>%
    # Assign DistType (Need to do after filtering FlySwim for W or F),
    # and add FlockID if there isn't one
    # NB: use seq_len(), not 1:nrow(). On 0 rows 1:nrow(.) is 1:0 == c(1, 0),
    # a length-2 vector, which case_when() rejects against a length-0
    # condition. all(is.na(x)) is TRUE for an empty vector, so this branch is
    # *guaranteed* to be taken when there are no observations.
    dplyr::mutate(DistType = assign.dist.type(.),
                  FlockID = dplyr::case_when(all(is.na(FlockID)) ~ seq_len(nrow(.)),
                                             TRUE ~ FlockID)) %>%
    dplyr::left_join(watches[, c("WatchID", "Sample.Label")], by = "WatchID") %>%
    dplyr::rename(object = FlockID,
                  size = Count,
                  distance = Distance) %>%
    dplyr::select(
      SurveyType,
      object,
      Sample.Label,
      Region.Label,
      ObserverName,
      DistMeth,
      DistType,
      Program,
      WatchID,
      Date,
      Year,
      Alpha,
      English,
      size,
      FlySwim,
      InTransect,
      distance,
      Visibility,
      Windspeed,
      Windforce,
      Glare,
      SeaState,
      TransFarEdge,
      Swell,
      LatStart,
      LongStart,
      ObsLat,
      ObsLong,
      ObsTime,
      ObsHeight,
      Association,
      Behaviour,
      Age,
      Plumage,
      Sex,
      DistanceCode
    ) %>%
    dplyr::mutate(
      Windspeed = dplyr::case_when(
        is.na(Windspeed) &
          !is.na(Windforce) ~ dplyr::left_join(., beaufort_conversion,
                                               by = c("Windforce" = "beaufort"))$speed.kts,
        TRUE ~ Windspeed
      ),
      weights = dplyr::case_when(
        DistanceCode %in% c("A", "B") ~ 2,
        DistanceCode %in% c("C", "D") ~ 1,
        TRUE ~ 0
      ),
      Visibility = dplyr::case_when(Visibility > 20 ~ 20,
                                    TRUE ~ Visibility),
      FlySwim = as.factor(FlySwim),
      distbegin = NA_integer_,
      distend = NA_integer_
    ) %>%
    droplevels

  # Make sure there's still some left
  if (nrow(obs) == 0)
    warning("No observations left after filtering!", immediate. = TRUE)

  # Clip to clip.area, exactly as the watches were above.
  #
  # NB: see the note on the watches clip above - ms_clip() errors rather than
  # returning an empty result when nothing overlaps, so check for overlap
  # first. Unlike the watches case this can fire even when the survey type does
  # cover the clip area, if it simply recorded no in-transect observations
  # inside it.
  obs_pts <- obs %>%
    sf::st_as_sf(coords = c("LongStart", "LatStart"), crs = sf::st_crs(inproj)) %>%
    sf::st_transform(sf::st_crs(4326)) %>% # for ms_clip below
    dplyr::select(object) # just keep object
  clip_area_4326 <- clip.area %>% sf::st_transform(sf::st_crs(4326))

  if (nrow(sf::st_filter(obs_pts, clip_area_4326)) == 0) {
    warning(sprintf(
      "create.survey.data: no %s observations overlap clip.area - returning no observations",
      dataset), immediate. = TRUE)
    obs <- obs_pts %>%
      dplyr::slice(0) %>%
      dplyr::left_join(obs, by = "object") %>%
      sf::st_transform(outproj)
  } else {
    obs <- obs_pts %>%
      rmapshaper::ms_clip(clip_area_4326) %>%   # do the clipping -
      dplyr::left_join(obs, by = "object") %>%  # add other cols back in
      sf::st_transform(outproj)
  }

  # Save as shapefile
  if (saveshp) {
    st.write.if.any(obs,
                    layer = paste0(file.prefix, "_obs.shp"),
                    what = "obs")
  }

  ##### After all that, now just create a single distdata for use in distance
  #sampling.This doesn't actually get used (unless we do analysis of all
  #seabirds). Rather, specific datasets for a given species/group are created by
  #create.dsm.data( in order to properly generate the samples where there were 0
  #observations of a given species.
  message("Creating distdata....")
  distdata <- obs %>%
    sf::st_drop_geometry() %>%
    droplevels

  # Set dataset attribute.
  #
  # NB: assigning a scalar into a column of a 0-row data frame errors with
  # "replacement has 1 row, data has 0", so size the value to the number of
  # rows. This keeps the column present (and correctly typed) on empty data so
  # downstream rbind()/bind_rows() still line up.
  distdata$Dataset <- rep(dataset, nrow(distdata))
  watches$Dataset <- rep(dataset, nrow(watches))

  the.data <-
    if (create_transects) {
      list(distdata = distdata,
           watches = sf::st_drop_geometry(watches),
           transects = transects)
    } else {
      list(
        distdata = distdata,
        watches = sf::st_drop_geometry(watches))
    }

  message("\nDone\n")
  the.data
}

#' Convert a watches data frame to a SpatialLinesDataFrame
#'
#' @param watches Data frame with columns \code{WatchID}, \code{LongStart},
#'   \code{LatStart}, \code{LongEnd}, \code{LatEnd}.
#' @return \code{SpatialLinesDataFrame} in WGS 84 (EPSG 4326).
#' @export
watches.to.lines <- function(watches) {

  lins <- watches %>%
    sp::split(.$WatchID) %>%
    lapply(function(x) {
      sp::Lines(list(sp::Line(
        matrix(c(x$LongStart, x$LatStart, x$LongEnd, x$LatEnd),
               nrow =  2, byrow = T))),
        x$WatchID[1L])
    }) %>%
    sp::SpatialLines()

  watches  <- as.data.frame(watches)
  rownames(watches) <- watches$WatchID
  l <- sp::SpatialLinesDataFrame(lins, watches)
  sp::proj4string(l) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  l
}


#' Custom GAM diagnostic plots and basis-dimension checks
#'
#' A wrapper around \code{\link[mgcv]{gam.check}} that plots Q-Q, residuals
#' vs linear predictor, a residual histogram, and response vs fitted values,
#' then runs \code{\link[mgcv]{k.check}}.
#'
#' @param b A fitted \code{gam} or \code{bam} object.
#' @param old.style If \code{TRUE}, use \code{qqnorm} instead of
#'   \code{\link[mgcv]{qq.gam}}.
#' @param type Residual type; one of \code{"deviance"}, \code{"pearson"}, or
#'   \code{"response"}.
#' @param k.sample Number of data points subsampled for the basis-dimension
#'   check.
#' @param k.rep Number of replicates for the basis-dimension check.
#' @param rep Number of simulations for \code{\link[mgcv]{qq.gam}}.
#' @param level Reference band coverage for \code{\link[mgcv]{qq.gam}}.
#' @param rl.col Colour for the reference line in \code{\link[mgcv]{qq.gam}}.
#' @param rep.col Colour for the simulation envelope.
#' @param ... Additional arguments passed to plot functions.
#' @return \code{invisible(NULL)}, called for its side-effects (plots and
#'   printed output).
#' @export
my.gam.check <- function(b, old.style = FALSE, type = c("deviance", "pearson",
                                                        "response"), k.sample = 5000, k.rep = 200, rep = 0, level = 0.9,
                         rl.col = 2, rep.col = "gray80", ...)
{
  type <- match.arg(type)
  resid <- residuals(b, type = type)
  linpred <- if (is.matrix(b$linear.predictors) && !is.matrix(resid))
    napredict(b$na.action, b$linear.predictors[, 1])
  else napredict(b$na.action, b$linear.predictors)

  if (old.style)  {
    message("qqnorm:")
    qqnorm(resid, ...)
  } else {
    message("qq.gam: ")
    mgcv::qq.gam(b, rep = rep, level = level, type = type, rl.col = rl.col,
                 rep.col = rep.col, ...)
  }
  plot(linpred, resid, main = "Resids vs. linear pred.", xlab = "linear predictor",
       ylab = "residuals", ...)
  hist(resid, xlab = "Residuals", main = "Histogram of residuals",
       ...)
  fv <- if (inherits(b$family, "extended.family"))
    predict(b, type = "response")
  else fitted(b)
  if (is.matrix(fv) && !is.matrix(b$y))
    fv <- fv[, 1]
  try(plot(fv, napredict(b$na.action, b$y), xlab = "Fitted Values",
           ylab = "Response", main = "Response vs. Fitted Values",
           ...))
  gamm <- !(b$method %in% c("GCV", "GACV", "UBRE", "REML",
                            "ML", "P-ML", "P-REML", "fREML"))
  if (gamm) {
    cat("\n'gamm' based fit - care required with interpretation.")
    cat("\nChecks based on working residuals may be misleading.")
  }
  cat("\n")
  kchck <- mgcv::k.check(b, subsample = k.sample, n.rep = k.rep)
  if (!is.null(kchck)) {
    cat("Basis dimension (k) checking results. Low p-value (k-index<1) may\n")
    cat("indicate that k is too low, especially if edf is close to k'.\n\n")
    printCoefmat(kchck, digits = 3)
  }
}

#' Observed vs expected density or abundance by covariate
#'
#' Aggregates observed and model-predicted counts (or densities) by levels of
#' \code{covar} and optionally plots them against each other.
#'
#' @param model Fitted \code{dsm} model object.
#' @param covar Character string naming the covariate column in
#'   \code{model$data} to aggregate by.
#' @param cut Optional numeric vector of breakpoints passed to \code{cut()}
#'   to bin a continuous covariate.
#' @param plotit If \code{TRUE}, produce a scatter plot of observed vs
#'   expected.
#' @param debug If \code{TRUE}, drop into \code{browser()} at the start.
#' @param ... Additional arguments passed to \code{plot()}.
#' @return 2-row matrix (Observed, Expected) with one column per level of
#'   \code{covar}.
#' @export
oe.dens <- function(model, covar, cut = NULL, plotit = FALSE, debug = FALSE, ...)
{
  if (debug) browser()
  #get data
  oe <- model$data

  # add in predictions as abundances
  resp <- as.character(model$formula[[2]])
  if (resp %in% c("D", "density", "Dhat", "density.est")) {
    oe$N <- predict(model) * oe$segment.area
  } else {
    oe$N <- predict(model)
  }
  if (!is.null(cut)) {
    oe[[covar]] <- cut(oe[[covar]], breaks = cut)
  }
  oe <- plyr::ddply(oe, covar, function(x) {
    data.frame(Observed = sum(x[resp]), Expected = sum(x$N), n = nrow(x))
  })
  cn <- oe[, 1]
  oe <- t(oe[, 2:3])
  colnames(oe) <- cn
  if (plotit) {
    maxlim <- max(oe[1,], oe[2,], na.rm = TRUE)
    plot(oe[1,], oe[2,], xlab = "obs", ylab = "exp",
         main = paste0("Obs vs exp at specific values of ", covar),
         xlim = c(0, maxlim), ylim = c(0, maxlim), ...)
  }
  return(oe)
}


#' Render a full preliminary DSM analysis report for one species
#'
#' Knits \code{Generic_full_test.Rmd} with the given species as a parameter
#' and writes the HTML output to \code{ResultsDir/<species>/}.
#'
#' Expects project globals \code{ResultsDir} and \code{RDir} in the calling
#' environment.
#'
#' @param species Character string species code (e.g. \code{"ATPU"}).
#' @return \code{invisible(NULL)}, called for its side-effect (rendered HTML).
#' @export
do.full.prelim <- function(species){

  if(!dir.exists(file.path(ResultsDir, species))) {
    dir.create(file.path(ResultsDir, species))
  }

  out.file <-
    file.path(ResultsDir,
              species,
              paste(species, "00.04_full_prelim.html", sep = "_"))
  message(sprintf("Rendering full prelim %s to %s",
                  species,
                  out.file))
  rmarkdown::render(
    file.path(RDir, "Generic_full_test.Rmd"),
    params = list(species = species),
    output_file = out.file
  )
}


#' Set default detection function spec values by platform class
#'
#' Populates \code{convert_units}, \code{cutpoints}, and
#' \code{distance.centers} in a DDF spec list based on whether the spec name
#' indicates a ship or aerial survey.
#'
#' @param ddf.def A single DDF spec list to be populated.
#' @param nm Character string DDF spec name of the form
#'   \code{"<dataset>_<Ship|Aerial>_<behav>_<D|N>"}.
#' @return \code{ddf.def} with platform-appropriate default values filled in.
#' @export
set.def.df.spec.values <- function(ddf.def, nm) {

  if(length(nm) != 1)
    stop("set.def.df.spec.values: length of ddf name is not 1.")


  # Note that the convert_units won't actually get used because we
  # only use ds() to compute the detection function and not
  # abundance.
  if(stringr::str_detect(nm, "Ship")) {
    # ECSAS and SOMEC Ship surveys have set cutpoints whereas aerial uses distbegin
    # and distend columns computed in 00.01_Extract_data.Rmd.
    ddf.def$convert_units = ecsas.ship.convert.units
    ddf.def$cutpoints = ecsas.ship.ctpoints
    ddf.def$distance.centers = ecsas.ship.distance.centers
  } else if (stringr::str_detect(nm, "Aerial")) {
    ddf.def$convert_units = ecsas.air.convert.units
    ddf.def$distance.centers = ecsas.air.distance.centers
  } else
    stop("set.def.df.spec.values: unrecognized ddf name.")
  ddf.def
}

#' Parse a DDF spec name into its four components
#'
#' Splits a name of the form \code{"<dataset>_<platform>_<behav>_<dist>"} and
#' returns the parts as a named list.
#'
#' @param nm Character string DDF spec name, e.g. \code{"ECSAS_Ship_F_D"}.
#' @return Named list with elements \code{dataset}, \code{platform_class},
#'   \code{behav}, and \code{dist_type}.
#' @export
parse.df.name <- function(nm){
  res <- stringr::str_split_1(nm, stringr::fixed("_"))
  list(dataset = res[1], platform_class = res[2], behav = res[3], dist_type = res[4])
}


#' Extract the dataset component from a DDF spec name
#'
#' @param nm Character string DDF spec name (e.g. \code{"ECSAS_Ship_F_D"}).
#' @return Character string dataset identifier (e.g. \code{"ECSAS"}).
#' @export
get.dataset <- function(nm) {
  parse.df.name(nm)$dataset
}

#' Extract the platform class from a DDF spec name
#'
#' @param nm Character string DDF spec name.
#' @return \code{"Ship"} or \code{"Aerial"}.
#' @export
get.platform.class <- function(nm) {
  parse.df.name(nm)$platform_class
}

#' Extract the behaviour component from a DDF spec name
#'
#' @param nm Character string DDF spec name.
#' @return \code{"F"} (flying) or \code{"W"} (water).
#' @export
get.behav <- function(nm) {
  parse.df.name(nm)$behav
}

#' Extract the distance-type component from a DDF spec name
#'
#' @param nm Character string DDF spec name.
#' @return \code{"D"} (distances available) or \code{"N"} (no distances).
#' @export
get.dist_type <- function(nm) {
  parse.df.name(nm)$dist_type
}

#' Drop DDF spec entries that no longer correspond to a live ddftype
#'
#' \code{01.02_Create_final_ddf_model_specs.Rmd} records the analyst's final
#' DDF choices by assigning into \code{df.mod.list[[species]]$<ddfname>}.  When
#' the ddftypes have been pruned to the survey types actually present (see
#' \code{\link{prune.ddf.globals}}), assigning to a name that is no longer in
#' the list silently \emph{appends} a new element rather than erroring, which
#' breaks the positional correspondence between \code{def.ddf.list} and
#' \code{ddftype_levels} that the rest of the pipeline relies on.
#'
#' This drops any such entries, reporting what was removed, so the hand-written
#' per-species blocks can be left untouched across SubProjects with different
#' survey coverage.
#'
#' @param df.mod.list Named list (by species) of named lists of DDF specs.
#' @param valid.names Character vector of the DDF spec names that should be
#'   retained; defaults to the project global \code{def.ddf.list}'s names.
#' @return \code{df.mod.list} with out-of-date entries removed.
#' @examples
#' \dontrun{
#' df.mod.list <- prune.df.mod.list(df.mod.list)
#' save(df.mod.list, file = dfModlistLoc)
#' }
#' @export
prune.df.mod.list <- function(df.mod.list, valid.names = names(def.ddf.list)) {

  checkmate::expect_list(df.mod.list, min.len = 1)
  checkmate::expect_character(valid.names, min.len = 1, any.missing = FALSE)

  dropped <- list()

  df.mod.list <- df.mod.list %>%
    purrr::imap(function(specs, species) {
      extra <- setdiff(names(specs), valid.names)
      if (length(extra) > 0) {
        dropped[[species]] <<- extra
        specs <- specs[names(specs) %in% valid.names]
      }
      specs
    })

  if (length(dropped) > 0) {
    warning(sprintf(
      paste0("prune.df.mod.list: dropped ddf spec(s) for ddftypes that are ",
             "not present in this study area:\n%s\nThis is expected when a ",
             "SubProject's study area has no coverage by one survey type - ",
             "the per-species blocks still set them, and they are ignored."),
      paste(sprintf("  %s: %s", names(dropped),
                    purrr::map_chr(dropped, paste, collapse = ", ")),
            collapse = "\n")),
      immediate. = TRUE)
  }

  df.mod.list
}


#' Determine which survey types a dataset actually contains
#'
#' Returns the distinct, sorted \code{SurveyType} values present in \code{dat}.
#' Most study areas are covered by both ship and aerial surveys, but some are
#' covered by only one (for example a study area no aerial survey has ever
#' flown), in which case the DDF types for the absent survey type have no data
#' at all and should be pruned - see \code{\link{prune.ddf.globals}}.
#'
#' @param dat Data frame (typically \code{the.data$watches}) with a
#'   \code{SurveyType} column.
#' @return Character vector of the survey types present, sorted.
#' @examples
#' \dontrun{
#' get.survey.types(the.data$watches)
#' #> [1] "Ship"
#' }
#' @export
get.survey.types <- function(dat) {

  checkmate::expect_data_frame(dat)
  checkmate::expect_names(names(dat), must.include = "SurveyType")

  types <- dat$SurveyType %>%
    as.character() %>%
    stats::na.omit() %>%
    unique() %>%
    sort()

  if (length(types) == 0)
    stop("get.survey.types: no non-NA SurveyType values found.")

  types
}


#' Prune the DDF globals to the survey types actually present
#'
#' \code{ddftype_levels}, \code{def.ddf.list} and \code{ddftype_to_platform}
#' each enumerate all combinations of survey type, behaviour and distance
#' availability, and are coupled \strong{positionally} - element \emph{i} of
#' each refers to the same DDF type (see the warning in
#' \code{analysis_settings.R} that they must be listed in the same order).
#' When a study area has no coverage by one survey type, that type's DDF types
#' have no observations and no segments, and carrying them forward creates
#' empty result folders, fails final DDF fitting, and trips the segment
#' bookkeeping in \code{\link{create.dsm.data}}.
#'
#' This prunes all three together, by position, so they cannot drift apart.
#' Pruning is a no-op when every survey type is present.
#'
#' Because the absent survey type contributes no segments, dropping its DDF
#' types yields the same \code{segdata} as keeping them as empty copies - the
#' change is structural, not statistical.
#'
#' @param survey.types Character vector of survey types present, as returned by
#'   \code{\link{get.survey.types}} (e.g. \code{c("Aerial", "Ship")}).
#' @param ddftype_levels Character vector of DDF type codes, whose first
#'   character is the survey type initial (e.g. \code{"SWD"}).
#' @param def.ddf.list Named list of default DDF specs, in the same order as
#'   \code{ddftype_levels}.
#' @param ddftype_to_platform Named character vector mapping DDF type to
#'   platform, named by \code{ddftype_levels}.
#' @param quiet If \code{TRUE}, suppress the message reporting what was pruned.
#' @return Named list with the pruned \code{ddftype_levels},
#'   \code{def.ddf.list} and \code{ddftype_to_platform}.
#' @examples
#' \dontrun{
#' pruned <- prune.ddf.globals(c("Ship"), ddftype_levels, def.ddf.list,
#'                             ddftype_to_platform)
#' pruned$ddftype_levels
#' #> [1] "SWD" "SFD" "SWN" "SFN"
#' }
#' @export
prune.ddf.globals <- function(survey.types,
                              ddftype_levels,
                              def.ddf.list,
                              ddftype_to_platform,
                              quiet = FALSE) {

  checkmate::expect_character(survey.types, min.len = 1, any.missing = FALSE)
  checkmate::expect_character(ddftype_levels, min.len = 1, any.missing = FALSE)
  checkmate::expect_list(def.ddf.list, min.len = 1)
  checkmate::expect_character(ddftype_to_platform, min.len = 1)
  checkmate::expect_flag(quiet)

  # The three globals are coupled positionally, so they must line up before we
  # prune by position.
  if (length(def.ddf.list) != length(ddftype_levels))
    stop(sprintf(
      paste0("prune.ddf.globals: def.ddf.list has %d entries but ",
             "ddftype_levels has %d. They must correspond one-to-one, in order."),
      length(def.ddf.list), length(ddftype_levels)))

  if (!identical(names(ddftype_to_platform), ddftype_levels))
    stop(paste0("prune.ddf.globals: names(ddftype_to_platform) must be exactly ",
                "ddftype_levels, in the same order."))

  # A ddftype belongs to a survey type if its first character matches that
  # survey type's initial ("S" for Ship, "A" for Aerial).
  known_initials <- unique(substr(ddftype_levels, 1, 1))
  wanted_initials <- unique(substr(survey.types, 1, 1))

  unknown <- setdiff(wanted_initials, known_initials)
  if (length(unknown) > 0)
    stop(sprintf(
      "prune.ddf.globals: survey type(s) %s do not match any ddftype in %s.",
      paste(sQuote(survey.types[substr(survey.types, 1, 1) %in% unknown]),
            collapse = ", "),
      paste(sQuote(ddftype_levels), collapse = ", ")))

  keep <- substr(ddftype_levels, 1, 1) %in% wanted_initials

  if (!any(keep))
    stop("prune.ddf.globals: pruning would remove every ddftype.")

  dropped <- ddftype_levels[!keep]

  if (!quiet) {
    if (length(dropped) == 0) {
      message(sprintf("Survey types present: %s. All %d ddftypes retained.",
                      paste(survey.types, collapse = ", "),
                      length(ddftype_levels)))
    } else {
      message(sprintf(
        paste0("Survey types present: %s. Dropping %d ddftype(s) with no ",
               "possible data: %s."),
        paste(survey.types, collapse = ", "),
        length(dropped),
        paste(dropped, collapse = ", ")))
    }
  }

  list(
    ddftype_levels = ddftype_levels[keep],
    def.ddf.list = def.ddf.list[keep],
    ddftype_to_platform = ddftype_to_platform[keep]
  )
}


#' Process all DDF specs for a given species
#'
#' Iterates over every element of \code{df.specs} calling
#' \code{\link{do.det.fcn.spec}} for each, then saves the updated spec list
#' with \code{\link{save.ddf.specs}}.  Called from
#' \code{Generic_1_ddf_fitting.Rmd}.
#'
#' @param df.specs Named list of DDF specification lists, one per
#'   dataset/platform/behaviour combination.
#' @param species Character string species code.
#' @param distdata Data frame of observation data for the species.
#' @param do.final.only If \code{TRUE}, fit only the pre-selected final model;
#'   otherwise fit all candidate models.
#' @param rerun If \code{TRUE}, re-fit models even when saved results exist.
#' @param parallel If \code{TRUE}, use parallel processing via
#'   \code{doSNOW}.
#' @param nCores Number of cores for parallel processing.
#' @param do.eda If \code{TRUE} (default), produce exploratory data plots.
#' @param ... Additional arguments passed to \code{\link{do.det.fcn.spec}}.
#' @return Updated \code{df.specs} list with fitted models attached.
#' @export
do.det.fcn.specs <-
  function(df.specs,
           species,
           distdata,
           do.final.only,
           rerun,
           parallel,
           nCores,
           do.eda = TRUE,
           ...) {

    ret <- df.specs %>%
      purrr::imap(
        do.det.fcn.spec,
        species = species,
        distdata = distdata,
        do.final.only = do.final.only,
        rerun = rerun,
        parallel = parallel,
        nCores = nCores,
        do.eda = do.eda,
        ...
      )

    # resave ddf spec for this species which may now have final fitted model
    #  and fitted distdata objects
    save.ddf.specs(ret, species)

    ret
  }

#' Save DDF specs for a single species to an RData file
#'
#' @param df.mod.specs Named list of DDF spec lists for the species.
#' @param species Character string species code; used to form the filename.
#' @return \code{invisible(NULL)}, called for its side-effect (file written).
#' @export
save.ddf.specs <- function(df.mod.specs, species) {
  filename <- file.path(RDataDir, paste0(species, dfModListSuffix))
  message(sprintf("Saving ddf specs for %s in '%s'", species, filename))
  save(df.mod.specs, file = filename)
}


#' Fit detection function(s) for one DDF spec and species
#'
#' Filters \code{distdata} for the appropriate dataset, platform class,
#' behaviour, and distance type; runs EDA plots; then calls
#' \code{\link{do.det.fcn}} to fit either all candidate detection functions or
#' just the pre-selected final model.  Stores the fitted model and augmented
#' distdata back into \code{df.spec} when \code{do.final.only = TRUE}.
#'
#' Expects project globals \code{spec.grps}, \code{seasons},
#' \code{ResultsDir}, and \code{ddftype_levels} in the calling environment.
#'
#' @param df.spec A single DDF specification list.
#' @param df.spec.nm Name of the spec, e.g. \code{"ECSAS_Ship_W_D"}.
#' @param species Character string species code.
#' @param distdata Data frame of observation data (all species; filtered
#'   internally).
#' @param do.final.only If \code{TRUE}, fit only the final model and store it
#'   in \code{df.spec}.
#' @param rerun If \code{TRUE}, re-fit even when saved results exist.
#' @param parallel If \code{TRUE}, use parallel processing.
#' @param nCores Number of cores for parallel processing.
#' @param do.eda If \code{TRUE}, produce EDA plots.
#' @param ... Additional arguments passed to \code{\link{do.det.fcn}}.
#' @return Updated \code{df.spec} with \code{fitted.model},
#'   \code{fitted.distdata}, and \code{ddftype} added when
#'   \code{do.final.only = TRUE}.
#' @export
do.det.fcn.spec <- function(df.spec,
                            df.spec.nm,
                            species,
                            distdata,
                            do.final.only,
                            rerun,
                            parallel,
                            nCores,
                            do.eda,
                            ...) {

  dataset <- get.dataset(df.spec.nm)
  platform_class <- get.platform.class(df.spec.nm)
  behav <- get.behav(df.spec.nm)
  dist_type <- get.dist_type(df.spec.nm)

  message(sprintf(
    "Doing %s, %s, for (%s, %s, %s, %s) ",
    species,
    ifelse(
      do.final.only,
      "final detection function",
      "all candidate det fcns"
    ),
    dataset,
    platform_class,
    behav,
    dist_type
  ))
  plot_title_suffix <- paste0("(",
                              paste(dataset, platform_class, behav, dist_type, sep = ", "),
                              ")")

  # Folder to save ddf fitting results when doing all candidate ddfs.
  folder <-
    file.path(ResultsDir,
              species,
              paste("DF Summaries", dataset, platform_class, behav, dist_type, sep = "_"))

  # Filter distdata based on whether this is a ddf for obs with no perp
  # distances, or a normal one. Record this fact in distdata for use by
  # create.dsm.data).
  if (dist_type == "N") {
    distdata <- dplyr::filter(
      distdata,
      dataset == dataset,
      SurveyType == platform_class,
      FlySwim == behav,
      Alpha %in% spec.grps[[species]],
      # only doing 300m strip transect, ignore others (only a few)
      !(DistanceCode %in% c("F", "G", "I", "K", "5")),
      # Weed out ones with NA distance but width of transect was more than 300m
      TransFarEdge <= 300 | is.na(TransFarEdge),
      (DistType != "Perp." | is.na(distance))
    )
  } else if (dist_type == "D"){
    # Normal ddf
    # Filter perp distances for this species, behav, platform_class, dataset
    distdata <- dplyr::filter(
      distdata,
      dataset == dataset,
      SurveyType == platform_class,
      FlySwim == behav,
      Alpha %in% spec.grps[[species]],
      DistType == "Perp.",
      # Only use obs with perp. distances.
      !is.na(distance)
    )
  } else
    stop("do.det.fcn.spec: unknown dist_type: ", dist_type)

  # assign season
  distdata <- assign.season(distdata, seasons[[species]])

  # # XXXX should make dynamic instead of relying on a static list.
  if (isTRUE(df.spec$remove.fishing))
    distdata %<>% dplyr::filter(!(CruiseID %in% fishingCruises))

  # EDA
  message("\nDistribution of distances:")
  print(table(distdata$distance, useNA = "always"))
  cat("\n")
  print(summary(distdata))
  # Look at covar plots - datasets with only 1 point will mess up geom_violin
  # so don't bother.
  if (isTRUE(do.eda) && nrow(distdata) > 1){
    ds.eda(
      distdata,
      species = species,
      do.final.only = do.final.only,
      df.spec = df.spec,
      suffix = plot_title_suffix)
  }


  # Only need to call do.det.fcn() if we're doing final ddfs (no matter if strip
  # trnasect or normal), OR we're doing all possible ddfs and dist_type is "D".
  # In the latter case, we are doing all possible ddfs as a side effect and
  # ignoring the return value from do.ds() called in do.det.fcn. We only bother
  # to call it to fit all possible models if dist_type is "D" since if it was
  # "N" then we already know what the final model will be (dummy_ddf) and there
  # is no need to fit "all possible" models since there aren't any.
  #
  # XXX Need to think about moving much of the above code into this conditional
  # since most of it doesn't need to be run in do.final == TRUE.
  if(do.final.only == TRUE || dist_type == "D"){
    # fit det fcn(s)
    res <- do.det.fcn(distdata = distdata,
                      dsetname = dataset,
                      species = species,
                      df.model = df.spec,
                      df.model.nm = df.spec.nm,
                      folder = folder,
                      do.final.only = do.final.only,
                      rerun = rerun,
                      parallel = parallel,
                      nCores = nCores,
                      ...)
    # When do.final.only is FALSE, res is either a list of models or a dataframe
    # of model specifications (if runModels was FALSE in do.ds) so no point in
    # trying to get the df_final and distdata elements.
    if (do.final.only) {
      df.spec$fitted.model <- res$df_final
      df.spec$fitted.distdata <- res$distdata
      # ddftype combines info on platform_class (aerial, ship) and behav(W, F, S).
      df.spec$ddftype <- substr(platform_class, 1, 1) %>% # "S" or "A"
        paste0(behav, dist_type) %>%
        factor(levels = ddftype_levels)
    }
  }

  # return the per-species df.model (with newly added fitted.model object if
  # do.final.only was true)
  df.spec
}

#' Create a strip-transect dummy detection function
#'
#' Wraps \code{\link[dsm]{dummy_ddf}} to create a strip-transect detection
#' function appropriate for the given dataset and platform class, using the
#' project-standard truncation distances.
#'
#' @param distdata Data frame of observation data for this DDF spec.
#' @param dataset Character string dataset name (currently \code{"ECSAS"}).
#' @param platform_class Character string; \code{"Ship"} or \code{"Aerial"}.
#' @return A \code{fake_ddf} object with \code{$data} set to \code{distdata}.
#' @export
create.strip.ddf <- function(distdata, dataset, platform_class) {
  # Set up appropriate one strip transect
  if (dataset == "ECSAS") {
    if (platform_class == "Aerial") {
      if (!is.numeric(distdata$object)) {
        warning(
          "create.strip.ddf: converting object id to numeric. This should be done in 00.01_Extract_data",
          immediate. = TRUE
        )
      }

      # XXX need to replace calculation of left with info from the aerial
      # transect for each observer??? Since it isn't really the smallest
      # distance if there are obs and it isn't really 0 if there are no obs.
      df_final <- dsm::dummy_ddf(
        distdata$object %>%
          stringr::str_replace("[a-zA-Z_]+", "") %>%
          as.numeric(),
        distdata$size,
        left = ifelse(nrow(distdata) == 0, 0, min(distdata$distbegin)),
        width = ecsas.air.right.trunc
      )
    } else if (platform_class == "Ship") {
      if (!is.numeric(distdata$object)) {
        warning(
          "create.strip.ddf: converting object id to numeric. This should be done in 00.01_Extract_data",
          immediate. = TRUE
        )
      }

      # For ECSAS ship we have defined cutpoints and no left trunc
      df_final <- dsm::dummy_ddf(
        distdata$object %>%
          stringr::str_replace("[a-zA-Z_]+", "") %>%
          as.numeric(),
        distdata$size,
        width = max(ecsas.ship.ctpoints)
      )
    } else
      # Platform_class == ...
      stop(sprintf(
        "create.strip.ddf: unrecognized survey platform class: %s",
        platform_class
      ))
  } else
    # dataset == ...
    # Add PQ stuff here later
    stop(sprintf(
      "create.strip.ddf: unrecognized survey dataset name: %s",
      dataset
    ))

  # replace the minimal data the dummy_ddf adds to df_final with our own with
  # more columns so that downstream processing can grab the distdata from df_final
  df_final$data <- distdata
  df_final
}

#' Fit detection function(s) for one dataset/species combination
#'
#' Either fits the single pre-selected final detection function
#' (\code{do.final.only = TRUE}) or runs all candidate models via
#' \code{\link{do.ds}} (\code{do.final.only = FALSE}).  When fitting the
#' final model, checks it with \code{\link{check.det.fcn}} and returns
#' augmented distdata.
#'
#' @param distdata Data frame of observation data (already filtered for the
#'   relevant species, behaviour, and platform).
#' @param species Character string species code.
#' @param df.model DDF specification list; used when \code{do.final.only =
#'   TRUE}.
#' @param df.model.nm Name of the DDF spec (e.g. \code{"ECSAS_Ship_W_D"}).
#' @param dsetname Character string dataset name.
#' @param do.final.only If \code{TRUE} (default), fit only the final model.
#' @param cleanFolder If \code{TRUE}, remove existing results from
#'   \code{folder} before fitting (only used when \code{do.final.only =
#'   FALSE}).
#' @param rerun If \code{TRUE}, re-fit even when saved results exist.
#' @param parallel If \code{TRUE}, use parallel processing.
#' @param folder Path to the folder where candidate model summaries are saved.
#' @param nCores Number of cores for parallel processing.
#' @param ... Additional arguments passed to \code{\link[Distance]{ds}}.
#' @return When \code{do.final.only = TRUE}: named list with elements
#'   \code{df_final} (fitted model) and \code{distdata} (augmented with
#'   \code{detProb} and \code{adjSize}).  When \code{do.final.only = FALSE}:
#'   return value of \code{\link{do.ds}}.
#' @export
do.det.fcn <- function(distdata,
                       species,
                       df.model = NULL,
                       df.model.nm = NULL,
                       dsetname = NULL,
                       do.final.only = TRUE,
                       cleanFolder = FALSE,
                       rerun = FALSE,
                       parallel = FALSE,
                       folder,
                       nCores = 1,
                       ...) {

  # set up species and spill specific pathnames
  # source(here::here("species settings.r"), echo = T, local = TRUE)

  dataset <- get.dataset(df.model.nm)
  platform_class <- get.platform.class(df.model.nm)

  # Set truncation distance if necessary for surveys without known set cutpoints
  if (dataset == "ECSAS") {
    if (platform_class == "Aerial") {
      truncation <-
        list(left = ifelse(nrow(distdata) == 0, 0, min(distdata$distbegin)),
             right = ecsas.air.right.trunc)
    } else if (platform_class == "Ship") {
      truncation <- max(df.model$cutpoints)
    } else
      stop(sprintf("do.det.fcn: unrecognized survey platform class: %s", platform_class))
  } else
    # Add PQ stuff here later
    stop(sprintf("do.det.fcn: unrecognized survey dataset name: %s", dataset))

  if (do.final.only) {
    # CleanFolder only applies when fitting all possible ddfs.
    if (cleanFolder)
      warning("do.det.fcn: do.final == TRUE, ignoring cleanFolder == TRUE.",
              immediate. = TRUE)

    # Doing final ddf
    if (df.model$strip) {
      df_final <- create.strip.ddf(distdata, dataset, platform_class)
    } else { # df.model$strip == FALSE

      if (rerun) {
        # (re-)fit the model - Doing a real ddf

        # Keep every column. This used to narrow distdata to a hard-coded list
        # to shrink the saved .rda files, but the strip-transect path does not
        # narrow (create.strip.ddf restores the full frame into df_final$data),
        # so the two paths produced fitted.distdata with different columns and
        # create.dsm.data's bind_rows padded the difference with NA. Any column
        # added upstream is now carried by both paths automatically.
        #
        # Still drop the redundant distance or distbegin/distend cols, which
        # ds() refuses to accept together.
        distdata <- check.distdata.cols(distdata)

        ###### use do.ds machinery to re-fit final model
        df_final <- do.ds(
          data = distdata,
          runModels = TRUE,
          incl.adj = FALSE,
          parallel = parallel,
          folder = folder,
          rerun = TRUE,
          models = with(
            df.model,
            data.frame(
              label =  create.model.name(final.key, final.formula, final.adj),
              key = final.key,
              # Can't insert a formula into a dataframe and as.character
              # separates the ~ from the rest, so paste it back together
              # then turn it back into a formulat in do.ds
              form =  paste(as.character(final.formula), collapse = " "),
              adj = ifelse(is.null(final.adj), "", final.adj)
            )
          ),
          max_adjustments = df.model$max.adjustments,
          nadj = df.model$final.nadj,
          cutpoints = df.model$cutpoints,
          convert_units = df.model$convert.units,
          truncation = truncation
        ) %>%
          # Note df_final will be a list with one element and we need to peel it
          # off.
          purrr::pluck(1)

      } else { # rerun = FALSE
        # Look for file with results of running this model already

        # Sometimes I forget to add as.formula() to final model specs, so
        # check it here.
        stopifnot(class(df.model$final.formula) == "formula")

        nam <- with(df.model,
                    create.model.name(final.key, final.formula, final.adj))
        pat <- paste0("AIC.*_", nam, ".RData")
        file <- list.files(folder, pattern = pat, full.names = TRUE,)
        if (length(file) != 1) {
          message("do.det.fcn: do.final == TRUE and rerun == FALSE. Looking for existing ddf model results file. Found ", length(file), " file(s), should be 1. Quitting.")
          stop()
        }

        # get the model
        load(file)
        df_final <- model
      } # End rerun = FALSE
    } # End strip == FALSE

    # Check model and augment distdata with detprob and adjSize
    distdata <- check.det.fcn(df_final, species, df.model.nm)

    return(list(df_final = df_final, distdata = distdata))
  } else {  # do.final.only == FALSE, doing all possible models
    # Do all candidate ddfs
    # folder <- ifelse(is.null(folder), summaryDir, folder) # folder is required now.

    # Remove redundant distance or distbegin/distend cols.
    distdata <- check.distdata.cols(distdata)
    do.ds(
      data = distdata,
      covars = df.model$candidate.covars,
      cutpoints = df.model$cutpoints,
      folder = folder,
      convert_units = df.model$convert.units,
      max_adjustments = df.model$max.adjustments,
      truncation = truncation,
      parallel = parallel,
      nCores = nCores,
      rerun = rerun,
      cleanFolder = cleanFolder,
      ...)
  }
}

#' Distance sampling exploratory data analysis plots
#'
#' Produces covariate-vs-distance plots for each variable listed in the DDF
#' spec, dispatching to \code{\link{plot_covar}}.
#'
#' @param distdata Data frame of observation data.
#' @param species Character string species code; used in plot titles.
#' @param do.final.only If \code{TRUE}, plot only the covariates in the final
#'   model formula; otherwise plot all candidate covariates.
#' @param df.spec DDF specification list supplying the formula and candidate
#'   covariates.
#' @param suffix Character string appended to plot titles to identify the
#'   dataset/platform/behaviour combination.
#' @return \code{invisible(NULL)}, called for its side-effect (plots).
#' @export
ds.eda <-
  function(distdata = NULL,
           species = NULL,
           do.final.only = T,
           df.spec = NULL,
           suffix = "") # to identify dataset, platform_class, and behaviour
  {

    if (nrow(distdata) == 0)
      return()

    if (do.final.only == TRUE) {
      vars <- all.vars(df.spec$final.formula)
    } else {
      vars <- df.spec$candidate.covars
    }

    # plot individual covars
    purrr::walk(
      vars,
      plot_covar,
      distdata = distdata,
      species = species,
      df.spec = df.spec,
      suffix = suffix
    )
  }


#' Plot a covariate against distance categories
#'
#' Produces a boxplot (for \code{"size"}), a violin plot (for factor
#' covariates), or a scatter plot with a linear smoother (for continuous
#' covariates).
#'
#' Named with an underscore rather than the project's usual dot separator:
#' \code{plot} is an S3 generic, so a function called \code{plot.covar} is
#' registered by roxygen as a \code{plot} method for the (non-existent) class
#' \code{"covar"} instead of being exported, and is then not callable by name.
#' This takes a covariate \emph{name}, not an object, so it is not a method.
#'
#' @param covar Character string name of the covariate column to plot.
#' @param distdata Data frame of observation data.
#' @param species Character string species code; used in the plot title.
#' @param df.spec DDF specification list supplying cutpoints and distance
#'   centres.
#' @param suffix Character string appended to the plot title.
#' @return \code{invisible(NULL)}, called for its side-effect (plot).
#' @export
plot_covar <-
  function(covar = NULL,
           distdata = NULL,
           species = NULL,
           df.spec = NULL,
           suffix = "") # to identify dataset, platform_class, and behaviour
  {
    title <- paste(species, suffix, sep = " ")
    if (covar == "size") {
      p <- ggplot2::ggplot(data = as.data.frame(distdata), ggplot2::aes(cut(distance, df.spec$cutpoints, right = FALSE), size))
      p <- p + ggplot2::geom_boxplot(varwidth = TRUE)
      p <- p + ggplot2::labs(x = "Distance Category", y = "Size", title = title)
      print(p)
    } else if (is.factor(distdata[, covar][[1]])) {
      p <- ggplot2::ggplot(data = as.data.frame(distdata), ggplot2::aes(get(covar), distance))
      p <- p + ggplot2::geom_violin(draw_quantiles = c(.25, .5, .75), scale = "count")
      p <- p + ggplot2::scale_y_continuous(labels = as.character(df.spec$distance.centers),
                                           breaks = df.spec$distance.centers)
      p <- p + ggplot2::labs(y = "Distance Category", x = covar, title = title)
      print(p)
    } else {
      p <- ggplot2::ggplot(data = as.data.frame(distdata), ggplot2::aes(get(covar), distance))
      p <- p + ggplot2::geom_point()
      p <- p + ggplot2::geom_smooth(method = "lm")
      p <-
        p + ggplot2::scale_y_continuous(
          labels = as.character(df.spec$distance.centers),
          breaks = df.spec$distance.centers
        )
      p <- p + ggplot2::labs(y = "Distance Category", x = covar, title = title)
      print(p)
    }
  }

#' Make an sf linestring from a two-row coordinate matrix
#'
#' @param xy2 Numeric matrix with 2 rows and 2 columns (lon, lat).
#' @return An \code{sfg} linestring object.
#' @export
make.line <- function(xy2){
  sf::st_linestring(matrix(xy2, nrow=2, byrow=TRUE))
}

#' Convert endpoint coordinate columns to an sfc linestring geometry
#'
#' @param df Data frame containing start/end coordinate columns.
#' @param names Character vector of four column names giving start lon, start
#'   lat, end lon, end lat.
#' @param crs CRS to assign to the returned geometry.
#' @return An \code{sfc_LINESTRING} object with one element per row of
#'   \code{df}.
#' @export
make.lines <- function(df, names=c("LongStart","LatStart","LongEnd","LatEnd"), crs){
  m = as.matrix(df[,names])
  lines = apply(m, 1, make.line, simplify=FALSE)
  sf::st_sfc(lines, crs = crs)
}

#' Convert a data frame of endpoint coordinates to an sf linestring object
#'
#' @param df Data frame with coordinate columns.
#' @param names Character vector of four column names (start lon, start lat,
#'   end lon, end lat).
#' @param crs CRS to assign to the output.
#' @return \code{df} as an \code{sf} object with a linestring geometry column.
#' @export
sf.pts.to.lines <- function(df, names=c("LongStart","LatStart","LongEnd","LatEnd"), crs){
  geom = make.lines(df, names, crs)
  df = sf::st_sf(df, geometry=geom)
  df
}

#' Assign season labels to observations based on date
#'
#' Matches each row's month-day to a season definition and adds a
#' \code{Season} factor column.  Works correctly across the new-year boundary
#' when a season definition spans December into January.
#'
#' @param dat Data frame (or sf object) to annotate.
#' @param season.def Named list of season definitions; each element is a
#'   named numeric vector with \code{"from"} and \code{"to"} entries encoded
#'   as \code{MMDD} integers.
#' @param datefield Character string giving the name of the date column in
#'   \code{dat}.
#' @param season.levels Factor levels for the returned \code{Season} column.
#'   \code{NULL} (the default) reproduces the historical behaviour exactly: the
#'   project global \code{season.names}, which this function used to reach for
#'   by lexical scope with no way for a caller to say otherwise, and which
#'   therefore made it unusable from any context that had not sourced
#'   \code{analysis_settings.R}. Falls back to \code{names(season.def)} when no
#'   such global exists.
#' @return \code{dat} with a \code{Season} factor column added.
#' @export
assign.season <- function(dat, season.def, datefield = "Date",
                          season.levels = NULL){

  if (is.null(season.levels))
    season.levels <- tryCatch(get("season.names", envir = globalenv()),
                              error = function(e) names(season.def))

  if (is.null(season.def)) {
    stop("assign.season: season.def is NULL. Could not find a season definition for this species!")
  }

  # Short circuit exit if no data
  if (nrow(dat) == 0)
    return(dat)

  # init
  dat$Season <- NA
  season.index <- rep(NA, nrow(dat))

  # convert dates to my format
  dates <- dplyr::pull(dat, datefield)
  dat$monthday <- lubridate::month(dates) * 100 + lubridate::day(dates)

  # step through each season with cheesy for loop
  for (i in seq_along(season.def)) {

    # Get logical index of matching rows
    if (season.def[[i]]["from"] > season.def[[i]]["to"]) {
      # wraps around new year
      criteria <-
        dat$monthday >= season.def[[i]]["from"] |
        dat$monthday <= season.def[[i]]["to"]
    } else {
      criteria <-
        dplyr::between(dat$monthday, season.def[[i]]["from"], season.def[[i]]["to"])
    }

    # assign the current season to the index of matching rows
    season.index[criteria] <- i
  }

  # Assign season
  dat$Season <- factor(
    names(season.def)[season.index],
    levels = season.levels,
  )

  if(any(is.na(dat$Season)))
    warning(sprintf("assign.season: %d rows failed to have season assigned",
                    sum(is.na(dat$Season))), immediate. = TRUE)
  dat
}


#' Write an sf object to a shapefile, skipping empty ones
#'
#' Wrapper around \code{\link[sf]{st_write}} that skips the write when
#' \code{dat} has no rows.  Writing an empty \code{sf} object fails for line
#' geometries because an empty \code{sfc} carries no geometry type for the ESRI
#' Shapefile driver to declare, and an empty layer is not useful in any case.
#' This arises whenever a study area has no coverage by one survey type.
#'
#' Note that when the write is skipped any pre-existing layer of the same name
#' is left untouched, so a stale file from an earlier run may remain on disk.
#' The skip is reported via \code{message()} so this is visible in the knitted
#' output.
#'
#' @param dat \code{sf} object to write.
#' @param dsn Directory path for the output shapefile.
#' @param layer Layer name (filename without extension) for the shapefile.
#' @param what Short description of the contents used in the messages
#'   (e.g. \code{"watches"}).
#' @param quiet If \code{TRUE}, suppress the "Saving ..." message.
#' @return \code{invisible(TRUE)} if written, \code{invisible(FALSE)} if
#'   skipped.
#' @export
st.write.if.any <- function(dat, dsn = ShapeDir, layer, what = "data",
                            quiet = FALSE) {

  checkmate::expect_class(dat, "sf")
  checkmate::expect_string(layer)
  checkmate::expect_string(what)

  if (nrow(dat) == 0) {
    message(sprintf(
      paste0("No %s to save - skipping shapefile '%s'. Any existing layer of ",
             "that name is now stale."),
      what, file.path(dsn, layer)))
    return(invisible(FALSE))
  }

  if (!quiet)
    message(sprintf("Saving %s to shapefile '%s'", what, file.path(dsn, layer)))

  # Will whine about discarded datum and abbreviated field names until these
  # warnings are removed from new rgdal.
  suppressWarnings(
    sf::st_write(
      dat,
      dsn = dsn,
      layer = layer,
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )
  )
  invisible(TRUE)
}


#' Draw a histogram, skipping empty data
#'
#' \code{\link[graphics]{hist}} errors with "invalid number of 'breaks'" on a
#' zero-length vector.  That happens for any diagnostic histogram of a dataset
#' or survey type with no coverage of the study area, which would otherwise
#' abort the knit partway through data extraction.
#'
#' Note the name deliberately does \strong{not} begin with \code{hist.}:
#' \code{hist} is an S3 generic, so roxygen registers any \code{hist.<x>}
#' function as a method for class \code{<x>} rather than exporting it, and it
#' would then not be callable by name.
#'
#' @param x Numeric vector to plot.
#' @param ... Further arguments passed to \code{\link[graphics]{hist}}.
#' @return The value of \code{hist()} if drawn, otherwise
#'   \code{invisible(NULL)}.
#' @examples
#' \dontrun{
#' draw.hist.if.any(numeric(0))  # no plot, no error
#' draw.hist.if.any(the.data$watches$WatchLenKm)
#' }
#' @export
draw.hist.if.any <- function(x, ...) {
  if (length(stats::na.omit(x)) == 0) {
    message("No data to plot - skipping histogram.")
    return(invisible(NULL))
  }
  graphics::hist(x, ...)
}


#' Convert a data frame to sf and write as a shapefile
#'
#' @param df Data frame to convert.
#' @param coords Character vector of two column names giving longitude and
#'   latitude.
#' @param crs Input CRS; default WGS 84 (EPSG 4326).
#' @param out.proj Target CRS to reproject to before writing.
#' @param dsn Directory path for the output shapefile.
#' @param layer Layer name (filename without extension) for the shapefile.
#' @return \code{invisible(NULL)}, called for its side-effect (file written).
#' @export
df.to.shapefile <- function(df,
                            coords = c("LongStart", "LatStart"),
                            crs = sf::st_crs(4326),
                            out.proj,
                            dsn = ShapeDir,
                            layer)
{
  sf::st_as_sf(df,
               coords = coords,
               crs = crs,
               remove = FALSE) %>%
    sf::st_transform(out.proj) %>%
    sf::st_write(
      dsn = dsn,
      layer = layer,
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )
}

#' Compute distances and their SD between successive GPS positions in a watch
#'
#' Adds \code{dists} and \code{dists.sd} columns to a single-row watch data
#' frame that contains a \code{posns} list column of GPS coordinates.
#'
#' @param watch Single-row data frame with a \code{posns} list column
#'   containing a matrix of GPS coordinates.
#' @return \code{watch} augmented with \code{dists} (list column of pairwise
#'   distances in metres) and \code{dists.sd} (SD of those distances, or
#'   \code{-1} if fewer than 3 positions).
#' @export
add.dist.sd <- function(watch) {
  if (nrow(watch$posns[[1]]) > 2) {
    watch$dists <-
      list(geodist::geodist(watch$posns[[1]], sequential = TRUE, measure = "geodesic"))
    watch$dists.sd <- round(sd(unlist(watch$dists)), 4)
  } else {
    watch$dists.sd <- -1
  }
  watch
}

#' Compute the GPS-track length of a watch in kilometres
#'
#' Filters \code{posns} to the time window of \code{watch} and sums geodesic
#' distances between consecutive positions.
#'
#' @param watch Single-row data frame with \code{WatchStartTime} and
#'   \code{WatchEndTime} columns.
#' @param posns Data frame of GPS positions with a \code{datetime} column.
#' @return Numeric scalar: total track length in km (0 if fewer than 2
#'   positions).
#' @export
get.gps.length <- function(watch, posns) {
  # get positions in this watch
  posns <- dplyr::filter(posns,
                         dplyr::between(posns$datetime, watch$WatchStartTime, watch$WatchEndTime))

  # sum distances between the points. Note that the units of val will
  # depend on the projection of the coords in posns, but we assume here
  # they are in meters
  if (nrow(posns) >= 2) {
    val <-
      list(geodist::geodist(posns, sequential = TRUE, measure = "geodesic")) %>%
      unlist %>%
      sum
  } else {
    val <- 0
  }

  val/1000
}


#' Rasterize one variable from a seasonal sf prediction grid
#'
#' Filters \code{obj} to \code{season}, converts to a \code{SpatVector}, and
#' rasterizes \code{variable} onto a grid whose resolution is
#' \code{predgridCellLength * 1000} metres (project global).
#'
#' @param season Character string season label used to filter \code{obj}.
#' @param obj \code{sf} data frame with a \code{Season} column and a column
#'   named \code{variable}.
#' @param variable Character string name of the column to rasterize.
#' @return A \code{SpatRaster} layer for the given season.
#' @export
make.season.raster <- function(season, obj, variable) {
  v <- dplyr::filter(obj, Season == season) %>%
    # Convert sf to SpatVector
    terra::vect()

  # Create a raster template with the same extent and resolution. Assumes
  # predgridCellLength is in km and raster projection units are metres.
  r <- terra::rast(v, resolution = predgridCellLength * 1000)

  # Rasterize, using an attribute field (e.g., "ID")
  ret <- terra::rasterize(v, r, field = variable)

  ret
}

#' Copy all objects from one environment to another
#'
#' Useful after reloading a background job result saved as an \code{.rda}
#' file: call \code{copy.env(loaded_env, .GlobalEnv)} to promote all objects
#' into the global environment.
#'
#' @param src Source environment.
#' @param dst Destination environment.
#' @return \code{invisible(NULL)}.
#' @export
copy.env <- function(src, dst) {
  for(n in ls(src, all.names=TRUE))
    assign(n, get(n, src), dst)
}

#' Build segment data with environmental covariates from watch data
#'
#' Creates the DSM segment data frame from \code{the.data$watches}, extracts
#' depth, depth gradient, SST, and SST gradient at each watch location from
#' pre-downloaded raster files, scales all covariates, and saves the result
#' as both an RData file and a shapefile.
#'
#' Expects project globals \code{predLayerStudyAreaDir}, \code{segdatloc},
#' and \code{ShapeDir} in the calling environment.
#'
#' \strong{This function does no spatial clipping.} It carried a
#' \code{study.area} argument until DSMHelper 0.15.0 and never referenced it,
#' while the documentation claimed the function clipped to the study area - it
#' does not, and the argument is gone. A segment's presence in the output is
#' decided entirely upstream, by the \code{clip.area} passed to
#' \code{\link{create.survey.data}}; its covariates come from the rasters in
#' \code{predLayerStudyAreaDir}, whose extent is set by \code{genpredrast} in
#' \code{00.02_Extract_env_rasters.rmd}. Those two must agree: a segment
#' outside the raster extent gets NA covariates here and is then dropped
#' outright by the NA filter in \code{Generic_2_dsm.Rmd}. That is the failure
#' mode to check first if a study-area buffer appears to have done nothing.
#'
#' @param the.data Named list with at least a \code{watches} element (as
#'   returned by \code{\link{create.survey.data}}).
#' @param inproj EPSG code or CRS for input watch coordinates.
#' @param outproj EPSG code or CRS for output segdata.
#' @param scale.factors Named list of mean and SD values used to standardise
#'   each covariate (e.g. \code{list(depth_mean = x, depth_sd = y, ...)}).
#' @param verbose If \code{TRUE}, print progress messages.
#' @return \code{sf} data frame of segment data with covariates attached.
#' @export
create.segdata <- function(the.data,
                           inproj,
                           outproj,
                           scale.factors,
                           verbose = FALSE) {

  # create initial segdata and reproject
  if (verbose) message("Creating initial segdata from watches")
  segdata <- the.data$watches %>%
    dplyr::rename(Effort = WatchLenKm) %>%
    dplyr::mutate(TransectID = as.character(TransectID),
                  year = lubridate::year(Date),
                  yday = lubridate::yday(Date),
                  MonthYear = format(Date, "%Y-%m"),
                  segment.area = Effort * TotalWidthKm) %>%
    sf::st_as_sf(coords = c("LongStart", "LatStart"), crs = inproj, remove = FALSE) %>%
    sf::st_transform(outproj) %>%
    cbind(sf::st_coordinates(.)) %>%
    dplyr::rename(x = X, y = Y)

  ###---------------------------------------------------------------------------
  #### Load rasters
  if (verbose) message("Loading environmental rasters")

  # get the dates of monthly rasters needed (ie. sst, sst.g, etc) so we can read
  # the needed files into a big SpatRaster
  dates.needed <- segdata$Date %>%
    as.character %>%
    stringr::str_sub(end = -4) %>%
    paste0("-16") %>%
    unique %>%
    stringr::str_sort()

  ### Depth and other static rasters
  # Depth
  depth <- terra::rast(file.path(predLayerStudyAreaDir, "depth.img"))

  # Depth gradient
  depth.g <- terra::rast(file.path(predLayerStudyAreaDir, "depth.g.img"))

  # SST
  #
  # A month with no raster is dropped here rather than being allowed to error.
  # terra::rast() on a missing file stops the whole run, which is a poor way to
  # find out that four watches fall outside the SST download - and the data can
  # widen without warning: a study-area buffer keeps effort that used to be
  # clipped away (issue #33), a new SOMEC delivery arrives with later dates, the
  # ECSAS year range moves. Those rows get NA sst and sst.g, and the NA filter in
  # Generic_2_dsm.Rmd then drops them exactly as it drops a segment with no
  # depth. Loud, because the right fix is usually to download the missing month.
  sst.available <- dates.needed[
    file.exists(file.path(predLayerStudyAreaDir, "sst",
                          paste0("sst.", dates.needed, ".img"))) &
    file.exists(file.path(predLayerStudyAreaDir, "sst",
                          paste0("sst.g.", dates.needed, ".img")))]

  sst.missing <- setdiff(dates.needed, sst.available)
  if (length(sst.missing)) {
    n.affected <- sum(paste0(segdata$MonthYear, "-16") %in% sst.missing)
    warning(sprintf(
      paste0("create.segdata: no SST raster for %d of %d months needed (%s). ",
             "%d of %d segments get NA sst/sst.g and will be dropped by the ",
             "NA filter before fitting. Download the missing month(s) if those ",
             "segments matter."),
      length(sst.missing), length(dates.needed),
      paste(sub("-16$", "", sst.missing), collapse = ", "),
      n.affected, nrow(segdata)), immediate. = TRUE)
  }

  if (!length(sst.available))
    stop("create.segdata: no SST raster exists for any month in the segment data.")

  sst <- terra::rast(file.path(predLayerStudyAreaDir, "sst",
                               paste0("sst.", sst.available, ".img")))

  # SST gradient
  sst.g <- terra::rast(file.path(predLayerStudyAreaDir, "sst",
                                 paste0("sst.g.", sst.available, ".img")))

  ###---------------------------------------------------------------------------
  ### Extract raster values at segdata locations

  ### Extract depth values
  if (verbose) message("Extracting depth at watch locations")
  segdata <-
    terra::extract(
      depth,
      terra::vect(segdata),
      bind = TRUE
    ) %>%
    sf::st_as_sf()

  ### Extract depth gradient values
  if (verbose) message("Extracting depth gradient at watch locations")
  segdata <-
    terra::extract(
      depth.g,
      terra::vect(segdata),
      bind = TRUE
    ) %>%
    sf::st_as_sf()

  ### Extract SST and SST gradient values
  #
  # Rows whose month has no raster are held out of the extract and given NA,
  # because terra::extract(layer = ) errors on a layer name it does not hold.
  # Both covariates use the same set of months, so one split serves both.
  # Ordering is restored explicitly rather than assumed.
  has.sst <- paste0(segdata$MonthYear, "-16") %in% sst.available

  extract.by.month <- function(dat, rast, prefix, out.name) {
    if (all(!has.sst)) {
      dat[[out.name]] <- NA_real_
      return(dat)
    }
    keep <- dat[has.sst, ]
    got <- terra::extract(
      rast,
      terra::vect(keep),
      layer = paste0(prefix, keep$MonthYear, "-16"),
      bind = TRUE
    ) %>%
      sf::st_as_sf() %>%
      dplyr::rename(!!out.name := value) %>%
      dplyr::select(-layer) # layer is off by one, though the value is correct

    if (all(has.sst)) return(got)

    drop <- dat[!has.sst, ]
    drop[[out.name]] <- NA_real_
    out <- rbind(got, drop[, names(got)])
    # Back into the incoming row order - rbind put the held-out rows last.
    out[order(c(which(has.sst), which(!has.sst))), ]
  }

  if (verbose) message("Extracting sst at watch locations")
  segdata <- extract.by.month(segdata, sst, "sst.", "sst")

  if (verbose) message("Extracting sst gradient at watch locations")
  has.sst <- paste0(segdata$MonthYear, "-16") %in% sst.available
  segdata <- extract.by.month(segdata, sst.g, "sst.g.", "sst.g")

  ###---------------------------------------------------------------------------
  ## Add scaled versions of all preds.
  if (verbose) message("Scaling covars")

  segdata %<>%
    dplyr::mutate(depth.sc = (depth - scale.factors$depth_mean)/scale.factors$depth_sd,
                  depth.g.sc = (depth.g - scale.factors$depth.g_mean)/scale.factors$depth.g_sd,
                  x.sc = (x - scale.factors$x_mean)/scale.factors$x_sd,
                  y.sc = (y - scale.factors$y_mean)/scale.factors$y_sd,
                  sst.sc = (sst - scale.factors$sst_mean)/scale.factors$sst_sd,
                  sst.g.sc = (sst.g - scale.factors$sst.g_mean)/scale.factors$sst.g_sd,
    )

  ###---------------------------------------------------------------------------
  ## Save segdata
  if (verbose) message("Saving results")
  save(segdata, file = segdatloc)

  st.write.if.any(segdata, layer = "segdata.shp", what = "segdata")

  if (verbose) print("Done.")
  segdata
}

#' Download one ERDDAP raster via rxtractogon
#'
#' \strong{Note: no longer used} due to timeout and lag issues; retained for
#' reference.  Calls \code{rerddapXtracto::rxtractogon()} for the given time
#' coordinate and returns a \code{RasterLayer}.
#'
#' @param tcoord Character string time coordinate (e.g.
#'   \code{"2020-06-15"}), or \code{NULL} for static datasets such as
#'   bathymetry.  Must be a single value; use with \code{purrr::map()} to
#'   iterate over dates.
#' @param dataset ERDDAP dataset identifier (e.g.
#'   \code{"jplMURSST41mday"}).
#' @param parameter Name of the variable within the dataset (e.g.
#'   \code{"sst"}).
#' @param xcoord Numeric vector of longitudes defining the bounding polygon.
#' @param ycoord Numeric vector of latitudes defining the bounding polygon.
#' @param plotit If \code{TRUE}, plot the downloaded raster.
#' @param saveit If \code{TRUE}, write the raster to \code{folder} in HFA
#'   format.
#' @param folder Directory path for saving the raster when \code{saveit =
#'   TRUE}.
#' @return A \code{RasterLayer} in the input projection of the ERDDAP source.
#' @export
rxtractogon.rast <-
  function(tcoord,
           dataset,
           parameter,
           xcoord,
           ycoord,
           plotit = FALSE,
           saveit = FALSE,
           folder) {

    # Some dataset (ie depth) don't require a tcoord
    if (is.null(tcoord)) {
      print(sprintf("Extracting %s from %s", parameter, dataset))
    } else {
      if (length(tcoord) != 1)
        stop(sprintf("tcoord should have length 1 but has length %d", length(tcoord)))

      print(sprintf("Extracting %s from %s for %s", parameter, dataset, tcoord))
    }

    dat <-
      rerddapXtracto::rxtractogon(
        rerddap::info(dataset),
        parameter =  parameter,
        xcoord = xcoord,
        ycoord = ycoord,
        tcoord = tcoord
      )
    layername <- names(dat)[1]
    dat <- purrr::pluck(dat, 1) # use 1 since it is always 1st element, but not always called same as value of parameter
    if (length(dim(dat)) > 2) # remove useless third dimension
      dat <- dat[,,1]

    rast <- dat %>%
      t %>%                   # rxtracto returns matrix in odd order with x and y transposed and south to north so: transpose
      .[nrow(.):1, ] %>%       # ... and reverse order of rows
      raster::raster(
        xmn = min(xcoord),
        xmx = max(xcoord),
        ymn = min(ycoord),
        ymx = max(ycoord),
        crs = latlongproj
      )

    if (plotit)
      # see note on terra::plot() in create.ncdf.rast - a bare plot() inside
      # this package resolves to graphics::plot and cannot draw a raster
      raster::plot(rast, main = sprintf("%s from %s for %s", parameter, dataset, tcoord))

    if (saveit) {

      if (is.null(tcoord))
        filename <- file.path(folder, paste(layername, "img", sep = "."))
      else
        filename <- file.path(folder, paste(layername, tcoord, "img", sep = "."))

      message("Saving downloaded ERDDAP raster to ", filename)
      if (!dir.exists(folder))
        dir.create(folder, recursive = TRUE)
      terra::writeRaster(rast, filename = filename, format = "HFA", overwrite = TRUE)
    }

    terra::rast()
  }

#' Apply a focal filter to every layer of a RasterStack
#'
#' \strong{Note: possibly no longer used.}
#'
#' @param x A \code{RasterStack} (or filename that can be read as one).
#' @param w Focal weight matrix passed to \code{\link[raster]{focal}}.
#' @param ... Additional arguments passed to \code{\link[raster]{focal}}.
#' @return A \code{RasterStack} with the focal operation applied to each
#'   layer.
#' @export
multi.focal <- function(x, w = matrix(1, nrow = 3, ncol = 3), ...) {

  if (is.character(x)) {
    x <- raster::brick(x)
  }
  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    terra::focal(x[[ind]], w = w, ...)
  }

  n <- seq(raster::nlayers(x))
  list <- lapply(X = n, FUN = fun, x = x, w = w, ...)

  out <- raster::stack(list)
  return(out)
}

#' Rebuild a raster stack from previously saved files
#'
#' Reads all files in \code{folder} whose names match \code{pattern} and
#' stacks them.  Useful for reloading cached rasters without re-downloading
#' from ERDDAP.
#'
#' @param folder Directory path to search.
#' @param pattern Regular expression matched against filenames (e.g.
#'   \code{"sst.+img$"}).
#' @return A \code{RasterStack} with layer names derived from the filenames.
#' @export
recreate.sst.mnth.from.files <- function(folder, pattern){
  files <- list.files(folder, pattern = pattern, full.names = T)
  r <- purrr::map(files, raster) %>%
    stack
  names(r) <- basename(files) %>%
    stringr::str_replace(stringr::fixed(".img"), "")
  r
}

#' Decode a CF-convention time axis to dates
#'
#' Handles the \code{"<unit> since <origin>"} form the CF conventions require,
#' e.g. \code{"seconds since 1970-01-01T00:00:00Z"}. Anything else returns
#' \code{NULL} rather than a guess - a wrong date here would be worse than no
#' date, because it would look like an answer.
#'
#' @param vals Numeric vector of time-axis values.
#' @param units The axis's \code{units} attribute.
#' @return A \code{Date} vector, or \code{NULL} if \code{units} is not a
#'   recognised CF time string.
#' @examples
#' decode.cf.time(1673827200, "seconds since 1970-01-01T00:00:00Z")
#' @export
decode.cf.time <- function(vals, units) {
  checkmate::assert_numeric(vals)
  checkmate::assert_string(units, na.ok = TRUE)

  if (is.na(units) || !grepl(" since ", units, fixed = TRUE)) return(NULL)

  parts  <- strsplit(units, " since ", fixed = TRUE)[[1]]
  unit   <- tolower(trimws(parts[1]))
  origin <- trimws(parts[2])

  # "1970-01-01T00:00:00Z" and "1970-01-01 00:00:00" are both legal.
  origin <- sub("T", " ", origin, fixed = TRUE)
  origin <- sub("Z$", "", origin)
  origin <- as.POSIXct(origin, tz = "UTC")
  if (is.na(origin)) return(NULL)

  mult <- switch(unit,
                 second = , seconds = , sec = , secs = 1,
                 minute = , minutes = , min = , mins = 60,
                 hour   = , hours   = , hr  = , hrs  = 3600,
                 day    = , days    = 86400,
                 NULL)
  if (is.null(mult)) return(NULL)

  as.Date(origin + vals * mult)
}


#' Spatial extent of an open netCDF handle
#'
#' Finds the longitude and latitude axes and returns the extent as cell
#' \emph{edges}, matching \code{terra::ext()} rather than the axis values,
#' which are cell centres. Half a cell is not much, but a coverage test that
#' is half a cell wrong is wrong in the direction that matters - it says a
#' file does not reach somewhere it does.
#'
#' Axes are found by their CF \code{units} (\code{degrees_east} /
#' \code{degrees_north}), falling back to the usual names only when a file
#' omits the attribute. Same reasoning as the time axis: the name is a
#' convention, the units are the standard.
#'
#' @param nc An open \code{ncdf4} handle.
#' @return Named list of \code{xmin}, \code{xmax}, \code{ymin}, \code{ymax},
#'   all \code{NA_real_} if the axes cannot be identified.
#' @keywords internal
nc.spatial.extent <- function(nc) {
  none <- list(xmin = NA_real_, xmax = NA_real_,
               ymin = NA_real_, ymax = NA_real_)

  units.of <- vapply(nc$dim, function(d)
    if (is.null(d$units)) NA_character_ else tolower(d$units), character(1))
  names.of <- tolower(names(nc$dim))

  pick <- function(unit, aliases) {
    i <- which(!is.na(units.of) & units.of == unit)
    if (!length(i)) i <- which(names.of %in% aliases)
    if (length(i)) nc$dim[[i[1]]]$vals else NULL
  }

  lon <- pick("degrees_east",  c("longitude", "lon", "x"))
  lat <- pick("degrees_north", c("latitude",  "lat", "y"))
  if (is.null(lon) || is.null(lat)) return(none)

  # Centres to edges. A single-valued axis has no spacing to halve.
  edges <- function(v) {
    v <- sort(v)
    if (length(v) < 2) return(c(v[1], v[1]))
    sp <- stats::median(diff(v))
    c(v[1] - sp / 2, v[length(v)] + sp / 2)
  }

  x <- edges(lon)
  y <- edges(lat)
  list(xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2])
}


#' Time coverage and spatial extent of a single netCDF file
#'
#' Reads the file's time axis - by CF units first, since the name is only a
#' convention - decodes it with \code{\link{decode.cf.time}}, and reads its
#' spatial extent with the same approach.
#'
#' Never throws: a file that will not open, has no time axis, or carries units
#' this cannot parse comes back as a row with \code{note} filled in and the
#' date columns \code{NA}. That is deliberate, because the point of scanning a
#' folder is to find the odd file out, and stopping on it would report the
#' problem by hiding every file after it. A file with no time axis still gets
#' its extent read - a bathymetry grid has no dates but very much has a box.
#'
#' @param f Path to a \code{.nc} file.
#' @return A one-row tibble: \code{file}, \code{n_times}, \code{first},
#'   \code{last}, \code{years}, \code{months}, \code{xmin}, \code{xmax},
#'   \code{ymin}, \code{ymax}, \code{note}. Extent is cell edges, as
#'   \code{terra::ext()} reports it.
#' @examples
#' \dontrun{
#' get.nc.file.times("GIS/Spatial covars/NetCDF/etopo180.nc")
#' }
#' @export
get.nc.file.times <- function(f) {
  checkmate::assert_file_exists(f)

  row <- function(note, ext = NULL, n_times = NA_integer_,
                  first = as.Date(NA), last = as.Date(NA),
                  years = NA_character_, months = NA_character_) {
    if (is.null(ext))
      ext <- list(xmin = NA_real_, xmax = NA_real_,
                  ymin = NA_real_, ymax = NA_real_)
    tibble::tibble(file = basename(f), n_times = n_times,
                   first = first, last = last,
                   years = years, months = months,
                   xmin = ext$xmin, xmax = ext$xmax,
                   ymin = ext$ymin, ymax = ext$ymax,
                   note = note)
  }

  nc <- try(ncdf4::nc_open(f), silent = TRUE)
  if (inherits(nc, "try-error"))
    return(row(paste("could not open:", trimws(attr(nc, "condition")$message))))
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # Extent first, so a file with no usable time axis still reports its box.
  ext <- nc.spatial.extent(nc)

  # By units, not by name: "time" is a convention, "since" is the standard.
  # Fall back to the usual names for a file that omits the units attribute.
  dim.units <- vapply(nc$dim, function(d)
    if (is.null(d$units)) NA_character_ else d$units, character(1))
  is.time <- !is.na(dim.units) & grepl(" since ", dim.units, fixed = TRUE)
  if (!any(is.time))
    is.time <- tolower(names(nc$dim)) %in% c("time", "t")
  if (!any(is.time)) return(row("no time dimension (static)", ext))

  d <- nc$dim[[which(is.time)[1]]]
  dates <- decode.cf.time(d$vals,
                          if (is.null(d$units)) NA_character_ else d$units)
  if (is.null(dates))
    return(row(paste0("time units not understood: ", d$units), ext))

  dates <- sort(dates)
  row("", ext,
      n_times = length(dates),
      first   = min(dates),
      last    = max(dates),
      years   = paste(sort(unique(format(dates, "%Y"))), collapse = ", "),
      months  = paste(sort(unique(format(dates, "%Y-%m"))), collapse = " "))
}


#' What years - and what area - is a folder of netCDF files actually from?
#'
#' Opens every netCDF in a folder and reports each file's time coverage and
#' spatial extent, then the folder as a whole: the overall span, months that
#' appear in more than one file, months missing from inside the span, and the
#' distinct bounding boxes present.
#'
#' It exists because the downloads are named by ERDDAP hash -
#' \code{jplMURSST41mday_14d0_c28b_14ce.nc} - so the filename says nothing
#' about what is inside, and a year that was never downloaded looks exactly
#' like a year that was. That matters because a segment whose month has no
#' raster is dropped by the NA filter before fitting, which is silent.
#'
#' \strong{Mixed extents are called out, because they are the quiet failure.}
#' These files are hand-downloaded a year at a time, so re-downloading part of
#' a series under a wider box is easy to do and leaves a folder whose covariate
#' silently changes footprint partway through the time series. One box is
#' reported as a single line; more than one is reported per box, with the files
#' listed.
#'
#' \strong{This reads the SOURCE netCDFs, not the processed rasters in
#' \code{predLayerStudyAreaDir}.} The two can disagree - the processed rasters
#' accumulate across runs and are not cleaned up, so a month whose source file
#' has since been removed can still have a raster. Comparing the two is the
#' point; do not read a clean report here as proof the pipeline has what it
#' needs, or vice versa.
#'
#' @param folder Folder to scan.
#' @param pattern Regex for the files to read.
#' @param report.gaps Also list months with no data between the earliest and
#'   latest timestep found anywhere in the folder. Assumes monthly coverage was
#'   wanted, which is what the SST downloads are; harmless otherwise, since a
#'   month counts as covered when any timestep falls inside it.
#' @return Invisibly, one row per file, as from \code{\link{get.nc.file.times}}.
#'   Prints a summary.
#' @examples
#' \dontrun{
#' report.nc.time.coverage("GIS/Spatial covars/NetCDF/sst")
#' }
#' @export
report.nc.time.coverage <- function(folder, pattern = "[.]nc$",
                                    report.gaps = TRUE) {
  checkmate::assert_directory_exists(folder)
  checkmate::assert_string(pattern)
  checkmate::assert_flag(report.gaps)

  files <- list.files(folder, pattern = pattern, full.names = TRUE,
                      ignore.case = TRUE)
  if (!length(files)) {
    message("No files matching ", pattern, " in ", folder)
    return(invisible(NULL))
  }

  res <- purrr::map_dfr(files, get.nc.file.times)

  cat(sprintf("\n%d netCDF file(s) in %s\n\n", nrow(res), folder))
  print(as.data.frame(res[, c("file", "n_times", "first", "last", "years",
                              "note")]),
        row.names = FALSE)

  ###--------------------------------------------------------------------------
  ### Spatial extent
  placed <- dplyr::filter(res, !is.na(xmin))
  if (!nrow(placed)) {
    cat("\nNo file has identifiable longitude/latitude axes.\n")
  } else {
    boxes <- placed %>%
      dplyr::count(xmin, xmax, ymin, ymax, name = "n_files") %>%
      dplyr::arrange(dplyr::desc(n_files))

    if (nrow(boxes) == 1) {
      cat(sprintf(
        "\nExtent: all %d file(s) share lon %.4f to %.4f, lat %.4f to %.4f\n",
        nrow(placed), boxes$xmin, boxes$xmax, boxes$ymin, boxes$ymax))
    } else {
      cat(sprintf("\nMIXED EXTENTS - %d different boxes in one folder:\n\n",
                  nrow(boxes)))
      print(as.data.frame(boxes), row.names = FALSE, digits = 8)
      cat("\nA covariate whose footprint changes partway through the series is",
          "\nalmost never what you want. Files by box:\n")
      for (i in seq_len(nrow(boxes))) {
        f <- placed$file[placed$xmin == boxes$xmin[i] &
                           placed$xmax == boxes$xmax[i] &
                           placed$ymin == boxes$ymin[i] &
                           placed$ymax == boxes$ymax[i]]
        cat(sprintf("  lon %.4f..%.4f lat %.4f..%.4f : %s\n",
                    boxes$xmin[i], boxes$xmax[i], boxes$ymin[i], boxes$ymax[i],
                    paste(f, collapse = ", ")))
      }
    }
    if (nrow(placed) < nrow(res))
      cat(sprintf("(%d file(s) had no identifiable lon/lat axes)\n",
                  nrow(res) - nrow(placed)))
  }

  ###--------------------------------------------------------------------------
  ### Time
  dated <- dplyr::filter(res, !is.na(first))
  if (!nrow(dated)) {
    cat("\nNo file in this folder carries a time axis.\n")
    return(invisible(res))
  }

  all.months <- sort(unique(unlist(strsplit(dated$months, " "))))
  cat(sprintf("\nOverall: %s to %s, %d distinct month(s), years %s\n",
              min(dated$first), max(dated$last), length(all.months),
              paste(sort(unique(unlist(strsplit(dated$years, ", ")))),
                    collapse = ", ")))

  # Overlap is worth knowing: these are hand-downloaded, so the same month
  # arriving in two files is easy to do and easy to miss.
  dup <- unlist(strsplit(dated$months, " "))
  dup <- names(which(table(dup) > 1))
  if (length(dup))
    cat(sprintf("Months present in more than one file: %s\n",
                paste(sort(dup), collapse = ", ")))

  if (report.gaps) {
    want <- format(seq(as.Date(paste0(min(all.months), "-01")),
                       as.Date(paste0(max(all.months), "-01")),
                       by = "month"), "%Y-%m")
    gaps <- setdiff(want, all.months)
    if (length(gaps))
      cat(sprintf("MISSING month(s) inside that span: %s\n",
                  paste(gaps, collapse = ", ")))
    else
      cat("No gaps: every month in that span is present.\n")
  }

  invisible(res)
}

#' Create a raster from one time slice of a NetCDF array
#'
#' Extracts the layer at position \code{index} from \code{dat}, flips the
#' y-dimension (which is inverted in NetCDF convention), optionally reprojects
#' or resamples to match \code{to}, and plots the result.
#'
#' @param the.date Date string labelling this layer (used in the plot title
#'   and passed as the raster name); use \code{""} for static datasets.
#' @param index Integer index of the time slice to extract from \code{dat}.
#' @param dat 3-D numeric array with dimensions (lon, lat, time), as returned
#'   by \code{ncdf4::ncvar_get()}.
#' @param datname Character string variable name (e.g. \code{"sst"}); used
#'   in the plot title.
#' @param x Numeric vector of longitude coordinates.
#' @param y Numeric vector of latitude coordinates.
#' @param inproj Projection string or EPSG code for the input data.
#' @param outproj Target projection for reprojection; ignored when \code{to}
#'   is supplied.
#' @param to Optional \code{SpatRaster} template; if supplied the result is
#'   reprojected, resampled, and masked to match it.
#' @return A \code{SpatRaster} for the requested time slice.
#' @export
create.ncdf.rast <-
  function(the.date,
           index,
           dat,
           datname,
           x,
           y,
           inproj,
           outproj = NULL,
           to = NULL) {

    message(sprintf(
      "Creating %s raster %s",
      datname,
      ifelse(as.character(the.date) == "", "", as.character(the.date))
    ))

    the.data <- switch(length(dim(dat)), NULL, dat, dat[,, index])

    # Make sure data dimensionality is sensible
    if (is.null(the.data))
      stop(paste0("create.ncdf.rast: illegal data dimension: ", dims ))

    res <-
      raster::raster(
        t(the.data),
        xmn = min(x),
        xmx = max(x),
        ymn = min(y),
        ymx = max(y),
        crs = sp::CRS(inproj)
      ) %>%
      raster::flip(direction = "y") %>%
      terra::rast()

    # is reprojection/resampling required
    if (!is.null(to)) {
      if (!is.null(outproj))
        warning("create.ncdf.rast: both 'outproj' and 'to' are provided, ignoring outproj",
                immediate. = TRUE)

      res <- terra::project(res, terra::rast(to), threads = TRUE) %>%
        terra::mask(study.area)
    } else if (!is.null(outproj))
      res <- terra::project(res, outproj, threads = TRUE) %>%
      terra::mask(study.area)

    # NB: terra::plot(), not bare plot(). This package's NAMESPACE has no
    # import() directives, so inside package code `plot` resolves up the
    # namespace chain to graphics::plot and never reaches terra's S4 generic
    # (which is only visible via the search path). A bare plot() on a
    # SpatRaster therefore dies with "invalid type passed to graphics
    # function", even though the identical call works from the console.
    terra::plot(res, main = paste(datname, the.date))
    res
  }



#' Convert a NetCDF file to a SpatRaster or RasterLayer
#'
#' Reads the named variable from a NetCDF file using \code{ncdf4} and returns
#' a single \code{SpatRaster} (multiple layers if the file contains multiple
#' time steps).
#'
#' @param filename Path to the NetCDF file.
#' @param dataset Character string variable name to extract (e.g.
#'   \code{"sst"}).
#' @param inproj Projection string or EPSG code for the NetCDF coordinates.
#' @param outproj Target projection; ignored when \code{to} is supplied.
#' @param to Optional \code{SpatRaster} template for reprojection and
#'   resampling.
#' @return A named \code{SpatRaster} (one layer per time step, or a single
#'   layer for static data).
#' @export
ncdf.to.raster <- function(filename,
                           dataset,
                           inproj,
                           outproj = NULL,
                           to = NULL) {

  message(sprintf("Extracting %s from netCDF file: %s", dataset, filename))

  x <- ncdf4::nc_open(filename)
  lon <- ncdf4::ncvar_get(x, "longitude")
  lat <- ncdf4::ncvar_get(x, "latitude")
  dat <- ncdf4::ncvar_get(x, dataset)

  # Get dimension descriptors
  dims <- tidync::tidync(filename) %>%
    tidync::activate(dataset) %>%
    tidync::hyper_dims()

  # netCDF may have multiple "layers" - one for each time step. Note there may
  # be only 1 time step.
  if ("time" %in% dims$name)
    dates <-
    (ncdf4::ncvar_get(x, "time") / 86400) %>%  # convert seconds → days since 1970-01-01
    lubridate::as_date(origin = lubridate::origin)
  else
    dates <- ""

  # get rasters for each date
  res <-
    purrr::map2(
      dates,
      seq_along(dates),
      create.ncdf.rast,
      dat = dat,
      datname = dataset,
      x = lon,
      y = lat,
      inproj = inproj,
      outproj = outproj,
      to = to
    )

  # If there was more than one date then stack 'em.
  # Otherwise, just peel off the single raster
  if(length(res) > 1){
    res <- terra::rast(res)
    names(res) <- paste(dataset, dates, sep = ".")
  } else {
    res <- res[[1]]
    names(res) <- dataset
  }
  res
}

#' Import and stack multiple NetCDF files as a SpatRaster
#'
#' Lists all files in \code{folder} matching \code{pattern}, calls
#' \code{\link{ncdf.to.raster}} on each, and stacks the results.
#'
#' @param folder Directory path containing the NetCDF files.
#' @param pattern Regular expression matched against filenames (e.g.
#'   \code{"^jplMURSST41mday.+nc$"}).
#' @param variable Character string variable name to extract from each file.
#' @param inproj Projection string or EPSG code for the NetCDF data.
#' @param outproj Target projection; ignored when \code{to} is supplied.
#' @param to Optional \code{SpatRaster} template for reprojection and
#'   resampling.  When both \code{outproj} and \code{to} are given,
#'   \code{outproj} is ignored with a warning.
#' @return A \code{SpatRaster} with all layers stacked (or a single
#'   \code{SpatRaster} when only one file is found).
#' @export
import.netCDF <-
  function(folder, pattern, variable, inproj, outproj = NULL, to = NULL) {

    files <- list.files(folder, pattern = pattern, full.names = T)

    if (length(files) == 0) {
      stop("import.netCDF: no matching filenames")
    } else {
      # reading one or more
      res <-
        purrr::map(
          files,
          ncdf.to.raster,
          dataset = variable,
          inproj = inproj,
          outproj = outproj,
          to = to
        )

      # If only 1 file, just return the first element of the list.
      # Note that if the single file had multiple layers (say years etc)
      # then the single element of res can still be a rasterbrick, but thats ok. I
      # just want to avoid returning a brick when a rasterlayer is expected.
      if (length(files) > 1) {
        res <- terra::rast(res)
      } else
        res <- res[[1]]

    }
    res
  }


#' Compute the mean raster for a given calendar month
#'
#' Selects all layers in \code{r} whose name encodes the given month and
#' returns their pixel-wise mean.  Used to create climatological monthly
#' averages from a multi-year monthly raster stack.
#'
#' @param mnth Integer calendar month (1 = January, …, 12 = December).
#' @param r \code{SpatRaster} with layer names containing \code{"-MM-"}
#'   date strings.
#' @return A single-layer \code{SpatRaster} representing the mean for
#'   \code{mnth}.
#' @export
rast.monthly.mean <- function(mnth, r){
  message("Getting monthly means for month ", mnth)
  r.mnths <- names(r) %>%
    stringr::str_split_fixed(stringr::fixed("."), n = Inf) %>%
    magrittr::extract(, 2) %>%
    stringr::str_split_fixed(stringr::fixed("-"), n = Inf) %>%
    magrittr::extract(, 2) %>%
    as.integer

  sel <- r[[which(r.mnths == mnth)]]

  # NB: terra::mean(), not bare mean(). mean() lives in base, so inside this
  # package it shadows terra's S4 method for SpatRaster and silently returns
  # NA ("argument is not numeric or logical") instead of the cell-wise mean
  # raster. The caller then fails with "none of the elements of x are a
  # SpatRaster". Same root cause as the terra::plot() note in create.ncdf.rast.
  terra::mean(sel)
}

#' Produce a dotchart for a single column of a data frame
#'
#' Silently does nothing if \code{varname} is not present in \code{dat}.
#'
#' @param varname Character string column name to plot.
#' @param dat Data frame containing the column.
#' @return \code{invisible(NULL)}, called for its side-effect (plot).
#' @export
do.dotchart <- function(varname, dat) {
  if (varname %in% names(dat))
    dotchart(dat[, varname],
             main = varname,
             xlab = "Values of variable",
             ylab = "Order of the data")
}

#' Compute seasonal mean for one dynamic variable across a month range
#'
#' Averages the per-row monthly columns \code{<var>.MM} (and their scaled
#' counterparts \code{<var>.MM_sc}) between \code{start} and \code{end}
#' months and returns the means as a named list.
#'
#' @param var Character string variable name (e.g. \code{"sst"}).
#' @param dat Data frame (normally a prediction grid) containing monthly
#'   value columns named \code{<var>.<month>}.
#' @param start Integer start month.
#' @param end Integer end month.
#' @return Named list with elements \code{<var>} and \code{<var>_sc}
#'   containing the row-wise means.
#' @export
get.seas.mean.var <- function(var, dat, start, end) {
  var.names <- paste(var, start:end, sep = ".")
  var.names.sc <- paste0(var.names, "_sc")
  ret = list(dat %>%
               dplyr::select(dplyr::all_of(var.names)) %>%
               rowMeans,
             dat %>%
               dplyr::select(dplyr::all_of(var.names.sc)) %>%
               rowMeans)
  names(ret) <- c(var, paste0(var, "_sc"))
  ret
}


#' Compute seasonal means for all dynamic variables for one season
#'
#' Filters \code{dat} to rows matching \code{seas}, then calls
#' \code{\link{get.seas.mean.var}} for each variable in \code{dyn.vars} to
#' compute the mean over the months defined in \code{season.spec}.
#'
#' @param seas Character string season label.
#' @param season.spec Named list of season boundary definitions for the
#'   species (from the project \code{seasons} list).
#' @param dat Data frame (normally a prediction grid) with monthly covariate
#'   columns and a \code{Season} column.
#' @param dyn.vars Character vector of dynamic variable names (e.g.
#'   \code{c("sst", "sst.g")}).
#' @return \code{dat} (filtered to \code{seas}) with seasonal mean columns
#'   appended.
#' @export
get.seas.mean <- function(seas, season.spec, dat, dyn.vars) {
  start <- season.spec[[seas]]["from"] %/% 100
  end <- season.spec[[seas]]["to"] %/% 100
  dat <- dplyr::filter(dat, as.character(Season) == seas)
  new.cols <- purrr::map(dyn.vars, get.seas.mean.var, dat = dat, start = start, end = end)
  cbind(dat, new.cols)
}


#' Copy prediction HTML reports for all species to a folder
#'
#' \strong{Note: may no longer be needed} now that results HTMLs are tracked
#' in git LFS.  Iterates over \code{spec.grps} (project global) and copies
#' each species prediction report to \code{folder}.
#'
#' @param folder Destination directory path (created if necessary).
#' @return \code{invisible(NULL)}, called for its side-effect (files copied).
#' @export
save.prediction.htmls <- function(folder) {
  if (!dir.exists(folder))
    dir.create(folder, recursive = TRUE)

  names(spec.grps) %>%
    purrr::map(function(species) {
      filename <- file.path(ResultsDir, species, paste0(species, "_3_prediction.html"))
      message(sprintf("Copying '%' in '%s'", filename, folder))
      file.copy(filename, folder)
    })
}


#' Reclassify, reproject, mask, and save a SpatRaster
#'
#' Passes \code{r} through \code{terra::classify()}, reprojects and resamples
#' to match \code{to}, masks to \code{to}, and writes the result to
#' \code{filename}.
#'
#' @param r Input \code{SpatRaster}.
#' @param class.arg Reclassification argument passed to
#'   \code{\link[terra]{classify}}: a 3-column (from, to, becomes), 2-column
#'   (is, becomes), or 1-column cut-point matrix.  Commonly
#'   \code{cbind(Inf, NA)} to convert infinite values to \code{NA}.
#' @param to \code{SpatRaster} template used for reprojection, resampling,
#'   and masking.
#' @param filename Output file path (overwritten if it exists).
#' @return \code{invisible(NULL)}, called for its side-effect (file written).
#' @export
reclassify.project.save <- function(r, class.arg, to, filename){
  outdir <- unique(dirname(filename))

  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE )

  r %>%
    terra::classify(rcl = class.arg) %>%
    terra::project(y = to,
                   method = "bilinear",
                   threads = TRUE) %>%
    terra::mask(to) %>%
    terra::writeRaster(filename = filename,
                       overwrite = TRUE)
}

#' Return names of dynamic environmental covariates (and their gradients)
#'
#' Reads the project-global \code{env_covar_spec} table and returns the names
#' of all dynamic variables, appending \code{".g"} suffixes for any that have
#' gradients enabled.
#'
#' @return Character vector of dynamic covariate names.
#' @export
dynamic.env.covar.names <- function(){
  dyn_vars <- dplyr::filter(env_covar_spec, var_type == "dynamic")
  grads <- dyn_vars$var_name[dyn_vars$do_gradient]
  if (length(grads) > 0)
    c(dyn_vars$var_name, paste0(grads, ".g"))
  else
    dyn_vars$var_name
}

#' Build a species-specific seasonal prediction grid
#'
#' Creates one seasonal copy of \code{predgrid} per entry in
#' \code{season.names}, computes seasonal mean values for all dynamic
#' covariates, drops the monthly columns, and then replicates the grid once per
#' level of the \code{platform} factor.  The combined grid is saved as a
#' shapefile.
#'
#' Expects project globals \code{seasons}, \code{season.names},
#' \code{ddftype_to_platform}, and \code{ShapeDir} in the calling
#' environment.
#'
#' @section Getting the platform levels right:
#' The number of prediction-grid copies must equal the number of segment-table
#' copies per segment, because [dsm.pred()] sums across them to get the combined
#' surface. Miller et al. (2021) hold that by construction - `segs2 <-
#' rbind(segs, segs)` and `pred2 <- rbind(pred, pred)` in their fulmar example -
#' and here it holds as long as each platform level is carried by exactly one
#' ddftype.
#'
#' The levels used to default to `unique(ddftype_to_platform)`, the declared
#' global. That is wrong whenever a species' fitted model has fewer levels than
#' the global declares - which happens routinely now that [create.dsm.data()]
#' drops ddftypes with no observations. Predicting on a grid with more levels
#' than the model has fails inside `mgcv` with "factor has new levels";
#' predicting with FEWER fails silently, summing a subset of the components and
#' under-reporting density with no error at all. So pass `platform.levels` from
#' the fitted model (`model$xlevels$platform`) rather than relying on the
#' default.
#'
#' The two side-effect arguments exist for consumers that want the seasonal
#' covariate values without the products built for prediction.
#' [assess.extrapolation()] is one: it needs the covariates the predictions
#' were made on, but must not rewrite the shapefile that
#' `Generic_3_prediction.Rmd` owns, and has no use for the platform copies
#' because `platform` is a model term rather than an environmental covariate.
#' Both default to `TRUE`, so the prediction path is unchanged.
#'
#' @param species Character string species code.
#' @param predgrid \code{sf} data frame representing the prediction grid (one
#'   row per cell, with monthly covariate columns).
#' @param write.shapefile Logical. Write `[species]_predgrid.shp` to
#'   \code{ShapeDir}? Defaults to \code{TRUE}.
#' @param replicate.platform Logical. Replicate the grid once per platform
#'   level, adding a \code{platform} column? Defaults to \code{TRUE}.
#' @param platform.levels Character vector of platform levels to replicate over,
#'   normally \code{model$xlevels$platform} for the model being predicted.
#'   Defaults to \code{unique(ddftype_to_platform)}; see *Getting the platform
#'   levels right*.
#' @return \code{sf} data frame with one row per (cell × season × platform)
#'   combination, containing seasonal mean covariates, with \code{platform} a
#'   factor whose levels are \code{platform.levels} in the given order. With
#'   \code{replicate.platform = FALSE}, one row per (cell × season) and no
#'   \code{platform} column.
#' @export
create.seasonal.predgrid <- function(species, predgrid,
                                     write.shapefile = TRUE,
                                     replicate.platform = TRUE,
                                     platform.levels = unique(ddftype_to_platform)) {
  checkmate::expect_string(species)
  checkmate::expect_class(predgrid, "sf")
  checkmate::expect_flag(write.shapefile)
  checkmate::expect_flag(replicate.platform)
  checkmate::expect_character(platform.levels, min.len = 1, any.missing = FALSE,
                              unique = TRUE)

  # Get species-specific season setting and create one predgrid copy per season.
  # TODO: what if there's no data for some seasons? Should there still be
  # a copy of every season?
  season.spec <- seasons[[species]]
  # do.call(rbind, ...) rather than purrr::list_rbind(): rbind dispatches to
  # rbind.sf and keeps the geometry column an sfc, which vctrs-based binding does
  # not guarantee. This is the literal generalisation of the
  # rbind(predgrid, predgrid, predgrid, predgrid) it replaces.
  ret <-
    do.call(rbind, replicate(length(season.names), predgrid, simplify = FALSE)) %>%
    dplyr::mutate(Season = as.factor(rep(season.names, each = nrow(predgrid))))

  # Get seasonal means for dynamic variables
  match.re <- c("[0-9]$", "[0-9]_sc$")
  p.geom <- sf::st_geometry(ret) # save geometry
  ret <- season.names %>%
    purrr::map_dfr(
      get.seas.mean,
      season.spec = season.spec,
      dat = sf::st_drop_geometry(ret),
      dyn.vars = dynamic.env.covar.names()
    ) %>%
    # Remove monthly values from predgrid to make it smaller now that we have
    # seasonal means computed, add area, and make back into sf object.
    dplyr::select(!dplyr::matches(match.re)) %>%
    cbind(p.geom) %>% # add geometry back in
    sf::st_sf()

  # Rename scaled columns to match what was in the model specs. This is harmless
  # if they are already named correctly with .sc.
  nms <- gsub( "_sc", ".sc", names(ret), fixed = TRUE)
  names(ret) <- nms

  # Save the seasonal prediction grid for GIS mapping. When doing multiple
  # species, there will be a separate predgrid for each species since it's
  # possible for the season boundaries to be species specific.
  if (write.shapefile)
    sf::st_write(ret,
                 dsn = ShapeDir,
                 layer = paste0(species, "_predgrid.shp"),
                 driver = "ESRI Shapefile",
                 delete_layer = TRUE
    )

  # Make a copy for each platform level in the model. Ordered platform-major,
  # season-major within platform, cell within season - dsm.pred() relies on that
  # ordering when it peels off one platform block as the template for the
  # combined surface.
  #
  # platform is made a factor with exactly platform.levels, in that order, so
  # predict() cannot silently coerce a character column against the model's own
  # level ordering.
  if (replicate.platform)
    ret <- replicate(length(platform.levels), ret, simplify = FALSE) %>%
      setNames(platform.levels) %>%
      purrr::list_rbind(names_to = "platform") %>%
      dplyr::mutate(platform = factor(platform, levels = platform.levels)) %>%
      sf::st_sf()


  # From multiddf paper code:
  # create an extra column to account for the variance propagation model
  # the variance propagation adds a random effect covariate named "XX"
  # which we can safely give the value 0 now the variance has been propagated.
  # pred$XX <- matrix(0, nrow=nrow(pred), ncol=3)
  ret
}

#' Produce leaflet prediction maps for all four seasons
#'
#' Iterates over \code{season.names} (project global), calling
#' \code{\link{do.pred.map}} for each season, and returns the resulting maps
#' as a named list.
#'
#' @param dat \code{sf} data frame of prediction grid cells with columns
#'   \code{subset}, \code{Season}, \code{Dens}, and geometry.
#' @param model Fitted \code{dsm} object; its \code{$data} element is used to
#'   add estimated abundance circles.
#' @param modname Character string model name; used in messages.
#' @param species Character string species code; used in messages and legend
#'   titles.
#' @param subs Which platform subset to map: \code{"Combined"} (the default and
#'   the usual choice), or any single platform level present in \code{dat}.
#'   The available levels come from the data rather than a fixed list, because
#'   they depend on \code{ddftype_to_platform} and on which ddftypes the species
#'   actually has - the old hardcoded \code{c("Combined", "F", "W")} silently
#'   stopped matching when either changed.
#' @param ... Additional arguments passed to \code{\link{do.pred.map}}.
#' @return Named list of leaflet map objects, one per season.
#' @export
do.pred.maps <-
  function(dat,
           model,
           modname,
           species,
           subs = "Combined",
           ...) {

    checkmate::expect_string(subs)
    available <- unique(as.character(dat$subset))
    if (!subs %in% available)
      stop(sprintf(
        "do.pred.maps: subset %s is not in these predictions. Available: %s.",
        sQuote(subs), paste(sQuote(available), collapse = ", ")))

    # Get data subset
    dat <- dplyr::filter(dat, subset == subs)
    segdata <- model$data %>%
      sf::st_as_sf()

    message(sprintf("%s, %s: Doing %s abundance prediction map for",
                    species, modname, subs))


    ret <- season.names %>%
      purrr::map(do.pred.map, dat, segdata, modname, species, subs, ...)

    names(ret) <- season.names
    ret
  }

#' Produce a single-season leaflet density prediction map
#'
#' Filters \code{dat} and \code{segdata} to \code{season}, simplifies
#' geometries, and builds a leaflet map with a filled-polygon density layer
#' and proportional-circle estimated abundance layer.
#'
#' @param season Character string season label.
#' @param dat \code{sf} prediction grid polygons with \code{Season} and
#'   \code{Dens} columns.
#' @param segdata \code{sf} segment data with \code{Season}, \code{LongStart},
#'   \code{LatStart}, \code{estAbund}, and \code{rawCount} columns.
#' @param modname Character string model name; used in legend.
#' @param species Character string species code; used in legend.
#' @param subset Character string platform subset label
#'   (\code{"Combined"}, \code{"F"}, or \code{"W"}).
#' @param samp_n Integer; if not \code{NA}, draw a random sample of this
#'   many prediction polygons (for performance testing).
#' @return A \code{leaflet} map object.
#' @export
do.pred.map <-
  function(season,
           dat,
           segdata,
           modname,
           species,
           subset,
           samp_n = NA) {


    message(sprintf("\t%s",season))

    # Filter by season, and create log Density for potential mapping - not currently
    # used.
    dat <- dat %>%
      dplyr::filter(Season == season) %>%
      sf::st_transform(latlongproj) %>%
      dplyr::mutate(lDens = dplyr::case_when(Dens == 0 ~ 0,
                                             TRUE ~ log(Dens))) %>%
      dplyr::select(Dens, geometry) %>%
      rmapshaper::ms_simplify()

    segdata <- segdata %>%
      dplyr::filter(Season == season) %>%
      get.combined.segdata() %>%
      sf::st_transform(latlongproj)

    # Plot only a sample of the polygons for efficiency? Typically used for
    # testing.
    if (!is.na(samp_n)) {
      index <- sample(1:nrow(dat), size = samp_n)
      dat <- dat[index,]
    }

    # Remove ridiculously large densities b/c they mess up the legend and swamp
    # everything else
    dat <- dplyr::mutate(dat,
                         Dens = dplyr::case_when(Dens > MAX_DENS_VALUE ~ NA,
                                                 TRUE ~ Dens))

    if ((n.na <- sum(is.na(dat$Dens))) > 0)
      message("Warning: ", n.na, " cells larger than ", MAX_DENS_VALUE, " were converted to NA")

    groups <- c("est abund", "Pred Dens")
    m <-
      leaflet::leaflet(
        data = dat,
        options = leaflet::leafletOptions(preferCanvas = TRUE)
      ) %>%
      # Options help to speed up rendering. NOTE - dont't use addProviderTiles
      # if you want to save the map and reload in a subsequent R session - it won't
      # work.
      leaflet::addTiles(options = leaflet::tileOptions(updateWhenZooming = FALSE,
                                                       updateWhenIdle = FALSE)) %>%
      leaflet::addMapPane("density", zIndex = 410) %>%
      leaflet::addMapPane("abund", zIndex = 420) %>%
      # Predicted density
      leaflet::addPolygons(
        fillColor = ~ pal_pred(Dens),
        color = ~ pal_pred(Dens),
        fillOpacity = 1.0,
        opacity = 1.0,
        weight = 1,
        group = "Pred Dens",
        options = leaflet::pathOptions(pane = "density")
      ) %>%
      # 0 Abund
      leaflet::addCircles(
        lng = ~ LongStart,
        lat = ~ LatStart,
        stroke = FALSE,
        radius = rep(1000, times = nrow(dplyr::filter(segdata, estAbund == 0))),
        color = "white",
        fillColor = "white",
        fillOpacity = 0.1,
        opacity = 0.9,
        group = "est abund",
        data = dplyr::filter(segdata, estAbund == 0),
        options = leaflet::pathOptions(pane = "abund")
      ) %>%
      # Est Abund
      leaflet::addCircles(
        lng = ~ LongStart,
        lat = ~ LatStart,
        radius = ~ estAbund * (max.circ.radius/max(estAbund)),#  scale so largest is 100km
        color = "black",
        weight = 1,
        popup = ~ htmlEscape(paste0("Abund ", round(estAbund, 3), ", raw ", round(rawCount, 3))),
        group = "est abund",
        data = dplyr::filter(segdata, estAbund != 0),
        options = leaflet::pathOptions(pane = "abund")
      ) %>%
      leaflet::addLegend(
        pal = pal_pred,
        opacity = 1,
        values = ~ Dens,
        # values = class_intervals$brks, # only for discrete color scale
        title = paste(species, season)
      ) %>%
      leafem::addMouseCoordinates() %>%
      leaflet::addScaleBar(position = "bottomright", options = leaflet::scaleBarOptions(imperial = FALSE)) %>%
      leaflet::addLayersControl(overlayGroups = groups,
                                options = leaflet::layersControlOptions(collapsed = FALSE)) %>%
      leaflet::hideGroup(c("est abund"))

    m
  }

#' Collapse per-platform segdata copies into a single combined segdata
#'
#' Species-specific segdata contains multiple copies of each segment (one per
#' ddftype).  This function strips the ddftype suffix from \code{Sample.Label}
#' and sums \code{estDens}, \code{estAbund}, and \code{rawCount} across the
#' copies, returning one row per original segment.
#'
#' @param segdata \code{sf} segment data frame in the format returned by
#'   \code{\link{create.dsm.data}}, with \code{Sample.Label} suffixes.
#' @return \code{sf} data frame with one row per segment and summed abundance
#'   estimates.
#' @export
get.combined.segdata <- function(segdata){
  # Note that the species-specific segdata has multiple copies for each original
  # segment - currently up to 4 per aerial and 4 per ship-based segments: one
  # for all combinations of [A]erial/[S]hip, W[ater]/F[lying], and
  # [L]ine/S[trip] transect with Sample.Label suffixes as specified by
  # ddftype_levels defined in analysis_settings.R, so we first need to chop these
  # suffixes to do the filter and then sum the estDen[sities].
  segdata <-
    dplyr::select(
      segdata,
      Sample.Label,
      estDens,
      estAbund,
      rawCount,
      SurveyType,
      x,
      y,
      LatStart,
      LongStart,
      geometry,
      Season,
      Date,
      platform
    ) %>%
    # Remove Sample_Label suffix (_x_x)
    dplyr::mutate(Sample.Label = stringr::str_sub(Sample.Label, 1, nchar(Sample.Label) -
                                                    nchar(ddftype_levels[1]) - 1)) %>%
    dplyr::arrange(Sample.Label)

  # There are now multiple consecutive rows for each segment.

  # Figure out how many copies of each segment there are. Complain if not all the
  # same. This relies on having a fixed symmetric setup where, for example,
  # the aerial segments are copied the same number of times as the ship ones
  # (ie once each for Water, Fly, and Strip)
  rl <- rle(segdata$Sample.Label)
  stopifnot(length(unique(rl$lengths)) == 1)
  ncopies <- rl$lengths[1]

  # Keep every nth row, summing values within each group of ncopies segdata rows,
  # and make it a column in  new data frame containing only every nth segdata row
  # in order to get one row per segment. Note that rollapply() uses the full
  # segdata (before selecting every nth row)

  keep <- rep(c(TRUE, rep(FALSE, times = ncopies - 1)),
              times = nrow(segdata) / ncopies)
  res <- segdata[keep,] %>%
    dplyr::mutate(
      estDens = zoo::rollapply(segdata$estDens, ncopies, by = ncopies, sum),
      estAbund = zoo::rollapply(segdata$estAbund, ncopies, by = ncopies, sum),
      rawCount = zoo::rollapply(segdata$rawCount, ncopies, by = ncopies, sum)
    ) %>%
    dplyr::select(-platform) # No longer makes any sense since it will have value of first row in group

  res
}

#' Export seasonal segdata and distdata shapefiles for a species
#'
#' Collapses per-platform segdata copies with
#' \code{\link{get.combined.segdata}}, optionally subsets to
#' \code{sample.labs}, and writes one shapefile per season plus a single
#' distdata shapefile.
#'
#' @param spec Character string species code.
#' @param sample.labs Optional character vector of \code{Sample.Label} values
#'   to restrict the output to.
#' @param folder Directory path for output shapefiles.
#' @param segdata \code{sf} segment data frame.
#' @param distdata Data frame of observation data.
#' @return \code{invisible(NULL)}, called for its side-effect (files written).
#' @export
create.species.shapefiles <-
  function(spec,
           sample.labs = NULL,
           folder,
           segdata,
           distdata) {

    # Get segdata for spec with fly/water combined.
    segdata <- get.combined.segdata(segdata)

    # Restrict to certain sample labels (useful for cropping to only those samples
    # in a certain spatial area)
    if (!is.null(sample.labs))
      segdata <- dplyr::filter(segdata, Sample.Label %in% sample.labs)

    # save as seasonal shapefiles if required
    for (seas in season.names) {
      layer <- paste(spec, seas, "segdata", sep = "_")

      # if shapefile doesn't exist or recreateSpecSegShapefiles is TRUE then
      # save the shapefile
      if (!file.exists(paste0(folder, "/", layer, ".shp")) ||
          recreateSpecSegShapefiles) {
        message(sprintf("Creating segdata shapefile for %s %s", spec, seas))

        dplyr::filter(segdata, Season == seas) %>%
          sf::st_write(
            dsn = folder,
            layer = layer,
            driver = "ESRI Shapefile",
            delete_layer = T
          )
      }
    }

    #### Now do same for distdata.
    # Distdata isn't seasonal b/c it is used in its entirety for the ddf (but
    # with season as covar if needed)
    layer <- paste(spec, "distdata", sep = "_")

    if (!file.exists(paste0(folder, "/", layer, ".shp")) ||
        recreateDistdataShapefile) {
      message(sprintf("Creating distdata shapefile for %s", spec))

      # Save distdata as a shapefile.
      distdata %>%
        # Needed since numeric ids can get too big for shapefile numbers
        dplyr::mutate(object = as.character(object)) %>%
        sf::st_as_sf(
          coords = c("LongStart", "LatStart"),
          crs = sf::st_crs(4326),
          remove = FALSE
        ) %>%
        sf::st_write(
          dsn = ShapeDir,
          layer = layer,
          driver = "ESRI Shapefile",
          delete_layer = TRUE
        )
    }
  }

#' Return the filename for the final model prediction raster
#'
#' Constructs the canonical filename for the GeoTIFF storing density
#' predictions for \code{spec} and \code{season}.
#'
#' @param spec Character string species code.
#' @param season Character string season label.
#' @return Character string filename (without directory path).
#' @export
get.final.prediction.name <- function(spec, season){

  # Get final model predictions raster filename
  modname <- final.dsm.models$dsm_final_name[final.dsm.models$species == spec]
  filename <- sprintf("%s.%s.%s.%d_sqkm.tif",
                      spec,
                      season,
                      modname,
                      predgridCellArea)
  filename
}

#' Return the filename for the final model CV raster
#'
#' Constructs the canonical filename for the GeoTIFF storing the coefficient
#' of variation of predictions for \code{spec} and \code{season}.
#'
#' @param spec Character string species code.
#' @param season Character string season label.
#' @return Character string filename (without directory path).
#' @export
get.final.variance.name <- function(spec, season){

  # Get final model variance raster filename
  modname <- final.dsm.models$dsm_final_name[final.dsm.models$species == spec]
  filename <- sprintf("%s.%s.%s.%d_sqkm_CV.tif",
                      spec,
                      season,
                      modname,
                      predgridCellArea)
  filename
}

# Parse a prediction raster filename into its components.
# Returns NULL if the name does not match the convention.
.parse_pred_filename <- function(filepath) {
  base  <- tools::file_path_sans_ext(basename(filepath))
  parts <- regmatches(
    base,
    regexec(
      "^([^.]+)\\.([^.]+)\\.(.+)\\.(\\d+_sqkm(?:_CV)?)$",
      base, perl = TRUE
    )
  )[[1]]
  if (length(parts) == 0) return(NULL)
  list(
    path     = filepath,
    species  = parts[2],
    season   = parts[3],
    model    = parts[4],
    cellsize = sub("_CV$", "", parts[5]),
    cv       = grepl("_CV$", parts[5])
  )
}

# Discover prediction files in folder matching the given filters.
# NULL filter values mean "no filter" (include all).
.discover_pred_files <- function(folder, species_filter, model_filter,
                                  cellsize_filter, season_filter, cv) {
  files   <- list.files(folder, full.names = TRUE)
  parsed  <- lapply(files, .parse_pred_filename)
  records <- Filter(Negate(is.null), parsed)
  records <- Filter(function(r) r$cv == cv, records)
  if (!is.null(species_filter))
    records <- Filter(function(r) r$species == species_filter, records)
  if (!is.null(model_filter))
    records <- Filter(function(r) r$model == model_filter, records)
  if (!is.null(cellsize_filter))
    records <- Filter(function(r) r$cellsize == cellsize_filter, records)
  if (!is.null(season_filter))
    records <- Filter(function(r) r$season %in% season_filter, records)
  records
}

# Build a leaflet map of a difference raster using a diverging blue-white-red
# palette centred on zero (blue = r2 > r1, red = r1 > r2).
.compare_pred_map <- function(diff_r, pair, compare_cv) {
  d_vals  <- as.vector(terra::values(diff_r, na.rm = TRUE))
  max_abs <- max(abs(d_vals))
  if (max_abs == 0) max_abs <- 1

  pal <- leaflet::colorNumeric(
    palette  = c("blue", "white", "red"),
    domain   = c(-max_abs, max_abs),
    na.color = "transparent"
  )

  r_wgs84 <- terra::project(diff_r, "EPSG:4326")
  r_rast  <- raster::raster(r_wgs84)

  title <- sprintf(
    "%s %s %s<br>vs %s %s %s%s",
    pair$species1, pair$season1, pair$model1,
    pair$species2, pair$season2, pair$model2,
    if (compare_cv) " (CV)" else ""
  )

  leaflet::leaflet() %>%
    leaflet::addTiles(
      options = leaflet::tileOptions(updateWhenZooming = FALSE,
                                     updateWhenIdle    = FALSE)
    ) %>%
    leaflet::addRasterImage(r_rast, colors = pal, opacity = 0.8) %>%
    leaflet::addLegend(
      pal     = pal,
      values  = raster::values(r_rast),
      title   = title,
      opacity = 1
    ) %>%
    leaflet::addScaleBar(
      position = "bottomright",
      options  = leaflet::scaleBarOptions(imperial = FALSE)
    )
}

# Find the first file in folder whose parsed components match exactly.
# Returns the file path, or NULL if not found.
.find_matching_file <- function(folder, species, season, model,
                                 cellsize, cv) {
  files <- list.files(folder, full.names = TRUE)
  for (f in files) {
    p <- .parse_pred_filename(f)
    if (is.null(p)) next
    if (p$species == species && p$season == season &&
        p$model == model && p$cellsize == cellsize && p$cv == cv)
      return(f)
  }
  NULL
}

#' Compare two sets of density surface model prediction rasters
#'
#' Finds matching raster files in two folders by parsing the standard DSMHelper
#' filename convention (\code{species.season.model.cellsize[_CV].ext}),
#' computes cell-by-cell differences, and saves a GeoTIFF difference raster
#' for each matched pair. Summary statistics are printed and returned.
#'
#' Two operating modes are available:
#' \itemize{
#'   \item \strong{Seasonal mode} (default, \code{season1 = NULL}): all
#'     seasons discovered in \code{folder1} are compared against the
#'     corresponding season in \code{folder2}. The \code{seasons} argument
#'     restricts which season labels are included.
#'   \item \strong{Single mode} (\code{season1} provided): exactly one pair
#'     is compared. \code{season2} defaults to \code{season1} when omitted,
#'     allowing same-season cross-model or cross-species comparisons.
#' }
#'
#' Both rasters in each pair must be geometrically identical (same extent,
#' resolution, and CRS) and must have \code{NA} in exactly the same cells.
#' Differences are computed as \code{r1 - r2}; positive values indicate that
#' the first raster predicts higher densities.
#'
#' Output difference rasters are always written as GeoTIFF regardless of the
#' input format.
#'
#' @param folder1 Character string. Directory containing the first raster set.
#' @param folder2 Character string. Directory containing the second raster
#'   set.
#' @param output_dir Character string or \code{NULL}. Directory for saved
#'   difference rasters; created if it does not exist. When \code{NULL}
#'   (default) difference rasters are not written to disk.
#' @param species1 Character string or \code{NULL}. Species code to match in
#'   \code{folder1} (e.g. \code{"NOGA"}). \code{NULL} includes all species.
#' @param model1 Character string or \code{NULL}. Model name to match in
#'   \code{folder1}. \code{NULL} includes all models.
#' @param species2 Character string or \code{NULL}. Species code to look up
#'   in \code{folder2}. \code{NULL} uses the same species as the matched
#'   \code{folder1} file.
#' @param model2 Character string or \code{NULL}. Model name to look up in
#'   \code{folder2}. \code{NULL} uses the same model as the matched
#'   \code{folder1} file.
#' @param cellsize Character string or \code{NULL}. Cell-size label (e.g.
#'   \code{"100_sqkm"}) to filter on. \code{NULL} includes all cell sizes
#'   found in \code{folder1}.
#' @param season1 Character string or \code{NULL}. Season label for the first
#'   raster in single mode. \code{NULL} (default) activates seasonal mode.
#' @param season2 Character string or \code{NULL}. Season for the second
#'   raster in single mode. Defaults to \code{season1} when \code{NULL}.
#'   Ignored in seasonal mode.
#' @param seasons Character vector or \code{NULL}. In seasonal mode, restrict
#'   comparisons to these season labels. \code{NULL} includes all seasons
#'   discovered in \code{folder1}.
#' @param compare_cv Logical. If \code{TRUE}, compare coefficient-of-variation
#'   rasters (\code{*_CV} files) instead of density rasters. Default
#'   \code{FALSE}.
#' @param plot Logical. If \code{TRUE}, produce an interactive leaflet map for
#'   each comparison using a diverging blue–white–red palette centred on zero
#'   (blue = second raster higher, red = first raster higher). Default
#'   \code{FALSE}.
#' @return A named list returned invisibly:
#'   \describe{
#'     \item{\code{$stats}}{Data frame with one row per comparison. Columns:
#'       \code{species1}, \code{season1}, \code{model1}, \code{species2},
#'       \code{season2}, \code{model2}, \code{mean_r1}, \code{min_r1},
#'       \code{max_r1}, \code{mean_r2}, \code{min_r2}, \code{max_r2},
#'       \code{mean_diff}, \code{min_diff}, \code{max_diff},
#'       \code{median_diff}, \code{sd_diff}, \code{pearson_r}.}
#'     \item{\code{$files}}{Character vector of saved difference raster
#'       paths.}
#'     \item{\code{$maps}}{Named list of leaflet map objects, one per
#'       comparison (empty when \code{plot = FALSE}).}
#'   }
#' @examples
#' \dontrun{
#' # Simplest case: compare all matching files across two folders
#' res <- compare_predictions("path/to/v1", "path/to/v2", "path/to/diffs")
#'
#' # Seasonal: compare all seasons for one species, two models
#' res <- compare_predictions(
#'   folder1    = "path/to/pred",
#'   folder2    = "path/to/pred",
#'   output_dir = "path/to/diffs",
#'   species1   = "NOGA",
#'   model1     = "dsm_nb_allpred",
#'   model2     = "dsm_nb_depth_only"
#' )
#'
#' # Single: compare Winter only
#' res <- compare_predictions(
#'   folder1    = "path/to/pred",
#'   folder2    = "path/to/pred",
#'   output_dir = "path/to/diffs",
#'   species1   = "NOGA",
#'   model1     = "dsm_nb_allpred",
#'   model2     = "dsm_nb_depth_only",
#'   season1    = "Winter"
#' )
#' }
#' @export
compare_predictions <- function(
  folder1,
  folder2,
  output_dir = NULL,
  species1   = NULL,
  model1     = NULL,
  species2   = NULL,
  model2     = NULL,
  cellsize   = NULL,
  season1    = NULL,
  season2    = NULL,
  seasons    = NULL,
  compare_cv = FALSE,
  plot       = FALSE
) {
  checkmate::expect_string(folder1, min.chars = 1)
  checkmate::expect_string(folder2, min.chars = 1)
  checkmate::expect_directory_exists(folder1)
  checkmate::expect_directory_exists(folder2)
  checkmate::expect_flag(compare_cv)
  if (!is.null(species1)) checkmate::expect_string(species1, min.chars = 1)
  if (!is.null(model1))   checkmate::expect_string(model1,   min.chars = 1)
  if (!is.null(species2)) checkmate::expect_string(species2, min.chars = 1)
  if (!is.null(model2))   checkmate::expect_string(model2,   min.chars = 1)
  if (!is.null(cellsize)) checkmate::expect_string(cellsize, min.chars = 1)
  if (!is.null(season1))  checkmate::expect_string(season1,  min.chars = 1)
  if (!is.null(season2))  checkmate::expect_string(season2,  min.chars = 1)
  if (!is.null(seasons))  checkmate::expect_character(seasons, min.len = 1)
  if (!is.null(output_dir))
    checkmate::expect_string(output_dir, min.chars = 1)
  checkmate::expect_flag(plot)

  seasonal_mode <- is.null(season1)
  season_filter <- if (seasonal_mode) seasons else season1

  candidates <- .discover_pred_files(
    folder          = folder1,
    species_filter  = species1,
    model_filter    = model1,
    cellsize_filter = cellsize,
    season_filter   = season_filter,
    cv              = compare_cv
  )

  if (length(candidates) == 0)
    stop("No matching prediction files found in folder1.")

  pairs <- list()
  for (cand in candidates) {
    s2   <- if (!is.null(species2)) species2 else cand$species
    m2   <- if (!is.null(model2))   model2   else cand$model
    sea2 <- if (!seasonal_mode) {
              if (!is.null(season2)) season2 else season1
            } else {
              cand$season
            }

    path2 <- .find_matching_file(
      folder   = folder2,
      species  = s2,
      season   = sea2,
      model    = m2,
      cellsize = cand$cellsize,
      cv       = compare_cv
    )

    if (is.null(path2)) {
      warning(sprintf(
        "No match in folder2 for: species=%s season=%s model=%s cellsize=%s",
        s2, sea2, m2, cand$cellsize
      ))
      next
    }

    pairs <- c(pairs, list(list(
      path1    = cand$path,
      path2    = path2,
      species1 = cand$species,
      season1  = cand$season,
      model1   = cand$model,
      cellsize = cand$cellsize,
      species2 = s2,
      season2  = sea2,
      model2   = m2
    )))
  }

  if (length(pairs) == 0)
    stop("No matching file pairs found between folder1 and folder2.")

  if (!is.null(output_dir)) create.dir.if.needed(output_dir)

  timestamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")
  stats_rows  <- list()
  saved_files <- character(0)
  maps        <- list()

  for (pair in pairs) {
    message(sprintf("Comparing:\n  %s\n  %s", pair$path1, pair$path2))

    r1 <- terra::rast(pair$path1)
    r2 <- terra::rast(pair$path2)

    terra::compareGeom(r1, r2, stopOnError = TRUE)

    na1 <- is.na(terra::values(r1))
    na2 <- is.na(terra::values(r2))
    if (!all(na1 == na2))
      stop(sprintf("NA masks differ between\n  %s\n  %s",
                   pair$path1, pair$path2))

    v1 <- as.vector(terra::values(r1, na.rm = TRUE))
    v2 <- as.vector(terra::values(r2, na.rm = TRUE))
    d  <- v1 - v2

    stats_rows <- c(stats_rows, list(data.frame(
      species1    = pair$species1,
      season1     = pair$season1,
      model1      = pair$model1,
      species2    = pair$species2,
      season2     = pair$season2,
      model2      = pair$model2,
      mean_r1     = mean(v1),
      min_r1      = min(v1),
      max_r1      = max(v1),
      mean_r2     = mean(v2),
      min_r2      = min(v2),
      max_r2      = max(v2),
      mean_diff   = mean(d),
      min_diff    = min(d),
      max_diff    = max(d),
      median_diff = median(d),
      sd_diff     = sd(d),
      pearson_r   = cor(v1, v2, method = "pearson")
    )))

    if (!is.null(output_dir) || plot) {
      diff_r <- r1 - r2
    }

    if (!is.null(output_dir)) {
      diff_path <- file.path(output_dir, sprintf(
        "diff.%s.%s.%s.vs.%s.%s.%s.%s.%s.tif",
        pair$species1, pair$season1, pair$model1,
        pair$species2, pair$season2, pair$model2,
        pair$cellsize, timestamp
      ))
      terra::writeRaster(diff_r, diff_path, datatype = "FLT4S",
                         overwrite = TRUE)
      message(sprintf("Saved: %s", diff_path))
      saved_files <- c(saved_files, diff_path)
    }

    if (plot) {
      map_key       <- sprintf(
        "%s.%s.%s.vs.%s.%s.%s",
        pair$species1, pair$season1, pair$model1,
        pair$species2, pair$season2, pair$model2
      )
      maps[[map_key]] <- .compare_pred_map(diff_r, pair, compare_cv)
    }
  }

  stats_df <- dplyr::bind_rows(stats_rows)
  message("\nComparison summary:")
  message(paste(capture.output(print(stats_df)), collapse = "\n"))

  invisible(list(stats = stats_df, files = saved_files, maps = maps))
}

#' How much of a predicted total comes from a handful of cells?
#'
#' A density surface can be finite, plausible cell by cell, and still have its
#' seasonal total decided by one grid cell. That is not a numerical failure and
#' nothing else in the pipeline flags it, because every value is well formed --
#' it is a property of the data. Highly zero-inflated effort with rare enormous
#' counts (a gull flock behind a trawler, say) gives a model one huge
#' observation to honour, and honouring it puts most of the predicted abundance
#' in a single place.
#'
#' Measured on \code{Atl IMRP}: the median species-season has 2.25 per cent of
#' its total in its ten largest cells, so the check is quiet almost everywhere.
#' Herring Gull in spring has \strong{94.5 per cent} -- 107,824 segments, only
#' 1.4 per cent of them holding any bird at all, and one segment holding 3,106
#' gulls. Its top cells run 44800, 637, 44.1, 32.3. Drop that one flock and the
#' seasonal total falls from about 48,000 birds to about 2,600.
#'
#' Concentration is also where the response family stops being a detail: where
#' one observation dominates, Tweedie and negative binomial disagree about how
#' likely a huge count is, and the totals diverge accordingly. The four
#' \code{Atl IMRP} species-seasons whose families disagree by more than 5x are
#' the four most concentrated.
#'
#' @param paths Character vector of raster files, one per species-season.
#' @param species,season Labels the same length as \code{paths}.
#' @param top_n How many of the largest cells to accumulate. Default 10.
#' @param map_limit Optional density above which the mapping code hides a cell,
#'   i.e. \code{MAX_DENS_VALUE}. Worth passing, because that limit is applied
#'   only when drawing leaflet maps -- "b/c they mess up the legend and swamp
#'   everything else" -- and never to the rasters that get copied to the
#'   versioned folder and shared. So the cells most likely to dominate a total
#'   are exactly the ones absent from the map you would check it against.
#' @return A tibble with one row per raster: \code{species}, \code{season},
#'   \code{cells} (finite cells), \code{total}, \code{max_cell}, and
#'   \code{pct_top_n}, sorted with the most concentrated first. With
#'   \code{map_limit}, also \code{cells_over_map_limit} and
#'   \code{pct_total_over_map_limit}: how much of the shipped total is invisible
#'   on the maps. A raster totalling zero gets \code{NA} concentration rather
#'   than a divide-by-zero. \code{max_cell_x} and \code{max_cell_y} give the
#'   coordinates of the largest cell, so a caller can ask what is special about
#'   it -- `03.70_Save_chosen_model_predictions.Rmd` uses them to look up its
#'   extrapolation status from [assess.extrapolation()].
#' @examples
#' \dontrun{
#' summarise.prediction.concentration(files, comb$species, comb$season)
#' }
#' @export
summarise.prediction.concentration <- function(paths, species, season, top_n = 10,
                                               map_limit = NULL) {
  checkmate::expect_character(paths, min.len = 1, any.missing = FALSE)
  checkmate::expect_atomic(species, len = length(paths))
  checkmate::expect_atomic(season, len = length(paths))
  checkmate::expect_count(top_n, positive = TRUE)
  checkmate::expect_number(map_limit, lower = 0, null.ok = TRUE)

  out <- purrr::pmap_dfr(list(paths, species, season), function(p, sp, se) {
    if (!file.exists(p))
      return(tibble::tibble(species = as.character(sp), season = as.character(se),
                            cells = NA_integer_, total = NA_real_,
                            max_cell = NA_real_, pct_top_n = NA_real_,
                            cells_over_map_limit = NA_integer_,
                            pct_total_over_map_limit = NA_real_,
                            max_cell_x = NA_real_, max_cell_y = NA_real_))
    r <- terra::rast(p)
    vall <- terra::values(r)[, 1]
    keep <- is.finite(vall)
    v <- vall[keep]
    tot <- sum(v)
    top <- sum(utils::head(sort(v, decreasing = TRUE), top_n))
    over <- if (is.null(map_limit)) v[0] else v[v > map_limit]
    # Coordinates of the largest cell, so callers can ask what is special about
    # it - 03.70 uses them to look up its extrapolation status.
    imax <- if (length(v)) which(keep)[which.max(v)] else NA_integer_
    xy <- if (is.na(imax)) c(NA_real_, NA_real_) else
      as.numeric(terra::xyFromCell(r, imax))
    tibble::tibble(species = as.character(sp), season = as.character(se),
                   cells = length(v), total = tot,
                   max_cell = if (length(v)) max(v) else NA_real_,
                   # A zero total is a legitimate answer for a species absent in
                   # a season; it just has no concentration to report.
                   pct_top_n = if (isTRUE(tot > 0)) 100 * top / tot else NA_real_,
                   cells_over_map_limit = length(over),
                   pct_total_over_map_limit =
                     if (isTRUE(tot > 0)) 100 * sum(over) / tot else NA_real_,
                   max_cell_x = xy[1], max_cell_y = xy[2])
  })
  dplyr::arrange(out, dplyr::desc(.data$pct_top_n))
}


# Gini coefficient of a non-negative vector.
#
# Used as one of the spikiness measures in compare.finalist.surfaces(). It is
# the scale-free one: multiplying every cell by a constant leaves it unchanged,
# so a Tweedie surface and a negative binomial surface predicting very different
# totals can still be asked which is the more unevenly distributed.
.gini <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2 || sum(v) <= 0) return(NA_real_)
  # Negative densities cannot come off a log link, but a raster read back from
  # disk is not a promise, and the formula below is only defined on
  # non-negative values.
  if (any(v < 0)) return(NA_real_)
  s <- sort(v)
  n <- length(s)
  (2 * sum(seq_len(n) * s) / (n * sum(s))) - (n + 1) / n
}

# cor() errors on a zero-variance argument rather than returning NA, and a
# season where a species is absent from the study area gives exactly that.
.safe.cor <- function(a, b, method) {
  if (length(a) < 3 || stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
  stats::cor(a, b, method = method)
}

# Spikiness of one surface, as a one-row tibble.
#
# Every measure here is deliberately scale-free apart from `total`, `max` and
# the over-limit pair, because the question is which family concentrates its
# prediction more, not which predicts more birds.
#
# `median_pos` is the median over cells holding any density at all: the median
# over every cell is often 0 in a season a species is largely absent from, and a
# ratio to it would then be undefined everywhere. It is RETURNED, not just used,
# because it is the denominator of `max_over_med`, `p99_over_med` and
# `n_spikes`, and on a near-empty surface it goes denormal - measured on
# NL_EXPL_DRL_RA RAZO Fall, a median positive cell of 2.0e-13 against a maximum
# of 0.31 birds/km2 gives `max_over_med` of 1.5e12, which reads as a
# catastrophic spike and is a vanishing denominator. Always read the ratio
# beside its denominator; `pct_top_n` and `gini` have no such failure mode.
#
# `n_over_limit` is the one that caught the real thing. Petrels Winter tw on
# NL_EXPL_DRL_RA has 2,267 cells above MAX_DENS_VALUE holding 99.9997% of the
# seasonal total and peaking at 1.7e10 birds/km2 - and its `pct_top_n` is only
# 14.9%, well under the level 03.70 warns on, because the blow-up is spread over
# thousands of cells rather than piled into one. A concentration measure looks
# for one big cell; this looks for many impossible ones.
.surface.spikes <- function(v, top_n, spike.ratio, map.limit, ref.dens,
                            ref.mult) {
  tot <- sum(v)
  srt <- sort(v, decreasing = TRUE)
  pos <- v[v > 0]
  med <- if (length(pos)) stats::median(pos) else NA_real_
  spikes <- if (isTRUE(med > 0)) v[v > spike.ratio * med] else v[0]
  p99 <- if (length(v)) unname(stats::quantile(v, 0.99)) else NA_real_
  over <- if (is.null(map.limit)) v[0] else v[v > map.limit]

  # A season the species was never recorded in has no ceiling to exceed. Zero
  # is the absence of a reference, not a reference of zero, and dividing by it
  # would make every such surface infinitely bad.
  has.ref <- !is.null(ref.dens) && !is.na(ref.dens) && ref.dens > 0
  ceiling <- if (has.ref) ref.mult * ref.dens else NA_real_
  over.ref <- if (has.ref) v[v > ceiling] else v[0]

  tibble::tibble(
    cells          = length(v),
    total          = tot,
    max            = if (length(srt)) srt[1] else NA_real_,
    median_pos     = med,
    p99            = p99,
    pct_top_n      = if (isTRUE(tot > 0))
                       100 * sum(utils::head(srt, top_n)) / tot else NA_real_,
    max_over_med   = if (isTRUE(med > 0)) srt[1] / med else NA_real_,
    p99_over_med   = if (isTRUE(med > 0)) p99 / med else NA_real_,
    n_spikes       = if (isTRUE(med > 0)) length(spikes) else NA_integer_,
    pct_in_spikes  = if (isTRUE(tot > 0) && isTRUE(med > 0))
                       100 * sum(spikes) / tot else NA_real_,
    n_over_limit   = if (is.null(map.limit)) NA_integer_ else length(over),
    pct_over_limit = if (is.null(map.limit) || !isTRUE(tot > 0)) NA_real_
                     else 100 * sum(over) / tot,
    obs_max        = if (has.ref) ref.dens else NA_real_,
    max_over_obs   = if (has.ref && length(srt)) srt[1] / ref.dens
                     else NA_real_,
    n_over_obs     = if (has.ref) length(over.ref) else NA_integer_,
    pct_over_obs   = if (has.ref && isTRUE(tot > 0))
                       100 * sum(over.ref) / tot else NA_real_,
    gini           = .gini(v))
}

#' Compare the prediction surfaces of two models, cell by cell
#'
#' Answers what the family cross-validation and the in-sample diagnostics both
#' leave open: when neither can separate the Tweedie and the negative binomial
#' finalist, do their predictions actually differ, and where?
#'
#' The comparison is deliberately two-sided. **Agreement** is the correlation
#' between the two surfaces and the ratio of their totals -- if those are tight
#' the choice of family does not matter for the product, and a tie can be left
#' as a tie. **Spikiness** is the rest of it, and is the half worth reading
#' first: two surfaces can correlate at 0.99 and still put wildly different
#' amounts of the total into a handful of cells, because that correlation is
#' carried by the tens of thousands of cells where both models predict almost
#' nothing.
#'
#' The spike measures are scale-free on purpose (`max_over_med`, `n_spikes`,
#' `gini`), so a family predicting twice the abundance is not thereby recorded
#' as twice as spiky. Two cautions come with them. The ones divided by
#' `median_pos` go meaningless on a near-empty surface, where that denominator
#' goes denormal -- `NL_EXPL_DRL_RA` RAZO Fall has a median positive cell of
#' 2.0e-13 and so a `max_over_med` of 1.5e12 off a maximum of 0.31 birds/km^2 --
#' which is why `median_pos` is returned beside them and should be read with
#' them. And none of them catches a surface that has blown up *everywhere*:
#' that is what `map.limit` is for.
#'
#' @section Calibrating the plausibility ceiling:
#' `max_over_obs` -- the largest predicted cell over the largest density ever
#' observed on a segment of the same species and season -- separates the
#' populations that a flat limit cannot. Measured over 88 `NL_EXPL_DRL_RA`
#' surfaces:
#'
#' \itemize{
#'   \item median **0.075**. The typical surface peaks thirteen times *below*
#'     anything ever seen, which is what should happen: a prediction cell is a
#'     seasonal mean over ten years and an observed segment density is one
#'     encounter on one day.
#'   \item 75th percentile 0.23, 90th 4.2, 95th 16.2. Twelve surfaces exceed 1.
#'   \item The legitimate tail tops out at **32.5** (LESP Winter under `nb`),
#'     with Shearwaters Spring at 19.0 and 11.1 -- concentrated, plausible,
#'     flock-driven surfaces.
#'   \item Then a gap of five orders of magnitude to **7.7e6**, Petrels Winter
#'     under `tw`.
#' }
#'
#' `ref.mult = 100` sits in that gap: three times above the worst legitimate
#' surface and four orders of magnitude below the pathology. It is a wide gap,
#' so the exact value is not delicate -- anything from 50 to 1,000 separates the
#' same two sets.
#'
#' Compare what the flat `MAX_DENS_VALUE` of 10,000 does on the same data. It
#' fires on 3 of 88 surfaces, and it cannot do better, because it is set near
#' the *top* of the observed range: real segments have held 7,183 birds/km2
#' (NOFU Spring), 4,091 (BLKI Fall) and 2,734 (BLKI Winter). For a razorbill in
#' spring, whose busiest segment ever held 16 birds/km2, a predicted cell of
#' 5,000 would be three hundred times anything ever recorded and the flat rule
#' would say nothing at all.
#'
#' @section What this found on NL_EXPL_DRL_RA:
#' Petrels Winter, fitted with `tw()` and the model `final.dsm.models` names, so
#' the surface `03.70` ships: **2,267 cells above `MAX_DENS_VALUE`**, together
#' holding 99.9997% of the seasonal total, peaking at 1.7e10 birds/km^2 for a
#' seasonal total of 3.95e12 birds. The `nb()` finalist on the same segments and
#' the same grid gives 6.3e6. The two surfaces correlate at 0.037.
#'
#' The existing concentration check does not see it. `pct_top_n` for that
#' surface is 14.9%, far below the level `03.70` warns on, because the blow-up
#' is spread over thousands of cells instead of piled into one -- and
#' `MAX_DENS_VALUE` is applied only when drawing maps, so every one of those
#' cells is absent from the picture and present in the file. A concentration
#' measure asks whether one cell dominates; `n_over_limit` asks whether many
#' cells are impossible, and only the second question had a useful answer here.
#'
#' Concentration is where the response family stops being a detail, which is why
#' this exists: where one observation dominates, Tweedie and negative binomial
#' disagree about how likely a huge count is, and the totals diverge
#' accordingly. See [summarise.prediction.concentration()] for the `Atl IMRP`
#' Herring Gull case that established the pattern.
#'
#' @param paths1,paths2 Character vectors of raster files to compare, one per
#'   season, in the same order. Both models must have been predicted on the same
#'   grid; the geometries are checked.
#' @param season Season labels, the same length as `paths1`.
#' @param species Character string species code, carried into the output.
#' @param model1,model2 Character strings naming the two models.
#' @param cell.area Cell area in km^2, i.e. `predgridCellArea`. The rasters hold
#'   density, so this is what turns a summed surface into an abundance.
#' @param top_n How many of the largest cells to accumulate and to list. Default
#'   10, matching `PRED_CONCENTRATION_TOP_N`.
#' @param spike.ratio A cell counts as a spike when it exceeds this multiple of
#'   the surface's median positive density. Default 10.
#' @param map.limit Optional flat density above which a cell is hidden from the
#'   maps, i.e. `MAX_DENS_VALUE`. Worth passing, but for a narrow reason: that
#'   limit is applied only when drawing leaflet maps and never to the rasters
#'   that get shared, so `n_over_limit` counts exactly the cells that are
#'   invisible on the map you would check a surface against and fully present in
#'   the file. It is a **visibility** measure, not a plausibility one -- see
#'   `ref.density` for that. `NULL` leaves `n_over_limit` and `pct_over_limit`
#'   `NA`.
#' @param ref.density The highest detection-corrected density ever observed on a
#'   segment, one value per season in the order of `season` (a scalar is
#'   recycled). From [observed.density.ceiling()]. This is the
#'   species-and-season-specific ceiling, and it is the right one: across
#'   `NL_EXPL_DRL_RA` the observed maximum spans 0 to 7,183 birds/km2 between
#'   species-seasons, so no single constant can be meaningful for all of them.
#'   A season with no observations at all gets no reference rather than a
#'   ceiling of zero. `NULL` leaves the `*_obs` columns `NA`.
#' @param ref.mult How many times `ref.density` a cell must exceed to be
#'   counted in `n_over_obs`. Default 100, measured rather than picked -- see
#'   the calibration section.
#' @return A list of two tibbles:
#'   \describe{
#'     \item{`$seasons`}{One row per season. Per-model columns are suffixed `1`
#'       and `2`; `total_ratio`, `max_ratio` and the correlations are the
#'       pairwise ones. `pct_cells_2x` is the share of cells where the two
#'       surfaces differ by more than a factor of two, counting only cells where
#'       at least one model predicts something -- the measure of *where* they
#'       disagree that a correlation cannot give. `top_n_shared` is how many of
#'       the two models' `top_n` cells are the same cells: whether they are
#'       spiky in the same places, as distinct from equally spiky.
#'       `n_over_limit` and `pct_over_limit` count cells hidden from the maps;
#'       `max_over_obs`, `n_over_obs` and `pct_over_obs` are the plausibility
#'       check against what was actually observed.}
#'     \item{`$cells`}{The `top_n` largest cells of each model, with the other
#'       model's value for the same cell, the coordinates and the ratio. This is
#'       what names the spots -- a caller can look each one up in the
#'       extrapolation assessment, exactly as `03.70` does for the chosen
#'       model.}
#'   }
#'   A season whose rasters are missing gets a row with `missing_raster` set
#'   rather than being dropped, so the caller can see which comparison did not
#'   happen.
#' @examples
#' \dontrun{
#' compare.finalist.surfaces(
#'   paths1 = file.path(predDir,
#'     sprintf("ATPU.%s.dsm_tw_allpred_abund_factor.4_sqkm.tif", season.names)),
#'   paths2 = file.path(predDir,
#'     sprintf("ATPU.%s.dsm_nb_allpred_abund_factor.4_sqkm.tif", season.names)),
#'   season = season.names, species = "ATPU",
#'   model1 = "dsm_tw_allpred_abund_factor",
#'   model2 = "dsm_nb_allpred_abund_factor",
#'   cell.area = predgridCellArea, map.limit = MAX_DENS_VALUE)
#' }
#' @export
compare.finalist.surfaces <- function(paths1, paths2, season, species,
                                      model1, model2, cell.area,
                                      top_n = 10, spike.ratio = 10,
                                      map.limit = NULL, ref.density = NULL,
                                      ref.mult = 100) {
  checkmate::expect_character(paths1, min.len = 1, any.missing = FALSE)
  checkmate::expect_character(paths2, len = length(paths1), any.missing = FALSE)
  checkmate::expect_atomic(season, len = length(paths1))
  checkmate::expect_string(species, min.chars = 1)
  checkmate::expect_string(model1, min.chars = 1)
  checkmate::expect_string(model2, min.chars = 1)
  checkmate::expect_number(cell.area, lower = 0)
  checkmate::expect_count(top_n, positive = TRUE)
  checkmate::expect_number(spike.ratio, lower = 1)
  checkmate::expect_number(map.limit, lower = 0, null.ok = TRUE)
  checkmate::expect_numeric(ref.density, lower = 0, null.ok = TRUE)
  checkmate::expect_number(ref.mult, lower = 1)

  # One reference per season, in the order the seasons were given. A scalar is
  # recycled so a caller with a single project-wide ceiling still works.
  if (!is.null(ref.density) && length(ref.density) == 1)
    ref.density <- rep(ref.density, length(paths1))
  if (!is.null(ref.density) && length(ref.density) != length(paths1))
    stop("compare.finalist.surfaces: ref.density must be length 1 or ",
         length(paths1), " (one per season), not ", length(ref.density))
  refs <- if (is.null(ref.density)) rep(NA_real_, length(paths1))
          else ref.density

  per.season <- purrr::pmap(list(paths1, paths2, season, refs),
                            function(p1, p2, se, ref) {
    se <- as.character(se)
    if (!file.exists(p1) || !file.exists(p2))
      return(list(
        seasons = tibble::tibble(species = species, season = se,
                                 model1 = model1, model2 = model2,
                                 missing_raster = TRUE),
        cells = NULL))

    r1 <- terra::rast(p1)
    r2 <- terra::rast(p2)
    terra::compareGeom(r1, r2, stopOnError = TRUE)

    v1all <- terra::values(r1)[, 1]
    v2all <- terra::values(r2)[, 1]

    # Each model's own finite cells for its own totals and spike measures, the
    # intersection for anything pairwise. These are the same set in every case
    # seen so far -- both surfaces come off one prediction grid -- but a family
    # that produced a non-finite value somewhere is itself the kind of finding
    # this function exists to surface, so it is counted rather than assumed
    # away.
    ok1 <- is.finite(v1all)
    ok2 <- is.finite(v2all)
    both <- ok1 & ok2

    s1 <- .surface.spikes(v1all[ok1], top_n, spike.ratio, map.limit,
                          ref, ref.mult)
    s2 <- .surface.spikes(v2all[ok2], top_n, spike.ratio, map.limit,
                          ref, ref.mult)

    a <- v1all[both]
    b <- v2all[both]

    # Cells where at least one model predicts something. A ratio between two
    # cells that both round to nothing is arithmetically enormous and
    # ecologically empty, so the disagreement measure is restricted to cells
    # carrying some density in at least one of the two.
    live <- (a > 0) | (b > 0)
    hi <- pmax(a[live], b[live])
    lo <- pmin(a[live], b[live])
    ratio <- hi / pmax(lo, .Machine$double.xmin)

    # The cells each model calls its largest, and what the other model says
    # about the same ground.
    top.of <- function(v, ok, mod, other) {
      idx <- which(ok)[utils::head(order(v[ok], decreasing = TRUE), top_n)]
      if (!length(idx)) return(NULL)
      xy <- terra::xyFromCell(r1, idx)
      tibble::tibble(
        species = species, season = se, top_of = mod, rank = seq_along(idx),
        cell = idx, x = xy[, 1], y = xy[, 2],
        dens = v[idx], other_dens = other[idx],
        ratio = v[idx] / other[idx])
    }

    seasons <- tibble::tibble(species = species, season = se,
                              model1 = model1, model2 = model2,
                              missing_raster = FALSE,
                              cells_compared = sum(both),
                              nonfinite1 = sum(!ok1),
                              nonfinite2 = sum(!ok2)) %>%
      dplyr::bind_cols(dplyr::rename_with(s1, ~ paste0(.x, "1")),
                       dplyr::rename_with(s2, ~ paste0(.x, "2"))) %>%
      dplyr::mutate(
        abund1       = .data$total1 * cell.area,
        abund2       = .data$total2 * cell.area,
        total_ratio  = if (isTRUE(s2$total > 0)) s1$total / s2$total
                       else NA_real_,
        max_ratio    = if (isTRUE(s2$max > 0)) s1$max / s2$max else NA_real_,
        pearson_r    = .safe.cor(a, b, "pearson"),
        spearman_r   = .safe.cor(a, b, "spearman"),
        # On the log scale as well, because the untransformed correlation
        # between two zero-inflated surfaces is set by their few largest cells
        # and reports near-perfect agreement almost regardless of the rest.
        log_r        = .safe.cor(log1p(a), log1p(b), "pearson"),
        live_cells   = sum(live),
        pct_cells_2x = if (sum(live)) 100 * sum(ratio > 2) / sum(live)
                       else NA_real_,
        median_abs_diff = if (length(a)) stats::median(abs(a - b))
                          else NA_real_,
        max_abs_diff    = if (length(a)) max(abs(a - b)) else NA_real_)

    cells <- dplyr::bind_rows(top.of(v1all, ok1, model1, v2all),
                              top.of(v2all, ok2, model2, v1all))
    seasons$top_n_shared <- if (!nrow(cells)) NA_integer_ else
      length(intersect(cells$cell[cells$top_of == model1],
                       cells$cell[cells$top_of == model2]))

    list(seasons = seasons, cells = cells)
  })

  list(seasons = dplyr::bind_rows(purrr::map(per.season, "seasons")),
       cells   = dplyr::bind_rows(purrr::map(per.season, "cells")))
}

#' Copy a species prediction HTML summary to the versioned predictions folder
#'
#' Copies the HTML report for the final model of \code{spec} into
#' \code{predVersionDir} (project global) and also copies the accompanying
#' \code{lib/} folder if needed.
#'
#' @param spec Character string species code.
#' @return \code{invisible(NULL)}, called for its side-effect (file copied).
#' @export
copy.prediction.summary <- function(spec){


  modname <- final.dsm.models$dsm_final_name[final.dsm.models$species == spec]
  source_path <- file.path(ResultsDir, spec, "Prediction summaries")
  create.dir.if.needed(predVersionDir)
  stopifnot(file.copy(file.path(source_path, paste0(modname, ".html")),
                      file.path(predVersionDir, paste0(spec,"_",  modname, ".html")),
                      overwrite = TRUE,
                      copy.date = TRUE))

  # Check if lib folder exists (contains needed .js files) and is not in
  # dest_path, and copy if needed.
  if (dir.exists(file.path(source_path, "lib")) &&
      !dir.exists(file.path(predVersionDir, "lib")))
    file.copy(
      file.path(source_path, "lib"),
      predVersionDir,
      recursive = TRUE,
      overwrite = TRUE
    )
}

#' Arrange four seasonal leaflet maps in a 2x2 HTML table
#'
#' \strong{Note: not currently used.}
#'
#' @param maps Named list of exactly four leaflet map objects, with names
#'   \code{"Spring"}, \code{"Summer"}, \code{"Fall"}, and \code{"Winter"}.
#' @return An HTML \code{tagList} containing a 2x2 table of maps.
#' @export
create.annual.map.grid <- function(maps) {
  if(length(maps) != 4)
    stop(sprintf("create_annual_map_grids: maps argument does contains %d maps - should be 4."),
         length(maps))

  res <-
    htmltools::tagList(tags$table(
      style = "width:100%",
      tags$tr(tags$td(htmltools::tagList(maps$Spring)),
              tags$td(htmltools::tagList(maps$Summer))),
      tags$tr(tags$td(htmltools::tagList(maps$Fall)),
              tags$td(htmltools::tagList(maps$Winter)))
    ))

  res
}

#' Save a list of leaflet maps to a timestamped HTML file
#'
#' Synchronises the maps with \code{leafsync::sync()} and saves the result to
#' \code{ResultsDir/<species>/Prediction summaries/} (project global).
#'
#' @param maps List of leaflet map objects (typically four seasonal maps).
#' @param modname Character string model name; used in the filename.
#' @param species Character string species code; used in the output path.
#' @return \code{invisible(NULL)}, called for its side-effect (file written).
#' @export
save.map <- function(maps, modname, species) {
  dirname <- here::here(ResultsDir, species, "Prediction summaries")
  if (!dir.exists(dirname))
    dir.create(dirname, recursive = TRUE)
  timestr <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- here::here(dirname, paste0(modname, "_", timestr, ".html"))
  message(sprintf("%s, %s: Saving map in %s.", species, modname, filename))
  list(htmltools::h2(paste0(species, "_", modname, "_", timestr)),
       leafsync::sync(maps)) %>%
    htmltools::tagList() %>%
    htmltools::save_html(file = filename)
}

#' Produce a patchwork of four seasonal ggplot prediction maps
#'
#' Iterates over seasons, calls \code{do.pred.map.ggplot} for each, and
#' combines the results with \code{patchwork::wrap_plots()}.
#'
#' @param dat \code{sf} prediction grid with \code{subset}, \code{Season},
#'   and \code{Dens} columns.
#' @param modname Character string model name; used in the plot title.
#' @param species Character string species code; used in messages.
#' @param subs Which platform subset to map: \code{"Combined"} (the default), or
#'   any single platform level present in \code{dat}. Validated against the data
#'   rather than a fixed list; see [do.pred.maps()].
#' @param ... Additional arguments passed to \code{do.pred.map.ggplot}.
#' @return A \code{patchwork} ggplot object.
#' @export
do.pred.maps.ggplot <-
  function(dat,
           modname,
           species,
           subs = "Combined",
           ...) {

    checkmate::expect_string(subs)
    available <- unique(as.character(dat$subset))
    if (!subs %in% available)
      stop(sprintf(
        paste0("do.pred.maps.ggplot: subset %s is not in these predictions. ",
               "Available: %s."),
        sQuote(subs), paste(sQuote(available), collapse = ", ")))

    # Get data subset
    dat <- dplyr::filter(dat, subset == subs)

    message(sprintf("%s, %s: Doing %s abundance prediction map for",
                    species, modname, subs))

    ret <- season.names %>%
      purrr::map(do.pred.map.ggplot, dat, modname, species, subs, ...)

    ret <- patchwork::wrap_plots(ret) + patchwork::plot_annotation(title = modname)
    ret
  }

#' Produce a single-season ggplot density prediction map
#'
#' @param season Character string season label.
#' @param dat \code{sf} prediction grid with \code{Season} and \code{Dens}
#'   columns.
#' @param modname Character string model name; used in messages.
#' @param species Character string species code; used in messages.
#' @param subs Platform subset label, used only for the plot title. Its caller
#'   [do.pred.maps.ggplot()] has already done the filtering and validation.
#' @param samp_n Integer; if not \code{NA}, plot a random sample of this many
#'   polygons.
#' @return A \code{ggplot} object.
#' @export
do.pred.map.ggplot <-
  function(season,
           dat,
           modname,
           species,
           subs = "Combined",
           samp_n = NA) {

    message(sprintf("\t%s",season))
    dat <- dat %>%
      dplyr::filter(Season == season) %>%
      dplyr::select(Dens, geometry) %>%
      rmapshaper::ms_simplify()


    # Plot only a sample of the polygons for efficiency? Typically used for
    # testing.
    if (!is.na(samp_n)) {
      index <- sample(1:nrow(dat), size = samp_n)
      dat <- dat[index, ]
    }

    # Remove ridiculously large densities b/c they mess up the legend and swamp
    # everything else
    dat <- dplyr::mutate(dat,
                         Dens = dplyr::case_when(Dens > MAX_DENS_VALUE ~ NA,
                                                 TRUE ~ Dens))

    ret <- ggplot2::ggplot(dat = dat) +
      ggplot2::geom_sf(ggplot2::aes(fill = Dens), color = NA)  +
      ggplot2::scale_fill_continuous(
        type = "viridis"
        # breaks = class_intervals$brks,
        # labels = round(class_intervals$brks, 4)
      ) +
      ggplot2::ggtitle(season)
    ret
  }


#' Quick leaflet map of watch start positions for debugging
#'
#' @param dat Data frame with \code{LongStart}, \code{LatStart}, and
#'   \code{WatchID} columns.
#' @return A \code{leaflet} map object.
#' @export
watch.map <- function(dat) {
  leaflet::leaflet(dat) %>%
    leaflet::addTiles() %>%
    leaflet::addCircles(
      lng = ~ LongStart,
      lat = ~ LatStart,
      radius = 1,
      popup = ~ htmlEscape(paste0("WatchID ", WatchID))
    )
}

#' Create a directory recursively if it does not already exist
#'
#' Returns \code{pathname} invisibly so the function can be used in a
#' pipeline.
#'
#' @param pathname Character string directory path to create.
#' @return \code{pathname} (invisibly).
#' @export
create.dir.if.needed <- function(pathname){
  if (!dir.exists(pathname)){
    dir.create(pathname, recursive = TRUE)
    message(sprintf("Creating needed folder %s", pathname))
  }
  return(invisible(pathname))

}


#' Initialise a default DDF model list structure for all species groups
#'
#' Creates one copy of \code{def.ddf.list} (project global) per entry in
#' \code{spec.grps} and fills in platform-appropriate default values via
#' \code{\link{set.def.df.spec.values}}.
#'
#' @param saveit If \code{TRUE}, save the list to \code{dfModlistLoc}
#'   (project global).
#' @return Named list of DDF spec lists, one per species group.
#' @export
init.df.mod.list <- function(saveit = FALSE) {

  # Create one copy of the basic list for each species group
  df.mod.list <- rep(list(def.ddf.list), length(spec.grps))
  names(df.mod.list) <- names(spec.grps)

  # Set default convert_units value based on whether survey is aerial or ship.
  # This may get updated dynamically by later processing.
  df.mod.list %<>%
    purrr::map( ~ purrr::imap(., set.def.df.spec.values))

  # Sometimes we want to save it (ie if it didn't already exist) but other times
  # we're just called to return a initial structure that can be modified (eg in
  # 01.02_Create_final_ddf_model_specs.Rmd)
  if (saveit)
    save(df.mod.list, file = dfModlistLoc)

  df.mod.list
}

#' Assign distance type by looking up DistMeth in the ECSAS database
#'
#' Joins \code{dat} against the \code{lkpDistMeth} table from the ECSAS
#' database and returns a factor of distance types (\code{"None"},
#' \code{"Perp."}, or \code{"Radial"}) based on \code{FlySwim} and
#' \code{DistMeth}.
#'
#' @param dat Data frame containing at minimum \code{DistMeth} and
#'   \code{FlySwim} columns.
#' @return Factor vector of distance types, the same length as \code{nrow(dat)}.
#' @export
assign.dist.type <- function(dat) {
  distmeth <- ECSASconnect::ECSAS.get.table(ecsas.path = ECSAS.Path, "lkpDistMeth")
  DistType <- dplyr::left_join(dat, distmeth, by = c("DistMeth" = "DistMethCode")) %>%
    dplyr::mutate(DistType = as.factor(
      dplyr::case_when(
        FlySwim == "F" ~ DistMethFlyHow,
        FlySwim == "W" ~ DistMethWaterHow,
        TRUE ~ NA_character_
      )
    )) %>%
    dplyr::pull(DistType)

  if (any(is.na(DistType)))
    warning("assign.dist.type: ",
            sum(is.na(DistType)),
            " rows could not be assigned a distance type",
            immediate. = TRUE)

  DistType
}


#' Impute a distance for observations that should have had one
#'
#' There are two ways an observation ends up without a perpendicular distance,
#' and they mean opposite things.
#'
#' The one that matters is **by design**: the watch's `DistMeth` says no
#' perpendicular distances are taken for that behaviour - `lkpDistMeth` gives
#' flying birds `None` under DistMeth 1 and 13 and `Radial` under 17 and 19, and
#' DistMeth 17 was standard from mid-2008 to the end of 2011. Those observations
#' are separate survey effort, they belong in the `_N` (strip/dummy) ddf, and
#' this function leaves them alone.
#'
#' The other is an **observer omission**: the DistMeth says perpendicular
#' distances are being recorded, and one simply was not. Those did not happen on
#' separate effort, so putting them in the `_N` ddf misattributes them. They used
#' to be swept up by a blanket `is.na(distance) ~ 0` fill, which is worse than
#' either option - it fabricates a detection exactly on the trackline, biasing
#' the detection function steeply toward zero distance.
#'
#' They are rare enough that the choice barely matters in aggregate - 8 rows of
#' 47,307 perpendicular observations on NL_EXPL_DRL_RA, 0.017% - but "rare" is
#' not "absent" and a silent zero is not a defensible value. Each gets the mean
#' observed distance for its own (`SurveyType`, `FlySwim`, `DistMeth`) group,
#' snapped to the nearest distance actually observed in that group.
#'
#' The snapping matters because these data are binned: ship observations carry
#' the bin centre in `distance` (0.025, 0.075, 0.15, 0.25 km) and have `distbegin`
#' / `distend` NA, since the binning happens at ddf-fitting time from the spec's
#' `cutpoints`. Assigning a raw group mean would put a value on the pile that is
#' not a bin centre; snapping keeps every distance one of the values the survey
#' can actually produce. Where `distbegin`/`distend` ARE populated (the aerial
#' path), the pair belonging to the chosen distance is copied across too.
#'
#' Deterministic - no RNG, so a rerun gives the same answer.
#'
#' @param distdata Data frame with `distance`, `DistType`, `SurveyType`,
#'   `FlySwim` and `DistMeth` columns. `distbegin`/`distend` are used if present.
#' @return `distdata` with imputable rows filled in. Rows whose `DistType` is not
#'   `"Perp."`, and rows in a group with no observed distances to average, are
#'   returned unchanged.
#' @export
impute.missing.perp.distances <- function(distdata) {
  checkmate::expect_data_frame(distdata)
  for (col in c("distance", "DistType", "SurveyType", "FlySwim", "DistMeth"))
    if (is.null(distdata[[col]]))
      stop("impute.missing.perp.distances: distdata has no '", col, "' column.")
  has.bins <- !is.null(distdata$distbegin) && !is.null(distdata$distend)

  target <- distdata$DistType == "Perp." & is.na(distdata$distance)
  target[is.na(target)] <- FALSE
  if (!any(target)) {
    message("impute.missing.perp.distances: nothing to impute.")
    return(distdata)
  }

  # Donors are the rows that DO have a distance, in the same group.
  donor <- which(distdata$DistType == "Perp." & !is.na(distdata$distance))
  key <- function(rows)
    paste(distdata$SurveyType[rows], distdata$FlySwim[rows],
          distdata$DistMeth[rows], sep = "\r")
  donor.key <- key(donor)

  n.done <- 0L
  for (i in which(target)) {
    pool <- donor[donor.key == key(i)]
    if (length(pool) == 0) next

    dists <- distdata$distance[pool]
    observed <- sort(unique(dists))
    chosen <- observed[which.min(abs(observed - mean(dists)))]

    distdata$distance[i] <- chosen
    if (has.bins) {
      # Copy the bin edges of a donor sitting at the chosen distance, if it has
      # any. Ship data has none; aerial does.
      src <- pool[distdata$distance[pool] == chosen &
                    !is.na(distdata$distbegin[pool])]
      if (length(src) > 0) {
        distdata$distbegin[i] <- distdata$distbegin[src[1]]
        distdata$distend[i] <- distdata$distend[src[1]]
      }
    }
    n.done <- n.done + 1L
  }

  message(sprintf(
    paste0("impute.missing.perp.distances: imputed %d of %d observation(s) that ",
           "should have had a perpendicular distance%s."),
    n.done, sum(target),
    if (n.done < sum(target))
      sprintf(" (%d had no same-group observations to average)",
              sum(target) - n.done) else ""))

  distdata
}

#' Generate the canonical detection function model name
#'
#' Combines the key function, formula, and adjustment term into the dot-
#' separated label used for output filenames and model lookup.
#'
#' @param key Character string key function (e.g. \code{"hn"}).
#' @param form Formula or character string formula; \code{~ 1} produces a
#'   key-only or adjustment-only name.
#' @param adj Character string adjustment term (\code{"cos"}, \code{"herm"},
#'   or \code{"poly"}), or \code{NULL}.
#' @return Character string model name.
#' @export
create.model.name <- function(key, form, adj) {
  # Was there a null formula?
  if (form == as.formula(~ 1)) {
    # Grab adjustment, converting names to shorthand used in filenames
    form <-
      if (is.null(adj))
        NULL
    else
      c(cos = "cos",
        herm = "hp",
        poly = "sp")[adj]
  } else {
    # Convert formula to dotted notation
    form <-
      stringr::str_replace(as.character(form)[2], stringr::fixed(" + "), ".")
  }
  # Deal with adjustment that might be NULL with ifelse()
  nam <- paste0(key, ifelse(is.null(form), "", paste0(".", form)))
  nam
}


#' Back up DSM summary RData files for one species
#'
#' Copies all \code{*.Rdata} files from
#' \code{ResultsDir/<species>/DSM Summaries/} into
#' \code{ResultsDir/Backups/<species>/DSM Summaries/}, appending the file
#' modification timestamp to each filename.
#'
#' @param species Character string species code.
#' @return \code{invisible(NULL)}, called for its side-effect (files copied).
#' @export
backup.dsm.summary <- function(species){
  folder <- here::here(ResultsDir, "Backups", species, "DSM Summaries")
  if (!dir.exists(folder))
    dir.create(folder, recursive = TRUE)

  message("Backing up ", species, " DSM summaries to ", folder)

  # src files
  src <- list.files(path = here::here(ResultsDir, species, "DSM Summaries"),
                    pattern = "*.Rdata", full.names = T)

  if (length(src) == 0){
    message("\tNo files to backup - quitting.")
    return
  }

  ## create dst filenames

  # add modification dates to files
  mtimes <- stringr::str_replace_all(file.info(src)$mtime, " ", "_") %>%
    stringr::str_replace_all(":", "") %>%
    stringr::str_replace("\\..*$", "") # remove trailing milliseconds
  stopifnot(length(src) == length(mtimes))
  dst <- src %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    here::here(folder, .) %>%
    paste0("_", mtimes, ".", tools::file_ext(src))

  file.copy(src, dst, overwrite = TRUE, copy.date = TRUE) %>%
    invisible
}

#' Back up DSM summary files for all species groups
#'
#' Calls \code{\link{backup.dsm.summary}} for every entry in the project
#' global \code{spec.grps}.
#'
#' @return \code{invisible(NULL)}.
#' @export
backup.dsm.summaries <- function(){
  names(spec.grps) %>%
    purrr::walk(backup.dsm.summary)
}


#' Delete all DSM summary RData files for all species groups
#'
#' @return \code{invisible(NULL)}, called for its side-effect (files deleted).
#' @export
remove.dsm.summaries <- function(){
  names(spec.grps) %>%
    purrr::walk(\(species){
      message("Removing DSM summaries for ", species)
      files <- list.files(path = here::here(ResultsDir, species, "DSM Summaries"),
                          pattern = "*.Rdata", full.names = T)
      file.remove(files)
    })
}


#' List dates available for a dynamic environmental covariate
#'
#' Scans the NetCDF files in \code{predLayerDir/NetCDF/<var.name>/} and
#' returns a data frame of available dates and their source filenames.
#'
#' @param var.name Character string variable name (e.g. \code{"sst"}).
#' @return Data frame with columns \code{date} (character) and
#'   \code{filename}, typically with 12 rows (one per month).
#' @export
get.covar.netCDF.dates <- function(var.name){
  # list files in the folder here(predLayerDir, "NetCDF", var.name)
  # get dates associated with files, convert to char and return
  #
  files <- list.files(here::here(predLayerDir, "NetCDF", var.name),
                      ".*\\.nc$",
                      full.names = TRUE)

  suppressWarnings(
    dates <- purrr::map_dfr(
      files,
      \(filenm) {
        data.frame(
          date = stars::read_ncdf(filenm, var = "time", proxy = TRUE) %>%
            stars::st_get_dimension_values("time") %>%
            stringr::str_replace(" UTC", ""),
          filename = filenm
        )
      }))

  dups <- duplicated(dates$date)
  if(any(dups)) {
    warning(sprintf(c("get.covar.netCDF.dates: variable %s: the following dates",
                      " were found in multiple files: "), var.name), immediate. = TRUE)
    dates[dups,]
  }

  dates
}

#' Retrieve, reproject, and clip one environmental covariate
#'
#' Reads the covariate described by \code{env_covar_spec} from local NetCDF
#' files, reprojects to \code{segProj} (project global), and clips to the
#' study area.  For dynamic covariates, checks that all required dates are
#' available and downloads missing ones if needed.
#'
#' Expects project globals \code{predLayerDir}, \code{segProj}, and
#' \code{study.area}.
#'
#' @param env_covar_spec Single-row data frame from the project
#'   \code{env_covar_spec} table with columns \code{var_name},
#'   \code{var_type}, \code{ERDDAP_dataset_name}, \code{netcdf_vars}, and
#'   \code{CRS}.
#' @param dates.needed Character vector of date strings required for dynamic
#'   covariates.
#' @param verbose If \code{TRUE} (default), print progress messages.
#' @return A \code{stars} object clipped to the study area.
#' @export
get.env.covar <- function(env_covar_spec,
                          dates.needed,
                          verbose = TRUE
) {


  var_name <- env_covar_spec$var_name

  message("get.env.covar: getting covariate '", var_name, "'")

  # Dynamic or static covar
  if (env_covar_spec$var_type == "static"){
    env_dat <- stars::read_ncdf(here::here(
      predLayerDir,
      "NetCDF",
      paste0(env_covar_spec$ERDDAP_dataset_name, ".nc")
    ), var = env_covar_spec$netcdf_vars)


  } else if (env_covar_spec$var_type == "dynamic") {
    # Dynamic involves making sure we have the dates we want.
    netcdf.dates.have  <- get.covar.netCDF.dates(var_name)
    dates.to.get <- setdiff(dates.needed, netcdf.dates.have$date)

    if (length(dates.to.get) > 0) {
      # Download needed files
      # XXXX TODO: need to use curl (or something ) to download the netcdf files
      # we need by constructing the correct URL.
    }

    # Now we should have all files we need.
    # Figure out which files to read, just in case we have files we don't need
    files <- dplyr::filter(netcdf.dates.have, date %in% dates.needed) %>%
      dplyr::pull(filename) %>%
      unique()

    # Read all needed netcdf files
    #
    # # This approach didn't work. Something to do with fact that some of the
    # .nc files have only one var in them ("sst") and some have both ("sst" and
    # "mask") due to different ways they were downloaded.
    #
    # res <- map(files, \(filenm) {
    #   stars::read_ncdf(filenm, var = var_name)
    # })
    #
    # # Collapse list of stars objects to single object with all dates combined.
    # env.dat <- Reduce(c, res)
    env_dat <- stars::read_stars(files, sub = var_name)
    if(is.na(sf::st_crs(env_dat)))
      sf::st_crs(env_dat) <- env_covar_spec$CRS
  } else
    stop("get.env.covar: variable ",
         var_name,
         ": illegal var_type '",
         env_covar_spec$var_type, "'")

  if (verbose)
    message("\tProjecting and clipping to study area...", appendLF = FALSE)

  # Re-project, clip, etc
  final <- sf::st_transform(env_dat, segProj) %>%
    `[`(study.area) %>%
    setNames(var_name) # b/c previous processing steps loose the name

  if (verbose)
    message("done.")

  final
}



#' Process one environmental covariate for both segdata and prediction grid
#'
#' Similar to \code{\link{get.env.covar}} but handles separate date sets for
#' the segment data and the prediction grid.
#'
#' Expects project globals \code{predLayerDir}, \code{segProj}, and
#' \code{study.area}.
#'
#' @param env_covar_spec Single-row data frame from the project
#'   \code{env_covar_spec} table.
#' @param segdata.dates.needed Character vector of dates needed for segdata
#'   extraction.
#' @param predgrid.dates.needed Character vector of dates needed for the
#'   prediction grid.
#' @return A \code{stars} object clipped to the study area.
#' @export
do.env.covar <- function(env_covar_spec,
                         segdata.dates.needed,
                         predgrid.dates.needed) {
  var_name <- env_covar_spec$var_name

  # Dynamic or static covar
  if (env_covar_spec$var_type == "static"){
    env.dat <- stars::read_ncdf(here::here(
      predLayerDir,
      "NetCDF",
      paste0(env_covar_spec$ERDDAP_dataset_name, ".nc")
    ), var = env_covar_spec$netcdf_vars)


  } else if (env_covar_spec$var_type == "dynamic") {
    # Dynamic involves making sure we have the dates we want.
    netcdf.dates.have  <- get.covar.netCDF.dates(var_name)
    all.dates.needed <- union(segdata.dates.needed, predgrid.dates.needed)
    dates.to.get <- setdiff(all.dates.needed, netcdf.dates.have)

    if (length(dates.to.get) > 0) {
      # Download needed files
      # XXXX TODO: need to use curl (or something ) to download the netcdf files
      # we need by constructing the correct URL.
    }

    # Now we should have all files we need.
    # Figure out which files to read
    files <- dplyr::filter(netcdf.dates.have, date %in% all.dates.needed) %>%
      dplyr::pull(filename) %>%
      unique()

    # Read all needed netcdf files
    res <- purrr::map(files, \(filenm) {
      stars::read_ncdf(filenm, var = var_name)
    })

    # Collapse list of stars objects to single object with all dates combined.
    # xxx this doesn't work, returned object seems to only have "attr" as it's
    # values
    env.dat <- Reduce(c, res)
  } else
    stop("do.env.covar: variable ",
         var_name,
         ": illegal var_type '",
         env_covar_spec$var_type, "'")

  # Re-project, clip, etc and save - xxx this line doesn't work
  final <- sf::st_transform(env.dat, segProj) %>%
    `[`(study.area)

  final
}


#' Create the standard folder structure for a new subproject
#'
#' Creates every per-subproject folder the pipeline writes into, so that a new
#' subproject can be run from \code{00.01} onwards without a step failing on a
#' missing directory.  Safe to re-run: existing folders are left alone.
#'
#' Note that every path is built from \code{subproj}, never from the
#' \code{SubProject} global.  Those are usually the same, but not while you are
#' setting up a new subproject - which is the only time this function is
#' called.
#'
#' Expects project globals \code{GenericRDataDir}, \code{GenericShapeDir},
#' \code{predLayerDir} and \code{GISDir}.  These are not subproject-specific,
#' unlike \code{predLayerStudyAreaDir} and \code{predLayerGridDir}, which is
#' why the subproject-level paths are assembled here rather than taken from
#' those globals.
#'
#' @param subproj Character string subproject identifier.
#' @return \code{invisible(NULL)}, called for its side-effect (directories
#'   created).
#' @examples
#' # create_subproject_folders("Atl IMRP")
#' @export
create_subproject_folders <- function(subproj) {
  checkmate::expect_string(subproj, min.chars = 1)

  create.dir.if.needed(file.path(GenericRDataDir, subproj))
  create.dir.if.needed(file.path(GenericShapeDir, subproj))
  create.dir.if.needed(here::here("Results", subproj))
  create.dir.if.needed(file.path(GISDir, "Predictions", subproj))
  create.dir.if.needed(file.path(GISDir, "Rasters", subproj))

  # Covariate rasters on the analysis grid, masked to the study area. 00.02
  # writes depth.img and depth.g.img at the top level, the individual monthly
  # sst/sst.g layers into sst/, and the monthly climatology into Predgrid/.
  # 00.03 and create.segdata() read them back from here.
  sa.dir <- file.path(predLayerDir, "Study area resolution & extent", subproj)
  create.dir.if.needed(sa.dir)
  create.dir.if.needed(file.path(sa.dir, "sst"))
  create.dir.if.needed(file.path(sa.dir, "Predgrid"))

  # The same covariates unmasked and buffered, which is what the gradients in
  # 00.02 are computed from and what use.env.raster.cache reads back. Same
  # three-way layout.
  grid.dir <- file.path(predLayerDir, "Analysis grid unmasked", subproj)
  create.dir.if.needed(grid.dir)
  create.dir.if.needed(file.path(grid.dir, "sst"))
  create.dir.if.needed(file.path(grid.dir, "Predgrid"))

  invisible(NULL)
}


#' Rasterize one variable from an sf object
#'
#' Creates a raster template from the extent of \code{sfobj} at resolution
#' \code{predgridCellLength * 1000} metres (project global) and rasterizes
#' \code{variable}.
#'
#' @param sfobj \code{sf} object containing the variable to rasterize.
#' @param variable Character string name of the column to rasterize.
#' @return A single-layer \code{SpatRaster}.
#' @export
make.raster <- function(sfobj, variable){
  # Convert sf to SpatVector
  v <- terra::vect(sfobj)

  # Create a raster template with the same extent and resolution. Assumes
  # predgridCellLength is in km and raster projection units are metres.
  r <- terra::rast(v, resolution = predgridCellLength * 1000)

  # Rasterize, using an attribute field (e.g., "ID")
  ret <- terra::rasterize(v, r, field = variable)

  ret
}

#' Check whether a gam/bam was fitted with discrete = TRUE
#'
#' @param model A \code{gam} or \code{bam} object from \code{mgcv}.
#' @return \code{TRUE} if the model was fitted with \code{discrete = TRUE},
#'   otherwise \code{FALSE}.
#' @export
is_discrete_gam <- function(model) {
  # Ensure it's a gam/bam object
  if (!inherits(model, "gam")) {
    stop("Model must be a gam/bam object from mgcv")
  }

  call_list <- as.list(model$call)

  # If "discrete" not supplied, default is FALSE
  if (!"discrete" %in% names(call_list)) {
    return(FALSE)
  }

  # Evaluate, in case it’s e.g. discrete = getOption("mgcv.discrete")
  isTRUE(eval(call_list$discrete, envir = parent.frame()))
}


#' Render a one-off RMarkdown file for a single species
#'
#' Sources \code{R/analysis_settings.R} (project-specific path) then renders
#' \code{rmdfile} with \code{species} as a parameter, timestamping the output
#' filename.
#'
#' @param rmdfile Character string filename (not full path) of the RMarkdown
#'   file to render, relative to \code{RDir}.
#' @param species Character string species code passed as a render parameter.
#' @return \code{invisible(NULL)}, called for its side-effect (HTML rendered).
#' @export
do_oneoff_render <- function(rmdfile, species) {
  source(here::here("R/analysis_settings.R"), echo = T)

  rmarkdown::render(
    here::here(RDir, rmdfile),
    params = list(species = species),
    output_file = here::here(
      ResultsDir,
      paste0(
        tools::file_path_sans_ext(rmdfile),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".html"
      )
    )
  )
}

# ============================================================================
# From ds_utils_0.6.R
# ============================================================================
# This file contains utils for distance sampling. Dave Fifield 2014.

#' Fit a single distance sampling model and save summary output
#'
#' Calls \code{\link[Distance]{ds}} with the specification in \code{mod},
#' writes a plain-text summary and an \code{.RData} file named
#' \code{AIC_<aic>_<label>} to \code{folder}, and returns the fitted model
#' object.  If fewer than \code{min_data_size_limit} observations are present
#' or fitting fails, a dummy \code{failed_model_*.txt} file is created and
#' \code{NA} is returned.
#'
#' @param mod Single-row data frame with columns \code{label}, \code{key},
#'   \code{adj}, and \code{form}.
#' @param data Observation data frame passed to \code{\link[Distance]{ds}}.
#' @param folder Directory path for output files; use \code{""} to skip
#'   file output.
#' @param logfileConn Log destination: \code{""} for stdout, or an open
#'   write-mode file connection.
#' @param verbose If \code{TRUE}, print extra debugging messages.
#' @param min_data_size_limit Minimum number of observations required to
#'   attempt model fitting.
#' @param ... Additional arguments passed to \code{\link[Distance]{ds}}.
#' @return Fitted \code{dsmodel} object, or \code{NA} on failure.
#' @export
run.ddf.model <-
  function(mod,
           data,
           folder = "",
           logfileConn = "",
           verbose,
           min_data_size_limit = 20,
           ...) {

    if (nrow(mod) > 1) stop("run.ddf.model: given more than 1 model to run!")

    cat(paste("Running", mod$label, "..."), file = logfileConn)

    # set up for adjustments
    adj = NULL
    if (mod$adj != "")
      adj <- mod$adj

    # Make sure we haven't specified both covars and adjustment terms. Need to
    # deal with case where formula is either a formula or a string.
    if (!is.null(adj) && mod$form != as.formula("~1") &&
        stringr::str_replace_all(mod$form, " ", "") != "~1") {
      mess <- sprintf("run.ddf.model: model (%s) has both covars and adjustment terms!",
                      mod$label)
      stop(mess)
    }

    monotonicity <- "none"
    if (mod$form == "~1")
      monotonicity <- "strict"

    model <- NULL
    res <- NULL
    if (nrow(data) >= min_data_size_limit) {
      # fit model - sometimes fitting fails
      res <- try({
        # create a new environment with only the things we want in attempt to keep
        # .rda file small when we eventually use save() which stores everything in the
        # environments in which any formula was created
        env <- new.env(parent = globalenv())
        env$data <- data
        env$mod <- mod
        env$adj <- adj
        env$monotonicity <- monotonicity
        env$dots <- list(...)

        if (logfileConn == "")
          model <- with(env, {
            model <- do.call('ds',
                             c(
                               list(
                                 data = data,
                                 key = mod$key,
                                 formula = as.formula(mod$form),
                                 monotonicity = monotonicity,
                                 adjustment = adj
                               ),
                               dots
                             ))
            return(model) # This just returns outside this block
          })
        else
          capture.output(model <- with(env, {
            model <- do.call('ds',
                             c(
                               list(
                                 data = data,
                                 key = mod$key,
                                 formula = as.formula(mod$form),
                                 monotonicity = monotonicity,
                                 adjustment = adj
                               ),
                               dots
                             ))
            return(model)
          }),
          file = logfileConn)
      })
    }

    # model fitting failed. Create dummy output file to stop model from being
    # attempted again.
    if (inherits(res, "try-error") || is.null(model)) {
      file.create(file.path(folder, paste0("failed_model_", mod$label, ".txt")))
      return(NA)
    }

    cat("finished!\n", file = logfileConn)

    # send summary output to file named "AIC_XXX.XXX_`mod$label.txt`"
    if (folder != "") {
      filename = paste(folder, paste0("AIC_", round(model$ddf$criterion, 3), "_",
                                      mod$label, ".txt"), sep = "/")

      if (verbose)
        cat(paste0("Opening model output file '", filename, "' \n"), file = logfileConn)
      modOutFileConn <- file(filename, open = "wt")

      # Create model output file
      if (exists("model")) {
        try.res <- try({
          # first spit out all the things that we want in the model summary
          ch <- mrds::ddf.gof(model$ddf, qq = FALSE)$chisquare$chi1
          summ <- summary(model)
          cat(paste0("Date: ", date(), "\n"), file = modOutFileConn)
          cat(paste0("Model name: ", mod$label, "\n"), file = modOutFileConn)
          cat(paste0("Key: ", summ$ds$key, "\n"), file = modOutFileConn)
          cat(paste0("Formula: ", model$ddf$ds$aux$ddfobj$scale$formula),
              "\n",
              file = modOutFileConn)
          cat(paste0("AIC: ",  round(model$ddf$criterion, 3), "\n"), file = modOutFileConn)
          #mrds::ddf.gof(model$ddf, qq=FALSE)$dsgof$CvM$p,
          cat(paste0("Chi_scores: ", paste(round((ch$observed - ch$expected) ^
                                                   2 / ch$expected, 3
          ), collapse = ", "), "\n"), file = modOutFileConn)
          cat(paste0("Chisquare p: ",  ch$p, "\n"), file = modOutFileConn)
          cat(paste0("Det prob: ", round(summ$ds$average.p, 3), "\n"), file = modOutFileConn)
          # May not exist for unif key model
          if ("average.p.se" %in% names(summ$ds)) {
            cat(paste0("SE(p): ", round(summ$ds$average.p.se, 5), "\n"), file = modOutFileConn)
            cat(paste0(
              "CV(p): ",
              round(summ$ds$average.p.se / summ$ds$average.p, 3),
              "\n"
            ), file = modOutFileConn)
          }
          cat("\n", file = modOutFileConn)
          capture.output(summary(model), file = modOutFileConn)
        })

        if (inherits(try.res, "try-error")) {
          cat("Creating model summary file failed.\n", file = logfileConn)
        }

        # Now save the model object as an .RData file
        save(model, file = sub(".txt", ".RData", filename))

      } else {
        cat("Model failed to fit.\n", file = logfileConn)
      }
      close(modOutFileConn)
    }

    if (verbose) cat("Returning model object\n", file = logfileConn)
    model
  }


#' Generate all n-way covariate combinations for a key function
#'
#' Produces a data frame of model specifications with all C(length(covars), n)
#' covariate combinations for the given key function.
#'
#' @param n Integer number of covariates to include per model (\code{0} for
#'   intercept-only).
#' @param key Character string key function (e.g. \code{"hn"}).
#' @param covars Character vector of candidate covariate names.
#' @return Data frame with columns \code{form}, \code{label}, \code{key}, and
#'   \code{adj}.
#' @export
gen.form.N <- function(n, key, covars){
  require(utils)

  # if n == 0 then no covars, just key function.
  if (n == 0) {
    return(data.frame(
      form = "~1",
      label = key,
      key = key,
      adj = "",
      stringsAsFactors = F
    ))
  }

  # get all combinations of chooseing n items from covars
  covarList <- combn(covars, n, simplify = F)

  plyr::ldply(covarList, function(item) {
    data.frame(form = paste("~", paste(item, collapse = " + ")), key = key,
               label = paste(key, paste(item, collapse = "."), sep = "."), adj = "", stringsAsFactors = F)
  }
  )
}

#' Fit all candidate distance sampling models
#'
#' Generates every combination of key function and covariates (up to
#' \code{length(covars)} covariates per model), optionally including
#' adjustment-term-only models, and fits each via
#' \code{\link{run.ddf.model}}.  Models with existing result files in
#' \code{folder} are skipped unless \code{rerun = TRUE}.
#'
#' @param data Observation data frame passed to \code{\link[Distance]{ds}}.
#' @param key Character vector of key functions to use; default
#'   \code{c("unif", "hn", "hr")}.
#' @param covars Character vector of candidate covariate names; \code{NULL}
#'   for no covariates.
#' @param runModels If \code{FALSE}, return the model specification data frame
#'   without fitting any models (useful for inspection).
#' @param parallel If \code{TRUE}, use \code{doSNOW} parallel processing.
#' @param nCores Number of parallel cluster nodes.
#' @param folder Directory path for output files.
#' @param logfile Path for a log file; use \code{""} to log to stdout.
#' @param rerun If \code{TRUE}, re-fit models whose output files already exist
#'   in \code{folder}.
#' @param verbose If \code{TRUE}, print extra messages.
#' @param models Optional data frame of pre-built model specifications
#'   (overrides automatic generation).
#' @param incl.adj If \code{TRUE} (default), prepend key+adjustment-only
#'   models from \code{adj.models} (project global).
#' @param cleanFolder If \code{TRUE}, delete all existing files from
#'   \code{folder} before fitting.
#' @param ... Additional arguments passed to \code{\link{run.ddf.model}}.
#' @return When \code{runModels = TRUE}: named list of fitted model objects.
#'   When \code{runModels = FALSE}: data frame of model specifications.
#' @export
do.ds <-
  function(data,
           key = c("unif", "hn", "hr"),
           covars = NULL,
           runModels = TRUE,
           parallel = FALSE,
           nCores = 2,
           folder = "",
           logfile = "",
           rerun = FALSE,
           verbose = FALSE,
           models = NULL,
           incl.adj = TRUE,
           cleanFolder = FALSE,
           ...) {

    # Models pre-supplied?
    if (!is.null(dim(models)[1]))
      modDat <- models
    else
      # generate models but don't add covars to "unif" key
      modDat <- purrr::map_dfr(key, function(k, covars) {
        if (k == "unif")
          gen.form.N(0, k, NULL)
        else
          purrr::map_dfr(0:length(covars), gen.form.N, k, covars)
      }, covars = covars)


    # append the adjustment only models
    if (incl.adj)
      modDat <- dplyr::bind_rows(adj.models, modDat)

    # setup logfile connection. If it's "" (the default) then logfileConn will be
    # "", which corresponds to stdout.
    if (logfile != "") {
      logfileConn <- file(logfile, "at") # open for append.
      cat(paste0("Logging info to ", logfile, "\n"))
    } else
      logfileConn <- logfile

    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

    cat(paste0("\n\ndoDS: ", date(), "\n"), file = logfileConn)

    if (parallel && logfile != "" && runModels)
      cat("WARNING: No logging to file for individual sub-processes when parallel == TRUE.\n", file = logfileConn)

    # only run models that currently have no results in folder. Useful when a
    # previous run has failed partway through and we want to avoid re-running a
    # bunch of models for which we already have results.
    if (folder != "") {
      if (runModels) cat(paste0("Saving output files in: ", folder, "\n"), file = logfileConn)

      # Clean the results folder?
      if (cleanFolder) {
        fl <- list.files(path = folder, full.names = TRUE)
        message(sprintf(
          "Cleaning results folder '%s'. Removing %d files.",
          folder,
          length(fl)
        ))
        file.remove(fl)
      }


      # Get list of models that already have result files
      if (!rerun) {
        fl <- sub(".txt", "",
                  list.files(path = folder,
                             pattern = ".*\\.txt",
                             full.names = T), fixed = T)
        existMods <- unlist(lapply(strsplit(basename(fl), "_"), function(x) x[3]))
        cat(paste0(sum(modDat$label %in% existMods), " of ", nrow(modDat), " models have already been run.\n"), file = logfileConn)
        modDat <- modDat[!(modDat$label %in% existMods),]
      }
    }

    if (runModels) {

      cat(paste0("Running ", nrow(modDat), " model(s).\n"), file = logfileConn)
      # Avoids setting up parallel cluster if not needed
      if(nrow(modDat) > 0) {
        if (parallel) {
          cat(paste0("Using parallel processing with ", nCores, " cores.\n"), file = logfileConn)
          cl <- parallel::makeCluster(min(nCores, nrow(modDat)), type = "SOCK")
          doSNOW::registerDoSNOW(cl)

          models <- plyr::dlply(
            modDat,
            plyr::.(label),
            .parallel = T,
            .paropts = list(
              .packages = "Distance",
              .verbose = TRUE,
              .errorhandling = "pass" # Don't think this makes any difference
            ),
            run.ddf.model,
            data = data,
            folder = folder,
            logfileConn = "",
            verbose = verbose,
            ...
          )
          parallel::stopCluster(cl)

          # Tried but had problems with some things not defined...
          # future::plan(multisession, workers = nCores)
          # models <- modDat %>%
          #   split(1:nrow(.)) %>%
          #   furrr::future_map(
          #     run.ddf.model,
          #     data = data,
          #     folder = folder,
          #     logfileConn = logfileConn,
          #     verbose = verbose,
          #     ...
          #   )
        } else { # parallel == FALSE
          models <- modDat %>%
            split(1:nrow(.)) %>%
            purrr::map(
              run.ddf.model,
              data = data,
              folder = folder,
              logfileConn = logfileConn,
              verbose = verbose,
              .progress = "Detection Function Progress",
              ...
            )

          # models <- plyr::dlply(modDat, .(label), .progress = progress_win(title = "Detection Function Progress"), run.ddf.model, data = data,
          #                       folder = folder,logfileConn = logfileConn, verbose = verbose, ...)
        }
      } else { # nrow(modDat) > 0
        models <- list()
      }

      cat(paste0("do.ds: finishing at ", date(), "\n"), file = logfileConn)

      if (logfile != "")
        close(logfileConn)

      return(models)
    } else { # runmodels == FALSE
      if (logfile != "")
        close(logfileConn)

      return(modDat)
    }
  }


#' Strip non-numeric characters and coerce to numeric
#'
#' Removes all characters not in \code{keep} and converts to numeric, without
#' generating the usual \code{NA}-coercion warning.
#'
#' @param x Character vector to convert.
#' @param keep Regular expression character class of characters to retain;
#'   default keeps digits, decimal point, and sign characters.
#' @return Numeric vector.
#' @export
destring <- function(x,keep="0-9.e+-") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep = ""),"",x)) )
}

#' Validate a fitted detection function and augment distdata
#'
#' Prints a model summary, plots the detection function, runs GOF tests, and
#' augments the distdata with \code{detProb} and \code{adjSize} columns.
#'
#' @param model Fitted \code{dsmodel} or \code{fake_ddf} object.
#' @param species Character string species code; used in plot titles.
#' @param modname Character string model name; used in plot titles.
#' @return Distdata data frame augmented with \code{detProb} and
#'   \code{adjSize} columns.
#' @export
check.det.fcn <-
  function(model = NULL,
           species = NULL,
           modname = NULL) {

    message("Running check.det.fcn...\n")
    print(summary(model))

    mod.dat <- df.get.data(model) # returns a list
    distdata <- mod.dat$distdata
    form <- mod.dat$form
    key <- mod.dat$key
    aic <- mod.dat$aic

    if (nrow(distdata) > 0) {
      # Note that fitted is in the right order since df.get.data uses:
      # fitted(dfobject$ddf)[as.character(distdata$object)] so we don't need to
      # do that here.
      distdata$detProb <- mod.dat$fitted
      distdata$adjSize <- distdata$size / distdata$detProb

      # make sure nothing went wrong
      if (any(is.na(distdata$detProb)))
        warning("Some detection probabilities are NA!", immediate. = T)

      cat("\n")
      cat(sprintf("Detection prob range: %s\n", paste(round(
        range(distdata$detProb), 4
      ), collapse = " - ")))
      cat(sprintf("Number of detection probs < 0.15: %d\n", nrow(dplyr::filter(
        distdata, detProb < .15
      ))))
      cat(sprintf("range of size: %s\n", paste(range(distdata$size), collapse = " - ")))
      cat(sprintf("range of adjusted size: %s\n", paste(round(
        range(distdata$adjSize), 2
      ), collapse = " - ")))

      # plot hist of det probs
      p <- ggplot2::ggplot(data = distdata, ggplot2::aes(x = detProb))
      p <- p + ggplot2::geom_histogram(binwidth = 0.1)
      p <- p + ggplot2::scale_x_continuous(breaks = seq(0, 1, .1))
      print(p)
    }

    # GOF testing
    message("GOF testing")
    if (!("fake_ddf" %in% class(model))){
      print(mrds::ddf.gof(model$ddf, asp = 1))
      # plot model
      plot(
        model,
        main = sprintf(
          "%s, %s, %s, %s AIC: %.3f - overall detection function",
          species,
          modname,
          key,
          form,
          aic
        )
      )
    } else
      message("Dummy ddf - no GOF or plot!")

    message("End of check.det.fcn\n")
    distdata
  }

#' Find segments with observations from more than one detection function type
#'
#' \strong{Note: needs updating} to use the current \code{ddftype} column
#' (currently references the legacy \code{det.fcn.type} column from the
#' subArctic analysis).
#'
#' @param distdata Observation data frame with a \code{det.fcn.type} column.
#' @param segdata Segment data frame with a \code{Sample.Label} column.
#' @return Data frame of segments with mixed detection function types.
#' @export
get.trouble <- function(distdata, segdata) {
  distdata %>%
    dplyr::group_by(Sample.Label) %>%
    dplyr::summarize(ndet = dplyr::n_distinct(det.fcn.type)) %>%
    dplyr::filter(ndet > 1) %>%
    dplyr::left_join(segdata, by = "Sample.Label") %>%
    dplyr::left_join(distdata, by = "Sample.Label") %>%
    dplyr::select(Sample.Label, ndet, det.fcn.type, df.type, distance, size, TransectType, Observer, object) %>%
    as.data.frame()
}


#' Extract distdata, formula, key, AIC, and fitted values from a DDF object
#'
#' Handles \code{dsmodel} objects and \code{fake_ddf} strip-transect objects.
#'
#' @param dfobject A \code{dsmodel} or \code{fake_ddf} object.
#' @return Named list with elements \code{distdata}, \code{form},
#'   \code{key}, \code{aic}, and \code{fitted}.
#' @export
df.get.data <- function(dfobject) {


  if ("dsmodel" %in% class(dfobject)) {
    distdata <- dfobject$ddf$data
    list(
      distdata = distdata,
      form = as.character(dfobject$ddf$call$dsmodel[[2]]$formula)[2],
      key = dfobject$ddf$call$dsmodel[[2]]$key,
      aic = dfobject$ddf$criterion,
      fitted = fitted(dfobject$ddf)[as.character(distdata$object)]
    )
  } else if ("fake_ddf" %in% class(dfobject)) {
    list(
      distdata = dfobject$data,
      form = "~1",
      key = "strip_transect",
      aic = NA,
      fitted = 1 # fitted det prob
    )
  } else if ("list" %in% class(dfobject)) {
    # Not sure what this case is for. Put stop() here to see if it's ever called
    dfobject$distdata
    stop("df.get.data: Don't know how to extract form, key, or aic from this object")
  }


  # distdata <- df.get.data(model)
  # form <- as.character(model$ddf$call$dsmodel[[2]]$formula)[2]
  # key <- model$ddf$call$dsmodel[[2]]$key
  # aic <- model$ddf$criterion

}

#' Load and check the top-ranked detection function model from a folder
#'
#' Reads the \code{.RData} file with the lowest AIC from \code{fold} and
#' calls \code{\link{check.det.fcn}} on it.  If \code{list.only = TRUE},
#' only prints the filename.
#'
#' \strong{Note: may need argument changes} — \code{behav} parameter may be
#' obsolete and \code{check.det.fcn} arguments may need updating.
#'
#' @param fold Directory path containing candidate model \code{.RData} files.
#' @param species Character string species code.
#' @param behav Character string behaviour label used in messages.
#' @param list.only If \code{TRUE}, only print the top model filename without
#'   loading or checking it.
#' @return \code{invisible(NULL)}.
#' @export
check.top.model <- function(fold, species, behav, list.only = FALSE) {
  files <- list.files(fold, ".+\\.RData$", full.names = TRUE)

  if (list.only) {
    message(sprintf(
      "Top model for %s %s %s is %s",
      basename(fold),
      species,
      behav,
      basename(files[1])
    ))
    return()
  }

  if (files[1] == "") {
    message(sprintf(
      "There were no model results for %s %s %s to load",
      basename(fold),
      species,
      behav
    ))
    return()
  }

  message(sprintf("Loading model results from %s", files[1]))
  load(files[1])
  try(check.det.fcn(model = model, paste(species, behav)))
}

#' Check all top detection function models for a species and behaviour
#'
#' \strong{Note: may need updating} — folder naming convention has changed.
#'
#' @param species Character string species code.
#' @param behav Character string behaviour label (\code{"F"} or \code{"W"}).
#' @param ... Additional arguments passed to \code{\link{check.top.model}}.
#' @return \code{invisible(NULL)}.
#' @export
check.top.models <- function(species, behav, ...){
  folder <- here::here("R", species, "DF Summaries", behav)

  if (!dir.exists(folder)) {
    message(sprintf("Folder %s does not exist!", folder))
    return()
  }

  dirs <- list.dirs(folder, recursive = FALSE)
  purrr::walk(dirs, check.top.model, species, behav, ...)
}


#' Produce a one-row summary data frame for a DSM model
#'
#' @param model A fitted \code{dsm} object, or a \code{try-error} if model
#'   fitting failed.
#' @return Single-row data frame with columns \code{response}, \code{terms},
#'   \code{AIC}, \code{REML}, \code{OverDisp}, and \code{Deviance_explained}.
#' @export
summarize.dsm <- function(model){

  # Check if model had error when it ran
  if (inherits(model, "try-error")) {
    data.frame(
      response = "Error",
      terms    = NA,
      AIC      = NA,
      REML     = NA,
      OverDisp = NA,
      "Deviance_explained" = NA
    )
  } else {
    summ <- summary(model)

    data.frame(
      response = model$family$family,
      terms    = paste(c(
        rownames(summ$pTerms.table), rownames(summ$s.table)
      ), collapse = ", "),
      AIC      = AIC(model),
      REML     = model$gcv.ubre,
      OverDisp = sum(resid(model, type = "pearson")) / model$df.residual,
      "Deviance_explained" = paste0(round(summ$dev.expl * 100, 2), "%")
    )
  }
}



# Response families the DSM model selection code knows how to rebuild and to
# simulate from. Keyed by get.family.key(). An environment rather than a literal
# so a SubProject can add a family without a package release - see
# register.dsm.family().
.dsm_families <- new.env(parent = emptyenv())

#' Register a response family for DSM model selection
#'
#' The selection code is agnostic to how many response families are in play:
#' \code{\link{get.family.finalists}} takes the AIC-best candidate in each family
#' it finds, and \code{\link{run.family.cv}} scores all of them against the same
#' folds. The only family-specific knowledge it needs is how to rebuild a family
#' object and how to draw from its predictive distribution, and that lives here.
#'
#' Adding a third or fourth family to \code{dsm.mod.specs} therefore takes one
#' call to this function -- from \code{analysis_settings.R} if the family is
#' specific to one SubProject -- rather than an edit inside the package.
#'
#' A constructor is required rather than reusing the family object supplied with
#' the model spec: mgcv's extended families carry mutable state in their
#' closures, so sharing one object across fits couples them and lets a failed fit
#' poison later ones. See \code{\link{fresh.family}}.
#'
#' @param key Family key, as \code{\link{get.family.key}} returns it, e.g.
#'   \code{"tw"}. Registering an existing key replaces it.
#' @param constructor Function of no arguments returning a newly built family
#'   object, e.g. \code{function() mgcv::tw()}.
#' @param simulate Function of \code{(fit, mu, n_sim)} returning
#'   \code{length(mu) * n_sim} draws from the fitted predictive distribution,
#'   varying fastest over \code{mu}. \code{fit} is the fitted model, so
#'   family parameters can be read from \code{fit$family} and \code{fit$sig2}.
#' @return Invisibly, the key.
#' @examples
#' \dontrun{
#' # a Poisson candidate added to dsm.mod.specs for one SubProject
#' register.dsm.family(
#'   "poisson",
#'   constructor = function() stats::poisson(),
#'   simulate = function(fit, mu, n_sim)
#'     stats::rpois(length(mu) * n_sim, lambda = rep(mu, n_sim)))
#' }
#' @seealso \code{\link{registered.dsm.families}}, \code{\link{fresh.family}},
#'   \code{\link{sim.response}}
#' @export
register.dsm.family <- function(key, constructor, simulate) {
  checkmate::expect_string(key, min.chars = 1)
  checkmate::expect_function(constructor, nargs = 0)
  checkmate::expect_function(simulate, args = c("fit", "mu", "n_sim"))

  assign(key, list(constructor = constructor, simulate = simulate),
         envir = .dsm_families)
  invisible(key)
}

#' Keys of the response families currently registered
#'
#' @return Character vector of keys, sorted.
#' @examples
#' registered.dsm.families()
#' @seealso \code{\link{register.dsm.family}}
#' @export
registered.dsm.families <- function() sort(ls(.dsm_families))

#' Look up a registered family, or fail with a useful message
#'
#' @param key Family key from \code{\link{get.family.key}}.
#' @param what Caller context, used to prefix the error.
#' @return The registry entry: a list of \code{constructor} and \code{simulate}.
#' @noRd
get.dsm.family <- function(key, what) {
  checkmate::expect_string(key)

  if (!exists(key, envir = .dsm_families, inherits = FALSE))
    stop(what, " is not registered. Registered: ",
         paste(registered.dsm.families(), collapse = ", "),
         ". Add it with register.dsm.family(\"", key,
         "\", constructor, simulate) - see ?register.dsm.family. Guessing a ",
         "constructor here would silently change the model, and guessing a ",
         "simulator would silently invalidate every cross-validation score.")

  get(key, envir = .dsm_families, inherits = FALSE)
}

.onLoad <- function(libname, pkgname) {
  # The two families dsm.mod.specs ships with. Anything else is registered by
  # the caller; nothing in the selection code assumes this pair.
  register.dsm.family(
    "tw",
    constructor = function() mgcv::tw(),
    # getTheta(TRUE) returns p for tw(); the scale lives in $sig2.
    simulate = function(fit, mu, n_sim)
      mgcv::rTweedie(rep(mu, n_sim), p = fit$family$getTheta(TRUE),
                     phi = fit$sig2))

  register.dsm.family(
    "nb",
    constructor = function() mgcv::nb(),
    simulate = function(fit, mu, n_sim)
      stats::rnbinom(length(mu) * n_sim, mu = rep(mu, n_sim),
                     size = fit$family$getTheta(TRUE)))
}


#' Short, stable key for a response family
#'
#' Normalises a family object's name so that fitted and unfitted objects agree
#' (\code{"Tweedie"} and \code{"Tweedie(p=1.167)"} both give \code{"tw"}), and so
#' that families beyond the two currently in \code{dsm.mod.specs} get a key of
#' their own rather than an error. Nothing downstream should assume the set of
#' keys is \code{c("tw", "nb")}; the families the selection code can rebuild and
#' simulate from are whatever \code{\link{register.dsm.family}} has been given.
#'
#' @param fam A family object, or the character name of one.
#' @return Length-one character key.
#' @examples
#' get.family.key(mgcv::tw())
#' get.family.key("Negative Binomial(0.017)")
#' @export
get.family.key <- function(fam) {
  nm <- if (is.character(fam)) fam else fam$family
  checkmate::expect_character(nm, len = 1, any.missing = FALSE)

  # Fitted families carry their estimated parameter, e.g. "Tweedie(p=1.167)".
  base <- tolower(trimws(sub("\\(.*$", "", nm)))

  switch(base,
         "tweedie"           = "tw",
         "negative binomial" = "nb",
         gsub("[^a-z0-9]+", "", base))
}


#' Pick the AIC-best candidate within each response family
#'
#' Step one of the two-step model selection. AIC cannot compare models across
#' response families when the response is non-integer -- a Tweedie likelihood is
#' a density and a negative binomial likelihood is a probability mass, so the two
#' are not on a common scale -- but it is entirely valid *within* a family, where
#' every candidate shares a response and a likelihood. Step two is
#' \code{\link{run.family.cv}}, which chooses between the finalists this returns.
#'
#' Works with whatever \code{mod.specs} contains: any number of candidates, any
#' number of families, any naming scheme. The family is read from each spec's
#' family object rather than parsed out of the model name.
#'
#' @param mod.res Named list of fitted DSMs as saved by \code{run.dsm.models()}.
#'   Entries that are \code{try-error}s, or that are missing from
#'   \code{mod.specs}, are dropped with a message.
#' @param mod.specs Data frame of candidate specifications with columns
#'   \code{modname}, \code{formula} and \code{family} -- i.e. \code{dsm.mod.specs}.
#' @return A tibble with one row per family present, columns \code{family_key},
#'   \code{modname}, \code{AIC}, \code{n_candidates} (how many candidates that
#'   family contributed), and list-columns \code{formula} (taken from the fitted
#'   model, so it carries the \code{offset(off.set)} term \code{dsm()} adds) and
#'   \code{family}. Ordered by \code{AIC}. Carries an \code{"all_candidates"}
#'   attribute holding every usable candidate in the same shape.
#' @examples
#' \dontrun{
#' load(here(RDataDir, "ATPU.dsm.RData"))
#' get.family.finalists(mod.res, dsm.mod.specs)
#' }
#' @export
get.family.finalists <- function(mod.res, mod.specs) {
  checkmate::expect_list(mod.res, min.len = 1)
  checkmate::expect_data_frame(mod.specs, min.rows = 1)
  if (!all(c("modname", "formula", "family") %in% names(mod.specs)))
    stop("get.family.finalists: mod.specs needs modname, formula and family columns")

  rows <- purrr::map(seq_len(nrow(mod.specs)), function(i) {
    modname <- mod.specs$modname[[i]]
    model   <- mod.res[[modname]]

    if (is.null(model)) {
      message(sprintf("get.family.finalists: no fitted model for %s, skipping", modname))
      return(NULL)
    }
    if (inherits(model, "try-error")) {
      message(sprintf("get.family.finalists: %s is a try-error, skipping", modname))
      return(NULL)
    }

    tibble::tibble(
      family_key = get.family.key(mod.specs$family[[i]]),
      modname    = modname,
      AIC        = summarize.dsm(model)$AIC,
      # The FITTED model's formula, not the spec's: dsm() appends
      # offset(off.set), and refitting without it drops the effort/area exposure
      # correction entirely - which makes nb() fail to converge on sparse data
      # and sends the predictions off by orders of magnitude.
      formula    = list(stats::formula(model)),
      family     = list(mod.specs$family[[i]]))
  })

  cand <- dplyr::bind_rows(rows)
  if (nrow(cand) == 0)
    stop("get.family.finalists: no usable fitted candidates")

  finalists <- cand %>%
    dplyr::group_by(.data$family_key) %>%
    dplyr::mutate(n_candidates = dplyr::n()) %>%
    dplyr::slice_min(.data$AIC, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$AIC)

  # Every candidate, so a matched-formulation comparison can be built later
  # without reloading the (multi-GB) model suite.
  attr(finalists, "all_candidates") <- cand
  finalists
}


#' Assign segments to spatial-block cross-validation folds
#'
#' A thin wrapper on \code{\link[blockCV]{cv_spatial}}. Squares of
#' \code{block_size} on the analysis projection are assigned to folds at random,
#' so each fold holds out whole regions rather than scattered segments. Random
#' hold-outs would leak neighbouring segments into training and flatter every
#' model equally, because the residuals are spatially structured.
#'
#' Blocks are square rather than blockCV's default hexagon so that
#' \code{block_size} keeps meaning an edge length, which is what
#' \code{dsm.options$family.cv.block.size.m} and everything written about it
#' assume.
#'
#' Fold membership depends on the coordinates, the block size, the seed and the
#' \strong{blockCV version} - the 4.0-0 release notes state that folds for a
#' given seed may differ from earlier versions. Record the version alongside the
#' seed in anything that has to be reproducible; \code{\link{run.family.cv}}'s
#' caller does this in the results file header.
#'
#' @param segdata Data frame with numeric \code{x} and \code{y} columns in
#'   projection units (metres).
#' @param n_folds Number of folds.
#' @param block_size Block edge length in projection units.
#' @param seed Random seed for block-to-fold assignment. The random number
#'   generator kind is pinned to Mersenne-Twister for the duration and restored
#'   afterwards, so the folds are the same whether or not the caller is running
#'   inside a \code{future} (which switches to L'Ecuyer-CMRG). blockCV calls
#'   \code{set.seed()} itself, and \code{set.seed()} resets whichever generator
#'   is current, so the pinning is still needed with the package doing the work.
#' @param crs Coordinate reference system of \code{x} and \code{y}, in any form
#'   \code{sf::st_crs()} accepts. Must be projected: blockCV interprets
#'   \code{block_size} as metres, and for a geographic CRS it divides by
#'   \code{deg_to_metre} instead, which silently produces blocks of the wrong
#'   size. Rejected rather than converted.
#' @param balance_column Optional column of \code{segdata} whose classes should
#'   be balanced across folds as well as the record count, e.g.
#'   \code{"SurveyType"}. \code{NULL} balances record counts only.
#'
#' @section Why blockCV rather than a grid of our own:
#' This was fifteen lines of \code{floor(x / block_size)} until issue #8. The
#' construction is the same either way; what the package adds is fold balancing
#' over repeated random draws, and a citable provenance (Valavi et al. 2019,
#' Methods in Ecology & Evolution 10:225-232) for a step a reader would otherwise
#' have to take on trust. Measured on \code{Atl IMRP} ATPU (468,424 segments,
#' 100 km blocks, 5 folds), balancing takes the largest-to-smallest fold ratio
#' from 1.36 to 1.11, and balancing on \code{SurveyType} takes the aerial share
#' of a fold from a 0.6-15.0\% spread to 3.2-8.4\%.
#'
#' @return Integer vector of fold membership, length \code{nrow(segdata)}, with
#'   an \code{"n_blocks"} attribute giving the number of occupied blocks. That
#'   count is worth reporting: it is the real sample size behind the folds, and a
#'   block size close to the study area extent can leave too few blocks to spread
#'   across \code{n_folds}.
#' @examples
#' \dontrun{
#' assign.blocks(data.frame(x = runif(100, 0, 5e5), y = runif(100, 0, 5e5)),
#'               n_folds = 5, block_size = 1e5, seed = 1, crs = 3347)
#' }
#' @export
assign.blocks <- function(segdata, n_folds, block_size, seed, crs,
                          balance_column = NULL) {
  checkmate::expect_data_frame(segdata, min.rows = 1)
  checkmate::expect_numeric(segdata$x, any.missing = FALSE)
  checkmate::expect_numeric(segdata$y, any.missing = FALSE)
  checkmate::expect_count(n_folds, positive = TRUE)
  checkmate::expect_number(block_size, lower = 0)
  checkmate::expect_string(balance_column, null.ok = TRUE)
  if (!is.null(balance_column))
    checkmate::expect_names(names(segdata), must.include = balance_column)

  # A geographic CRS does not fail, it quietly rescales block_size by
  # deg_to_metre and returns blocks of the wrong size. Refuse it.
  crs <- sf::st_crs(crs)
  if (is.na(crs))
    stop("assign.blocks: crs is missing or unrecognised. blockCV needs a ",
         "projected CRS in metres to interpret block_size.")
  if (isTRUE(sf::st_is_longlat(crs)))
    stop("assign.blocks: crs is geographic (long/lat). block_size is metres, ",
         "and blockCV would rescale it by deg_to_metre. Pass the projected ",
         "analysis CRS (segProj) instead.")

  cols <- c("x", "y", balance_column)
  pts  <- sf::st_as_sf(segdata[, cols, drop = FALSE], coords = c("x", "y"),
                       crs = crs)

  # Pin the generator, not just the seed. Under furrr's seed = TRUE the workers
  # run L'Ecuyer-CMRG, and set.seed() on a different generator yields a
  # different sequence - which silently produced different folds inside workers
  # than in the main process, from an identical seed. blockCV seeds itself, but
  # set.seed() resets the CURRENT generator, so this still applies.
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  blocks <- blockCV::cv_spatial(
    x         = pts,
    column    = balance_column,
    size      = block_size,
    k         = n_folds,
    hexagon   = FALSE,          # squares: block_size stays an edge length
    selection = "random",
    # 100 attempts at an even split, against one draw before. Costs ~90 s at
    # 468,424 points against ~4 s unbalanced (Atl IMRP, 80 cores) - which is
    # why run.family.cv() takes a precomputed folds argument.
    balance   = TRUE,
    iteration = 100L,
    seed      = seed,
    biomod2   = FALSE,          # we do not use the biomod2 fold matrix
    plot      = FALSE,
    report    = FALSE,
    progress  = FALSE)

  # blockCV builds the grid to the data extent, so the block count can differ by
  # one or two from a grid anchored at the projection origin (196 against 197 on
  # Atl IMRP). Neither anchoring is more correct; the count is reported so the
  # difference is visible rather than surprising.
  ids <- blocks$folds_ids
  if (anyNA(ids))
    stop("assign.blocks: blockCV left ", sum(is.na(ids)), " segment(s) outside ",
         "every block. Increase cv_spatial()'s extend argument.")

  structure(as.integer(ids), n_blocks = nrow(blocks$blocks))
}



#' Continuous ranked probability score from predictive samples
#'
#' CRPS is a proper scoring rule evaluated on the response scale, so unlike a
#' log-score it carries no dominating-measure baggage and is defined for discrete
#' and continuous predictive distributions alike. That is what makes it usable to
#' compare a negative binomial against a Tweedie.
#'
#' \code{CRPS(F, y) = E|X - y| - 0.5 * E|X - X'|}, with the second expectation
#' evaluated by the sorted-sample identity
#' \code{E|X - X'| = (2 / m^2) * sum_i (2i - m - 1) * x_(i)}, which is
#' O(m log m) per observation rather than the O(m^2) all-pairs form.
#'
#' @param sim_mat Numeric matrix of predictive draws, one row per observation and
#'   one column per draw.
#' @param obs Numeric vector of observed values, length \code{nrow(sim_mat)}.
#' @return Numeric vector of per-observation CRPS. Lower is better.
#' @examples
#' calc.crps(matrix(rnorm(200), nrow = 2), c(0, 1))
#' @export
calc.crps <- function(sim_mat, obs) {
  checkmate::expect_matrix(sim_mat, mode = "numeric")
  checkmate::expect_numeric(obs, len = nrow(sim_mat), any.missing = FALSE)

  n_draw <- ncol(sim_mat)
  term_1 <- rowMeans(abs(sim_mat - obs))
  sorted <- t(apply(sim_mat, 1, sort))
  term_2 <- as.vector(sorted %*% (2 * seq_len(n_draw) - n_draw - 1)) * (2 / n_draw^2)

  term_1 - 0.5 * term_2
}


#' Randomised probability integral transform from predictive samples
#'
#' Uniform on (0, 1) when the predictive distribution is correct. The random
#' tie-break makes this valid for discrete predictive distributions and, unlike
#' DHARMa's quantile residuals, applies identical treatment to every family --
#' DHARMa decides whether to jitter from the declared family, which makes its
#' statistic a poor basis for choosing *between* families.
#'
#' @param sim_mat Numeric matrix of predictive draws, one row per observation.
#' @param obs Numeric vector of observed values, length \code{nrow(sim_mat)}.
#' @return Numeric vector of PIT values in [0, 1].
#' @examples
#' calc.pit(matrix(rpois(200, 3), nrow = 2), c(2, 4))
#' @export
calc.pit <- function(sim_mat, obs) {
  checkmate::expect_matrix(sim_mat, mode = "numeric")
  checkmate::expect_numeric(obs, len = nrow(sim_mat), any.missing = FALSE)

  n_below <- rowSums(sim_mat <  obs)
  n_equal <- rowSums(sim_mat == obs)

  (n_below + stats::runif(length(obs)) * n_equal) / ncol(sim_mat)
}


#' Draw from a fitted model's predictive distribution
#'
#' Simulates new responses at supplied fitted means using the fitted family's own
#' parameters. Add a branch here when a new family is added to
#' \code{dsm.mod.specs}; the error is deliberate rather than a silent fallback,
#' because a wrong predictive distribution would quietly corrupt every score.
#'
#' @param fit Fitted \code{bam}/\code{gam} object.
#' @param mu Numeric vector of fitted means on the response scale.
#' @param n_sim Number of draws per element of \code{mu}.
#' @return Numeric matrix, \code{length(mu)} rows by \code{n_sim} columns.
#' @examples
#' \dontrun{
#' sim.response(fit, predict(fit, type = "response"), 200)
#' }
#' @export
sim.response <- function(fit, mu, n_sim) {
  checkmate::expect_class(fit, "gam")
  checkmate::expect_numeric(mu, any.missing = FALSE)
  checkmate::expect_count(n_sim, positive = TRUE)

  key  <- get.family.key(fit$family)
  spec <- get.dsm.family(key,
                         what = paste0("sim.response: family '",
                                       fit$family$family, "'"))

  draws <- spec$simulate(fit, mu, n_sim)
  if (length(draws) != length(mu) * n_sim)
    stop("sim.response: the simulator registered for '", key, "' returned ",
         length(draws), " draws, expected ", length(mu) * n_sim)

  matrix(draws, nrow = length(mu))
}


#' A fresh, unused copy of a response family
#'
#' mgcv's extended families (\code{tw()}, \code{nb()}) carry mutable state in
#' their closures' environment: fitting writes the estimated parameter back into
#' the object it was handed. Reusing one object across fits therefore couples
#' them, and a *failed* fit can leave it in a state that makes every later fit
#' using it fail too.
#'
#' That is not hypothetical. \code{dsm.mod.specs} holds one family object per
#' candidate, and \code{\link{get.family.finalists}} hands back a reference to
#' it, so every species and every fold in a run share the same \code{nb()}. In a
#' sequential run over eleven species the first species failed three folds and
#' every subsequent species then failed all five, 0 for 50.
#'
#' Rebuilds from the family key using the constructor registered for it, so it
#' covers whatever has been passed to \code{\link{register.dsm.family}} rather
#' than a fixed pair. A family carrying non-default arguments will not survive
#' the round trip - hence the explicit error for anything unregistered, rather
#' than a silent fallback that would quietly change the model.
#'
#' @param fam A family object.
#' @return A newly constructed family object of the same kind.
#' @examples
#' f <- mgcv::nb()
#' identical(environment(f$getTheta), environment(fresh.family(f)$getTheta))
#' @export
fresh.family <- function(fam) {
  key <- get.family.key(fam)
  spec <- get.dsm.family(key,
                         what = paste0("fresh.family: family '",
                                       if (is.character(fam)) fam else fam$family,
                                       "'"))
  spec$constructor()
}


#' Spatial-block cross-validation of one species' family finalists
#'
#' Step two of the two-step model selection. Refits each finalist from
#' \code{\link{get.family.finalists}} on every spatial-block training set and
#' scores it on the held-out block, using criteria that compare across response
#' families because they live on the response scale.
#'
#' Every finalist sees identical folds, so the comparison is paired -- but a fold
#' can still be dropped for one family and not another when its fit fails or
#' diverges, which breaks that pairing. The caller is responsible for excluding
#' such folds from every family before averaging. Models are
#' refitted with \code{mgcv::bam()} rather than through \code{dsm()} because the
#' stored model already carries everything needed: \code{$data} holds the segment
#' data including the \code{off.set} offset column, and \code{formula()} carries
#' the offset term, so there is no detection function to reattach.
#'
#' @param species Character species/group label, copied into the output.
#' @param segdata Segment data from a fitted model's \code{$data}. Needs
#'   \code{x}, \code{y} and the response named in the finalists' formulae.
#' @param finalists Tibble from \code{\link{get.family.finalists}}.
#' @param n_folds,block_size,n_sim Cross-validation geometry: number of folds,
#'   block edge length in projection units, predictive draws per held-out
#'   segment.
#' @param seed Seed for block assignment; predictive draws use \code{seed +
#'   fold}, so every finalist sees the same draws within a fold.
#' @param crs Projected CRS of \code{segdata}'s coordinates, passed to
#'   \code{\link{assign.blocks}}. Required unless \code{folds} is supplied.
#' @param balance_column Optional \code{segdata} column whose classes are
#'   balanced across folds, passed to \code{\link{assign.blocks}}.
#' @param folds Optional precomputed fold vector, one entry per row of
#'   \code{segdata} and matched to it by position. Supply it to avoid repeating
#'   the balancing search across calls that share segment data; \code{NULL}
#'   computes it here.
#' @param nthreads Passed to \code{mgcv::bam()}.
#' @param max.pred.ratio How many times the largest observed value a fold's
#'   largest prediction may reach before the fit is called diverged and the fold
#'   is dropped. Guards the simulation step: \code{mgcv::rTweedie} draws a
#'   Poisson count per observation and then that many gamma variates, so absurd
#'   fitted values ask for an absurd vector and the error is an allocation
#'   failure rather than anything catchable-looking. Complements
#'   \code{dsm.options$family.cv.diverged.ratio}, which screens scores after the
#'   fact; this screens predictions before a score exists at all.
#' @param allow.slow.refit If a fold still fails after a single-threaded retry,
#'   whether to refit it with \code{discrete = FALSE}. Correct but very
#'   expensive: roughly 280x the discrete path (474 min against 1.7 min per fit
#'   on \code{NL_EXPL_DRL_RA}, 153,344 segments, 80 cores). Treat the ratio as
#'   the transferable number; either way a handful of folds can dominate an
#'   entire run. Defaults to \code{FALSE}, in which case the fold is dropped and
#'   reported.
#' @param calib.covar Column in \code{segdata} to aggregate the
#'   observed-vs-expected calibration check over, or \code{NULL} to skip it.
#'   Defaults to \code{"platform"}, which is present in the segment data whether
#'   or not the formula uses it -- factor-smooth and no-factor formulations still
#'   get a calibration score.
#' @param progress If \code{TRUE}, report each fold as it completes.
#' @return A tibble, one row per finalist per fold: \code{species},
#'   \code{family_key}, \code{modname}, \code{fold}, \code{n_test},
#'   \code{n_blocks}, \code{CRPS},
#'   \code{MAE}, \code{RMSE}, \code{PIT_KS}, \code{cover90}, \code{calib} (mean
#'   |1 - observed/expected| over the levels of \code{calib.covar}), one
#'   \code{OE_<level>} column per level, and \code{mins}.
#' @examples
#' \dontrun{
#' run.family.cv("ATPU", dsm_data, finalists, n_folds = 5, seed = 20260820)
#' }
#' @export
run.family.cv <- function(species, segdata, finalists,
                          n_folds = 5, block_size = 1e5, n_sim = 200,
                          seed = 1, crs = NULL, balance_column = NULL,
                          folds = NULL, nthreads = 1, calib.covar = "platform",
                          allow.slow.refit = FALSE, max.pred.ratio = 1e4,
                          progress = TRUE) {
  checkmate::expect_character(species, len = 1, any.missing = FALSE)
  checkmate::expect_data_frame(segdata, min.rows = 1)
  checkmate::expect_data_frame(finalists, min.rows = 1)
  checkmate::expect_number(max.pred.ratio, lower = 1)
  if (!all(c("family_key", "modname", "formula", "family") %in% names(finalists)))
    stop("run.family.cv: finalists must come from get.family.finalists()")

  # Reuse a precomputed assignment when the caller has one. Balancing costs
  # ~90 s at 468,424 points (Atl IMRP, 80 cores), and the caller runs this
  # function twice per species - once on the finalists, once on the matched
  # formulation - over identical segdata.
  fold_id <- if (is.null(folds)) {
    assign.blocks(segdata, n_folds, block_size, seed, crs, balance_column)
  } else {
    # A fold vector is matched to segdata BY POSITION, so a caller passing one
    # built from different rows would silently score the wrong segments.
    checkmate::expect_integerish(folds, len = nrow(segdata), any.missing = FALSE)
    folds
  }
  if (isTRUE(progress))
    message(sprintf("%s: %d segments, %d occupied blocks, %d folds, sizes %s",
                    species, nrow(segdata), attr(fold_id, "n_blocks"),
                    length(unique(fold_id)),
                    paste(table(fold_id), collapse = "/")))

  resp <- as.character(finalists$formula[[1]][[2]])

  out <- NULL
  for (i in seq_len(nrow(finalists))) {
    for (k in seq_len(n_folds)) {
      train_idx <- which(fold_id != k)
      test_idx  <- which(fold_id == k)
      started   <- Sys.time()

      # discrete = TRUE is much faster but is the fragile path: on sparse data the
      # nb() fits sit close to the edge of convergence and fail intermittently
      # under parallel load, with "Error in if (sum(uconv))". Fall back to the
      # stable path rather than dropping the fold, which would silently turn a
      # two-family comparison into a one-family one. Mirrors
      # dsm.options$refit.without.discrete in the final-fit step.
      # A fresh family per fit. The object in dsm.mod.specs is shared by every
      # species and every fold, and mgcv mutates it, so without this a single
      # failed fit poisons all the ones after it. See fresh.family().
      fit_args <- list(finalists$formula[[i]],
                       data   = segdata[train_idx, ],
                       family = fresh.family(finalists$family[[i]]),
                       method = "fREML")

      fit <- try(do.call(mgcv::bam, c(fit_args, list(discrete = TRUE,
                                                     nthreads = nthreads))),
                 silent = TRUE)
      fit_path <- "discrete"

      # The failure is intermittent and load-related, so retry single-threaded
      # first: it removes the non-determinism in the threaded accumulation and
      # costs roughly the same. Only then consider discrete = FALSE, which is
      # correct but wildly expensive - roughly 280x the discrete path (474 min
      # against 1.7 min per fit on NL_EXPL_DRL_RA, 153,344 segments, 80 cores),
      # so it is off unless asked for.
      if (inherits(fit, "try-error") && nthreads != 1) {
        message(sprintf("run.family.cv: %s %s fold %d failed (%s); retrying single-threaded",
                        species, finalists$modname[i], k,
                        sub("
.*", "", fit[1])))
        fit_args$family <- fresh.family(finalists$family[[i]])
        fit <- try(do.call(mgcv::bam, c(fit_args, list(discrete = TRUE,
                                                       nthreads = 1))),
                   silent = TRUE)
        fit_path <- "discrete_1thread"
      }

      if (inherits(fit, "try-error") && isTRUE(allow.slow.refit)) {
        message(sprintf("run.family.cv: %s %s fold %d still failing; refitting with discrete = FALSE (expect hours)",
                        species, finalists$modname[i], k))
        fit_args$family <- fresh.family(finalists$family[[i]])
        fit <- try(do.call(mgcv::bam, c(fit_args, list(discrete = FALSE))),
                   silent = TRUE)
        fit_path <- "exact"
      }

      if (inherits(fit, "try-error")) {
        message(sprintf("run.family.cv: %s %s fold %d failed both ways: %s",
                        species, finalists$modname[i], k,
                        sub("
.*", "", fit[1])))
        next
      }

      # discrete = FALSE is essential here, not an optimisation choice. A model
      # fitted with discrete = TRUE re-discretises whatever newdata it is given,
      # so the SAME held-out segment gets a different prediction depending on
      # which other segments happen to share its fold - about 0.3% on this data.
      # Several species are separated by less than that, so composition-dependent
      # predictions could decide a verdict. discrete = FALSE predicts exactly and
      # is invariant to the composition of newdata; it also puts folds that fell
      # back to a non-discrete fit on the same footing as the rest.
      pred_mu <- try(stats::predict(fit, newdata = segdata[test_idx, ],
                                    type = "response", discrete = FALSE),
                     silent = TRUE)
      if (inherits(pred_mu, "try-error")) {
        message(sprintf("run.family.cv: %s %s fold %d could not predict: %s",
                        species, finalists$modname[i], k,
                        sub("\n.*", "", pred_mu[1])))
        next
      }
      obs_y <- segdata[[resp]][test_idx]

      # A fit can converge to something mgcv returns without complaint and that
      # is still useless: theta driven to an extreme, or training blocks that
      # omit part of the covariate range. The tell is in the predictions, before
      # any score exists to judge - so judge them here, against the data they are
      # meant to predict.
      #
      # This is not defensive padding. Both known pathologies on Atl IMRP
      # Shearwaters are caught by it: the nb fit that scored 3e17 (predictions of
      # 8e18 birds per segment), and the tw fit that never produced a score at
      # all, because mgcv::rTweedie draws a Poisson count per observation and
      # then that many gamma variates - so absurd fitted values ask for an absurd
      # vector. It asked for 1198.6 Gb and took a three-hour run down with it.
      scale_ref <- max(c(obs_y[is.finite(obs_y)], 1))
      if (!all(is.finite(pred_mu)) || max(pred_mu) > max.pred.ratio * scale_ref) {
        message(sprintf("run.family.cv: %s %s fold %d diverged - largest prediction %.3g against a largest observed %.3g; dropping the fold",
                        species, finalists$modname[i], k,
                        suppressWarnings(max(pred_mu)), scale_ref))
        next
      }

      # Everything from here to the score is wrapped, for the same reason the
      # fit is: a bad model should cost one fold, not the whole run.
      scored <- try({
        # Same reasoning as assign.blocks(): fix the generator as well as the
        # seed, so scores are reproducible whether or not the caller is inside a
        # future.
        old_kind <- RNGkind()
        on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
        suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
        set.seed(seed + k)
        sim_mat <- sim.response(fit, pred_mu, n_sim)
        RNGkind(old_kind[1], old_kind[2], old_kind[3])

        # 90% predictive interval coverage: the direct check on whether a
        # family's tails are too light to carry prediction variance downstream.
        list(sim_mat = sim_mat,
             lower   = apply(sim_mat, 1, stats::quantile, 0.05),
             upper   = apply(sim_mat, 1, stats::quantile, 0.95))
      }, silent = TRUE)

      if (inherits(scored, "try-error")) {
        message(sprintf("run.family.cv: %s %s fold %d could not be scored: %s",
                        species, finalists$modname[i], k,
                        sub("\n.*", "", scored[1])))
        next
      }
      sim_mat <- scored$sim_mat
      lower   <- scored$lower
      upper   <- scored$upper

      row <- tibble::tibble(
        species    = species,
        family_key = finalists$family_key[i],
        modname    = finalists$modname[i],
        fold       = k,
        n_test     = length(test_idx),
        # Constant within a species, but carried per row so the block count -
        # the real sample size behind the folds - survives into the results file
        # and the report, which otherwise only see the fitting messages.
        n_blocks   = attr(fold_id, "n_blocks"),
        CRPS       = mean(calc.crps(sim_mat, obs_y)),
        MAE        = mean(abs(obs_y - pred_mu)),
        RMSE       = sqrt(mean((obs_y - pred_mu)^2)),
        PIT_KS     = as.numeric(suppressWarnings(
                       stats::ks.test(calc.pit(sim_mat, obs_y), "punif")$statistic)),
        cover90    = mean(obs_y >= lower & obs_y <= upper),
        fit_path   = fit_path)

      # Observed/expected by the calibration covariate, if there is one. Kept
      # generic: one column per level, plus the mean distance from 1 so the
      # verdict logic downstream never has to know the levels.
      if (!is.null(calib.covar) && calib.covar %in% names(segdata)) {
        grp <- segdata[[calib.covar]][test_idx]
        oe  <- tapply(obs_y, grp, sum) / tapply(pred_mu, grp, sum)
        row <- dplyr::bind_cols(
          row,
          tibble::as_tibble(stats::setNames(as.list(oe), paste0("OE_", names(oe)))),
          tibble::tibble(calib = mean(abs(1 - oe), na.rm = TRUE)))
      }

      row$mins <- as.numeric(difftime(Sys.time(), started, units = "mins"))
      out <- dplyr::bind_rows(out, row)

      if (isTRUE(progress))
        message(sprintf("  %s %-34s fold %d done (%.1f min) CRPS=%.4f",
                        species, finalists$modname[i], k,
                        row$mins, row$CRPS))
    }
  }
  out
}


#' Apply dsm_var_gam to one chunk of prediction data
#'
#' Calls \code{\link[dsm]{dsm_var_gam}} on \code{dat} using the offset stored
#' in \code{dat$.my.off.set}.  Used internally by \code{\link{get.per.cell.var}}.
#'
#' @param dat Data frame (single prediction grid chunk) with a
#'   \code{.my.off.set} column.
#' @param this.dsm Fitted \code{dsm} object.
#' @return Named list with elements \code{pred.var} (per-cell variance) and
#'   \code{pred} (numeric prediction vector), carrying a \code{chunk.stats}
#'   attribute: \code{cells}, \code{mins} and \code{peak.gb} for this chunk.
#'   \code{\link{get.per.cell.var}} summarises those across chunks. They exist
#'   because one ATPU chunk of 220,000 rows ran 5x slower than nine others of
#'   the same size and reached a 31.7 GB working set, which nothing in the code
#'   explains (issue #46); without per-chunk numbers the next long run would
#'   again only be diagnosable from outside, with \code{tasklist}.
#' @export
apply.dsm.var <- function(dat, this.dsm){
  t0 <- proc.time()[["elapsed"]]

  res <- dsm::dsm_var_gam(this.dsm, dat, purrr::map(dat, ".my.off.set"))

  out <- list(pred.var = res$pred.var, pred = unlist(res$pred))

  # gc() here both reports memory and returns the chunk's, which the serial
  # path used to do with a bare gc() call - see get.per.cell.var().
  #
  # Column 6 of the gc() matrix is "max used (Mb)", summed over Ncells/Vcells.
  # It is the high-water mark of this R process, not of this chunk alone, so
  # across chunks on one worker it only ever rises. That is the number wanted:
  # what made issue #46 visible from outside was one worker process holding
  # 31.7 GB while nine held 172 MB.
  attr(out, "chunk.stats") <- c(
    cells   = length(dat),
    mins    = (proc.time()[["elapsed"]] - t0) / 60,
    peak.gb = sum(gc()[, 6]) / 1024
  )

  out
}

#' Compute per-cell variance for a DSM in memory-safe chunks
#'
#' Splits the prediction grid into \code{nchunks} pieces, applies
#' \code{\link{apply.dsm.var}} to each chunk (optionally in parallel), and
#' returns the results.  Chunking avoids the out-of-memory errors that arise
#' when \code{dsm_var_gam} is called on a very large prediction grid.
#'
#' @param this.dsm Fitted \code{dsm} object.
#' @param df Prediction grid data frame.
#' @param nchunks Integer number of chunks to split \code{df} into. Set it
#'   well above \code{nodes} when running in parallel - the chunks are
#'   dispatched with load balancing, so \code{nchunks == nodes} leaves nothing
#'   to balance and the step runs as long as its slowest chunk. Smaller chunks
#'   also cap per-worker memory. \code{nchunks} does not change the answer when
#'   \code{exact.predict} is \code{TRUE}.
#' @param exact.predict If \code{TRUE} (default) and the model was fitted with
#'   \code{discrete = TRUE}, drop its \code{$dinfo} so that \code{dsm_var_gam()}
#'   predicts exactly. \code{dsm_var_gam()} takes no \code{discrete} argument, so
#'   this is the only way to reach it: without it each chunk is discretised
#'   separately and a cell's variance depends on which chunk it fell in, meaning
#'   \code{nchunks} silently changes the answer. Only the evaluation changes -
#'   coefficients and \code{Vp} are untouched.
#' @param off.set Numeric offset (cell area): either a scalar or a vector the
#'   same length as \code{nrow(df)}.
#' @param parallel If \code{TRUE}, process chunks in parallel with
#'   \code{parLapplyLB}.
#' @param nodes Integer number of cluster nodes for parallel processing.
#' @return List of per-chunk results from \code{\link{apply.dsm.var}}, in the
#'   order of the chunks, each carrying a \code{chunk.stats} attribute.
#'   \code{\link{report.chunk.stats}} is called on the way out and prints them.
#' @export
get.per.cell.var <- function(this.dsm,
                             df,
                             nchunks,
                             off.set = 1,
                             parallel = FALSE,
                             nodes = 1,
                             exact.predict = TRUE
) {

  # If off.set is a vector add it to the df so it will be split properly as well.
  if (length(off.set) > 1 && length(off.set) != nrow(df)) {
    stop("get.per.cell.var: off.set vector is not same length as df")
  }

  # dsm_var_gam() predicts internally and takes no discrete argument, so there
  # is no way to ask it for exact prediction. But a model fitted with
  # discrete = TRUE re-discretises whatever newdata it is handed, so with
  # nchunks > 1 each chunk is discretised on its own and a cell's variance
  # depends on which chunk it landed in - change nchunks and the answer changes.
  #
  # Dropping $dinfo makes predict.bam() take the exact path, which is precisely
  # what predict(discrete = FALSE) does. Verified identical to that call, and
  # composition-independent, while leaving the coefficients and Vp untouched -
  # this changes how the fitted model is evaluated, not the fit itself.
  if (isTRUE(exact.predict) && !is.null(this.dsm$dinfo)) {
    message(sprintf(paste0(
      "get.per.cell.var: model was fitted with discrete = TRUE; predicting ",
      "exactly instead%s. Pass exact.predict = FALSE to keep the discretised ",
      "path."),
      if (nchunks > 1)
        sprintf(", so per-cell variance does not depend on which of the %d chunks a cell falls in",
                nchunks) else ""))
    this.dsm$dinfo <- NULL
  } else if (!is.null(this.dsm$dinfo) && nchunks > 1) {
    warning(sprintf(paste0(
      "get.per.cell.var: model was fitted with discrete = TRUE and ",
      "exact.predict is FALSE, so each of the %d chunks is discretised ",
      "separately and per-cell variance depends on chunk boundaries."), nchunks),
      call. = FALSE)
  }

  df$.my.off.set <- off.set

  # split data into chunks for processing
  if (nchunks > 1) {
    dat.split <- split(df, cut(1:nrow(df), nchunks, FALSE))
  } else{
    # Process all of df in one chunk, make it a list so it can be processed
    # by either map() of parallel::parLapply() below.
    dat.split <- list(df)
  }

  # split each chunk in dat.split into sublists with 1 cell per element.
  dat.split <- purrr::map(dat.split, ~ split(.x, 1:nrow(.x)))

  if (parallel) {

    cl <- parallel::makeCluster(nodes)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Could just execute the things we need instead.
    # parallel::clusterEvalQ(cl, source(here::here("R/analysis_settings.R")))
    parallel::clusterEvalQ(cl, {
      library(dsm)
      library(purrr)
    })

    # Export the model and the worker function - not dat.split, which used to
    # be here too. That sent every chunk to every worker on top of the chunk
    # the dispatcher already delivers: on NL_EXPL_DRL_RA, 2.2 million single-row
    # data frames per worker, each needing 220,000 of them (issue #46).
    #
    # The model goes by export rather than as a parLapplyLB argument because an
    # argument is serialised once per task: with nchunks > nodes the model would
    # be shipped repeatedly, which is the same mistake in a different place.
    #
    # apply.dsm.var is exported by name rather than reached as
    # DSMHelper::apply.dsm.var so the workers run whatever the driver is
    # running. That matters when the package is being developed by sourcing
    # functions.R, where the installed version is a different, older one.
    parallel::clusterExport(cl, c("apply.dsm.var", "this.dsm"),
                            envir = environment())

    # environment(worker) <- globalenv() is load-bearing, not tidiness. A
    # closure defined here carries this frame, and this frame holds df and
    # dat.split - so leaving it attached would ship the whole prediction grid
    # again, by the back door, and undo the fix above. Detached, both names
    # resolve to the copies clusterExport put in each worker's global env.
    worker <- function(chunk) apply.dsm.var(chunk, this.dsm)
    environment(worker) <- globalenv()

    # parLapplyLB with chunk.size = 1 rather than parLapply: static scheduling
    # made the step as slow as its slowest chunk, and on ATPU nine workers sat
    # idle for over five hours while the tenth finished. Load balancing needs
    # nchunks > nodes to have anything to balance - see the argument docs.
    # Results still come back in the order of dat.split.
    print(system.time(
      res <- parallel::parLapplyLB(cl, dat.split, worker, chunk.size = 1)
    ))

  } else {  # non-parallel version
    # apply the function to the chunks serially with map.
    #
    # This used to end each iteration with a bare gc(), which is what map()
    # then returned - so the serial path handed back gc() matrices instead of
    # results and map(res, "pred.var") came out empty. apply.dsm.var() now
    # calls gc() itself, for the chunk stats, so the collection still happens.
    print(system.time(
      res <- purrr::map(dat.split, apply.dsm.var, this.dsm = this.dsm)
    ))
  }

  report.chunk.stats(res)

  res
}

#' Summarise the per-chunk cost of a get.per.cell.var() run
#'
#' Prints the \code{chunk.stats} attribute that \code{\link{apply.dsm.var}}
#' attaches to each chunk result: cells, minutes and peak memory, plus the
#' spread across chunks. A single slow chunk is the failure mode this exists to
#' make visible (issue #46), and it is invisible in a wall-clock total.
#'
#' @param res List of chunk results from \code{\link{get.per.cell.var}}.
#' @return \code{invisible(NULL)}; called for the printed summary.
#' @examples
#' \dontrun{
#' res <- get.per.cell.var(m, predgrid, nchunks = 40, parallel = TRUE, nodes = 10)
#' report.chunk.stats(res)
#' }
#' @export
report.chunk.stats <- function(res) {
  checkmate::expect_list(res, min.len = 1)

  stats <- purrr::map(res, ~ attr(.x, "chunk.stats"))
  if (any(purrr::map_lgl(stats, is.null))) return(invisible(NULL))

  st <- as.data.frame(do.call(rbind, stats))
  st$chunk <- seq_len(nrow(st))

  message(sprintf(
    "get.per.cell.var: %d chunks, %s cells, %.1f min of chunk time",
    nrow(st), format(sum(st$cells), big.mark = ","), sum(st$mins)))
  message(sprintf(
    "  slowest chunk %d at %.1f min vs median %.1f (%.1fx), peak %.1f GB",
    st$chunk[which.max(st$mins)], max(st$mins), stats::median(st$mins),
    max(st$mins) / stats::median(st$mins), max(st$peak.gb)))

  # The ATPU run that prompted this was 5x. Anything near that is the same
  # unexplained thing, not a rough edge in the chunking.
  if (max(st$mins) > 2 * stats::median(st$mins))
    warning(sprintf(paste0(
      "get.per.cell.var: chunk %d took %.1fx the median. Chunks are equal ",
      "sized, so this is not the split. See issue #46."),
      st$chunk[which.max(st$mins)], max(st$mins) / stats::median(st$mins)),
      call. = FALSE)

  invisible(NULL)
}


#' Compute a density estimate with uncertainty from a DSM
#'
#' \strong{Note: not currently used.}  \code{Generic_4_variance.Rmd} does this
#' job, and does it per cell as well as in total.  Kept because the lognormal
#' CI below is not available anywhere else.
#'
#' Both the offset and the divisor come from \code{predgrid$area}, not from the
#' \code{predgridCellArea} global.  Those are not the same thing: cells clipped
#' by the study area boundary are smaller than a whole cell (measured on
#' \code{NL_EXPL_DRL_RA} at 2 km, 1,039 of 183,197 cells run down to 1.908 of
#' 3.996 sq km), and the nominal constant is itself a rounding of the real cell
#' size, because \code{rast(resolution = )} adjusts the cell to fit the study
#' area extent in whole cells.  Using the constant inflated the offset on every
#' coastal cell and then divided the total by an area the grid does not cover.
#'
#' @param dsm_final Fitted \code{dsm} object.
#' @param predgrid Prediction grid with an \code{area} column in square km, as
#'   produced by \code{00.03_Create_prediction_grids.Rmd}.
#' @return Named list with elements \code{pred.est}, \code{CV}, \code{SE},
#'   and \code{CI} (a three-element vector giving the 5%, mean, and 95%
#'   lognormal confidence interval).  \code{pred.est} is a density, in
#'   individuals per square km.
#' @export
get.dens.est <- function(dsm_final, predgrid) {
  checkmate::expect_multi_class(predgrid, c("sf", "data.frame"))
  if (!"area" %in% names(predgrid))
    stop("get.dens.est: predgrid has no 'area' column. It is added by ",
         "00.03_Create_prediction_grids.Rmd and is the per-cell area in sq km.")
  checkmate::expect_numeric(predgrid$area, lower = 0, any.missing = FALSE,
                            min.len = 1)

  # use dsm.var.gam to get estimated abundance params
  densEst <- summary(dsm::dsm.var.gam(dsm_final, predgrid,
                                      off.set = predgrid$area))

  # The estimates in dsm.var.gam are summed over the entire predgrid, so divide
  # by the total area the grid actually covers to get a density.
  total.area <- sum(predgrid$area)
  densEst %<>%
    purrr::map_at(c("pred.est", "se"), ~ .x / total.area)

  #calculate a lognormal CI for the density est.
  cv.square <- densEst$cv^2
  asymp.ci.c.term <- exp(1.96*sqrt(log(1 + cv.square))) # stolen from print.summary.dsm.var.R. Also in buckland et al 2001.
  list(pred.est = densEst$pred.est, CV = densEst$cv, SE = densEst$se, CI = c("%5" = densEst$pred.est/asymp.ci.c.term,
                                                                             "Mean" = densEst$pred.est, "95%" = densEst$pred.est * asymp.ci.c.term))
}

#' Pretty-print the output of get.dens.est
#'
#' \strong{Note: not currently used.}
#'
#' Named with an underscore rather than the project's usual dot separator:
#' \code{print} is an S3 generic, so a function called \code{print.dens.est}
#' is registered by roxygen as a \code{print} method for the (non-existent)
#' class \code{"dens.est"} instead of being exported. Its argument is a plain
#' named list, so it is not a method.
#'
#' @param densEst Named list as returned by \code{\link{get.dens.est}}.
#' @return \code{invisible(NULL)}, called for its side-effect (printed output).
#' @export
print_dens_est <- function(densEst){
  cat("Density estimate:\n\n")

  cat("Approximate asymptotic confidence interval:\n")
  print(densEst$CI)

  cat("\n")
  cat("Point estimate                 :", densEst$pred.est,"\n")
  cat("Standard error                 :", densEst$SE,"\n")
  cat("Coefficient of Variation       :", densEst$CV,"\n")
}


#' Render a Generic analysis RMarkdown file for one species
#'
#' Knits one of the \code{Generic_x_xxx.Rmd} pipeline files with
#' \code{species} as a parameter and writes the HTML output to
#' \code{ResultsDir/<species>/}.  Intermediate files are written to
#' \code{tempdir()} to allow concurrent renders.
#'
#' Expects project globals \code{ResultsDir} and \code{RDir}.
#'
#' @param species Character string species code.
#' @param file Character string filename of the RMarkdown file (not full
#'   path), relative to \code{RDir} (e.g. \code{"Generic_2_dsm.Rmd"}).
#' @return Character string \code{species}, invisibly.
#' @export
do.generic.render <- function(species, file){

  suffix <- stringr::str_replace(file, "^Generic", "") %>%
    stringr::str_replace("Rmd$", "html")
  out.file <- file.path(ResultsDir, species, paste0(species, suffix))

  # Make sure output dir exists
  if (!dir.exists(dirname(out.file)))
    dir.create(dirname(out.file), recursive = TRUE)
  render_file <- file.path(RDir, file)

  # Output message in render pane if called from a knitted document.
  if (isTRUE(getOption('knitr.in.progress'))) {
    cat(sprintf("\nRendering generic file %s for %s to %s", file, species, out.file),
        file = stderr())
  }

  message(sprintf("Rendering generic file %s for %s to %s", file, species, out.file))
  rmarkdown::render(
    render_file,
    params = list(species = species),
    intermediates_dir = tempdir(),
    output_file = out.file
  )

  return(species)
}


#' Strip the ddftype suffix from a segment label
#'
#' `Sample.Label` carries the ddftype as a trailing `_SWD` / `_SFD` / `_SWN` /
#' `_SFN` / `_AWD` ... suffix, so it is unique per **row**, not per segment:
#' `segdata` holds one row per (segment × ddftype) and every label differs.
#' Deduplicating on `Sample.Label` is therefore a silent no-op — on Atl IMRP it
#' leaves all 468,424 rows in place rather than the 117,106 segments they
#' represent.
#'
#' @param sample.label Character vector of `Sample.Label` values.
#' @return Character vector of segment identifiers, one per physical segment.
#' @examples
#' segment.id.from.label(c("ECSAS_shp_967841155_SWD", "ECSAS_shp_967841155_SFD"))
#' @export
segment.id.from.label <- function(sample.label) {
  checkmate::expect_character(sample.label, any.missing = FALSE)
  sub("_[A-Z]{3}$", "", sample.label)
}


#' Assess extrapolation for one species and season
#'
#' Runs [dsmextra::extrapolation_analysis()] to identify prediction-grid cells
#' whose covariate values are novel relative to the segments the model was
#' fitted on, and maps the result back onto the geometry of the prediction
#' raster.
#'
#' ExDet (Mesgaran et al. 2014) is negative under **univariate** extrapolation
#' (a covariate outside its sampled range), lies in 0-1 where conditions are
#' **analogue** to the sample, and exceeds 1 under **combinatorial**
#' extrapolation (each covariate in range, but the combination unobserved). The
#' most influential covariate (MIC) names the covariate responsible.
#'
#' @section Why the result is resampled rather than joined:
#' The prediction grid is not a complete lattice — it is the ocean-only subset
#' of one, and 2.77 per cent of Atl IMRP's cells do not sit on the 10 km lattice
#' at all. `dsmextra` therefore rasterises it onto a lattice of its own, whose
#' origin is offset from the prediction raster's by half a cell. Joining on
#' coordinates silently matches nothing. This function resamples
#' (`method = "near"`) onto `pred.raster` instead, so every returned row is a
#' cell of the surface being assessed. Cells `dsmextra` could not place come
#' back `NA` and are counted in the return value rather than guessed at.
#'
#' @section Running before the predictions exist:
#' The assessment itself needs only the segments and the prediction grid, so it
#' is run **before** `03.00_Do_all_prediction.Rmd` - its findings are an input to
#' masking, not a commentary on predictions already made. `pred.raster` is
#' therefore optional. When it is `NULL` the target geometry is rasterized from
#' `predgrid` by [make.raster()], which is the same template
#' [make.season.raster()] builds when the predictions are written, so the `cell`
#' indices in the returned `cells` tibble address the prediction rasters that
#' `03.00` will later write. Verified on Atl IMRP ATPU Fall: identical
#' resolution, extent, CRS and 222 x 154 dimensions, 12,970 non-`NA` cells on
#' each, and no cell non-`NA` in one but not the other.
#'
#' Passing a `pred.raster` explicitly is still supported and gives the same
#' answer; it is worth doing only when assessing a surface that is not on the
#' prediction grid's geometry.
#'
#' @param species Character string species code.
#' @param season Character string season name; must be one of `season.names`.
#' @param segdata Segment data frame (or `sf`) for `species`, containing
#'   `Season`, `Sample.Label` and every name in `covariate.names`.
#' @param predgrid `sf` prediction grid with monthly covariate columns, as
#'   loaded from `prediction_grids.rda`.
#' @param pred.raster `SpatRaster` defining the target geometry, or `NULL`
#'   (the default) to build one from `predgrid`. See *Running before the
#'   predictions exist*.
#' @param covariate.names Character vector of covariates to assess. Defaults to
#'   the project global `extrap.covars`.
#' @param crs Projected coordinate system. Defaults to the project global
#'   `segProj`.
#' @param resolution Target raster resolution in map units, passed to
#'   `dsmextra`. Required because the grid is irregular; defaults to
#'   `predgridCellLength * 1000`.
#' @param compute.nearby Logical. Also compute the Gower's-distance %N surface?
#'   Costs roughly 30x the ExDet computation: on Atl IMRP HERG Spring (26,956
#'   segments, 12,970 grid cells, 5 covariates, single-threaded) ExDet took
#'   1.9 s and nearby 63 s, so the full 14 species x 4 seasons is about an hour.
#' @param verbose Logical, passed through to `dsmextra`.
#'
#' @section Maps:
#' Maps are **not** generated here. `extrapolation_analysis(map.generate =
#' TRUE)` calls `print()` on each leaflet widget, which is the wrong behaviour
#' inside a function and drops the widgets on the floor when it is called from
#' a loop. The report calls [dsmextra::map_extrapolation()] itself on
#' `$extrap$extrapolation` and `$extrap$nearby`, where it can pass the sightings
#' and tracks overlays and control display.
#' @return A list with `extrap` (a list holding the `dsmextra` objects:
#'   `$extrapolation`, `$compare`, and `$nearby` when computed), `cells` (a tibble
#'   with one row per non-`NA` cell of `pred.raster`: `cell`, `x`, `y`,
#'   `ExDet`, `type`, `mic_name`), and `info` (a one-row tibble of counts:
#'   samples used, cells assessed, cells unplaced, and the NA drops).
#' @references Mesgaran MB, Cousens RD, Webber BL (2014). Here be dragons.
#'   Diversity & Distributions 20:1147-1159.
#' @examples
#' \dontrun{
#' res <- assess.extrapolation("HERG", "Spring", segdata, predgrid, pred.raster)
#' subset(res$cells, type == "univariate" & mic_name == "depth.g")
#' }
#' @export
assess.extrapolation <- function(species, season, segdata, predgrid,
                                 pred.raster = NULL,
                                 covariate.names = extrap.covars,
                                 crs = segProj,
                                 resolution = predgridCellLength * 1000,
                                 compute.nearby = TRUE,
                                 verbose = FALSE) {
  checkmate::expect_string(species)
  checkmate::expect_choice(season, season.names)
  checkmate::expect_data_frame(segdata)
  checkmate::expect_class(predgrid, "sf")
  if (!is.null(pred.raster))
    checkmate::expect_class(pred.raster, "SpatRaster")
  checkmate::expect_character(covariate.names, min.len = 1, any.missing = FALSE)
  checkmate::expect_number(resolution, lower = 0)
  checkmate::expect_flag(compute.nearby)

  seg <- sf::st_drop_geometry(segdata)
  missing.cols <- setdiff(c(covariate.names, "Season", "Sample.Label"), names(seg))
  if (length(missing.cols))
    stop("assess.extrapolation: segdata is missing column(s): ",
         paste(missing.cols, collapse = ", "))

  # ---- samples: this season only, one row per segment -----------------------
  seg <- seg[as.character(seg$Season) == season, , drop = FALSE]
  if (!nrow(seg))
    stop("assess.extrapolation: no ", species, " segments in season ", season)
  seg$.seg_id <- segment.id.from.label(seg$Sample.Label)
  samples <- seg[!duplicated(seg$.seg_id), covariate.names, drop = FALSE]
  n.samples.raw <- nrow(samples)
  samples <- samples[stats::complete.cases(samples), , drop = FALSE]
  if (!nrow(samples))
    stop("assess.extrapolation: every ", species, " ", season,
         " segment has an NA in ", paste(covariate.names, collapse = "/"))

  # ---- prediction grid: the same seasonal covariate values the predictions
  # were made on. Neither side effect of create.seasonal.predgrid() is wanted
  # here - see its documentation.
  pgrid <- create.seasonal.predgrid(species, predgrid,
                                    write.shapefile = FALSE,
                                    replicate.platform = FALSE) %>%
    sf::st_drop_geometry()
  pgrid <- pgrid[as.character(pgrid$Season) == season, , drop = FALSE]
  missing.cols <- setdiff(covariate.names, names(pgrid))
  if (length(missing.cols))
    stop("assess.extrapolation: seasonal prediction grid is missing column(s): ",
         paste(missing.cols, collapse = ", "),
         ". Check that they are covariates create.seasonal.predgrid() retains.")
  pgrid <- pgrid[, covariate.names, drop = FALSE]
  n.grid.raw <- nrow(pgrid)
  # dsmextra drops NA rows silently; do it here so the count can be reported.
  pgrid <- pgrid[stats::complete.cases(pgrid), , drop = FALSE]

  message(sprintf(
    "assess.extrapolation: %s %s - %d segments (%d before NA drop), %d grid cells (%d before)",
    species, season, nrow(samples), n.samples.raw, nrow(pgrid), n.grid.raw))

  # ---- dsmextra -------------------------------------------------------------
  # compute_extrapolation/compute_nearby are called directly rather than through
  # extrapolation_analysis(), which has no `resolution` argument and so cannot
  # pass one down. This grid always needs one: it is the ocean-only subset of a
  # lattice and 2.77% of Atl IMRP's cells are off that lattice entirely, so
  # dsmextra always rasterises and errors out without a resolution.
  extrap <- list()
  extrap$extrapolation <- dsmextra::compute_extrapolation(
    samples           = samples,
    covariate.names   = covariate.names,
    prediction.grid   = pgrid,
    coordinate.system = crs,
    resolution        = resolution,
    verbose           = verbose)

  extrap$compare <- dsmextra::compare_covariates(
    extrapolation.type   = "both",
    extrapolation.object = extrap$extrapolation,
    n.covariates         = NULL,
    create.plots         = FALSE,
    display.percent      = TRUE,
    verbose              = verbose)

  if (compute.nearby)
    extrap$nearby <- dsmextra::compute_nearby(
      samples           = samples,
      covariate.names   = covariate.names,
      prediction.grid   = pgrid,
      coordinate.system = crs,
      nearby            = 1,
      resolution        = resolution,
      verbose           = verbose)

  # No prediction raster to resample onto - this normally runs before 03.00 has
  # written any. Rasterizing the grid itself gives the same geometry those
  # rasters will have; see "Running before the predictions exist" above.
  if (is.null(pred.raster)) {
    template.field <- intersect(c("area", "x", "depth"), names(predgrid))[1]
    if (is.na(template.field))
      stop("assess.extrapolation: predgrid has no column to build a raster ",
           "template from (looked for area, x, depth). Pass pred.raster.")
    pred.raster <- make.raster(predgrid, template.field)
  }

  cells <- extrapolation.cells.on.raster(extrap$extrapolation, pred.raster,
                                         covariate.names)

  info <- tibble::tibble(
    species          = species,
    season           = season,
    n_samples        = nrow(samples),
    n_samples_nadrop = n.samples.raw - nrow(samples),
    n_grid           = nrow(pgrid),
    n_grid_nadrop    = n.grid.raw - nrow(pgrid),
    n_cells          = nrow(cells),
    n_cells_unplaced = sum(is.na(cells$ExDet)))

  list(extrap = extrap, cells = cells, info = info)
}


#' Map a dsmextra result onto the geometry of a prediction raster
#'
#' Resamples the ExDet and MIC rasters returned by
#' [dsmextra::compute_extrapolation()] onto `pred.raster` and returns one row
#' per non-`NA` cell of that raster. See the *Why the result is resampled*
#' section of [assess.extrapolation()] for why a coordinate join does not work.
#'
#' @param ex A `dsmextra` extrapolation object with a `$rasters` element.
#' @param pred.raster `SpatRaster` defining the target geometry.
#' @param covariate.names Character vector used to resolve the integer MIC
#'   index to a covariate name.
#' @return A tibble with `cell`, `x`, `y`, `ExDet`, `type`, `mic_name`. `type`
#'   is `"univariate"`, `"analogue"`, `"combinatorial"`, or `NA` where
#'   `dsmextra` placed no value.
#' @export
extrapolation.cells.on.raster <- function(ex, pred.raster, covariate.names) {
  checkmate::expect_list(ex)
  checkmate::expect_class(pred.raster, "SpatRaster")
  checkmate::expect_character(covariate.names, min.len = 1, any.missing = FALSE)
  if (is.null(ex$rasters$ExDet$all))
    stop("extrapolation.cells.on.raster: no $rasters$ExDet$all in the dsmextra object.")

  ExDet <- terra::resample(terra::rast(ex$rasters$ExDet$all), pred.raster,
                           method = "near")
  MIC   <- terra::resample(terra::rast(ex$rasters$mic$all), pred.raster,
                           method = "near")

  keep <- !is.na(terra::values(pred.raster)[, 1])
  ev   <- terra::values(ExDet)[, 1]
  mv   <- terra::values(MIC)[, 1]

  # dsmextra codes "no MIC" as 0. Indexing a vector with 0 silently drops the
  # element rather than returning NA, so map it out before subsetting.
  mv[!is.na(mv) & mv == 0] <- NA_integer_

  xy <- terra::xyFromCell(pred.raster, seq_along(ev))

  tibble::tibble(
    cell     = seq_along(ev),
    x        = xy[, 1],
    y        = xy[, 2],
    ExDet    = ev,
    type     = dplyr::case_when(is.na(ev) ~ NA_character_,
                                ev < 0    ~ "univariate",
                                ev > 1    ~ "combinatorial",
                                TRUE      ~ "analogue"),
    mic_name = ifelse(is.na(mv), NA_character_, covariate.names[mv]))[keep, ]
}


#' Summarise how much predicted abundance sits in extrapolated cells
#'
#' Joins an [assess.extrapolation()] result to a prediction surface and reports,
#' per extrapolation type and per most-influential covariate, how many cells are
#' involved and what share of the seasonal total they carry.
#'
#' ExDet is blind to the response: a cell can be wildly novel and hold no
#' predicted birds, or be perfectly analogue and hold most of them. This is the
#' function that tells the two apart, and it is the bridge between the
#' extrapolation assessment and the concentration check in
#' `03.70_Save_chosen_model_predictions.Rmd`.
#'
#' @param cells The `cells` tibble from [assess.extrapolation()].
#' @param pred.raster `SpatRaster` of predicted abundance, on the same geometry
#'   `cells` was built against.
#' @return A list with `by_type` and `by_mic` tibbles (`cells`, `abundance`,
#'   `pct_of_total`, most-concentrated first) and `total`, the summed surface.
#' @export
summarise.extrapolation.abundance <- function(cells, pred.raster) {
  checkmate::expect_data_frame(cells)
  checkmate::expect_class(pred.raster, "SpatRaster")

  cells$pred <- terra::values(pred.raster)[, 1][cells$cell]
  total <- sum(cells$pred, na.rm = TRUE)
  pct   <- function(x) if (total > 0) 100 * sum(x, na.rm = TRUE) / total else NA_real_

  by_type <- cells %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(cells = dplyr::n(),
                     abundance = sum(.data$pred, na.rm = TRUE),
                     pct_of_total = pct(.data$pred),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$pct_of_total))

  by_mic <- cells %>%
    dplyr::filter(!is.na(.data$mic_name)) %>%
    dplyr::group_by(type, mic_name) %>%
    dplyr::summarise(cells = dplyr::n(),
                     abundance = sum(.data$pred, na.rm = TRUE),
                     pct_of_total = pct(.data$pred),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$pct_of_total))

  list(by_type = by_type, by_mic = by_mic, total = total)
}


# do.extrapolation() used to live here. It was the Hibernia version, taking
# spill/dataset and rendering a Generic_0_extrapolation.Rmd that does not exist
# in this project, so it errored on any call. The replacement would have been a
# near-copy of do.generic.render(), which already renders Generic_*.Rmd per
# species to Results/[SubProject]/[species]/. 02.60_Assess_Extrapolation.Rmd
# uses that instead.


#' Comprehensive diagnostics for a fitted DSM
#'
#' Runs a full suite of checks including smooth plots, gratia appraise,
#' mgcv \code{\link[mgcv]{gam.check}}, DHARMa residual tests, spatial
#' autocorrelation tests, observed vs expected plots, and (optionally)
#' concurvity checks and variograms.
#'
#' @param dsm_final Fitted \code{dsm} object.
#' @param modname Character string model name; used in messages and titles.
#' @param segdata Segment data frame (used for spatial autocorrelation tests
#'   and residual-vs-term plots).
#' @param smoother.plots If \code{TRUE}, produce \code{gratia::draw()} smoother
#'   plots.
#' @param brief If \code{TRUE} (default), skip the more expensive diagnostics
#'   (concurvity, autocorrelogram, variograms).
#' @return Invisibly, a named list of the statistics computed along the way -
#'   \code{k.check}, \code{residuals}, \code{zeroinflation}, \code{spatial},
#'   \code{overdispersion}, \code{oe.platform}, \code{oe.depth}, and with
#'   \code{brief = FALSE} also \code{concurvity} and \code{variogram}. Entries
#'   that could not be computed hold the \code{try-error}. Pass it to
#'   \code{\link{interpret.dsm.checks}} for plain-English verdicts. Still called
#'   mainly for its side-effects: the printing and plotting are unchanged, and
#'   the statistics are captured as they are produced rather than recomputed,
#'   because \code{DHARMa::simulateResiduals()} is the slowest thing here.
#' @section Reproducibility under parallelism:
#'   Safe to call from a \code{future}/\code{furrr} worker: the generator is
#'   pinned to Mersenne-Twister for the duration and restored afterwards, so the
#'   seeded steps give the same numbers as a serial run. Without that pin they
#'   do not - measured on ATPU, an unpinned worker returned a k-index of 0.833
#'   against 0.896 serial and a Moran's I of -0.00216 against 0.000964, from
#'   identical seeds, because furrr's workers run L'Ecuyer-CMRG and
#'   \code{set.seed()} resets whichever generator is current.
#' @export
check.dsm <- function(dsm_final,
                      modname,
                      segdata,
                      smoother.plots = FALSE,
                      brief = TRUE # Only do partial diagnostics to speed it up
) {
  # accessing columns with segdata[, termlab] below doesn't work with tbls or
  # sf objects.
  segdata <- as.data.frame(segdata)

  # Pin the generator, not just the seed. Three things in here are stochastic
  # and seeded - k.check(), the DHARMa simulations that follow it, and the
  # location subsample for the spatial test - and set.seed() resets whichever
  # generator is CURRENT. Under furrr's future.seed = TRUE the workers run
  # L'Ecuyer-CMRG, so the same seed yields a different sequence there than in
  # the main process, and the k-index and Moran's I would silently differ
  # between a serial and a parallel run of identical code. This is the same
  # trap assign.blocks() documents, where it produced different CV folds.
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  # Several dsm-package functions do partial matching so turn off and re-enable
  # at end.
  options(warnPartialMatchDollar = FALSE)

  message("\n===========================================================================\n\nDSM model summary for model ", modname)
  print(summary(dsm_final))
  plot(dsm_final, scale = 0) # Plot smooths w/ different y-axis scale for each term

  if(isTRUE(smoother.plots)) {
    # Smoother plots
    message("Smoother plots")
    par(mfrow = c(1,1))
    # XXXX Residual plotting in gratia is not working b/c of some mess-up with the name of the
    # offset column in dsm. need to debug some more
    # gratia::draw(dsm_final, wrap = FALSE, residuals = dsm.options$do.dsm.residuals) %>%
    gratia::draw(dsm_final, wrap = FALSE) %>%
      print
  }

  if (any(grepl("s(x.sc, y.sc", as.character(dsm_final$formula), fixed = T))) {
    mgcv::vis.gam(dsm_final,  view = c("x.sc","y.sc"), main = "s(x.sc,y.sc) (response scale)",
                  type = "response", asp = 1, plot.type = "contour")
    mgcv::vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = 0, phi = 45,
                  main = "s(x.sc,y.sc) (response scale)", type = "response",
                  asp = 1, ticktype = "detailed")

    mgcv::vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = 60, phi = 45,
                  main = "s(x.sc,y.sc) (response scale)", type = "response",
                  asp = 1, ticktype = "detailed")

    mgcv::vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = -60, phi = 45,
                  main = "s(x.sc,y.sc) (response scale)", type = "response",
                  asp = 1, ticktype = "detailed")
  }

  # Shouldn't these be seasonal?
  if (any(grepl("s(x, y", as.character(dsm_final$formula), fixed = T))) {
    mgcv::vis.gam(dsm_final,  view = c("x","y"), main = "s(x, y) (response scale)",
                  type = "response", asp = 1, plot.type = "contour")
    mgcv::vis.gam(dsm_final,  view = c("x","y"), theta = 0, phi = 45,
                  main = "s(x,y) (response scale)", type = "response",
                  asp = 1, ticktype = "detailed")

    mgcv::vis.gam(dsm_final,  view = c("x","y"), theta = 60, phi = 45,
                  main = "s(x,y) (response scale)", type = "response",
                  asp = 1, ticktype = "detailed")

    mgcv::vis.gam(dsm_final,  view = c("x","y"), theta = -60, phi = 45,
                  main = "s(x,y) (response scale)", type = "response",
                  asp = 1, ticktype = "detailed")
  }


  # Sometimes whines about S3 methods.
  message("Gratia checks")
  try(print(suppressWarnings(gratia::appraise(dsm_final))))

  # Gam checks from MGCV
  par(mfrow = c(1,1))
  message("MGCV checks")
  try(my.gam.check(dsm_final))
  message("dsm::rqgam_check():")
  dsm::rqgam_check(dsm_final)

  # Everything worth interpreting is collected into `checks` as it is computed,
  # and returned. Recomputing it afterwards would mean a second
  # DHARMa::simulateResiduals(), which is the slowest thing in here.
  checks <- list(modname = modname, n = nrow(segdata),
                 family = dsm_final$family$family)

  # Basis dimension adequacy. my.gam.check() prints this but does not return it;
  # k.check() is what it calls internally, so take it from there.
  #
  # Seeded, because k.check() is stochastic - it compares the residual
  # autocorrelation against randomly resampled neighbours, so successive calls on
  # the SAME model give different answers. Three calls on ATPU's
  # dsm_tw_allpred_abund_factor returned k-index 0.884, 0.789 and 0.897 for
  # s(x,y), which straddles the threshold this check reads. Without a seed the
  # verdict would not reproduce from one run to the next.
  set.seed(get0("K_CHECK_SEED", ifnotfound = 20260828))
  checks$k.check <- try(mgcv::k.check(dsm_final), silent = TRUE)

  # Remove "dsm" class to make DHARMa happy
  message("DHARMa checks")
  simmod <- dsm_final
  class(simmod) <- class(simmod)[-1] # Simulate residuals doesn't like dsm class
  sims <- DHARMa::simulateResiduals(fittedModel = simmod)
  plot(sims)
  print(checks$residuals <- DHARMa::testResiduals(sims))
  print(checks$zeroinflation <- DHARMa::testZeroInflation(sims))

  # Spatial autocorrelation.
  #
  # This has never once run. DHARMa needs unique coordinates and the multi-ddf
  # design guarantees the opposite: every segment appears once per ddftype at
  # the same x,y, so the call failed with "requires unique x,y values" on all 22
  # models of the last NL_EXPL_DRL_RA run - and, being wrapped in try(), failed
  # quietly enough that the report simply had no spatial test in it.
  #
  # recalculateResiduals() aggregates to one residual per group, which is what
  # the commented-out block below was reaching for.
  #
  # The group has to be the LOCATION, not the segment. Aggregating by segment is
  # the obvious reading of "one residual per sample unit" and it is not enough:
  # transects revisit places, so ATPU's 38,443 segments sit on only 38,147
  # distinct coordinates, with 277 locations shared by up to 7 segments each.
  # That leaves duplicate x,y and the test refuses exactly as before - tested,
  # and it does.
  #
  # The earlier attempt grouped by location correctly and died allocating
  # memory. It would: the test builds a dense n x n distance matrix, so 38,147
  # locations need ~11.6 GB and Atl IMRP's would need ~110 GB. So cap it - a
  # random sample tests the same hypothesis with slightly less power, at ~200 MB
  # for the default 5,000.
  checks$spatial <- try({
    seg.id <- segment.id.from.label(segdata$Sample.Label)
    loc.id <- paste(segdata$x, segdata$y, sep = "_")
    recal <- DHARMa::recalculateResiduals(sims, group = loc.id)
    # recalculateResiduals returns groups in the order of factor(group) levels.
    locs <- data.frame(loc = loc.id, x = segdata$x, y = segdata$y) %>%
      dplyr::group_by(loc) %>%
      dplyr::summarise(x = dplyr::first(x), y = dplyr::first(y),
                       .groups = "drop") %>%
      dplyr::arrange(match(loc, levels(factor(loc.id))))
    resids <- recal$scaledResiduals
    stopifnot(length(resids) == nrow(locs))
    cap <- get0("SPATIAL_AUTOCORR_MAX_N", ifnotfound = 5000)
    n.locs <- nrow(locs)
    if (n.locs > cap) {
      set.seed(get0("SPATIAL_AUTOCORR_SEED", ifnotfound = 20260828))
      keep <- sort(sample(n.locs, cap))
      message(sprintf(
        "Spatial autocorrelation: %d locations, testing a random %d of them.",
        n.locs, cap))
      locs <- locs[keep, ]; resids <- resids[keep]
    }
    out <- DHARMa::testSpatialAutocorrelation(resids, locs$x, locs$y,
                                              plot = FALSE)
    out$n.tested <- nrow(locs)
    out$n.locations <- n.locs
    out$n.segments <- dplyr::n_distinct(seg.id)
    out
  }, silent = TRUE)
  if (!inherits(checks$spatial, "try-error"))
    print(checks$spatial)
  else
    message("Spatial autocorrelation test failed: ",
            conditionMessage(attr(checks$spatial, "condition")))

  # If more than one resid at a given location.
  # Note still use try() since this may fail to allocate enough memory if
  # size of data is too big.
  # if(inherits(res, "try-error")) {
  #   message("testSpatialAUtocorrelation failed: aggregating resids spatially")
  #
  # This always fails trying to allocate more than 256GB (the RAM I have so
  # let's not bother)
  # segdata$loc <- paste(as.character(segdata$x), as.character(segdata$y), sep = "_")
  # recal <- recalculateResiduals(sims, group = segdata$loc)
  # locs <- segdata %>%
  #   distinct(loc, .keep_all = TRUE)
  # try(testSpatialAutocorrelation(recal, locs$x, locs$y))
  #
  # # May fail due to not enough memory for big jobs
  # if(inherits(res, "try-error")) {
  #   message("Spatial autocorrecation test failed for aggregated data. Garbage collecting...")
  #   gc()
  # }
  # }

  message("Doing resids vs model terms")
  # plot resids vs each term in model
  try(attr(terms(simmod), "term.labels") %>%
        purrr::walk(function(termlab, sims, segdata) {
          DHARMa::plotResiduals(sims, segdata[, termlab], xlab = termlab)
        }, sims = sims, segdata = segdata))

  if (!brief) {
    # Concurvity
    message("Concurvity checks\nEach term with whole of rest of model")
    checks$concurvity <- try(mgcv::concurvity(dsm_final), silent = TRUE)
    if (!inherits(checks$concurvity, "try-error"))
      print(round(checks$concurvity, digits = 3))
    message(
      "Concurvity of pairwise terms ('estimate' measure presented)\nEach row shows how terms in columns depend on the term in that row."
    )
    checks$concurvity.pairwise <-
      try(mgcv::concurvity(dsm_final, full = FALSE)[["estimate"]], silent = TRUE)
    if (!inherits(checks$concurvity.pairwise, "try-error"))
      print(round(checks$concurvity.pairwise, digits = 3))
    message("Plot is non-symmetric, showing how terms on y-axis depend on terms on the x-axis")
    try(dsm::vis_concurvity(dsm_final))
  }

  par(mfrow = c(1,1))
  # check observed vs expected. See Miller et al 2021 pg 11 (of 18) for
  # reccommendation to use the "platform" variable to aggregate
  # by.
  message("Observed vs expected plot")
  checks$oe.platform <- try(oe.dens(dsm_final, covar = "platform", plotit = T),
                            silent = TRUE)
  if (!inherits(checks$oe.platform, "try-error")) print(checks$oe.platform)
  # depth is continuous, so bin it - without cut, oe.dens aggregates by every
  # unique depth value, giving one point (and one table column) per segment.
  checks$oe.depth <- try(oe.dens(dsm_final, covar = "depth", cut = 10, plotit = T),
                         silent = TRUE)
  if (!inherits(checks$oe.depth, "try-error")) print(checks$oe.depth)
  # oe.dens(dsm_final, covar = "depth.g", plotit = T)
  # oe.dens(dsm_final, covar = "sst", plotit = T)
  # oe.dens(dsm_final, covar = "sst.g", plotit = T)
  # # oe.dens(dsm_final, covar = "year", plotit = T)
  # # oe.dens(dsm_final, covar = "yday", plotit = T)

  # Overdispersion: Pearson chi-square / residual df. ~1 indicates the
  # mean-variance relationship is adequate; >1 overdispersed, <1 underdispersed.
  # Note the scale parameter is estimated (not fixed at 1) for the Tweedie and
  # negative binomial families used here, so treat this as a rough guide.
  message("Overdispersion statistic (Pearson chi-sq / resid df)")
  print(checks$overdispersion <- OD_dsm_final <-
          sum(resid(dsm_final, type = "pearson")^2)/dsm_final$df.res)

  # Bubble plot
  message("Doing bubbleplot")
  mydata <-
    data.frame(
      E = resid(dsm_final, type = "pearson"),
      x = segdata$x / 1000,
      y = segdata$y / 1000
    )
  sp::coordinates(mydata) <- ~ x + y

  print(sp::bubble(
    mydata,
    "E",
    col = c("black", "red"),
    main = "Residuals",
    xlab = "X-coords",
    ylab = "y-coords"
  ))

  if (!brief) {
    # check autocorellogram
    # Build the transect label from cruiseid + date + platform so the
    # duplicated segments (one for fly, one for swim) do not land in the same
    # "Transect". This previously used FlySwim. FlySwim is not gone - it is an
    # observation-level column, still present in the.data$distdata - but it has
    # never been part of the DSM segment data that check.dsm() sees, so it
    # silently contributed nothing to the label and every fly/swim pair shared
    # one. platform is the segment-level equivalent. Note the older
    # Generic_2_dsm_bam_test.R does build a segdata carrying FlySwim, which is
    # likely where the original reference came from.
    # create segment label as %H:%M:%S
    message("Doing autocorellogram")
    dsm_final$data <- dsm_final$data %>%
      dplyr::mutate(
        tr.lab = paste(
          dsm_final$data$CruiseID,
          dsm_final$data$Date,
          dsm_final$data$platform,
          sep = "_"
        ),
        seg.lab = format(lubridate::as_datetime(dsm_final$data$StartTime), "%H:%M:%S")
      )
    par(mfrow = c(1, 1))
    dsm::dsm_cor(
      dsm_final,
      Transect.Label = "tr.lab",
      Segment.Label = "seg.lab",
      max.lag = 20
    )

    message("Doing variograms")
    V <- (gstat::variogram(E ~ 1, mydata))
    plot(
      x = V$dist,
      y = V$gamma,
      xlab = "Distance (km)",
      ylab = "Semi-variance",
      pch = 16,
      cex = 2 * V$np / max(V$np)
    )

    # Fit the full-extent variogram so there is something to interpret. The five
    # plots below are the same empirical variogram at different cutoffs and V is
    # overwritten by each, so nothing used to survive this block. nugget/sill is
    # the number that matters: near 1 means the residuals are spatially
    # structureless, which is what a well-specified spatial model should leave
    # behind.
    checks$variogram <- try({
      # Capture the fit's own warnings rather than letting them scroll past.
      # fit.variogram() reports non-convergence as a warning and returns a
      # value anyway, so without this a failed fit is indistinguishable from a
      # good one - and on NL_EXPL_DRL_RA that mattered: nugget/sill piled up on
      # exactly 0.00 and exactly 1.00 while Moran's I said 0.001 on the same
      # residuals. Those are boundary solutions, not measurements.
      warns <- character(0)
      fit <- withCallingHandlers(
        gstat::fit.variogram(V, gstat::vgm("Exp"), warn.if.neg = FALSE),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      nug <- if ("Nug" %in% fit$model) fit$psill[fit$model == "Nug"] else 0
      part <- sum(fit$psill[fit$model != "Nug"])
      list(nugget = nug, partial.sill = part, sill = nug + part,
           range = max(fit$range), nugget.ratio = nug / (nug + part),
           model = as.character(fit$model[fit$model != "Nug"])[1],
           singular = isTRUE(attr(fit, "singular")),
           converged = !any(grepl("convergence", warns, ignore.case = TRUE)),
           warnings = warns,
           # The lag geometry, so the fitted range can be judged against the
           # distances the empirical variogram actually saw. A range outside
           # them is extrapolation, not a measurement.
           first.lag = min(V$dist), last.lag = max(V$dist), n.lags = nrow(V))
    }, silent = TRUE)
    if (!inherits(checks$variogram, "try-error"))
      message(sprintf(
        "  fitted %s variogram: nugget %.3g, sill %.3g, range %.3g km, nugget/sill %.2f",
        checks$variogram$model, checks$variogram$nugget, checks$variogram$sill,
        checks$variogram$range, checks$variogram$nugget.ratio))
    else
      message("  variogram fit failed: ",
              conditionMessage(attr(checks$variogram, "condition")))

    V <- (gstat::variogram(E ~ 1, mydata, cutoff = 100))
    plot(
      x = V$dist,
      y = V$gamma,
      xlab = "Distance (km)",
      ylab = "Semi-variance",
      pch = 16,
      cex = 2 * V$np / max(V$np)
    )


    V <- (gstat::variogram(E ~ 1, mydata, cutoff = 10))
    plot(
      x = V$dist,
      y = V$gamma,
      xlab = "Distance (km)",
      ylab = "Semi-variance",
      pch = 16,
      cex = 2 * V$np / max(V$np)
    )

    V <- (gstat::variogram(E ~ 1, mydata, cutoff = 5))
    plot(
      x = V$dist,
      y = V$gamma,
      xlab = "Distance (km)",
      ylab = "Semi-variance",
      pch = 16,
      cex = 2 * V$np / max(V$np)
    )

    V <- (gstat::variogram(E ~ 1, mydata, cutoff = 2))
    plot(
      x = V$dist,
      y = V$gamma,
      xlab = "Distance (km)",
      ylab = "Semi-variance",
      pch = 16,
      cex = 2 * V$np / max(V$np)
    )
  }
  # Turn partial match warnings back on
  options(warnPartialMatchDollar = TRUE)

  checks$brief <- brief
  invisible(checks)
}


#' Default thresholds for interpreting DSM checking results
#'
#' Every boundary the interpretation uses, in one place so it can be seen and
#' overridden. Pass a partial list to \code{\link{interpret.dsm.checks}} to
#' change individual entries.
#'
#' Each entry is \code{c(watch, problem)} and is read as "ok below the first,
#' watch between, problem above" unless the check says otherwise.
#'
#' @return Named list of thresholds.
#' @export
default.check.thresholds <- function() {
  list(
    dispersion.hi   = c(1.1, 1.5),   # DHARMa ratio, above 1
    dispersion.lo   = c(0.9, 0.7),   # and below 1 - under-dispersion
    uniformity.d    = c(0.02, 0.05), # KS statistic
    zeroinfl        = c(1.05, 1.2),  # ratio of observed to simulated zeros
    outliers        = c(1.5, 3),     # observed / expected outlier frequency
    kindex          = c(0.9, 0.8),   # k-index, lower is worse
    edf.ratio       = c(0.5, 0.8),   # edf/k', only damning together with kindex
    concurvity      = c(0.5, 0.8),   # "worst" measure
    morans.i        = c(0.05, 0.15)  # residual spatial autocorrelation
    # No nugget.ratio: the variogram is reported, never scored - see the
    # variogram block of interpret.dsm.checks() for the measurements.
  )
}


#' Turn DSM checking results into plain-English verdicts
#'
#' Reads the list returned by \code{\link{check.dsm}} and says, for each
#' diagnostic, what the number means and whether it is a problem. One row per
#' check, with a verdict of \code{"ok"}, \code{"watch"} or \code{"problem"}.
#'
#' @section Why p-values are reported but never scored:
#' These models are fitted to 100,000+ segment rows, and at that size a
#' significance test detects departures far too small to matter - ATPU's
#' dispersion is 1.22 with p < 2.2e-16, which is a 22% effect reported as
#' overwhelming evidence. Every verdict below is therefore taken from the effect
#' size. The p-values are printed alongside because they are what the underlying
#' tests report, not because they carry the decision.
#'
#' @section Basis size needs two numbers, not one:
#' \code{gam.check}'s advice is "low p-value (k-index < 1) may indicate that k is
#' too low, **especially if edf is close to k'**", and the second half is the half
#' that matters. ATPU's \code{s(x,y)} smooths have k-index 0.76 with p < 2e-16,
#' which reads as a clear problem until you notice edf is 16-27 against a k' of
#' 99 - the basis is nowhere near saturated and there is nothing to fix. So a
#' smooth is only flagged when the k-index is low *and* the basis is being used
#' up. That is the specific wrong answer this function exists to avoid.
#'
#' @param checks List returned by \code{\link{check.dsm}}.
#' @param thresholds Named list overriding \code{\link{default.check.thresholds}};
#'   only the entries given are replaced.
#' @return A tibble with columns \code{check}, \code{statistic}, \code{value},
#'   \code{verdict} and \code{comment}, with the model name in the
#'   \code{"modname"} attribute.
#' @export
interpret.dsm.checks <- function(checks, thresholds = list()) {
  checkmate::expect_list(checks)
  checkmate::expect_list(thresholds)
  th <- utils::modifyList(default.check.thresholds(), thresholds)

  rows <- list()
  add <- function(check, statistic, value, verdict, comment)
    rows[[length(rows) + 1]] <<- dplyr::tibble(
      check = check, statistic = statistic, value = value,
      verdict = verdict, comment = comment)

  ok <- function(x) !is.null(x) && !inherits(x, "try-error")
  # "ok below watch, problem above problem", for a statistic where bigger is worse
  band.hi <- function(v, t) if (v >= t[2]) "problem" else if (v >= t[1]) "watch" else "ok"
  # and where smaller is worse
  band.lo <- function(v, t) if (v <= t[2]) "problem" else if (v <= t[1]) "watch" else "ok"

  ## ---- dispersion ----------------------------------------------------------
  if (ok(checks$residuals$dispersion)) {
    d <- unname(checks$residuals$dispersion$statistic)
    p <- checks$residuals$dispersion$p.value
    v <- if (d >= 1) band.hi(d, th$dispersion.hi) else band.lo(d, th$dispersion.lo)
    add("Dispersion", "DHARMa ratio", sprintf("%.3f", d), v, sprintf(
      "Residual spread is %.0f%% %s than the fitted model implies%s. (p = %s, not used - see note on n.)",
      abs(d - 1) * 100, if (d >= 1) "wider" else "narrower",
      switch(v, ok = ", which is within tolerance",
             watch = " - worth noting but not disqualifying",
             problem = " - the mean-variance relationship is wrong"),
      format.pval(p, digits = 2)))
  }

  ## ---- uniformity ----------------------------------------------------------
  if (ok(checks$residuals$uniformity)) {
    D <- unname(checks$residuals$uniformity$statistic)
    v <- band.hi(D, th$uniformity.d)
    add("Uniformity", "KS D", sprintf("%.4f", D), v, sprintf(
      "The largest gap between the residual distribution and uniform is %.1f%%%s.",
      D * 100,
      switch(v, ok = ", so the residuals are distributed as they should be",
             watch = ", a mild departure from uniform",
             problem = ", so the distributional assumption is not holding")))
  }

  ## ---- outliers ------------------------------------------------------------
  if (ok(checks$residuals$outliers)) {
    o <- checks$residuals$outliers
    ratio <- unname(o$estimate) / unname(o$null.value)
    v <- band.hi(ratio, th$outliers)
    add("Outliers", "observed / expected", sprintf("%.2f", ratio), v, sprintf(
      "%.2f%% of points fall outside the simulated range against %.2f%% expected%s.",
      unname(o$estimate) * 100, unname(o$null.value) * 100,
      switch(v, ok = " - unremarkable", watch = " - a mild excess",
             problem = " - a real excess of extreme values")))
  }

  ## ---- zero inflation ------------------------------------------------------
  if (ok(checks$zeroinflation)) {
    z <- unname(checks$zeroinflation$statistic)
    v <- band.hi(z, th$zeroinfl)
    add("Zero inflation", "obs / sim zeros", sprintf("%.3f", z), v, sprintf(
      "The data hold %.1f%% %s zeros than the model simulates%s.",
      abs(z - 1) * 100, if (z >= 1) "more" else "fewer",
      switch(v, ok = " - the zero behaviour is captured",
             watch = " - a mild excess of zeros",
             problem = " - consider a zero-inflated or hurdle formulation")))
  }

  ## ---- basis size ----------------------------------------------------------
  if (ok(checks$k.check) && is.matrix(checks$k.check)) {
    kc <- as.data.frame(checks$k.check)
    names(kc) <- c("kprime", "edf", "kindex", "pvalue")[seq_len(ncol(kc))]
    kc$ratio <- kc$edf / kc$kprime
    # Only a problem when the basis is BOTH poorly fitting and nearly used up.
    kc$v <- mapply(function(ki, r) {
      if (ki <= th$kindex[2] && r >= th$edf.ratio[2]) "problem"
      else if (ki <= th$kindex[1] && r >= th$edf.ratio[1]) "watch"
      else "ok"
    }, kc$kindex, kc$ratio)
    # Worst verdict first, then the most-used basis within it. Sorting by lowest
    # k-index instead would report a smooth penalised to 0 edf as the worst case,
    # which is exactly the row where a low k-index cannot matter.
    worst <- kc[order(match(kc$v, c("problem", "watch", "ok")), -kc$ratio), ][1, ]
    v <- worst$v
    add("Basis size (k)", "worst k-index [edf/k']",
        sprintf("%.2f [%.2f]", worst$kindex, worst$ratio), v, sprintf(
      "Lowest k-index is %.2f on %s, using %.0f of %.0f available df (%.0f%%)%s.",
      worst$kindex, rownames(worst), worst$edf, worst$kprime, worst$ratio * 100,
      switch(v,
        ok = ". A low k-index alone is not evidence of too-small k when the basis is barely used, which is the case here",
        watch = ". Worth watching - the basis is moderately used and fitting imperfectly",
        problem = ". The basis is nearly exhausted and fitting badly, so k is genuinely too low")))
  }

  ## ---- concurvity ----------------------------------------------------------
  #
  # High, real, and reported rather than scored - but not for the reason it
  # first looks like.
  #
  # The tempting explanation is that a by = Season smooth is zero outside its own
  # level, so the other levels can mimic it and inflate the measure structurally.
  # That is wrong, and measurably so: the pairwise concurvity between a smooth
  # and its own siblings in other seasons is exactly 0.000 for every term on
  # ATPU, which it must be, since they are supported on disjoint rows.
  #
  # What is actually entangled is different covariates within the SAME season.
  # Every term's worst partner is one of those, and the reason is geography:
  # depth is a smooth function of position, so s(depth):SeasonX and
  # s(x,y):SeasonX compete for one signal. Regressing the covariates on location
  # for ATPU gives depth 98.4% of deviance explained, depth.g 50.2%, sst 35.2% -
  # and the concurvity ranks the same way, s(x,y) vs s(depth) at 0.978 down to
  # s(depth.g):SeasonWinter at 0.524.
  #
  # So this is a genuine caveat: the individual shapes of the depth, depth.g and
  # sst smooths cannot be read as separate effects. It is not scored because it
  # is a property of the covariate set rather than of a fitted model - it says
  # nothing about whether this model is adequate, and it cannot separate the two
  # family finalists, which share a formula and therefore a model matrix.
  if (ok(checks$concurvity) && is.matrix(checks$concurvity)) {
    keep <- colnames(checks$concurvity) != "para"
    if (any(keep)) {
      est <- checks$concurvity["estimate", keep]
      wrst <- checks$concurvity["worst", keep]
      add("Concurvity", "worst estimate", sprintf("%.3f", max(est, na.rm = TRUE)),
          "reference", sprintf(paste(
            "Most entangled term is %s (estimate %.2f, worst %.2f). Real, and",
            "driven by depth being close to a function of position, so the depth",
            "and location smooths compete for one signal - read their individual",
            "shapes with that in mind. Not scored: it describes the covariate",
            "set, not this model's adequacy, and it cannot separate two finalists",
            "that share a formula."),
            names(est)[which.max(est)], max(est, na.rm = TRUE),
            wrst[which.max(est)]))
    }
  }

  ## ---- spatial autocorrelation ---------------------------------------------
  if (ok(checks$spatial)) {
    I <- unname(checks$spatial$statistic["observed"])
    v <- band.hi(abs(I), th$morans.i)
    tested <- checks$spatial$n.tested
    total <- checks$spatial$n.locations
    add("Spatial autocorrelation", "Moran's I", sprintf("%.4f", I), v, sprintf(
      "Residual spatial correlation is %.3f%s%s.", I,
      switch(v, ok = ", effectively none - the spatial smooth has absorbed the structure",
             watch = " - some structure remains unmodelled",
             problem = " - substantial structure remains, so the spatial term is not capturing it"),
      if (!is.null(tested) && !is.null(total) && tested < total)
        sprintf(" (a random %s of %s survey locations)",
                format(tested, big.mark = ","), format(total, big.mark = ",")) else ""))
  } else {
    add("Spatial autocorrelation", "Moran's I", NA_character_, "not run",
        "The test could not be computed for this model.")
  }

  ## ---- variogram (reported, never scored) ---------------------------------
  # Reference-only, for the same reason as concurvity: it describes the data
  # geometry rather than this model's adequacy, and at this design it is not
  # identified. Measured on NL_EXPL_DRL_RA, 22 models over 38,147 locations:
  #
  #   - the empirical variogram is flat and noisy (gamma 5.08, 5.88, 6.43,
  #     6.01, 4.99 ... 3.63 over 15 lags), which is what Moran's I ~ 0.001 on
  #     the same residuals says: there is no residual structure to fit
  #   - so fit.variogram() wanders to the parameter boundaries. 7 of 22 were
  #     singular, 13 of 22 did not converge, 18 of 22 were one or the other
  #   - 7 fitted a range of 5-11 km when the FIRST lag centre is 19.4 km and
  #     bins are 30 km wide, so nugget and structure cannot be separated at
  #     all; 5 fitted 2,155-3,507 km against a last lag of 465 km
  #   - 21 of 22 were degenerate or out of range. Scoring them produced 10
  #     "problem" verdicts that contradicted the spatial test on the same
  #     residuals, and the single survivor scored "watch" on the luck of a fit
  #
  # A check that resolves 1 of 22 is not a check. Moran's I answers the same
  # question directly, with a permutation reference, and is scored. The five
  # variogram plots stay, because reading them by eye is still worth doing.
  if (ok(checks$variogram)) {
    vg  <- checks$variogram
    nr  <- vg$nugget.ratio
    # first.lag/last.lag are absent from assessments saved before this check.
    have.lags <- !is.null(vg$first.lag)
    flags <- c(
      if (isTRUE(vg$singular)) "the fit is singular",
      if (identical(vg$converged, FALSE)) "the fit did not converge",
      if (have.lags && vg$range < vg$first.lag)
        sprintf("the fitted range (%.3g km) is below the first lag (%.3g km)",
                vg$range, vg$first.lag),
      if (have.lags && vg$range > vg$last.lag)
        sprintf("the fitted range (%.3g km) is beyond the last lag (%.3g km)",
                vg$range, vg$last.lag))

    add("Variogram", "nugget / sill", sprintf("%.2f", nr), "reference", paste(
      sprintf("%.0f%% of residual variance is at zero distance, range %.3g km.",
              nr * 100, vg$range),
      if (length(flags))
        sprintf("Not identified here - %s.", paste(flags, collapse = "; ")) else
        "This fit is identified, which is the exception.",
      "Reported, not scored: the exponential fit is unidentified for almost",
      "every model at this design. Read the Moran's I row above for the",
      "residual spatial structure question, and the plots by eye."))
  }

  ## ---- Pearson overdispersion (reported, never scored) ---------------------
  if (ok(checks$overdispersion))
    add("Pearson overdispersion", "chi-sq / df",
        sprintf("%.2f", checks$overdispersion), "reference", paste(
      "Reported for continuity, not scored. The scale parameter is estimated",
      "rather than fixed at 1 for the tw() and nb() families used here, so this",
      "is not on the same footing as DHARMa's dispersion ratio above and the two",
      "routinely disagree."))

  out <- dplyr::bind_rows(rows)
  attr(out, "modname") <- checks$modname
  attr(out, "n") <- checks$n
  out
}


#' Compare two family finalists on their diagnostics
#'
#' The model selection step (\code{02.05}) chooses the family by spatial-block
#' cross-validation, and where the families are indistinguishable it keeps
#' whatever \code{final.dsm.models} already names. This looks at the two
#' finalists' diagnostics and says whether the one that was not chosen is
#' materially better - evidence that could overturn a tie.
#'
#' @section What it will and will not overturn:
#' A decisive cross-validation verdict is an out-of-sample result and is not
#' overturned by in-sample diagnostics: where \code{decisive} is \code{TRUE} this
#' reports the disagreement and recommends nothing. Only \code{[tie-break]} and
#' \code{[kept]} entries - where CV could not separate the families - are open to
#' being changed.
#'
#' The rule is deliberately a stated comparison rather than a weighted score, so
#' that a recommendation can be argued with: the challenger wins only if it has
#' strictly fewer \code{problem} verdicts and no more \code{watch} verdicts, or
#' the same number of \code{problem} and at least \code{watch.margin} fewer
#' \code{watch}. Checks marked \code{reference} never count - see
#' \code{\link{interpret.dsm.checks}} for why concurvity and Pearson
#' overdispersion are among them.
#'
#' @param incumbent Tibble from \code{\link{interpret.dsm.checks}} for the model
#'   \code{final.dsm.models} currently names.
#' @param challenger Tibble from \code{\link{interpret.dsm.checks}} for the other
#'   family finalist.
#' @param decisive Logical; was the CV verdict decisive for this species?
#' @param watch.margin Integer; how many fewer \code{watch} verdicts the
#'   challenger needs when the \code{problem} counts are equal. Default 2.
#' @return A one-row tibble: \code{recommend} (logical), \code{winner},
#'   \code{reason}, and the two verdict tallies.
#' @export
compare.dsm.finalists <- function(incumbent, challenger, decisive,
                                  watch.margin = 2) {
  checkmate::expect_data_frame(incumbent)
  checkmate::expect_data_frame(challenger)
  checkmate::expect_flag(decisive)
  checkmate::expect_count(watch.margin)

  tally <- function(x, what) sum(x$verdict == what, na.rm = TRUE)
  inc <- c(problem = tally(incumbent, "problem"), watch = tally(incumbent, "watch"))
  cha <- c(problem = tally(challenger, "problem"), watch = tally(challenger, "watch"))
  inc.name <- attr(incumbent, "modname")
  cha.name <- attr(challenger, "modname")

  better <- (cha["problem"] < inc["problem"] && cha["watch"] <= inc["watch"]) ||
            (cha["problem"] == inc["problem"] &&
               cha["watch"] <= inc["watch"] - watch.margin)

  tallies <- sprintf("%s has %d problem / %d watch, %s has %d problem / %d watch",
                     inc.name, inc["problem"], inc["watch"],
                     cha.name, cha["problem"], cha["watch"])

  if (decisive) {
    rec <- FALSE
    reason <- if (better)
      sprintf(paste("CV was decisive for %s, so the model stands. Note the",
                    "diagnostics disagree: %s. Worth a look, but an out-of-sample",
                    "result is not overturned by in-sample diagnostics."),
              inc.name, tallies)
    else
      sprintf("CV was decisive and the diagnostics agree: %s.", tallies)
  } else if (better) {
    rec <- TRUE
    reason <- sprintf(paste("%s [chosen]: CV could not separate the families, and",
                            "%s has the better diagnostics - %s."),
                      get0("DSM_DIAGNOSTIC_TAG", ifnotfound = "diag-override"),
                      cha.name, tallies)
  } else {
    rec <- FALSE
    reason <- sprintf(paste("CV could not separate the families and the",
                            "diagnostics do not either: %s. Existing entry stands."),
                      tallies)
  }

  dplyr::tibble(
    recommend = rec,
    winner = if (rec) cha.name else inc.name,
    incumbent = inc.name, challenger = cha.name,
    inc.problem = unname(inc["problem"]), inc.watch = unname(inc["watch"]),
    cha.problem = unname(cha["problem"]), cha.watch = unname(cha["watch"]),
    decisive = decisive,
    reason = reason)
}



#' Add observation counts and estimated abundance to segment data
#'
#' Joins detection-probability-corrected observation counts from
#' \code{distdata} onto \code{segdata}, filling segments with no observations
#' with zeros, and adds \code{rawCount}, \code{estAbund}, and \code{estDens}
#' columns.
#'
#' @param segdata \code{sf} segment data frame.
#' @param distdata Observation data frame already filtered to the species of
#'   interest and augmented with \code{adjSize}.
#' @return \code{segdata} augmented with \code{rawCount}, \code{estAbund},
#'   and \code{estDens} columns.
#' @export
augment.segdata <- function(segdata, distdata) {
  newsegdata <- distdata %>%
    dplyr::group_by(Sample.Label) %>%
    dplyr::summarize(estAbund = sum(adjSize), rawCount = sum(size)) %>%
    dplyr::right_join(segdata, by = "Sample.Label") %>%
    sf::st_as_sf()
  newsegdata$estAbund[is.na(newsegdata$estAbund)] <- 0
  newsegdata$rawCount[is.na(newsegdata$rawCount)] <- 0
  newsegdata$estDens <- newsegdata$estAbund/newsegdata$segment.area
  newsegdata
}


#' Adjust smooth k for temporal covariates when data are sparse
#'
#' Inspects \code{form} for \code{s(year, ...)} and \code{s(yday, ...)}
#' terms.  If the number of unique values in \code{segdata} is less than
#' \code{k}, the term is either removed (fewer than 2 unique values) or
#' replaced with a version using a reduced \code{k}.
#'
#' @param form A model formula.
#' @param segdata Segment data frame containing \code{year} and/or
#'   \code{yday} columns.
#' @param k Default basis dimension; smooth terms with fewer than \code{k}
#'   unique values are adjusted.
#' @return Updated formula.
#' @export
adjust.time.covars <- function(form, segdata, k = 10) {
  # check if form contains year or yday
  # check unique number of values for each
  # if either is less than default 10 then

  # check for year
  new.k <- length(unique(segdata$year))
  if (grepl("s(year", as.character(form)[3], fixed = TRUE) &&
      (new.k < k)) {
    message("Removing year from formula for lack of data")
    form <- update(form, . ~ . - mgcv::s(year, bs = "ts"))

    # enough data to have term at all?
    if (new.k >= 2) {
      message("\t... and replacing it with smooth where k = ", new.k)
      form %<>% as.character
      form <- sprintf('%s ~ %s + s(year, bs = "ts", k = %d)', form[2],
                      form[3], new.k) %>%
        as.formula
    }
  }

  # check for yday
  new.k <- length(unique(segdata$yday))
  if (grepl("s(year", as.character(form)[3], fixed = TRUE) &&
      (new.k < k)) {
    message("Removing yday from formula for lack of data")
    form <- update(form, . ~ . - mgcv::s(yday, bs = "ts"))

    # enough data to have term at all?
    if (new.k >= 2) {
      message("\t... and replacing it with smooth where k = ", new.k)
      form %<>% as.character
      form <- sprintf('%s ~ %s + s(yday, bs = "ts", k = %d)', form[2],
                      form[3], new.k) %>%
        as.formula
    }
  }

  form
}

#' Fit or reload a single DSM model
#'
#' Calls \code{\link[dsm]{dsm}} with the specification in \code{mod.def} and
#' saves the result to \code{folder/<modname>.Rdata}.  If
#' \code{rerun.dsms = FALSE} and the file already exists, the saved result is
#' loaded instead.
#'
#' @param mod.def Single-row data frame with list columns \code{formula} and
#'   \code{family} and a character column \code{modname}.
#' @param ddf.obj List of detection function objects in the order expected by
#'   \code{\link[dsm]{dsm}}.
#' @param segment.data \code{sf} or data frame of segment data.
#' @param observation.data Data frame of observation data.
#' @param method Smoothing parameter estimation method; default
#'   \code{"REML"}.
#' @param convert.units Numeric conversion factor passed to
#'   \code{\link[dsm]{dsm}}.
#' @param control List of control options passed to \code{\link[dsm]{dsm}};
#'   default keeps data with \code{keepData = TRUE}.
#' @param folder Directory path for saving the model \code{.Rdata} file.
#' @param rerun.dsms If \code{TRUE} (default), always refit; if \code{FALSE},
#'   load a previously saved result when available.
#' @param ... Additional arguments passed to \code{\link[dsm]{dsm}}.
#' @return Fitted \code{dsm} object, or a \code{try-error} on failure.
#' @export
run.dsm.model <- function(mod.def,
                          ddf.obj,
                          segment.data,
                          observation.data,
                          method = "REML",
                          convert.units = 1,
                          control = list(keepData = T),
                          folder,
                          rerun.dsms = TRUE,
                          ...) {
  # logfileConn <- file(description = paste0("E:/test_", Sys.getpid(), ".txt"), open = "at")
  # source(here::here("R/analysis_settings.R"))
  # logfileConn <- file(here::here(ResultsDir,
  #                         paste0("dsm_logfile_", Sys.getpid(), ".txt")), "at")
  #
  # needed if being called from future_map() on a worker process. If not,
  # it messes up trying to access any sf object (segment.data) when using
  # the s2 spherical geometry package.
  sf::sf_use_s2(FALSE)

  filename <- file.path(folder, paste0(mod.def$modname, ".Rdata"))

  # Can we reuse a previously saved fit? A file holding a try-error records a
  # failure, not a result, so treat it as if it were absent. Reusing one would
  # make a single bad run permanent: every later rerun.dsms = FALSE call would
  # reload the same error in seconds and never refit, with nothing in the
  # return value to distinguish that from a fit that had just failed.
  model <- NULL
  if (rerun.dsms == FALSE && file.exists(filename)) {
    # Load saved model result from file
    message(sprintf(
      "Loading previously saved model results for %s.",
      mod.def$modname
    ))
    load(filename)
    if (inherits(model, "try-error")) {
      message(sprintf(
        paste0("run.dsm.model: %s holds a failed fit from an earlier run. ",
               "Refitting rather than reusing it."),
        basename(filename)
      ))
      model <- NULL
    }
  }

  if (is.null(model)) {
    # Rerun the dsm
    message("Running dsm model ", mod.def$modname)

    # Adjust formula
    form <- adjust.time.covars(mod.def$formula[[1]], segment.data)

    # A factor needs at least two levels to be used as a model term. mgcv would
    # otherwise fail with the opaque "contrasts can be applied only to factors
    # with 2 or more levels", so check here where the cause can be named.
    #
    # Note this does NOT fire merely because a study area has only one survey
    # type: platform is bird behaviour (W/F), not survey type, so it keeps both
    # levels in that case. It fires when a species was never recorded in one
    # behaviour class, which would break every candidate model that has a
    # platform term (currently all of them).
    if ("platform" %in% all.vars(form) && !is.null(segment.data$platform)) {
      plat_levels <- unique(stats::na.omit(as.character(segment.data$platform)))
      if (length(plat_levels) < 2)
        stop(sprintf(
          paste0("run.dsm.model: model '%s' has a platform term but segdata ",
                 "has only the platform level %s. A factor needs >= 2 levels. ",
                 "This species has no observations in the other behaviour ",
                 "class; use a model without a platform term (the _nofactor ",
                 "variants in dsm.mod.specs) or drop this species."),
          mod.def$modname, paste(sQuote(plat_levels), collapse = ", ")))
    }

    # Call dsm(). Note the list indexing for formula and family since these are
    # list columns in mod.def. Temporarily disable warnings about partial
    # matches of list/dataframe names since mgcv is full of such things.
    options(warnPartialMatchDollar = FALSE)
    print(system.time(model <- try(dsm::dsm(
      formula = form,
      ddf.obj = ddf.obj,
      segment.data = segment.data,
      observation.data = observation.data,
      family = mod.def$family[[1]],
      method = method,
      control = control,
      convert.units = convert.units,
      segment.area = segment.data$segment.area,
      ...
    ),
    outFile = stdout())))
    options(warnPartialMatchDollar = TRUE)


    if(!inherits(model, "try-error")){
      # if we used "bam" then the data is not kept even if keepData == TRUE, so
      # add it back in. Sometimes "data" is in model but it is NA.
      #
      # The test must stay scalar: model$data is normally a data frame, so
      # is.na() on it returns a matrix and if() then fails with "the condition
      # has length > 1". Testing for a usable data frame covers all three bad
      # cases (absent, NULL, or a bare NA) in one scalar expression.
      kd <- control[["keepData"]]
      data_ok <- is.data.frame(model$data) && nrow(model$data) > 0
      if (isTRUE(kd) && !data_ok) {
        model$data <-
          dsm:::make.data(
            response = as.character(mod.def$formula[[1]])[2],
            segdata = segment.data,
            obsdata = observation.data,
            ddfobject = ddf.obj,
            family = mod.def$family[[1]],
            group = FALSE,
            convert.units = 1,
            availability = 1,
            segment.area = segment.data$segment.area
          )
      }

      # check if my computed response is equal to that from dsm, just to make sure
      # we understand how response is being computed. Note the "as.data.frame"
      # is used to remove the atribs from segment.data b/c it's an sf object.
      stopifnot(
        all.equal(
          model$data %>% dplyr::arrange(Sample.Label) %>% magrittr::extract("estAbund"),
          segment.data %>% dplyr::arrange(Sample.Label) %>% as.data.frame %>% magrittr::extract("estAbund")
        )
      )
    } else { # Model failed to fit.
      message("run.dsm.model: dsm() failed: ", model)
    }

    # Only a real fit is worth caching. Saving a try-error here is what made
    # the failure permanent, so leave the file alone instead: if an earlier
    # good fit is on disk it stays usable, and if nothing is there the next
    # run retries from scratch.
    if (inherits(model, "try-error")) {
      message(sprintf(
        "run.dsm.model: not saving failed fit for %s to %s",
        mod.def$modname, filename
      ))
    } else {
      message(sprintf("Saving model result to %s", filename))
      save(model, file = filename, compress = FALSE)
    }
  } # End rerun dsm

  model
}


#' Check that summing the platform levels reconstitutes the total density
#'
#' The multi-ddf design replicates the segment table once per ddftype, each copy
#' carrying the segment's full area, and recovers the total by summing the
#' per-platform predictions ([dsm.pred()]). That only works when each platform
#' level is carried by exactly one segdata copy per segment. If several ddftypes
#' share a level, the GAM has nothing to tell those copies apart, fits their mean
#' rather than letting them sum, and every density comes out divided by the
#' number of copies - silently.
#'
#' This compares the model-implied combined density against a design-based
#' estimate over the surveyed area, using only the fitted model. It needs no
#' predictions and is cheap, so it can run before the expensive steps.
#'
#' The two are not expected to agree exactly: the GAM is penalised and its fitted
#' values need not sum to the observed total. A ratio near 1 is the pass; the
#' failure this exists to catch is a ratio near \code{1 / copies.per.platform},
#' which was 0.494 on NL_EXPL_DRL_RA ATPU before \code{ddftype_to_platform} gave
#' each ddftype its own level.
#'
#' @section Projects with both survey types:
#' Ship and aerial contribute different numbers of segments, and a ddftype can
#' survive for one survey type and not the other, so the levels need not all rest
#' on the same segments. Two consequences:
#'
#' \code{copies.per.platform} is unaffected and remains the structural test. It
#' is a maximum over per-(segment, platform) counts, and a segment belongs to one
#' survey type, so it can only exceed 1 when two ddftypes of the *same* survey
#' type share a level - the collapse being tested for.
#'
#' \code{ratio} is weaker when the levels cover different segments, because
#' \code{design.density} spreads every bird over the whole surveyed area while a
#' level present for only one survey type is fitted over that survey type's part
#' of it. The message says so when it happens; read
#' \code{copies.per.platform} rather than the ratio in that case.
#'
#' @param model Fitted \code{dsm} object with \code{control = list(keepData =
#'   TRUE)}, so \code{model$data} carries \code{abundance.est},
#'   \code{segment.area}, \code{platform}, \code{ddftype_orig} and
#'   \code{Sample.Label}.
#' @param tolerance Numeric; how far the ratio may sit from 1 before a warning is
#'   issued. Defaults to 0.15.
#' @return Invisibly, a one-row tibble with \code{copies.per.segment},
#'   \code{copies.per.platform}, \code{surveyed.area.sqkm},
#'   \code{design.density}, \code{model.density} and \code{ratio}.
#' @export
check.platform.combination <- function(model, tolerance = 0.15) {
  checkmate::expect_class(model, "dsm")
  checkmate::expect_number(tolerance, lower = 0)

  dat <- model$data
  if (is.null(dat) || is.null(dat$abundance.est))
    stop("check.platform.combination: model$data has no abundance.est. Fit with ",
         "control = list(keepData = TRUE).")
  dat <- sf::st_drop_geometry(dat)
  dat$.fit <- stats::fitted(model)
  dat$.seg <- segment.id.from.label(dat$Sample.Label)

  n.segments <- dplyr::n_distinct(dat$.seg)
  copies.per.segment <- nrow(dat) / n.segments

  # The count that matters. Anything above 1 means several ddftypes share a
  # platform level and the sum below will under-report by that factor.
  #
  # This one is safe in a project with both survey types even though they may
  # contribute different numbers of segments: it is a maximum over per-(segment,
  # platform) counts, so how many segments each survey type has cannot affect it.
  # A ship segment appears in the S ddftypes only and an aerial segment in the A
  # ddftypes only, so neither can reach 2 unless two ddftypes of the SAME survey
  # type share a platform level - which is exactly the collapse being tested for.
  copies.per.platform <- dat %>%
    dplyr::count(.seg, platform) %>%
    dplyr::pull(n) %>%
    max()

  # Surveyed area is the area of the DISTINCT segments.
  #
  # It used to be sum(segment.area) / copies.per.segment, which is right only
  # when every segment has the same number of copies. That fails in a project
  # with both survey types where the two keep different numbers of ddftypes -
  # say ship keeps all four and aerial loses AWN for want of observations. Then
  # copies.per.segment is a count-weighted average while the numerator is an
  # area-weighted sum, and the two do not cancel. With 200 ship segments of
  # 0.30 km2 and 50 aerial of 0.80, copies.per.segment is 3.8 and the old
  # expression returns 94.74 km2 for a true 100 - and sprintf("%d", 3.8) then
  # errors outright. Taking the distinct segments' area is exact however the
  # copies fall.
  seg.area <- dat %>%
    dplyr::group_by(.seg) %>%
    dplyr::summarise(area = dplyr::first(segment.area),
                     n_areas = dplyr::n_distinct(segment.area),
                     .groups = "drop")
  if (any(seg.area$n_areas > 1))
    stop("check.platform.combination: segment.area differs between copies of ",
         "the same segment, so there is no single surveyed area to compare ",
         "against. create.segdata.copy() carries it through unchanged, so this ",
         "means something downstream has altered it.")
  surveyed.area <- sum(seg.area$area)
  design.density <- sum(dat$abundance.est) / surveyed.area

  # Each platform level's density is its fitted total over the area of ITS
  # copies; summing across levels is what dsm.pred() does on the grid.
  by.platform <- dat %>%
    dplyr::group_by(platform) %>%
    dplyr::summarise(fit = sum(.fit), area = sum(segment.area),
                     .groups = "drop") %>%
    dplyr::mutate(dens = fit / area)
  model.density <- sum(by.platform$dens)

  ratio <- model.density / design.density

  # Whether the levels rest on the same ground. They do whenever every segment
  # appears in every level, which is the usual case. They do not when a ddftype
  # survives for one survey type but not the other: that level is then estimated
  # over a subset of the study area while design.density spreads its birds over
  # all of it, so the ratio mixes two spatial supports and loses precision as a
  # calibration measure. copies.per.platform is unaffected and remains the
  # structural test.
  mixed.support <- !isTRUE(all.equal(range(by.platform$area)[1],
                                     range(by.platform$area)[2]))

  copies.txt <- if (isTRUE(all.equal(copies.per.segment,
                                     round(copies.per.segment))))
    sprintf("%g", copies.per.segment)
  else
    sprintf("%.2f on average (levels do not all cover the same segments)",
            copies.per.segment)

  message(sprintf(
    paste0("Platform combination check: %s segment copies over %d platform ",
           "level(s) (%s), %d copies per platform."),
    copies.txt, nrow(by.platform),
    paste(by.platform$platform, collapse = ", "), copies.per.platform))
  message(sprintf(
    "  design-based %.4f /sqkm, model-implied %.4f /sqkm, ratio %.3f",
    design.density, model.density, ratio))

  if (mixed.support)
    message("  NB: platform levels cover different segments (areas ",
            paste(sprintf("%s %.0f", by.platform$platform, by.platform$area),
                  collapse = ", "),
            " sqkm), so the ratio mixes spatial supports - read ",
            "copies.per.platform, not the ratio.")

  if (copies.per.platform > 1)
    warning(sprintf(
      paste0("check.platform.combination: %d segdata copies share each platform ",
             "level, so the combined prediction is roughly 1/%d of the true ",
             "density. Give each ddftype its own level in ddftype_to_platform."),
      copies.per.platform, copies.per.platform), immediate. = TRUE)
  else if (abs(ratio - 1) > tolerance && !mixed.support)
    warning(sprintf(
      paste0("check.platform.combination: model-implied density is %.3f x the ",
             "design-based estimate (tolerance %.2f). Copies per platform is 1, ",
             "so this is not the platform-collapse failure - look at the fit."),
      ratio, tolerance), immediate. = TRUE)

  invisible(dplyr::tibble(
    copies.per.segment = copies.per.segment,
    copies.per.platform = copies.per.platform,
    mixed.support = mixed.support,
    surveyed.area.sqkm = surveyed.area,
    design.density = design.density,
    model.density = model.density,
    ratio = ratio))
}


#' Generate density predictions from a DSM model
#'
#' Calls \code{predict} on the model named \code{modname} from \code{mod.res}
#' for all rows of \code{predgrid} (all seasons × platform levels), computes
#' a \code{"Combined"} subset by summing across platform levels, writes a
#' combined shapefile, and saves four-season rasters as GeoTIFFs.
#'
#' Predictions are made with \code{discrete = FALSE} so that a cell's predicted
#' density depends only on that cell's covariates and the model, not on the
#' composition of the grid it happens to be predicted alongside. Note this does
#' not reach \code{dsm::dsm_var_gam()} in the variance step, which takes no
#' \code{discrete} argument - only refitting without discrete does.
#'
#' Expects project globals \code{ShapeDir}, \code{predDir}, \code{season.names},
#' and \code{predgridCellArea}.
#'
#' @param modname Character string model name; used to look up the model in
#'   \code{mod.res}.
#' @param species Character string species code; used in output filenames.
#' @param mod.res Named list of fitted \code{dsm} objects.
#' @param predgrid \code{sf} prediction grid with one row per
#'   (cell × season × platform) combination and all model predictors present.
#'   The \code{area} column is used as the prediction offset.
#' @return \code{predgrid} augmented with \code{NHat}, \code{Dens}, and
#'   \code{subset} columns, with rows for each platform level plus a
#'   \code{"Combined"} set.
#' @export
dsm.pred <-
  function(modname,
           species,
           mod.res,
           predgrid) {

    message(sprintf("Predicting %s for dsm model %s", species,
                    modname))

    # Extract model and create a subset column in predgrid from platform
    # Might have been just able to use platform.
    model <- mod.res[[modname]]

    # Reconcile predgrid platform levels against the ones the model was
    # actually fitted with.
    #
    # These are derived independently: segdata$platform levels come from the
    # data (as.factor() in create.dsm.data), whereas the predgrid copies come
    # from whatever create.seasonal.predgrid() was given - which defaults to the
    # *declared* global ddftype_to_platform. They agree only if every declared
    # platform survives into the fitted data, and create.dsm.data() now drops
    # ddftypes with no observations, so routinely they do not.
    #
    # Both directions are errors, for different reasons:
    #
    #   predgrid has a level the model lacks - predict() would fail deep inside
    #     mgcv with "factor has new levels". Loud, but obscure.
    #
    #   model has a level the predgrid lacks - nothing fails. The combined
    #     surface below would simply sum fewer components than the model has and
    #     under-report density, silently. This is the dangerous one, and it is
    #     exactly the failure that halved the predictions before ddftype_to_platform
    #     was given one level per ddftype.
    #
    # Fail here for both, where the cause can be named.
    if (is.null(predgrid$platform))
      stop("dsm.pred: predgrid has no 'platform' column. It is added by ",
           "create.seasonal.predgrid(replicate.platform = TRUE), which is the ",
           "default; assess.extrapolation() is the only caller that turns it off.")
    plat_levels <- levels(factor(predgrid$platform))
    mod_levels <- model$xlevels$platform
    if (!is.null(mod_levels)) {
      extra <- setdiff(plat_levels, mod_levels)
      if (length(extra) > 0)
        stop(sprintf(
          paste0("dsm.pred: predgrid has platform level(s) %s that model '%s' ",
                 "was not fitted with (model has %s). Pass the model's own ",
                 "levels to create.seasonal.predgrid(platform.levels = ...) ",
                 "rather than relying on the global ddftype_to_platform."),
          paste(sQuote(extra), collapse = ", "), modname,
          paste(sQuote(mod_levels), collapse = ", ")))

      missing <- setdiff(mod_levels, plat_levels)
      if (length(missing) > 0)
        stop(sprintf(
          paste0("dsm.pred: model '%s' was fitted with platform level(s) %s ",
                 "that the predgrid does not have (predgrid has %s). The ",
                 "combined surface sums across platform levels, so predicting ",
                 "without these would under-report density by roughly %d%% ",
                 "with no error. Pass the model's own levels to ",
                 "create.seasonal.predgrid(platform.levels = ...)."),
          modname, paste(sQuote(missing), collapse = ", "),
          paste(sQuote(plat_levels), collapse = ", "),
          round(100 * length(missing) / length(mod_levels))))
    }

    # One row per cell x season x platform, in that nesting. If this does not
    # hold, the head()/split() arithmetic building the combined surface below is
    # meaningless - it would silently mix cells across blocks.
    n.cells <- nrow(predgrid) / (length(season.names) * length(plat_levels))
    if (n.cells != round(n.cells))
      stop(sprintf(
        paste0("dsm.pred: predgrid has %d rows, which is not a whole number of ",
               "cells x %d seasons x %d platform levels. It is not the grid ",
               "create.seasonal.predgrid() builds."),
        nrow(predgrid), length(season.names), length(plat_levels)))
    stopifnot(all(table(predgrid$platform, predgrid$Season) == n.cells))

    ret <- predgrid %>%
      dplyr::mutate(subset = platform)

    # Predict from model for all seasons and values of platform, which
    # may or may not be a factor in the given model. So, in a no_factor model
    # we will end up with n identical copies of predictions: 1 for each value
    # of platform. In a factor model these predictions for each level of platform
    # are different.
    # discrete = FALSE forces exact prediction. A model fitted with
    # discrete = TRUE re-discretises whatever newdata it is handed, so a cell's
    # prediction depends on the composition of the grid it is predicted with -
    # change the grid extent, or predict in batches, and unchanged cells move.
    # predict.dsm() strips the "dsm" class and forwards ... to predict.bam(), so
    # this reaches the right place. It is a no-op for a model that was fitted
    # without discrete (e.g. via dsm.options$refit.without.discrete), which is
    # the other, broader way to get the same guarantee.
    ret <- ret %>%
      dplyr::mutate(NHat = predict(model, newdata=ret, off.set=ret$area,
                                   discrete = FALSE),
                    Dens = NHat/area)

    # Create summed (ie Combined Nhat) across all levels of platform.
    ret <- ret %>%
      # Peel off one copy of all seasons predgrid template. ie nrow(ret) divided
      # by the number of platforms.
      head(nrow(.) / length(unique(.$platform))) %>%
      # Sum NHats across platforms splitting the dataframe into a list with one
      # element per platform, extracting the NHat columns and summing with reduce
      dplyr::mutate(subset = "Combined",
                    NHat =  split(ret, ~ platform) %>%
                      purrr::map(~ .x$NHat) %>%
                      purrr::reduce(`+`),
                    Dens = NHat/area,
                    platform = NA) %>%
      rbind(ret)   # tack on original rows for each platform level



    # if (do.plots) {
    #   # Plot combined result
    #   combined_plot <- ggplot2::ggplot() +
    #     ggplot2::geom_sf(
    #       data = filter(ret, subset == "Combined"),
    #       mapping = ggplot2::aes(colour = Dens, fill = Dens)
    #     ) +
    #     ggplot2::labs(x = "", y = "", fill = "Dens") +
    #     ggtitle("Fly+Swim Combined") +
    #     theme_minimal() +
    #     scale_colour_viridis_c(option = "E") +
    #     scale_fill_viridis_c(option = "E")
    #   print(combined_plot)
    #
    #   # Plot individual fly and swim results
    #   ind_plot <- ggplot2::ggplot() +
    #     ggplot2::geom_sf(
    #       data = filter(ret, subset != "Combined"),
    #       mapping = ggplot2::aes(colour = Dens, fill = Dens)
    #     ) +
    #     ggplot2::labs(x = "", y = "", fill = "Dens") +
    #     ggplot2::facet_wrap(vars(FlySwim)) +
    #     theme_minimal() +
    #     scale_colour_viridis_c(option = "E") +
    #     scale_fill_viridis_c(option = "E")
    #   print(ind_plot)
    # }

    # Save as a shapefile
    ret %>%
      dplyr::filter(subset == "Combined") %>%
      sf::st_write(
        dsn = ShapeDir,
        layer = paste(species, modname, "predictions", sep = "_"),
        driver = "ESRI Shapefile",
        delete_layer = TRUE
      )



    # create 4-lyr seasonal raster of combined values for plotting and
    # saving
    message("Creating and saving seasonal rasters.")
    pp_raster <- season.names %>%
      purrr::map(make.season.raster,
                 obj = dplyr::filter(ret, subset == "Combined"),
                 variable = "Dens") %>%
      terra::rast()
    names(pp_raster) <- season.names

    # produce 4 panel plot. Not necessary since this is done at end of Generic_3_prediction.rmd
    # plot(pp_raster, main = paste(season.names, paste(species, modname, sep = "_")))

    # Save as a raster
    terra::writeRaster(pp_raster,
                       file.path(
                         predDir,
                         sprintf(
                           "%s.%s.%s.%d_sqkm.tif",
                           species,
                           season.names,
                           modname,
                           predgridCellArea
                         )
                       ),
                       datatype = "FLT4S",
                       overwrite = T)
    ret
  }



#' Render the DDF fitting RMarkdown report for one species
#'
#' @param species Character string species code.
#' @param do.final If \code{TRUE}, render the final DDF report; otherwise
#'   render the candidate DDF report.
#' @param rerun If \code{TRUE}, re-fit models even when saved results exist.
#' @param parallel If \code{TRUE}, use parallel processing.
#' @param nCores Number of cores for parallel processing.
#' @param do.eda If \code{TRUE}, produce EDA plots during rendering.
#' @param cleanFolder If \code{TRUE}, clean the results folder before fitting.
#' @return \code{invisible(NULL)}, called for its side-effect (HTML rendered).
#' @export
do.det.fcn.render <- function(species,
                              do.final,
                              rerun,
                              parallel,
                              nCores,
                              do.eda,
                              cleanFolder) {


  if (do.final)
    suffix <- "01_final_ddf.html"
  else
    suffix <- "01_candidate_ddf.html"

  out.file <- file.path(ResultsDir, species, paste(species, suffix, sep = "_"))
  message(sprintf("Rendering generic ddf fitting for %s to %s", species, out.file))

  # Output message in render pane if called from a knitted document.
  if (isTRUE(getOption('knitr.in.progress'))) {
    cat(sprintf("\nRendering generic ddf fitting for %s to %s", species, out.file),
        file = stderr())
  }

  # Make sure output dir exists
  if (!dir.exists(dirname(out.file)))
    dir.create(dirname(out.file), recursive = TRUE)

  rmarkdown::render(
    file.path(here::here("R"), "Generic_1_ddf_fitting.rmd"),
    params = list(
      species = species,
      do.final = do.final,
      rerun = rerun,
      parallel = parallel,
      nCores = nCores,
      do.eda = do.eda,
      cleanFolder = cleanFolder
    ),
    output_file = out.file
  )
}

#' Produce a summary table of candidate detection functions for one DDF spec
#'
#' @param df.spec.name Character string DDF spec name (used to locate the
#'   results folder and label the output).
#' @param species Character string species code.
#' @return Data frame of candidate model statistics sorted by AIC.
#' @export
candidate.detfcn.summary <- function(df.spec.name, species) {

  folder <- file.path(ResultsDir, species, paste("DF Summaries", df.spec.name, sep = "_"))
  do.ds.det.fcn.checks(folder, species = species, dsetname = df.spec.name)

  get.ds.res(folder) %>%
    dplyr::mutate(dataset = df.spec.name, species = species) %>%
    dplyr::relocate(species, dataset, Model, Key, Formula, DetProb, AIC, deltaAIC, )
}


#' Collect detection function results from summary text files
#'
#' Reads all \code{AIC*.txt} files in \code{folder} with
#' \code{\link{get.stats.file}} and returns a data frame sorted by AIC with a
#' \code{deltaAIC} column.
#'
#' @param folder Directory path containing \code{AIC_*.txt} summary files.
#' @return Data frame with one row per fitted model and columns \code{Model},
#'   \code{Key}, \code{Formula}, \code{AIC}, \code{deltaAIC}, and several
#'   GOF statistics.
#' @export
get.ds.res <- function(folder) {
  fl <- list.files(path = folder, pattern = "AIC.*\\.txt", full.names = T)

  # If all ddf failed there will be no results so we can't use "order"
  # to sort by aic, so just return the empty dataframe
  if (length(fl) == 0) {
    data.frame(
      Model = NA,
      Key  = NA,
      Formula = NA,
      AIC  = NA,
      ChiScores  = NA,
      ChisquareP  = NA,
      DetProb = NA,
      DetSE = NA,
      DetCV = NA,
      NObs = NA,
      NCov = NA,
      NCovSE = NA,
      NCovCV = NA,
      deltaAIC = NA
    )
  } else {
    # XXXX replace with dplyr
    res <- plyr::ldply(fl, get.stats.file)
    res$deltaAIC <- res$AIC - min(res$AIC, na.rm = T)
    res[order(res$AIC), ]
  }
}

#' Load and plot each fitted detection function in a results folder
#'
#' Iterates over all \code{AIC*.RData} files in \code{folder} and calls
#' \code{\link{check.det.fcn}} on each loaded model.
#'
#' @param folder Directory path containing \code{AIC_*.RData} model files.
#' @param species Character string species code; passed to
#'   \code{\link{check.det.fcn}}.
#' @param dsetname Character string dataset name; used in plot titles.
#' @return \code{invisible(NULL)}, called for its side-effects (diagnostic
#'   plots).
#' @export
do.ds.det.fcn.checks <- function(folder, species, dsetname) {
  fl <- list.files(path = folder, pattern = "AIC.*\\.RData", full.names = T)

  purrr::walk(fl, function(filename) {
    message(sprintf("Checking detection function: %s", filename))
    load(filename)
    check.det.fcn(model,
                  species = species,
                  stringr::str_replace(basename(folder), stringr::fixed("DF Summaries_"), ""))
  })
}


#' Parse a detection-function summary text file into a one-row data frame
#'
#' Reads the plain-text summary file written by \code{\link{do.ds}} and
#' extracts model name, key function, formula, AIC, chi-square GOF statistics,
#' detection probability, and abundance estimates.
#'
#' @param path Full path to a \code{.txt} summary file.
#' @return One-row data frame with columns \code{Model}, \code{Key},
#'   \code{Formula}, \code{AIC}, \code{ChiScores}, \code{ChisquareP},
#'   \code{DetProb}, \code{DetSE}, \code{DetCV}, \code{NObs}, \code{NCov},
#'   \code{NCovSE}, \code{NCovCV}.
#' @export
get.stats.file <- function(path) {
  message("Getting stats file ", path)

  lines <- readLines(path)

  getStat <- function(lines, string){
    if (length(line <- grep(string, lines)) == 0) {
      warning(path, ": No line containing ", string, "!", immediate. = T)
      return("Not available")
    }
    strsplit(lines[line], string)[[1]][2]
  }

  res <- data.frame(Model = getStat(lines, "^Model name: "),
                    Key  = getStat(lines, "Key: "),
                    Formula  = getStat(lines, "Formula: "),
                    AIC  = as.numeric(getStat(lines, "AIC: ")),
                    ChiScores  = getStat(lines, "Chi_scores: "),
                    ChisquareP  = as.numeric(destring(getStat(lines, "Chisquare p: "))),
                    DetProb = as.numeric(destring(getStat(lines, "Det prob: "))),
                    DetSE = as.numeric(destring(getStat(lines, "^SE\\(p\\): "))),
                    DetCV = as.numeric(destring(getStat(lines, "CV\\(p\\): "))),
                    NObs = as.numeric(destring(getStat(lines, "Number of observations :"))),
                    NCov = as.numeric(destring(strsplit(getStat(lines, "N in covered region "), " ")[[1]][1])),
                    NCovSE = as.numeric(destring(strsplit(getStat(lines, "N in covered region "), " ")[[1]][2])),
                    NCovCV = as.numeric(destring(strsplit(getStat(lines, "N in covered region "), " ")[[1]][3])),
                    stringsAsFactors = F)
  res
}


#' Add ddftype suffix and ddfobj.orig to fitted.distdata in a DDF spec
#'
#' Appends the \code{ddftype} value as a suffix to \code{Sample.Label} in
#' \code{df.mod.spec$fitted.distdata} and adds an integer \code{ddfobj.orig}
#' column derived from the \code{ddftype} factor.  Used by
#' \code{\link{create.dsm.data}} to ensure unique \code{Sample.Label} values
#' across the combined \code{segdata} and to record each observation's original
#' detection-function index before sequential renumbering.
#'
#' @param df.mod.spec Single element of a DDF model spec list, as created by
#'   \code{\link{do.det.fcn.spec}}.
#' @return \code{df.mod.spec} with \code{fitted.distdata} augmented by a
#'   \code{ddftype}-suffixed \code{Sample.Label} and an integer
#'   \code{ddfobj.orig} column.
#' @export
augment.distdata <- function(df.mod.spec) {
  df.mod.spec$fitted.distdata <- df.mod.spec$fitted.distdata %>%
    dplyr::mutate(Sample.Label = paste(Sample.Label, df.mod.spec$ddftype, sep = "_"),
                  ddfobj.orig = as.integer(df.mod.spec$ddftype))
  df.mod.spec
}


#' Create a platform-specific segdata copy for one DDF spec
#'
#' Filters \code{init_segdata} to the survey type (ship or aerial) implied by
#' \code{df.mod.spec$ddftype}, appends the ddftype as a \code{Sample.Label}
#' suffix, and adds \code{ddftype_orig} and \code{platform} columns.  The
#' \code{platform} column is derived via the \code{ddftype_to_platform} lookup
#' vector which must exist in the calling environment.
#'
#' @param df.mod.spec Single element of a DDF model spec list; must have a
#'   scalar \code{ddftype} field.
#' @param init_segdata Segment data frame containing a \code{SurveyType}
#'   column with values \code{"Ship"} and \code{"Aerial"}.
#' @return Filtered and augmented segment data frame for the given ddftype.
#' @export
create.segdata.copy <- function(df.mod.spec, init_segdata) {

  # This all only works if df.mod.spec$ddftype is a single value (which it
  # should be), but I'm paranoid.
  stopifnot(length(unique(df.mod.spec$ddftype)) == 1)

  segdata <- switch(
    substr(df.mod.spec$ddftype, 1, 1),
    S = dplyr::filter(init_segdata, SurveyType == "Ship"),
    A = dplyr::filter(init_segdata, SurveyType == "Aerial")
  ) %>%
    dplyr::mutate(
      Sample.Label = paste(Sample.Label, df.mod.spec$ddftype, sep = "_"),
      ddftype_orig = df.mod.spec$ddftype,
      # Lookup the platform type for this ddftype. May or may not have a 1:1
      # correspondence and ddftype_to_platform gives flexibility.
      platform = ddftype_to_platform[as.character(df.mod.spec$ddftype)]
    )
  segdata
}

#' Assemble distdata, segdata, and DDF list for DSM fitting
#'
#' Builds the three data structures required by \code{\link[dsm]{dsm}} for a
#' given species: (1) combined observation data (\code{distdata}) with
#' sequentially renumbered \code{ddfobj} values; (2) replicated segment data
#' (\code{segdata}) — one copy per \code{ddftype} — augmented with
#' \code{Season}, \code{segment.area}, \code{platform}, and \code{ddfobj}; and
#' (3) the list of fitted detection-function objects (\code{ddfs}), with any
#' ddftype that had no observations excluded.
#'
#' Expects project globals \code{seasons}, \code{ddftype_levels},
#' \code{def.ddf.list}, and \code{ddftype_to_platform} in the calling
#' environment.
#'
#' @param species Character string species code; used to select season
#'   definitions from \code{seasons[[species]]}.
#' @param df.mod.specs List of DDF model spec objects, as returned by
#'   \code{\link{do.det.fcn.specs}}.
#' @param init_segdata Segment data frame with a \code{SurveyType} column
#'   (\code{"Ship"} or \code{"Aerial"}), as produced by
#'   \code{\link{create.segdata}}.
#' @return Named list with elements \code{distdata}, \code{segdata}, and
#'   \code{ddfs}.
#' @export
create.dsm.data <- function(species, df.mod.specs, init_segdata) {

  # The ddfobj bookkeeping below turns each spec's ddftype factor into an
  # integer (augment.distdata) and indexes conv/def.ddf.list with it, so the
  # specs' factor levels must be exactly the current ddftype_levels.
  #
  # They diverge whenever ddftype_levels changes - most easily by pruning to
  # the survey types present (prune.ddf.globals) - while the fitted ddfs saved
  # in <species>_dfModList.rda still carry the old levels. The symptom is an
  # inscrutable "Something went wrong assigning ddfobjs" further down, so check
  # up front and say what to do about it.
  spec_levels <- unique(unlist(lapply(df.mod.specs, function(s) levels(s$ddftype))))
  if (!is.null(spec_levels) && !identical(spec_levels, ddftype_levels))
    stop(sprintf(
      paste0("create.dsm.data: the fitted ddf specs use ddftype levels (%s) ",
             "that differ from the current ddftype_levels (%s). The saved ddfs ",
             "are stale - re-run 01.03_Do_final_ddf_fitting.Rmd for '%s' so the ",
             "ddfs are refitted against the current ddftypes."),
      paste(spec_levels, collapse = ", "),
      paste(ddftype_levels, collapse = ", "),
      species))

  #### Create distdata ---------------------

  # Create the observation data that will be passed to dsm(). Note that some of
  # these distdata components specified by df.mod.specs may have no obs and will
  # thus not be included
  distdata <- df.mod.specs %>%
    purrr::map(augment.distdata) %>%
    # Extract list of fitted.distdatas
    purrr::map("fitted.distdata") %>%
    # Remove datasets with no obs - note use of base R Filter()
    Filter(function(x) nrow(x) > 0, .) %>%
    purrr::map_dfr( ~ .) %>%  # Convert list to single dataframe
    dplyr::mutate(
      # Need to renumber the ddfobj values sequentially from 1 to the number
      # of ddfs, whereas they may currently have gaps in the numeric sequence.
      # For example, if there were no obs in aerial_F_D (ddfobj.orig == 4) or
      # aerial_F_N (ddfobj.orig == 8) then the
      # ddfobj.orig would be numbered 1, 2, 3, 5, 6, 7 but should be renumbered
      # as 1:6.
      ddfobj = dplyr::dense_rank(ddfobj.orig)
    )

  # Make sure we still have some obs in distdata and all rows  are assigned a
  # ddfobj.
  stopifnot(nrow(distdata) > 0)
  stopifnot(sum(table(distdata$ddfobj)) == nrow(distdata))

  #### Create segdata --------------------

  # Keep a segdata copy only for ddftypes that actually have observations.
  #
  # This used to emit a copy for every ddftype regardless, pointing the ones
  # with no obs at min(newddfs) so dsm() would not trip over a ddfobj with no
  # observations. Their abundance.est came out 0 (no Sample.Label matches), so
  # they contributed a full study area of structural zeros. The long-standing
  # "XXXX Does having all these 0's bias the gam????" note asked what that did.
  #
  # It halved the density, and it did so through ddftype_to_platform. While D
  # and N shared a platform level, an all-zero copy sat beside a populated one
  # with nothing in the model to tell them apart, so the GAM fitted their mean:
  # on NL_EXPL_DRL_RA ATPU, fitted SWD = fitted SWN = 1109 against observed 2147
  # and 0. Now that each ddftype has its own platform level (see
  # ddftype_to_platform in analysis_settings.R) an empty ddftype would instead
  # become an all-zero LEVEL, whose intercept a log link cannot fit sensibly.
  #
  # Either way the copy carries no information, so drop it. The same note
  # worried that this makes the number of copies vary per species and that
  # "downstream code (ie. in prediction step) would know how to parse out the
  # rows properly" - which is why create.seasonal.predgrid() now takes its
  # platform levels from the fitted model rather than from the global.
  #
  # Note this drops only ddftypes with NO observations at all. A ddftype with
  # few observations is kept and warned about below: the _N types are
  # strip/dummy ddfs (generic.strip.ddf.spec) so nothing is estimated from those
  # observations, and the exposure is the DSM intercept, not the ddf.
  origddfs <- sort(unique(distdata$ddfobj.orig))
  kept <- purrr::map_int(df.mod.specs, ~ as.integer(.x$ddftype)) %in% origddfs
  if (any(!kept))
    message(sprintf(
      "create.dsm.data: dropping segdata %s for ddftype %s - no observations.",
      ifelse(sum(!kept) == 1, "copy", "copies"),
      paste(purrr::map_chr(df.mod.specs[!kept], ~ as.character(.x$ddftype)),
            collapse = ", ")))
  df.mod.specs <- df.mod.specs[kept]

  # Every surviving ddftype has observations, so conv just renumbers the
  # original ddfobj values sequentially to match the renumbering done in
  # distdata above. conv is indexed by the original ddftype integer.
  nddfs <- length(def.ddf.list)
  newddfs <- sort(unique(distdata$ddfobj))
  conv <- rep(NA_integer_, nddfs)
  conv[origddfs] <- newddfs

  # Cylce through the list of ddf model specs adding a copy of segdata for each.
  segdata <- df.mod.specs %>%
    purrr::map(create.segdata.copy, init_segdata = init_segdata) %>%
    purrr::map_dfr(~ .) %>% # Convert to one large dataframe
    assign.season(seasons[[species]], datefield = "Date") %>%
    dplyr::select(
      SurveyType,
      Sample.Label,
      CruiseID,
      Program,
      ObserverName,
      DistMeth,
      CalcDurMin,
      ObsHeight,
      TransectID,
      TransectSides,
      Dataset,
      year,
      yday,
      Date,
      StartTime,
      EndTime,
      LatStart,
      LatEnd,
      LongStart,
      LongEnd,
      x,
      y,
      depth,
      depth.g,
      sst,
      sst.g,
      Effort,
      segment.area,
      platform,
      ddftype_orig,
      Season,
      dplyr::ends_with(".sc")
    ) %>%
    # Now that platform has been created, turn it into a factor, and assign
    # ddfobj based on conversion vector created above
    dplyr::mutate(platform = as.factor(platform),
                  ddfobj = conv[ddftype_orig]) %>%
    # dsm() uses observations in segdata and associated detection function to
    # calculate the gam response internally. Calculate it here so I can check it's
    # done right and I understand what's going on internally in dsm().
    augment.segdata(distdata)

  ##### Make sure assignment of ddfobj went properly
  #
  # Every ddftype that reaches here has observations, so each one contributes
  # exactly one copy of its survey type's segments and every segment gets a real
  # ddfobj. That is a stronger and much simpler statement than the arithmetic
  # this replaced, which existed only to account for the orphan copies that are
  # no longer created.
  stopifnot(!anyNA(segdata$ddfobj))

  # Get survey type of each ddf. This requires that the order of ddftype_levels
  # and def.ddf.list have the same order.
  survey_type_index <- substr(ddftype_levels, 1, 1)

  # get aerial and ship initial segdata sizes and name them "A" and "S".
  #
  # NB: table() only returns entries for survey types actually present in the
  # data, so a study area covered by only one survey type (eg. one with no
  # aerial coverage at all) leaves sizes["A"] as NA. Build the vector over every
  # survey type implied by ddftype_levels so absent ones are 0, not NA.
  expected_types <- unique(survey_type_index)
  sizes <- table(init_segdata$SurveyType) %>%
    setNames(names(.) %>% substr(1, 1))
  sizes <- setNames(as.vector(sizes[expected_types]), expected_types)
  sizes[is.na(sizes)] <- 0

  expected_rows <- sizes[survey_type_index[origddfs]]
  actual_rows <- as.vector(table(segdata$ddfobj)[as.character(newddfs)])
  if (!isTRUE(all.equal(as.vector(expected_rows), actual_rows)))
    stop(sprintf(
      paste0("create.dsm.data: segdata copy sizes are wrong for '%s'. ddftype ",
             "%s should contribute %s segments respectively but contributed %s."),
      species, paste(ddftype_levels[origddfs], collapse = ", "),
      paste(expected_rows, collapse = ", "),
      paste(actual_rows, collapse = ", ")))

  # Warn about a platform level resting on very few observations. It is not an
  # error and nothing is dropped: the _N ddftypes are strip/dummy ddfs so their
  # observation count does not threaten a detection function. What it threatens
  # is the DSM, where each ddftype now carries its own platform intercept. On
  # NL_EXPL_DRL_RA, SWN is empty for nine of eleven species and rests on 1 (BLKI)
  # and 3 (NOFU) observations for the other two, because lkpDistMeth gives water
  # birds perpendicular distances under almost every DistMeth in use.
  sparse <- table(distdata$ddfobj.orig)
  sparse <- sparse[sparse < SPARSE_PLATFORM_LEVEL_OBS]
  if (length(sparse) > 0)
    warning(sprintf(
      paste0("create.dsm.data: '%s' has platform level(s) resting on very few ",
             "observations: %s. Kept, but check that the level's fitted ",
             "intercept and its contribution to the combined prediction are ",
             "sensible."),
      species,
      paste(sprintf("%s (n=%d)",
                    ddftype_to_platform[as.integer(names(sparse))],
                    as.vector(sparse)), collapse = ", ")),
      immediate. = TRUE)

  ### Create ddfs ----------------------------------------

  # df.mod.specs is already filtered to the ddftypes with observations, in
  # ddftype_levels order, and distdata$ddfobj was dense-ranked over the same
  # ascending order - so element i here is ddfobj i. Indexing by
  # unique(distdata$ddfobj.orig) would be wrong now, since those are positions in
  # the UNFILTERED list.
  stopifnot(identical(unname(purrr::map_int(df.mod.specs,
                                            ~ as.integer(.x$ddftype))),
                      as.integer(origddfs)))
  ddfs <- purrr::map(df.mod.specs, "fitted.model")

  ### TODO: need to deal with segments where WhatCount was something funky
  ### like only counting gannets. In that case, if we are doing a different species
  ### then these segments will be zeros but should not be included at all. For
  ### now I think this NOGA exception is the only case.

  list(distdata = distdata, segdata = segdata, ddfs = ddfs)
}


#' List detection function models with no results file in a folder
#'
#' Calls \code{\link{do.ds}} with \code{runModels = FALSE} to enumerate all
#' candidate model combinations, then returns the subset that has no
#' corresponding results file in \code{folder}.  Useful for diagnosing stalled
#' or hung model runs.
#'
#' @param folder Directory path to search for existing results files.
#' @param covars Character vector of covariate names for model enumeration;
#'   defaults to the project global \code{dfCovars}.
#' @return Data frame of unfinished models (same structure as the model table
#'   produced by \code{\link{do.ds}}).
#' @export
list.unfinished.ddfs <- function(folder, covars = dfCovars){
  do.ds(
    data = NULL,
    folder = folder,
    runModels = FALSE,
    covars = covars
  )
}

#' Create a sentinel file marking a detection function model as failed
#'
#' Writes an empty \code{failed_model_<label>.txt} file into \code{folder} so
#' that downstream code treats this model as having failed without waiting for
#' it to finish.  Use when a model run hangs and must be marked failed by hand.
#'
#' @param mod Model descriptor with a \code{label} field used to construct the
#'   filename.
#' @param folder Directory path in which to create the sentinel file.
#' @return \code{TRUE} (invisibly) on success, as returned by
#'   \code{\link[base]{file.create}}.
#' @export
create.failed.ddf.file <- function(mod, folder) {
  file.create(file.path(folder, paste0("failed_model_", mod$label, ".txt")))
}


#' Integer disaggregation factor for a covariate's gradient grid
#'
#' Chooses how many fine cells to divide each analysis cell into, so the fine
#' cell is as close as possible to - but never finer than - the source raster's
#' own cell size. Returns 1 when the analysis grid is already at or below
#' source resolution, in which case no disaggregation happens and the gradient
#' is computed on the analysis grid itself.
#'
#' The source spacing is the mean of the x and y cell size in km. For a lon/lat
#' raster those differ - 0.01 deg is 1.112 km N-S everywhere but 0.758 km E-W
#' at 47N - so the mean is a compromise between discarding real E-W detail and
#' inventing N-S detail that is not there.
#'
#' @param src SpatRaster of the covariate at its source resolution.
#' @param coarse.km numeric(1). Analysis grid cell length in km.
#' @param mid.lat numeric(1). Latitude at which to convert degrees of longitude
#'   to km; ignored when `src` is already projected.
#' @return integer(1), at least 1.
#' @examples
#' # MUR sst, 0.01 deg, at 47N: 10 km grid -> 10, 2 km grid -> 2 (1 km cells)
#' # ETOPO depth, 1 arc-min:    10 km grid -> 6,  2 km grid -> 1 (no disagg)
#' @export
gradient.disagg.factor <- function(src, coarse.km, mid.lat) {
  checkmate::expect_class(src, "SpatRaster")
  checkmate::expect_number(coarse.km, lower = 0, finite = TRUE)
  checkmate::expect_number(mid.lat, lower = -90, upper = 90)

  r <- terra::res(src)
  src.km <- if (terra::is.lonlat(src)) {
    mean(c(r[2] * 111.195, r[1] * 111.195 * cos(mid.lat * pi / 180)))
  } else {
    mean(r) / 1000
  }
  max(1L, as.integer(floor(coarse.km / src.km)))
}

#' Covariate and gradient on the analysis grid, differentiated at fine scale
#'
#' Projects `src` onto a fine grid nested inside `coarse`, optionally fills NA
#' holes so the focal operator does not erode the coast, computes the Belkin &
#' O'Reilly (2009) gradient there in units per km, and reduces both the
#' covariate and the gradient back onto `coarse` with `GRADIENT_AGG_FUN`.
#'
#' @param src SpatRaster at source resolution, any CRS.
#' @param coarse SpatRaster template - the buffered analysis grid (gradrast).
#' @param coarse.km numeric(1). Cell length of `coarse` in km.
#' @param mid.lat numeric(1). Study area centre latitude.
#' @param clamp.lower numeric(1) or NULL. Values below this are clamped up to
#'   it before differentiating. Used for depth, where land is clamped to the
#'   waterline rather than removed.
#' @param fill logical(1). Whether to NA-fill before differentiating. FALSE for
#'   covariates with no gaps.
#' @return list with `covar` and `grad`, both SpatRasters on `coarse`.
#'
#' @section Project globals:
#' Reads \code{GRADIENT_NA_FILL_PASSES}, \code{SOBEL_KERNEL_GAIN} and
#' \code{GRADIENT_AGG_FUN}, all defined in \code{analysis_settings.R} with
#' the measurements that set them.  They are globals rather than arguments to
#' match the rest of this package, and because they are properties of the
#' analysis rather than of a single call.
#' @export
fine.gradient <- function(src, coarse, coarse.km, mid.lat,
                          clamp.lower = NULL, fill = TRUE) {
  checkmate::expect_class(src, "SpatRaster")
  checkmate::expect_class(coarse, "SpatRaster")
  checkmate::expect_number(coarse.km, lower = 0, finite = TRUE)
  checkmate::expect_flag(fill)

  n <- gradient.disagg.factor(src, coarse.km, mid.lat)
  fine <- if (n > 1) terra::disagg(coarse, fact = n) else coarse
  message(sprintf("  fine grid: factor %d, cell %.3f km", n, coarse.km / n))

  x <- terra::project(src, fine, method = "average", threads = TRUE)
  names(x) <- names(src)

  g.in <- x
  if (!is.null(clamp.lower))
    g.in <- terra::clamp(g.in, lower = clamp.lower, values = TRUE)
  if (fill)
    for (i in seq_len(GRADIENT_NA_FILL_PASSES))
      g.in <- terra::focal(g.in, w = 3, fun = mean, na.rm = TRUE,
                           na.policy = "only")

  g <- grec::getGradients(g.in, method = "BelkinOReilly2009") /
    (SOBEL_KERNEL_GAIN * (coarse.km / n))

  if (n > 1) {
    x <- terra::aggregate(x, fact = n, fun = GRADIENT_AGG_FUN, na.rm = TRUE)
    g <- terra::aggregate(g, fact = n, fun = GRADIENT_AGG_FUN, na.rm = TRUE)
  }

  # Confine the gradient to the covariate's own footprint. Without this the
  # NA fill can leave a gradient in an analysis cell that is entirely land, so
  # sst.g would be defined where sst is not - measured as exactly one cell on
  # Atl IMRP, harmless because the dsm.covars.ns filter would drop that segment
  # on sst anyway, but a covariate defined outside its own data is not
  # something to leave lying around.
  g <- terra::mask(g, x[[1]])

  names(x) <- names(src)
  list(covar = x, grad = g)
}


#' Refuse to delete a non-empty prediction version folder
#'
#' \code{03.70_Save_chosen_model_predictions.Rmd} and
#' \code{04.70_Save_chosen_model_prediction_variance.Rmd} both begin by
#' \code{unlink()}ing \code{predVersionDir} so that copied files get fresh
#' dates.  That is correct when you are correcting a version in place, and
#' destructive when you are not: those folders are shared outside the project
#' and are not in git, so a version deleted here is gone.
#'
#' Keeping the previous version means setting \code{predVersionDir} to a NEW
#' folder name in \code{analysis_settings.R} first.  Nothing enforced that -
#' it was a manual step documented in CLAUDE.md and Notes.docx and nowhere in
#' the code, which is to say it depended on remembering.  This function makes
#' the code stop instead.
#'
#' Call it immediately before the \code{unlink()}.  It is a no-op when the
#' folder does not exist or is empty, which is the normal case for a new
#' version.
#'
#' @param dir Character path to the version folder, normally
#'   \code{predVersionDir}.
#' @param allow.overwrite Logical.  \code{TRUE} permits the delete and warns
#'   instead of stopping.  Normally \code{ALLOW_PRED_VERSION_OVERWRITE} from
#'   \code{analysis_settings.R}, which is tracked, so an override left switched
#'   on shows up in \code{git diff}.
#' @return \code{invisible(TRUE)} if an existing version is about to be
#'   overwritten, \code{invisible(FALSE)} if there is nothing there.  Called
#'   for its side effect of stopping.
#' @examples
#' # guard.version.dir(predVersionDir, ALLOW_PRED_VERSION_OVERWRITE)
#' @export
guard.version.dir <- function(dir, allow.overwrite) {
  checkmate::expect_string(dir, min.chars = 1)
  checkmate::expect_flag(allow.overwrite)

  if (!dir.exists(dir))
    return(invisible(FALSE))

  n <- length(list.files(dir, all.files = TRUE, no.. = TRUE, recursive = TRUE))
  if (n == 0)
    return(invisible(FALSE))

  if (allow.overwrite) {
    warning(sprintf(
      paste0("guard.version.dir: overwriting %d file(s) in an existing version ",
             "folder because ALLOW_PRED_VERSION_OVERWRITE is TRUE:\n  %s"),
      n, dir), immediate. = TRUE)
    return(invisible(TRUE))
  }

  stop(sprintf(
    paste0(
      "predVersionDir already holds %d file(s) and this step deletes it:\n",
      "  %s\n\n",
      "These folders are shared outside the project and are not in git, so ",
      "deleting one loses it.\n\n",
      "To KEEP that version: set predVersionDir in analysis_settings.R to a new ",
      "folder name (e.g. bump Version_N and the date), then re-run.\n",
      "To REPLACE it deliberately: set ALLOW_PRED_VERSION_OVERWRITE <- TRUE in ",
      "analysis_settings.R."),
    n, dir), call. = FALSE)
}

#' The highest density this species was ever actually seen at, by season
#'
#' A prediction-surface sanity check needs a ceiling, and a flat one cannot
#' serve every species. Measured across `NL_EXPL_DRL_RA`, the largest
#' detection-corrected density on any segment ranges from **0 to 7,183
#' birds/km2** depending on species and season -- a factor of seven thousand. A
#' single constant is therefore far too permissive for a razorbill in spring
#' (whose busiest segment ever held 16 birds/km2) and barely above the
#' legitimate range for a fulmar in spring (7,183). This returns the
#' species-and-season-specific reference instead.
#'
#' What it measures is deliberately the observed **maximum**, not a quantile:
#' the statement it supports is "the model predicts a seasonal-mean density in
#' one cell higher than the single busiest instantaneous encounter ever
#' recorded", which is interpretable without further calibration. `p999` and
#' `p99` come back too, because the maximum is one flock and worth seeing beside
#' a less brittle number.
#'
#' Predictions should ordinarily sit **well below** this. A prediction cell is a
#' seasonal mean over `predgrid.years`; an observed segment density is one
#' encounter on one day. Measured over 88 NL_EXPL_DRL_RA surfaces, the median
#' ratio of largest predicted cell to observed maximum is **0.075** -- the
#' typical surface peaks thirteen times below anything ever seen.
#'
#' @param segdata The `segdata` saved by the DSM step, with `estAbund`,
#'   `segment.area`, `Season` and `Sample.Label` columns. An `sf` object is
#'   fine; the geometry is dropped.
#' @return A tibble with one row per season: `season`, `n_seg` (physical
#'   segments, not rows), `n_pos`, `obs_max`, `obs_p999`, `obs_p99`,
#'   `obs_mean`. A season the species was never recorded in gets `obs_max` 0,
#'   which callers must treat as "no reference", not as a ceiling of zero.
#' @examples
#' \dontrun{
#' load(here(RDataDir, "ATPU_distdata_segdata_ddfs.Rda"))
#' observed.density.ceiling(segdata)
#' }
#' @export
observed.density.ceiling <- function(segdata) {
  checkmate::expect_data_frame(segdata, min.rows = 1)
  for (col in c("estAbund", "segment.area", "Season", "Sample.Label"))
    if (is.null(segdata[[col]]))
      stop("observed.density.ceiling: segdata has no '", col, "' column")

  sd <- segdata
  if (inherits(sd, "sf")) sd <- sf::st_drop_geometry(sd)

  # segdata holds one row per segment x ddftype, so a raw row density is only
  # part of the segment's birds. The surfaces this is compared against are the
  # summed "Combined" ones, so sum the ddftype rows back to the segment first.
  # Sample.Label carries the ddftype suffix and is unique per row - see
  # segment.id.from.label().
  sd %>%
    dplyr::mutate(
      .seg = segment.id.from.label(as.character(.data$Sample.Label))) %>%
    dplyr::group_by(.data$Season, .data$.seg) %>%
    dplyr::summarise(abund = sum(.data$estAbund, na.rm = TRUE),
                     area = dplyr::first(.data$segment.area),
                     .groups = "drop") %>%
    dplyr::filter(.data$area > 0) %>%
    dplyr::mutate(dens = .data$abund / .data$area) %>%
    dplyr::group_by(.data$Season) %>%
    dplyr::summarise(
      n_seg    = dplyr::n(),
      n_pos    = sum(.data$dens > 0),
      obs_max  = max(.data$dens),
      obs_p999 = unname(stats::quantile(.data$dens, 0.999)),
      obs_p99  = unname(stats::quantile(.data$dens, 0.99)),
      obs_mean = mean(.data$dens),
      .groups  = "drop") %>%
    dplyr::rename(season = "Season")
}

# The document scaffold and stylesheet for write.finalist.summary.page().
#
# Kept in its own function purely so the page builder above stays readable; it
# is one long string and has no logic in it. Everything the page draws comes
# from the tokens defined here, in all three theme states the page can be read
# in: an explicit light choice, an explicit dark choice, and the default where
# neither is stamped and only the OS preference distinguishes them.
.finalist.page.head <- function() {
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Two Families, One Question</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,400;0,600;1,400&amp;family=Source+Sans+3:ital,wght@0,400;0,600&amp;family=IBM+Plex+Mono:wght@400;500;600&amp;display=swap">
<style>
*,*::before,*::after{box-sizing:border-box}
body,h1,h2,h3,p,dl,dd,figure,table{margin:0}
:root{
  --ground:#F5F6F7; --surface:#FFFFFF; --sunk:#ECEFF1; --ink:#131A20;
  --muted:#5C6975; --faint:#8794A0; --rule:#DCE1E5; --rule-firm:#C3CBD2;
  --a:#A8681F; --a-wash:#F6EDE0; --b:#1A626B; --b-wash:#E3EEEF;
  --crit:#A02D20; --crit-wash:#F8E8E5; --ok:#1D6640;
  --shadow:0 1px 2px rgba(19,26,32,.06),0 6px 20px rgba(19,26,32,.05);
  --measure:65ch; --wide:62rem;
  --s1:.35rem; --s2:.7rem; --s3:1.1rem; --s4:1.75rem; --s5:2.75rem; --s6:4.25rem;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --ground:#0F1418; --surface:#161C22; --sunk:#1B232A; --ink:#E4EAEF;
  --muted:#97A5B1; --faint:#6F7E8B; --rule:#27313A; --rule-firm:#3A4650;
  --a:#DD9E52; --a-wash:#2B2216; --b:#55AEB7; --b-wash:#12262A;
  --crit:#E57764; --crit-wash:#2C1815; --ok:#5FBF8A;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 6px 20px rgba(0,0,0,.3);
}}
:root[data-theme="dark"]{
  --ground:#0F1418; --surface:#161C22; --sunk:#1B232A; --ink:#E4EAEF;
  --muted:#97A5B1; --faint:#6F7E8B; --rule:#27313A; --rule-firm:#3A4650;
  --a:#DD9E52; --a-wash:#2B2216; --b:#55AEB7; --b-wash:#12262A;
  --crit:#E57764; --crit-wash:#2C1815; --ok:#5FBF8A;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 6px 20px rgba(0,0,0,.3);
}
body{background:var(--ground);color:var(--ink);
  font-family:"Source Sans 3",ui-sans-serif,system-ui,sans-serif;
  font-size:17px;line-height:1.62;-webkit-font-smoothing:antialiased}
.page{max-width:var(--wide);margin:0 auto;padding:var(--s6) var(--s4);
  display:flex;flex-direction:column;gap:var(--s6)}
.prose{max-width:var(--measure)}
.prose>*+*{margin-top:var(--s3)}
h1,h2,h3{font-family:Spectral,Georgia,"Times New Roman",serif;font-weight:600;
  text-wrap:balance;line-height:1.18;letter-spacing:-.008em}
h1{font-size:clamp(2.1rem,5vw,3rem)}
h2{font-size:clamp(1.5rem,3vw,1.9rem)}
h3{font-size:1.16rem}
.eyebrow{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.74rem;
  font-weight:500;letter-spacing:.13em;text-transform:uppercase;color:var(--faint)}
code{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.87em;
  background:var(--sunk);padding:.1em .34em;border-radius:3px}
.masthead{display:flex;flex-direction:column;gap:var(--s3)}
.masthead .lede{font-family:Spectral,Georgia,serif;font-size:1.22rem;
  line-height:1.5;color:var(--muted);max-width:58ch}
.meta{display:flex;flex-wrap:wrap;gap:var(--s1) var(--s3);padding-top:var(--s3);
  border-top:1px solid var(--rule);font-family:"IBM Plex Mono",monospace;
  font-size:.78rem;color:var(--faint)}
.fam{font-family:"IBM Plex Mono",monospace;font-weight:600;font-size:.82em;
  padding:.12em .45em;border-radius:3px;white-space:nowrap}
.fam-a{color:var(--a);background:var(--a-wash)}
.fam-b{color:var(--b);background:var(--b-wash)}
.finding{background:var(--surface);border:1px solid var(--rule);
  border-left:4px solid var(--crit);border-radius:4px;box-shadow:var(--shadow);
  padding:var(--s4);display:flex;flex-direction:column;gap:var(--s3)}
.finding.is-ok{border-left-color:var(--ok)}
.finding .eyebrow{color:var(--crit)}
.finding.is-ok .eyebrow{color:var(--ok)}
.finding h2{font-size:1.42rem}
.finding p{max-width:62ch}
.decades{display:flex;flex-direction:column;gap:var(--s2);margin-top:var(--s2)}
.decade-row{display:grid;grid-template-columns:5.5rem 1fr;gap:var(--s3);
  align-items:center}
.decade-label{font-family:"IBM Plex Mono",monospace;font-size:.8rem;text-align:right}
.decade-track{position:relative;height:1.55rem;background:var(--sunk);
  border-radius:2px;overflow:hidden}
.decade-fill{position:absolute;inset:0 auto 0 0;border-radius:2px}
.decade-fill.a{background:var(--a)}
.decade-fill.b{background:var(--crit)}
.decade-val{position:absolute;top:50%;transform:translateY(-50%);
  font-family:"IBM Plex Mono",monospace;font-size:.78rem;font-weight:600;
  padding:0 .5rem;white-space:nowrap;color:var(--surface)}
.decade-scale{display:flex;justify-content:space-between;
  font-family:"IBM Plex Mono",monospace;font-size:.68rem;color:var(--faint);
  padding-left:calc(5.5rem + var(--s3))}
.caption{font-size:.84rem;color:var(--muted);max-width:62ch}
.stats{display:grid;gap:1px;grid-template-columns:repeat(auto-fit,minmax(11rem,1fr));
  background:var(--rule);border:1px solid var(--rule);border-radius:4px;overflow:hidden}
.stat{background:var(--surface);padding:var(--s3)}
.stat dt{font-family:"IBM Plex Mono",monospace;font-size:.69rem;letter-spacing:.08em;
  text-transform:uppercase;color:var(--faint);margin-bottom:var(--s1)}
.stat dd{font-family:"IBM Plex Mono",monospace;font-variant-numeric:tabular-nums;
  font-size:1.3rem;font-weight:600;line-height:1.15}
.stat dd .unit{font-size:.78rem;font-weight:400;color:var(--muted);display:block;
  margin-top:.2rem}
.stat.is-crit dd{color:var(--crit)}
.table-wrap{overflow-x:auto;border:1px solid var(--rule);border-radius:4px;
  background:var(--surface)}
table{width:100%;border-collapse:collapse;font-size:.88rem}
caption{caption-side:top;text-align:left;padding:var(--s3) var(--s3) var(--s2);
  color:var(--muted);font-size:.86rem}
th,td{padding:.5rem .8rem;text-align:left;white-space:nowrap}
thead th{font-family:"IBM Plex Mono",monospace;font-size:.7rem;letter-spacing:.06em;
  text-transform:uppercase;color:var(--faint);font-weight:500;
  border-bottom:1px solid var(--rule-firm)}
tbody tr+tr td{border-top:1px solid var(--rule)}
td.n,th.n{font-family:"IBM Plex Mono",monospace;font-variant-numeric:tabular-nums;
  text-align:right}
tbody tr.is-crit td{background:var(--crit-wash)}
.win-a{color:var(--a);font-weight:600}
.win-b{color:var(--b);font-weight:600}
.tally{display:grid;gap:var(--s3);grid-template-columns:repeat(auto-fit,minmax(15rem,1fr))}
.tally-item{background:var(--surface);border:1px solid var(--rule);border-radius:4px;
  padding:var(--s3)}
.tally-item h3{font-family:"Source Sans 3",sans-serif;font-size:.88rem;
  font-weight:600;margin-bottom:var(--s2)}
.bar{display:flex;height:1.5rem;border-radius:2px;overflow:hidden;background:var(--sunk)}
.bar span{display:grid;place-items:center;font-family:"IBM Plex Mono",monospace;
  font-size:.74rem;font-weight:600;color:var(--surface)}
.bar .b-a{background:var(--a)}
.bar .b-b{background:var(--b)}
.bar .b-tie{background:var(--faint)}
.tally-item .caption{margin-top:var(--s2);font-size:.78rem}
section{display:flex;flex-direction:column;gap:var(--s4)}
section>h2{max-width:var(--measure)}
footer{border-top:1px solid var(--rule);padding-top:var(--s4);color:var(--muted);
  font-size:.88rem;display:flex;flex-direction:column;gap:var(--s2)}
@media (max-width:34rem){
  .decade-row{grid-template-columns:4rem 1fr}
  .decade-scale{padding-left:calc(4rem + var(--s3))}
  .page{padding:var(--s5) var(--s3)}
}
</style>
</head>
<body>'
}

# Minimal HTML escaping for values that reach the summary page. Species and
# model names come from the pipeline rather than from a user, so this is
# belt-and-braces - but a page that gets shared should not be able to be
# derailed by an ampersand in a species group name.
.esc <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# Format a number for display: plain with thousands separators when it is of a
# human size, scientific when it is not. A summary page that prints
# 3951058971034.7 has told the reader nothing they can hold in their head.
.hnum <- function(x, sig = 3) {
  if (length(x) != 1 || is.na(x)) return("&mdash;")
  if (x != 0 && (abs(x) >= 1e5 || abs(x) < 1e-3)) {
    e <- floor(log10(abs(x)))
    m <- signif(x / 10^e, sig)
    sprintf("%s &times; 10<sup>%d</sup>", format(m, trim = TRUE), e)
  } else {
    format(signif(x, sig), big.mark = ",", scientific = FALSE, trim = TRUE)
  }
}

.pct <- function(x, dp = 1) if (is.na(x)) "&mdash;" else
  sprintf(paste0("%.", dp, "f%%"), x)

#' Write the shareable summary page for the finalist prediction comparison
#'
#' The per-species reports and the roll-up answer the question in detail; this
#' is the one page to hand someone who was not in the room. Everything on it is
#' computed from the arguments -- no finding, figure or species name is written
#' into the template -- so it stays true after a refit, and a SubProject where
#' nothing is wrong gets a page that says so rather than an empty version of
#' the page where something was.
#'
#' The output is a standalone HTML document that opens on its own. The body is
#' delimited by `<!-- ARTIFACT:BEGIN -->` and `<!-- ARTIFACT:END -->` sentinels,
#' because publishing it as an Artifact needs the content without the
#' `<!doctype>`/`<head>`/`<body>` scaffold, which that publisher supplies
#' itself. Extract between the sentinels rather than re-deriving the page.
#'
#' @param seasons The roll-up tibble: one row per species-season, as
#'   [compare.finalist.surfaces()] returns them, bound across species, with an
#'   `incumbent` column added.
#' @param cells The matching spike-cell tibble.
#' @param path Where to write the file.
#' @param subproject SubProject name, for the masthead.
#' @param cell.km2 Prediction cell area, for the masthead.
#' @param top_n,ref.mult,map.limit The settings the numbers were computed under,
#'   so the page states its own thresholds instead of leaving them implicit.
#' @param warn.pct The concentration level the shipping step warns at.
#' @return `path`, invisibly.
#' @examples
#' \dontrun{
#' write.finalist.summary.page(seasons, cells,
#'   file.path(ResultsDir, paste0(SubProject, "_finalist_comparison.html")),
#'   subproject = SubProject, cell.km2 = predgridCellArea)
#' }
#' @export
write.finalist.summary.page <- function(seasons, cells, path, subproject,
                                        cell.km2, top_n = 10, ref.mult = 100,
                                        map.limit = NULL, warn.pct = 50) {
  checkmate::expect_data_frame(seasons, min.rows = 1)
  checkmate::expect_data_frame(cells)
  checkmate::expect_string(path, min.chars = 1)
  checkmate::expect_string(subproject, min.chars = 1)
  checkmate::expect_number(cell.km2, lower = 0)

  fam1 <- seasons$model1[1]; fam2 <- seasons$model2[1]
  k1 <- if (!is.null(seasons$fam1)) seasons$fam1[1] else "1"
  k2 <- if (!is.null(seasons$fam2)) seasons$fam2[1] else "2"
  n_sp <- length(unique(seasons$species))
  n_se <- length(unique(seasons$season))

  # ---- the headline -------------------------------------------------------
  # Worst surface by how far it exceeds what was actually observed. That is the
  # measure with a defensible zero point; if no reference was available it falls
  # back to the flat limit, and if neither fired the page says nothing is wrong.
  imp <- seasons %>%
    dplyr::select(dplyr::any_of(c("species", "season", "max_over_obs1",
                                  "max_over_obs2"))) %>%
    tidyr::pivot_longer(dplyr::starts_with("max_over_obs"),
                        names_to = "which", values_to = "r") %>%
    dplyr::filter(!is.na(.data$r)) %>%
    dplyr::arrange(dplyr::desc(.data$r))

  worst <- if (nrow(imp) && imp$r[1] > ref.mult) {
    w <- seasons %>% dplyr::filter(.data$species == imp$species[1],
                                   .data$season == imp$season[1])
    list(row = w, fam = if (imp$which[1] == "max_over_obs1") k1 else k2,
         ratio = imp$r[1])
  } else NULL

  # ---- agreement ----------------------------------------------------------
  med <- function(v) stats::median(v, na.rm = TRUE)
  n_big <- sum(seasons$total_ratio > 1.5 | seasons$total_ratio < 1 / 1.5,
               na.rm = TRUE)

  spikier <- function(a, b) {
    ok <- !is.na(a) & !is.na(b)
    c(k1 = sum(a[ok] > b[ok]), tie = sum(a[ok] == b[ok]), k2 = sum(b[ok] > a[ok]))
  }
  t_gini <- spikier(seasons$gini1, seasons$gini2)
  t_mom  <- spikier(seasons$max_over_med1, seasons$max_over_med2)
  t_tot  <- c(k1 = sum(seasons$total_ratio > 1, na.rm = TRUE), tie = 0,
              k2 = sum(seasons$total_ratio < 1, na.rm = TRUE))

  bar <- function(t, label, note) {
    parts <- c(
      if (t[["k1"]] > 0) sprintf('<span class="b-a" style="flex:%d;">%s %d</span>',
                                 t[["k1"]], .esc(k1), t[["k1"]]),
      if (t[["tie"]] > 0) sprintf('<span class="b-tie" style="flex:%d;"></span>',
                                  t[["tie"]]),
      if (t[["k2"]] > 0) sprintf('<span class="b-b" style="flex:%d;">%s %d</span>',
                                 t[["k2"]], .esc(k2), t[["k2"]]))
    sprintf(paste0('<div class="tally-item"><h3>%s</h3>',
                   '<div class="bar" role="img" aria-label="%s %d, tie %d, %s %d">',
                   '%s</div><p class="caption">%s</p></div>'),
            .esc(label), .esc(k1), t[["k1"]], t[["tie"]], .esc(k2), t[["k2"]],
            paste(parts, collapse = ""), .esc(note))
  }

  # ---- tables -------------------------------------------------------------
  row.html <- function(cells.v, crit = FALSE)
    sprintf('<tr%s>%s</tr>', if (crit) ' class="is-crit"' else "",
            paste(cells.v, collapse = ""))
  td <- function(x, n = FALSE) sprintf('<td%s>%s</td>',
                                       if (n) ' class="n"' else "", x)

  flagged <- seasons %>%
    dplyr::mutate(
      .td = !is.na(.data$total_ratio) &
        (.data$total_ratio > 1.5 | .data$total_ratio < 1 / 1.5),
      .cd = xor(.data$pct_top_n1 >= warn.pct,
                .data$pct_top_n2 >= warn.pct) %in% TRUE,
      .rd = if (is.null(seasons$max_over_obs1)) FALSE else
        (.data$max_over_obs1 > ref.mult | .data$max_over_obs2 > ref.mult) %in% TRUE) %>%
    dplyr::filter(.data$.td | .data$.cd | .data$.rd) %>%
    dplyr::arrange(dplyr::desc(.data$.rd),
                   dplyr::desc(abs(log(.data$total_ratio))))

  flag.rows <- if (!nrow(flagged)) "" else paste(vapply(seq_len(nrow(flagged)),
    function(i) {
      r <- flagged[i, ]
      inc <- if (is.null(r$incumbent) || is.na(r$incumbent)) "&mdash;" else
        if (grepl(paste0("_", k1, "_"), r$incumbent)) k1 else k2
      row.html(c(
        td(.esc(r$species)), td(.esc(r$season)),
        td(.hnum(r$abund1), TRUE), td(.hnum(r$abund2), TRUE),
        td(.hnum(r$total_ratio), TRUE),
        td(if (is.null(r$max_over_obs1)) "&mdash;" else
             sprintf('<span class="win-a">%s</span> / <span class="win-b">%s</span>',
                     .hnum(r$max_over_obs1, 2), .hnum(r$max_over_obs2, 2)), TRUE),
        td(sprintf('<span class="fam fam-%s">%s</span>',
                   if (inc == k1) "a" else "b", .esc(inc)))),
        crit = isTRUE(r$.rd))
    }, character(1)), collapse = "\n")

  # ---- headline block -----------------------------------------------------
  if (is.null(worst)) {
    head.block <- sprintf(paste0(
      '<div class="finding is-ok"><p class="eyebrow">What the comparison found</p>',
      '<h2>No surface exceeds what was observed.</h2>',
      '<p>Across all %d species-seasons, neither family predicts a cell more ',
      'than %s times the highest density ever recorded on a segment of the same ',
      'species and season. Where the two families differ they differ in degree, ',
      'not in kind, and the table below says where.</p></div>'),
      nrow(seasons), format(ref.mult, big.mark = ","))
  } else {
    w <- worst$row
    other <- if (worst$fam == k1) k2 else k1
    a.tot <- w$abund1; b.tot <- w$abund2
    big <- if (worst$fam == k1) a.tot else b.tot
    sml <- if (worst$fam == k1) b.tot else a.tot
    n.over <- if (worst$fam == k1) w$n_over_obs1 else w$n_over_obs2
    p.over <- if (worst$fam == k1) w$pct_over_obs1 else w$pct_over_obs2
    mx <- if (worst$fam == k1) w$max1 else w$max2
    obs <- if (!is.null(w$obs_max1)) w$obs_max1 else NA_real_

    # A log axis spanning both totals, rounded out to whole decades, so the two
    # bars are honestly placed rather than merely ordered.
    lo <- floor(log10(max(min(sml, big), .Machine$double.xmin))) - 1
    hi <- ceiling(log10(max(sml, big))) + 1
    # Pick the number of tick labels so the decades divide evenly, rather than
    # stretching the axis until five of them fit: five labels over a nine-decade
    # span come out 10^5, 10^7, 10^10, 10^12, 10^14, and widening the span to
    # fix that pushes the bars away from the edge they should reach.
    gaps <- which((hi - lo) %% (3:5) == 0)
    n.tick <- if (length(gaps)) (3:5)[gaps[length(gaps)]] + 1 else 5
    # A span of one or two decades has fewer decades than labels, and asking
    # for five gives repeats (10^0, 10^0, 10^1, 10^2, 10^2).
    n.tick <- min(n.tick, hi - lo + 1)
    posn <- function(v) max(0, min(100, (log10(max(v, 10^lo)) - lo) / (hi - lo) * 100))
    ticks <- paste(sprintf("<span>10<sup>%d</sup></span>",
                           round(seq(lo, hi, length.out = n.tick))),
                   collapse = "")

    head.block <- sprintf(paste0(
      '<div class="finding"><p class="eyebrow">What the comparison found</p>',
      '<h2>%s in %s is not a tie. It is a runaway.</h2>',
      '<p>Under <span class="fam fam-%s">%s</span>%s the surface puts ',
      '<strong>%s cells</strong> more than %s&times; the highest density ever ',
      'recorded on a segment of this species and season (%s birds/km<sup>2</sup>), ',
      'peaking at %s. Those cells hold %s of the seasonal total. The ',
      '<span class="fam fam-%s">%s</span> finalist, fitted to the same segments ',
      'on the same grid, gives a total %s. The two surfaces correlate at %s.</p>',
      '<div class="decades" role="img" aria-label="Predicted totals on a log scale">',
      '<div class="decade-row"><span class="decade-label">',
      '<span class="fam fam-%s">%s</span></span><span class="decade-track">',
      '<span class="decade-fill a" style="width:%.1f%%;"></span>',
      '<span class="decade-val" style="left:%.1f%%;color:var(--ink);">&nbsp;%s birds</span>',
      '</span></div>',
      '<div class="decade-row"><span class="decade-label">',
      '<span class="fam fam-%s">%s</span></span><span class="decade-track">',
      '<span class="decade-fill b" style="width:%.1f%%;"></span>',
      '<span class="decade-val" style="right:%.1f%%;">%s birds&nbsp;</span>',
      '</span></div>',
      '<div class="decade-scale">%s</div></div>',
      '<p class="caption">Seasonal totals on a log scale &mdash; the bar length is ',
      'the exponent. On a linear scale the smaller bar would be invisible.</p>',
      '<dl class="stats">',
      '<div class="stat is-crit"><dt>Cells past the observed ceiling</dt>',
      '<dd>%s<span class="unit">of %s, under <span class="fam fam-%s">%s</span></span></dd></div>',
      '<div class="stat is-crit"><dt>Their share of the total</dt>',
      '<dd>%s<span class="unit">the surface <em>is</em> the outliers</span></dd></div>',
      '<div class="stat"><dt>Largest cell</dt><dd>%s<span class="unit">birds per km<sup>2</sup></span></dd></div>',
      '<div class="stat"><dt>Correlation of the two</dt><dd>%s<span class="unit">Pearson, over live cells</span></dd></div>',
      '</dl></div>'),
      .esc(w$species), tolower(.esc(w$season)),
      if (worst$fam == k1) "a" else "b", .esc(worst$fam),
      if (!is.null(w$incumbent) && !is.na(w$incumbent) &&
          grepl(paste0("_", worst$fam, "_"), w$incumbent))
        " &mdash; the model this species currently ships &mdash; " else ", ",
      format(n.over, big.mark = ","), format(ref.mult, big.mark = ","),
      .hnum(obs), .hnum(mx), .pct(p.over, 4),
      if (worst$fam == k1) "b" else "a", .esc(other),
      sprintf("%s times smaller", .hnum(big / sml)),
      .hnum(w$pearson_r, 2),
      "a", .esc(k1), posn(a.tot), posn(a.tot), .hnum(a.tot),
      "b", .esc(k2), posn(b.tot), 100 - posn(b.tot), .hnum(b.tot),
      ticks,
      format(n.over, big.mark = ","), format(w$cells1, big.mark = ","),
      "a", .esc(worst$fam),
      .pct(p.over, 4), .hnum(mx), .hnum(w$pearson_r, 2))
  }

  body <- sprintf(paste0(
'<!-- ARTIFACT:BEGIN -->
<div class="page">
<header class="masthead">
  <p class="eyebrow">%s &nbsp;&middot;&nbsp; pipeline step 02.15</p>
  <h1>Two Families, One Question</h1>
  <p class="lede">Cross-validation could not choose between the two response
  families for these species. Neither could the model diagnostics. So we asked
  the question neither of them asks &mdash; whether the tie matters &mdash; by
  predicting both surfaces and comparing them cell by cell.</p>
  <p class="meta"><span>%d species &times; %d seasons = %d comparisons</span>
  <span>%s prediction cells each</span><span>%s km&sup2; cells</span>
  <span>%s</span></p>
</header>

%s

<section>
  <h2>Everywhere else, the two families mostly agree</h2>
  <div class="prose">
    <p>Across all %d species-seasons the median ratio of seasonal totals is
    <strong>%s</strong>, the median rank correlation between the surfaces is
    <strong>%s</strong>, and on the log scale &mdash; the one that reflects the
    surface as it would be mapped &mdash; the median correlation is
    <strong>%s</strong>. For most species in most seasons the unresolved tie
    genuinely does not change the product.</p>
    <p>%s of the %d are exceptions, listed below: the seasonal totals differ by
    more than a factor of 1.5, the concentration crosses the %s%% level the
    shipping step warns at, or a surface runs past what was ever observed.</p>
  </div>
  <div class="table-wrap"><table>
    <caption>Species-seasons where the family choice changes the answer.</caption>
    <thead><tr><th>Species</th><th>Season</th><th class="n">%s total</th>
    <th class="n">%s total</th><th class="n">ratio</th>
    <th class="n">max &divide; observed max</th><th>currently ships</th></tr></thead>
    <tbody>%s</tbody>
  </table></div>
  <p class="caption">Totals are the summed density surface in birds. Ratio is
  %s &divide; %s. The last numeric column is each family&rsquo;s largest cell as a
  multiple of the highest density ever recorded on a segment of that species and
  season &mdash; below 1 is normal, since a prediction is a seasonal mean and an
  observation is one encounter.</p>
</section>

<section>
  <h2>Which family is spikier?</h2>
  <div class="prose">
    <p>Measured three ways over the %d species-seasons. The margins are not
    dramatic, and the exception is the one that matters most.</p>
  </div>
  <div class="tally">%s%s%s</div>
  <div class="prose">
    <p>The two families are also rarely spiky in the <em>same places</em>. Of the
    %d largest cells in each surface, the median overlap is <strong>%s</strong>
    &mdash; %d of the %d species-seasons share none at all. Agreement on the
    shape of a surface does not imply agreement on where the birds pile up.</p>
  </div>
</section>

<section>
  <h2>How the ceiling is set</h2>
  <div class="prose">
    <p>A flat density limit cannot serve every species. On this SubProject the
    highest density ever recorded on a segment ranges from %s to %s
    birds/km&sup2; between species-seasons, so one constant is far too permissive
    for the small species and barely above the legitimate range for the large
    ones.</p>
    <p>The ceiling used here is therefore <strong>species- and
    season-specific</strong>: each surface is judged against the busiest segment
    ever recorded for that species in that season, and a cell is counted when it
    exceeds %s&times; it. The median surface peaks <strong>%s&times;</strong> that
    reference &mdash; below it, as it should be, since a prediction cell is a
    seasonal mean over ten years and an observed segment density is a single
    encounter on a single day.</p>
  </div>
</section>

<footer>
  <p>Generated by <code>02.15_Compare_finalist_predictions.Rmd</code> from
  <code>DSMHelper::compare.finalist.surfaces()</code>. Per-species reports, the
  full %d-row table and the spike-cell coordinates are in the results folder
  beside this page.</p>
</footer>
</div>
<!-- ARTIFACT:END -->'),
    .esc(subproject), n_sp, n_se, nrow(seasons),
    format(seasons$cells1[1], big.mark = ","),
    format(cell.km2, big.mark = ","), format(Sys.Date(), "%d %B %Y"),
    head.block,
    nrow(seasons), .hnum(med(seasons$total_ratio)),
    .hnum(med(seasons$spearman_r)), .hnum(med(seasons$log_r)),
    nrow(flagged), nrow(seasons), format(warn.pct),
    .esc(k1), .esc(k2), flag.rows, .esc(k1), .esc(k2),
    nrow(seasons),
    bar(t_gini, "Higher Gini across the whole surface",
        "Inequality over every cell: 0 flat, 1 all in one place."),
    bar(t_mom, "Higher largest cell divided by median cell",
        "Scale-free, so a bigger total is not counted as spikier."),
    bar(t_tot, "Larger seasonal total",
        "Which family predicts more birds, season by season."),
    top_n, .hnum(med(seasons$top_n_shared)),
    sum(seasons$top_n_shared == 0, na.rm = TRUE), nrow(seasons),
    .hnum(min(c(seasons$obs_max1, seasons$obs_max2), na.rm = TRUE)),
    .hnum(max(c(seasons$obs_max1, seasons$obs_max2), na.rm = TRUE)),
    format(ref.mult, big.mark = ","),
    .hnum(med(c(seasons$max_over_obs1, seasons$max_over_obs2)), 2),
    nrow(seasons))

  writeLines(c(.finalist.page.head(), body, "</body>", "</html>"), path)
  invisible(path)
}
