
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



#' Build observation and watch tables from raw ECSAS data
#'
#' Takes raw ECSAS data, filters observations, clips both watches and
#' observations to the study area, and optionally assembles watches into
#' transects.  Returns a list with elements \code{distdata}, \code{watches},
#' and (if \code{create_transects = TRUE}) \code{transects}.
#'
#' Expects \code{study.area} (an \code{sf} polygon) to exist in the calling
#' environment.
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
                               intransect.only = TRUE) {


  coll = checkmate::makeAssertCollection()
  checkmate::assert_data_frame(raw.dat, add = coll)
  checkmate::assert(
    checkmate::check_class(study.area, "sf"),
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

  # clip to study area
  #
  # NB: rmapshaper::ms_clip() throws "Not compatible with STRSXP: [type=list]"
  # rather than returning an empty result when nothing overlaps, so check for
  # overlap first. This happens for a whole dataset when the study area has no
  # coverage by this survey type at all.
  watch_pts <- watches %>%
    sf::st_as_sf(coords = c("LongStart", "LatStart"), crs = sf::st_crs(inproj)) %>%
    sf::st_transform(sf::st_crs(4326)) %>% # for ms_clip below
    dplyr::select(WatchID) # just keep WatchID
  study_area_4326 <- study.area %>% sf::st_transform(sf::st_crs(4326))

  if (nrow(sf::st_filter(watch_pts, study_area_4326)) == 0) {
    warning(sprintf(
      "create.survey.data: no %s watches overlap study.area - returning no watches",
      dataset), immediate. = TRUE)
    watches <- watch_pts %>%
      dplyr::slice(0) %>%
      dplyr::left_join(watches, by = "WatchID") %>%
      sf::st_transform(outproj)
  } else {
    watches <- watch_pts %>%
      rmapshaper::ms_clip(study_area_4326) %>%   # do the clipping -
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

  # Clip to study area
  #
  # NB: see the note on the watches clip above - ms_clip() errors rather than
  # returning an empty result when nothing overlaps, so check for overlap
  # first. Unlike the watches case this can fire even when the survey type does
  # cover the study area, if it simply recorded no in-transect observations
  # inside it.
  obs_pts <- obs %>%
    sf::st_as_sf(coords = c("LongStart", "LatStart"), crs = sf::st_crs(inproj)) %>%
    sf::st_transform(sf::st_crs(4326)) %>% # for ms_clip below
    dplyr::select(object) # just keep object
  study_area_4326 <- study.area %>% sf::st_transform(sf::st_crs(4326))

  if (nrow(sf::st_filter(obs_pts, study_area_4326)) == 0) {
    warning(sprintf(
      "create.survey.data: no %s observations overlap study.area - returning no observations",
      dataset), immediate. = TRUE)
    obs <- obs_pts %>%
      dplyr::slice(0) %>%
      dplyr::left_join(obs, by = "object") %>%
      sf::st_transform(outproj)
  } else {
    obs <- obs_pts %>%
      rmapshaper::ms_clip(study_area_4326) %>%   # do the clipping -
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
#' @return \code{dat} with a \code{Season} factor column added.
#' @export
assign.season <- function(dat, season.def, datefield = "Date"){

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
    levels = season.names,
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
#' @param the.data Named list with at least a \code{watches} element (as
#'   returned by \code{\link{create.survey.data}}).
#' @param study.area \code{sf} polygon defining the study area.
#' @param inproj EPSG code or CRS for input watch coordinates.
#' @param outproj EPSG code or CRS for output segdata.
#' @param scale.factors Named list of mean and SD values used to standardise
#'   each covariate (e.g. \code{list(depth_mean = x, depth_sd = y, ...)}).
#' @param verbose If \code{TRUE}, print progress messages.
#' @return \code{sf} data frame of segment data with covariates attached.
#' @export
create.segdata <- function(the.data,
                           study.area,
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
  files <- file.path(predLayerStudyAreaDir, "sst", paste0("sst.", dates.needed, ".img"))
  sst <- terra::rast(files)

  # SST gradient
  files <- file.path(predLayerStudyAreaDir, "sst", paste0("sst.g.", dates.needed, ".img"))
  sst.g <- terra::rast(files)

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

  ### Extract SST values
  if (verbose) message("Extracting sst at watch locations")
  needed.layers <- paste0("sst.", segdata$MonthYear, "-16")
  segdata <-
    terra::extract(
      sst,
      terra::vect(segdata),
      layer = needed.layers,
      bind = TRUE
    ) %>%
    sf::st_as_sf() %>%
    dplyr::rename(sst = value) %>%
    dplyr::select(-layer) # Note that layer is off by one even though the sst values is correct

  ### Extract SST gradient values
  if (verbose) message("Extracting sst gradient at watch locations")
  needed.layers <- paste0("sst.g.", segdata$MonthYear, "-16")
  segdata <-
    terra::extract(
      sst.g,
      terra::vect(segdata),
      layer = needed.layers,
      bind = TRUE
    ) %>%
    sf::st_as_sf() %>%
    dplyr::rename(sst.g = value) %>%
    dplyr::select(-layer) # Note that layer is off by one even though the sst values is correct

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
#' Creates four seasonal copies of \code{predgrid}, computes seasonal mean
#' values for all dynamic covariates, drops the monthly columns, and then
#' replicates the grid once per level of the \code{platform} factor (one copy
#' per level of \code{ddftype_to_platform}).  The combined grid is saved as a
#' shapefile.
#'
#' Expects project globals \code{seasons}, \code{season.names},
#' \code{ddftype_to_platform}, and \code{ShapeDir} in the calling
#' environment.
#'
#' @param species Character string species code.
#' @param predgrid \code{sf} data frame representing the prediction grid (one
#'   row per cell, with monthly covariate columns).
#' @return \code{sf} data frame with one row per (cell × season × platform)
#'   combination, containing seasonal mean covariates.
#' @export
create.seasonal.predgrid <- function(species, predgrid) {
  # Get species-specific season setting and create predgrid.all.seas with 4 seasons
  season.spec <- seasons[[species]]
  ret <-
    rbind(predgrid, predgrid, predgrid, predgrid) %>%
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
  sf::st_write(ret,
               dsn = ShapeDir,
               layer = paste0(species, "_predgrid.shp"),
               driver = "ESRI Shapefile",
               delete_layer = TRUE
  )

  # Make a copy for each value of platform in model (technically only needed by
  # factor and fs models but make copy anyways and run.dsm.pred will handle it
  # correctly for nofactor model).
  ret <- replicate(length(unique(ddftype_to_platform)), ret, simplify = FALSE) %>%
    setNames(unique(ddftype_to_platform)) %>%
    purrr::list_rbind(names_to = "platform") %>%
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
#' @param subs Which platform subset to map: \code{"Combined"} (default),
#'   \code{"F"} (flying), or \code{"W"} (water).
#' @param ... Additional arguments passed to \code{\link{do.pred.map}}.
#' @return Named list of leaflet map objects, one per season.
#' @export
do.pred.maps <-
  function(dat,
           model,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           ...) {


    # Get data subset
    subs <- match.arg(subs)
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
#' @param subs Platform subset: \code{"Combined"} (default), \code{"F"}, or
#'   \code{"W"}.
#' @param ... Additional arguments passed to \code{do.pred.map.ggplot}.
#' @return A \code{patchwork} ggplot object.
#' @export
do.pred.maps.ggplot <-
  function(dat,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           ...) {


    # Get data subset
    subs <- match.arg(subs)
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
#' @param subs Platform subset label.
#' @param samp_n Integer; if not \code{NA}, plot a random sample of this many
#'   polygons.
#' @return A \code{ggplot} object.
#' @export
do.pred.map.ggplot <-
  function(season,
           dat,
           modname,
           species,
           subs = c("Combined", "F", "W"),
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
#' Expects project globals \code{GenericRDataDir}, \code{GenericShapeDir},
#' \code{predLayerDir}, \code{GISDir}, and \code{SubProject}.
#'
#' @param subproj Character string subproject identifier.
#' @return \code{invisible(NULL)}, called for its side-effect (directories
#'   created).
#' @export
create_subproject_folders <- function(subproj) {
  create.dir.if.needed(file.path(GenericRDataDir, subproj))
  create.dir.if.needed(file.path(GenericShapeDir, subproj))
  create.dir.if.needed(here::here("Results", subproj))
  create.dir.if.needed(file.path(predLayerDir, "Study area resolution & extent", subproj))
  create.dir.if.needed(file.path(GISDir, "Predictions", SubProject))
  create.dir.if.needed(file.path(GISDir, "Rasters", SubProject))

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



#' Short, stable key for a response family
#'
#' Normalises a family object's name so that fitted and unfitted objects agree
#' (\code{"Tweedie"} and \code{"Tweedie(p=1.167)"} both give \code{"tw"}), and so
#' that families beyond the two currently in \code{dsm.mod.specs} get a key of
#' their own rather than an error. Nothing downstream should assume the set of
#' keys is \code{c("tw", "nb")}.
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
#' Squares of \code{block_size} on the analysis projection are assigned to folds
#' at random, so each fold holds out whole regions rather than scattered
#' segments. Random hold-outs would leak neighbouring segments into training and
#' flatter every model equally, because the residuals are spatially structured.
#'
#' @param segdata Data frame with numeric \code{x} and \code{y} columns in
#'   projection units (metres).
#' @param n_folds Number of folds.
#' @param block_size Block edge length in projection units.
#' @param seed Random seed for block-to-fold assignment. The random number
#'   generator kind is pinned to Mersenne-Twister for the duration and restored
#'   afterwards, so the folds are the same whether or not the caller is running
#'   inside a \code{future} (which switches to L'Ecuyer-CMRG).
#' @return Integer vector of fold membership, length \code{nrow(segdata)}, with
#'   an \code{"n_blocks"} attribute giving the number of occupied blocks. That
#'   count is worth reporting: it is the real sample size behind the folds, and a
#'   block size close to the study area extent can leave too few blocks to spread
#'   across \code{n_folds}.
#' @examples
#' assign.blocks(data.frame(x = runif(100, 0, 5e5), y = runif(100, 0, 5e5)),
#'               5, 1e5, 1)
#' @export
assign.blocks <- function(segdata, n_folds, block_size, seed) {
  checkmate::expect_data_frame(segdata, min.rows = 1)
  checkmate::expect_numeric(segdata$x, any.missing = FALSE)
  checkmate::expect_numeric(segdata$y, any.missing = FALSE)
  checkmate::expect_count(n_folds, positive = TRUE)
  checkmate::expect_number(block_size, lower = 0)

  # Pin the generator, not just the seed. Under furrr's seed = TRUE the workers
  # run L'Ecuyer-CMRG, and set.seed() on a different generator yields a
  # different sequence - which silently produced different folds inside workers
  # than in the main process, from an identical seed.
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  set.seed(seed)
  block_id     <- paste(floor(segdata$x / block_size),
                        floor(segdata$y / block_size), sep = "_")
  unique_block <- unique(block_id)
  # rep_len then sample spreads folds as evenly as the block count allows.
  block_fold   <- stats::setNames(
    sample(rep_len(seq_len(n_folds), length(unique_block))), unique_block)

  structure(unname(block_fold[block_id]), n_blocks = length(unique_block))
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

  key <- get.family.key(fit$family)
  n   <- length(mu)

  if (key == "tw") {
    # getTheta(TRUE) returns p for tw(); the scale lives in $sig2.
    matrix(mgcv::rTweedie(rep(mu, n_sim),
                          p = fit$family$getTheta(TRUE), phi = fit$sig2),
           nrow = n)
  } else if (key == "nb") {
    matrix(stats::rnbinom(n * n_sim, mu = rep(mu, n_sim),
                          size = fit$family$getTheta(TRUE)),
           nrow = n)
  } else {
    stop("sim.response: no predictive simulation defined for family '",
         fit$family$family, "'")
  }
}


#' Spatial-block cross-validation of one species' family finalists
#'
#' Step two of the two-step model selection. Refits each finalist from
#' \code{\link{get.family.finalists}} on every spatial-block training set and
#' scores it on the held-out block, using criteria that compare across response
#' families because they live on the response scale.
#'
#' Every finalist sees identical folds, so the comparison is paired. Models are
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
#' @param nthreads Passed to \code{mgcv::bam()}.
#' @param calib.covar Column in \code{segdata} to aggregate the
#'   observed-vs-expected calibration check over, or \code{NULL} to skip it.
#'   Defaults to \code{"platform"}, which is present in the segment data whether
#'   or not the formula uses it -- factor-smooth and no-factor formulations still
#'   get a calibration score.
#' @param progress If \code{TRUE}, report each fold as it completes.
#' @return A tibble, one row per finalist per fold: \code{species},
#'   \code{family_key}, \code{modname}, \code{fold}, \code{n_test}, \code{CRPS},
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
                          seed = 1, nthreads = 1, calib.covar = "platform",
                          progress = TRUE) {
  checkmate::expect_character(species, len = 1, any.missing = FALSE)
  checkmate::expect_data_frame(segdata, min.rows = 1)
  checkmate::expect_data_frame(finalists, min.rows = 1)
  if (!all(c("family_key", "modname", "formula", "family") %in% names(finalists)))
    stop("run.family.cv: finalists must come from get.family.finalists()")

  fold_id <- assign.blocks(segdata, n_folds, block_size, seed)
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
      fit_args <- list(finalists$formula[[i]],
                       data   = segdata[train_idx, ],
                       family = finalists$family[[i]],
                       method = "fREML")

      fit <- try(do.call(mgcv::bam, c(fit_args, list(discrete = TRUE,
                                                     nthreads = nthreads))),
                 silent = TRUE)
      discrete_used <- TRUE

      if (inherits(fit, "try-error")) {
        message(sprintf("run.family.cv: %s %s fold %d failed under discrete = TRUE (%s); retrying without",
                        species, finalists$modname[i], k,
                        sub("
.*", "", fit[1])))
        fit <- try(do.call(mgcv::bam, c(fit_args, list(discrete = FALSE))),
                   silent = TRUE)
        discrete_used <- FALSE
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
      pred_mu <- stats::predict(fit, newdata = segdata[test_idx, ],
                                type = "response", discrete = FALSE)
      obs_y   <- segdata[[resp]][test_idx]

      # Same reasoning as assign.blocks(): fix the generator as well as the seed,
      # so scores are reproducible whether or not the caller is inside a future.
      old_kind <- RNGkind()
      suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
      set.seed(seed + k)
      sim_mat <- sim.response(fit, pred_mu, n_sim)
      RNGkind(old_kind[1], old_kind[2], old_kind[3])

      # 90% predictive interval coverage: the direct check on whether a family's
      # tails are too light to carry prediction variance downstream.
      lower <- apply(sim_mat, 1, stats::quantile, 0.05)
      upper <- apply(sim_mat, 1, stats::quantile, 0.95)

      row <- tibble::tibble(
        species    = species,
        family_key = finalists$family_key[i],
        modname    = finalists$modname[i],
        fold       = k,
        n_test     = length(test_idx),
        CRPS       = mean(calc.crps(sim_mat, obs_y)),
        MAE        = mean(abs(obs_y - pred_mu)),
        RMSE       = sqrt(mean((obs_y - pred_mu)^2)),
        PIT_KS     = as.numeric(suppressWarnings(
                       stats::ks.test(calc.pit(sim_mat, obs_y), "punif")$statistic)),
        cover90    = mean(obs_y >= lower & obs_y <= upper),
        discrete   = discrete_used)

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
#'   \code{pred} (numeric prediction vector).
#' @export
apply.dsm.var <- function(dat, this.dsm){
  res <- dsm::dsm_var_gam(this.dsm, dat, purrr::map(dat, ".my.off.set"))

  list(pred.var = res$pred.var, pred = unlist(res$pred))
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
#' @param nchunks Integer number of chunks to split \code{df} into.
#' @param off.set Numeric offset (cell area): either a scalar or a vector the
#'   same length as \code{nrow(df)}.
#' @param parallel If \code{TRUE}, process chunks in parallel with
#'   \code{parLapply}.
#' @param nodes Integer number of cluster nodes for parallel processing.
#' @return List of per-chunk results from \code{\link{apply.dsm.var}}.
#' @export
get.per.cell.var <- function(this.dsm,
                             df,
                             nchunks,
                             off.set = 1,
                             parallel = FALSE,
                             nodes = 1
) {

  # If off.set is a vector add it to the df so it will be split properly as well.
  if (length(off.set) > 1 && length(off.set) != nrow(df)) {
    stop("get.per.cell.var: off.set vector is not same length as df")
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

    # Could just execute the things we need instead.
    # parallel::clusterEvalQ(cl, source(here::here("R/analysis_settings.R")))
    parallel::clusterEvalQ(cl, {
      library(dsm)
      library(purrr)
    })

    # Need envir arg or else it won't find data objects when being rendered.
    # takes about 1 min
    parallel::clusterExport(
      cl,
      c(
        "apply.dsm.var",
        "dat.split",
        "this.dsm"
      ),
      envir = environment()
    )

    # Run the function
    system.time(res <- parallel::parLapply(
      cl,
      dat.split,
      apply.dsm.var,
      this.dsm = this.dsm
    )
    )

    parallel::stopCluster(cl)
  } else {  # non-parallel version
    # apply the function to the chunks serially with map
    print(system.time(
      res <- purrr::map(dat.split, \(chunk) {
        apply.dsm.var(chunk, this.dsm)
        gc()
      }
      )
    ))
  }

  res
}


#' Compute a density estimate with uncertainty from a DSM
#'
#' \strong{Note: not currently used.}
#'
#' @param dsm_final Fitted \code{dsm} object.
#' @param predgrid Prediction grid data frame.
#' @return Named list with elements \code{pred.est}, \code{CV}, \code{SE},
#'   and \code{CI} (a three-element vector giving the 5%, mean, and 95%
#'   lognormal confidence interval).
#' @export
get.dens.est <- function(dsm_final, predgrid) {

  # use dsm.var.gam to get estimated abundance params
  densEst <- summary(dsm::dsm.var.gam(dsm_final, predgrid, off.set = predgridCellArea))

  # The estimates in dsm.var.gam are summed for the entire predgrid study area
  # so we need to divide by the number of predgrid cells * cellarea.
  densEst %<>%
    purrr::map_at(c("pred.est", "se"), ~ .x/(predgridCellArea * nrow(predgrid)))

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


#' Render the extrapolation analysis report for one dataset
#'
#' @param spill Character string spill/project identifier used in the output
#'   filename.
#' @param dataset Character string dataset name passed as a render parameter.
#' @param debug If \code{TRUE}, drop into \code{browser()} at the start.
#' @return \code{invisible(NULL)}, called for its side-effect (HTML rendered).
#' @export
do.extrapolation <- function(spill, dataset, debug = FALSE){
  browser(expr = debug)

  out.file <- file.path(here::here(), paste(spill, dataset, "0_extrapolation.html", sep = "_"))
  message(sprintf("Rendering generic extrapolation  for %s to %s", dataset, out.file))
  rmarkdown::render(file.path(here::here(), "Generic_0_extrapolation.Rmd"),
                    params = list(spill = spill, dataset = dataset),
                    output_file = out.file)
}


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
#' @return \code{invisible(NULL)}, called for its side-effects (plots and
#'   printed output).
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

  # Remove "dsm" class to make DHARMa happy
  message("DHARMa checks")
  simmod <- dsm_final
  class(simmod) <- class(simmod)[-1] # Simulate residuals doesn't like dsm class
  sims <- DHARMa::simulateResiduals(fittedModel = simmod)
  plot(sims)
  print(DHARMa::testResiduals(sims))
  print(DHARMa::testZeroInflation(sims))
  # May fail if x,y locations are not unique
  res <- try(DHARMa::testSpatialAutocorrelation(sims, segdata$x, segdata$y))
  if (!inherits(res, "try-error"))
    print(res)

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
    try(print(mgcv::concurvity(dsm_final) %>% round(digits = 3)))
    message(
      "Concurvity of pairwise terms ('estimate' measure presented)\nEach row shows how terms in columns depend on the term in that row."
    )
    try(print(mgcv::concurvity(dsm_final, full = FALSE)[["estimate"]] %>% round(digits = 3)))
    message("Plot is non-symmetric, showing how terms on y-axis depend on terms on the x-axis")
    try(dsm::vis_concurvity(dsm_final))
  }

  par(mfrow = c(1,1))
  # check observed vs expected. See Miller et al 2021 pg 11 (of 18) for
  # reccommendation to use the "platform" variable to aggregate
  # by.
  message("Observed vs expected plot")
  try(print(oe.dens(dsm_final, covar = "platform", plotit = T)))
  # depth is continuous, so bin it - without cut, oe.dens aggregates by every
  # unique depth value, giving one point (and one table column) per segment.
  try(print(oe.dens(dsm_final, covar = "depth", cut = 10, plotit = T)))
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
  print(OD_dsm_final <- sum(resid(dsm_final, type = "pearson")^2)/dsm_final$df.res)

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


#' Generate density predictions from a DSM model
#'
#' Calls \code{predict} on the model named \code{modname} from \code{mod.res}
#' for all rows of \code{predgrid} (all seasons × platform levels), computes
#' a \code{"Combined"} subset by summing across platform levels, writes a
#' combined shapefile, and saves four-season rasters as GeoTIFFs.
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
    # from the *declared* global ddftype_to_platform in
    # create.seasonal.predgrid(). They agree as long as every declared platform
    # occurs in the data, but a study area covered by only one survey type (or
    # a species never recorded in one behaviour class) can break that, and
    # predict() would then fail deep inside mgcv with "factor has new levels".
    # Fail here instead, where the cause is obvious.
    mod_levels <- model$xlevels$platform
    if (!is.null(mod_levels)) {
      extra <- setdiff(unique(as.character(predgrid$platform)), mod_levels)
      if (length(extra) > 0)
        stop(sprintf(
          paste0("dsm.pred: predgrid has platform level(s) %s that model '%s' ",
                 "was not fitted with (model has %s). The predgrid is built ",
                 "from the global ddftype_to_platform, so prune it to the ",
                 "survey types actually present."),
          paste(sQuote(extra), collapse = ", "), modname,
          paste(sQuote(mod_levels), collapse = ", ")))
    }

    ret <- predgrid %>%
      dplyr::mutate(subset = platform)

    # Predict from model for all seasons and values of platform, which
    # may or may not be a factor in the given model. So, in a no_factor model
    # we will end up with n identical copies of predictions: 1 for each value
    # of platform. In a factor model these predictions for each level of platform
    # are different.
    ret <- ret %>%
      dplyr::mutate(NHat = predict(model, newdata=ret, off.set=ret$area),
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

  # Deal with ddfobj numbering when there are ddfs that had no observations
  # Create a conversion vector, conv, which will be used to set the ddfobj variable
  # below. There are two cases:
  #
  # 1) for segdata copies whose corresponding ddf actually had observations, the
  # original ddfobj number (set from ddftype) may not be right if there were any
  # preceding ddfs in df.mod.specs with no obs. Thus, we set the ddfobj to the
  # matching renumbered ddfobj created in distdata processing above
  #
  # 2) for segdata copies whose corresponding ddf had no observations, we need
  # to change the ddfobj. Otherwise, dsm will get upset when it tries to find
  # the observations forthese segments whose ddfobj points to a ddf (dummy or
  # otherwise???) with no observations due to a sanity check in dsm which
  # probably should be modified. To get around this, for segments whose ddfobj is
  # currently a ddf with no obs, we change the ddfobj for those segments to
  # point to any valid ddfobj which actually does have observations (for
  # exammple, below I use the first good ddfobj min(ddfobj)). Then, the abundance.est
  # in these segments (as computed by dsm:::make.data()) will still be 0 since
  # none of the Sample.Labels in the fitted.distdata for that substituted ddfobj
  # will match the Sample.Label in these segments (ie in this copy of segdata).
  # (XXXX Does having all these 0's bias the gam???? See Notes.docx May 23, 2025).
  #
  # Note it might be better just to these segdata copies with no obs altogether
  # but I'm not sure this is valied (nor am I sure keeping them is valid). Also,
  # that would mean that the number of init.segdata copies in the final segdata
  # would vary dynamically and this would need to be kept track of so that
  # downstream code (ie. in prediction step) would know how to parse out the
  # rows properly. Not sure what would be best
  #
  # conv can then be indexed by the original ddftype  to get the correct ddfobj
  # for each segment.
  nddfs <- length(def.ddf.list)
  origddfs <- sort(unique(distdata$ddfobj.orig))
  newddfs <- sort(unique(distdata$ddfobj))
  conv <- rep(NA, nddfs)
  conv[origddfs] <- newddfs
  conv[is.na(conv)] <- min(newddfs)

  # Cylce through the list of ddf model specs adding a copy of segdata for each.
  segdata <- df.mod.specs %>%
    # Insert a copy of the right segdata (ship vs aerial) for each ddf spec.
    # regardless of whether there were any obs in that ddf. If there weren't
    # then all segments for this ddf spec will end up with zero obs in the
    # response variable created by dsm:::make_data()
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
  # of each segdata copy (either aerial of ship) and making sure
  # get ddfs with no obs
  mingood <- min(newddfs)
  conv[mingood] <- NA # ignore the actuall good ddf that the others point to
  no_obs_ddfs <- which(conv == mingood)

  # Get survey type of each ddf. This requires that the order of ddftype_levels
  # and def.ddf.list have the same order.
  survey_type_index <- substr(ddftype_levels, 1,1)

  # get aerial and ship initial segdata sizes and name them "A" and "S".
  #
  # NB: table() only returns entries for survey types actually present in the
  # data, so a study area covered by only one survey type (eg. one with no
  # aerial coverage at all) leaves sizes["A"] as NA. That NA propagates through
  # sum() into tot_no_obs_ddf_segs, making the check below evaluate to NA and
  # abort with "missing value where TRUE/FALSE needed". Build the vector over
  # every survey type implied by ddftype_levels so absent ones are 0, not NA.
  expected_types <- unique(survey_type_index)
  sizes <- table(init_segdata$SurveyType) %>%
    setNames(names(.) %>% substr(1,1))
  sizes <- setNames(as.vector(sizes[expected_types]), expected_types)
  sizes[is.na(sizes)] <- 0

  # Get the number of segments that had a ddf with no obs that have used
  # min(newddfs) as their ddfobj.
  tot_no_obs_ddf_segs<- sum(sizes[survey_type_index[no_obs_ddfs]])

  # The number of segments with ddfobj equal to mingodd should equal the number
  # it would normally have had if no others from no_obs ddfs pointed to it
  # plus the number no_obs_ddf segments that are pointing to it. If not, something
  # went wrong.
  if (table(segdata$ddfobj)[mingood] !=
      sizes[survey_type_index[mingood]] + tot_no_obs_ddf_segs) {
    message(
      "Something went wrong assigning ddfobjs. Number of segments with ",
      mingood,
      " for ddfobj is ",
      table(segdata$ddfobj)[mingood],
      " but should be ",
      sizes[survey_type_index[mingood]] + tot_no_obs_ddf_segs
    )
    stop()
  }

  ### Create ddfs ----------------------------------------

  # Extract ddfs from df.mod.specs and drop any which don't have any observations
  # or else dsm() will get upset. Note we use the original ddfobj numbering
  # before it was recoded since this extracts the correct ddfs.
  ddfs <- purrr::map(df.mod.specs, "fitted.model") %>%
    magrittr::extract(unique(distdata$ddfobj.orig))

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
