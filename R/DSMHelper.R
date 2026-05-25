# DSMHelper - Distance Sampling and Density Surface Modelling Utilities
# Combined from functions.R and ds_utils_0.6.R (Atlantic DSM project)
#
# TODO - functions that reference project-specific globals and need updating
# before they will work outside the Atlantic DSM project:
#   create.survey.data()        - references segProj, study.area
#   create.segdata()            - references segProj, study.area, the.data
#   do.generic.render()         - references ResultsDir, RDir
#   do_oneoff_render()          - sources 'analysis settings.r' directly
#   do.full.prelim()            - references spec.grps and other globals
#   create_subproject_folders() - references GenericRDataDir, GenericShapeDir,
#                                 predLayerDir, GISDir, SubProject
#   get.seas.mean()             - references predgrid
#   various env covar fns       - reference segProj, study.area, RasterDir
# ============================================================================

# ============================================================================
# From functions.R
# ============================================================================
# Make sure both distance and distbegin/distend columns are not both in data
# otherwise ds() (as of Distance 1.0.9) will get mad. If distbegin and distend
# are both all NA (as for ECSAS ship data) then just remove them. If they're not
# NA then (e.g., in the case of aerial data) then get rid of distance after 
# making sure distance == (distbegin + distend)/2
check.distdata.cols <- function(data) {
  if (!is.null(data$distance) && !is.null(data$distbegin) && !is.null(data$distend)) {
    if (all(is.na(data$distbegin)) && all(is.na(data$distend))) {
      message("check.distdata.cols: removing distbegin and distend columns that are all NA from data")
      data <- select(data, -distbegin, -distend)
    } else if (all.equal(data$distance, (data$distbegin + data$distend) / 2)) {
      message(
        "check.distdata.cols: removing distance column because data has distbegin and distend columns"
      )
      data <- select(data, -distance)
    } else {
      stop(
        paste0(
          "run.ddf.model: data has both distance and non-NA distbegin/distend ",
          "columns but data$distance != (data$distbegin + data$distend)/2"
        )
      )
    }
  }
  data
}    


# rel.folder = path (relative to here()) to place the output file in
check.problem.data = function(alldat,
                              rel.folder,
                              leave = NULL,
                              fixed_db = NULL,
                              quiet = TRUE,
                              filename = "suspicious_posns.csv",
                              ...) {
  
  # Look for watches with problematic watches
  probs <- ECSAS.find.suspicious.posn(alldat, leave = leave, ...) %>%
    mutate(prev_fixed = WatchID %in% fixed_db)

  if (nrow(probs) > 1) {
    message (sprintf("%d problem watches found. %d of these were flagged as previously fixed in db",
                     nrow(probs), sum(probs$prev_fixed)))

    if (!quiet) {
      message("Problem watchIDs:")
      print(probs$WatchID)

      message("Problem watchIDs that were flagged as previously fixed in db:")
      print(filter(probs, prev_fixed == TRUE)$WatchID)

      hist(probs$dist_diff_km)
      hist(probs$pct_diff)
    }

    probs %>%
      select(CruiseID, WatchID, ObserverName, PlatformName, Date, StartTime, EndTime,
             LatStart, LongStart, LatEnd, LongEnd, WatchLenKm, PlatformSpeed, CalcDurMin,
             prev_fixed, dist_dr_km, dist_geo_km, dist_diff_km, pct_diff) %T>%
      write_csv(file = here(rel.folder, filename))
  } else {
    probs <- NULL

    # None found - clean up old files
    message("No problem watches found.")

    # create empty file or truncate if it exists
    file.create(file.path(here(rel.folder, filename)))

    #remove shapefile
    # file.remove(list.files(ShapeDir, pattern = paste0(filename, "\\..*"), full.names = T))
  }

  probs
}



# Creates data structure needed for distance sampling and dsm analysis.
#
# Takes a set of raw data and creates separate obs, and watches clipped to the
# study area and stores them along with the distdata objects in a list called
# "the.data" and returns it. file.prefix is a string that is used to form the
# filenames for various shapefiles. Typically either "ECSAS.ship" etc.
#
# dataset = character string identifying the dataset. Good for multi dataset
# analyses.
# inproj, outproj - projections of input and output datasets respectively
# file.prefix - filename prefix for obs and watches shapefiles that gets
#     prepended to "_obs.shp" and "_watches.shp" to form the output shapefile
#     names.
# saveshp - whether to save shapefiles or not.
# intransect.only - only include observations where InTransect was TRUE?
#
# Note that obs with no count, no distmeth, non-birds, FlySwim != "W" or "F" or
#   association == 18 are filtered out. 
create.survey.data <- function(raw.dat = NULL, 
                               dataset = NULL,
                               file.prefix = NULL,
                               inproj = 4326,
                               outproj = segProj,
                               saveshp = TRUE,
                               intransect.only = TRUE) {

  
  coll = makeAssertCollection()
  assert_data_frame(raw.dat, add = coll)
  assert(
    check_class(study.area, "sf"),
    add = coll
  )
  assert(
    check_string(file.prefix),
    add = coll
  )
  reportAssertions(coll)

  
  # create watches
  message("Creating watches...")
  watches <- raw.dat %>%
    select(
      SurveyType,
      TransectID,
      Program,
      CruiseID,
      WatchID,
      ObserverName,
      TransFarEdge,
      DistMeth,
      Date,
      StartTime,
      EndTime,
      LatStart,
      LongStart,
      LatEnd,
      LongEnd,
      WatchLenKm,
      ObsHeight,
      CalcDurMin,
      TotalWidthKm,
      WhatCount
    ) %>%
    mutate(
      Sample.Label = WatchID,
      # For ecsas, TransectSides is 1 for ship, and assigned elsewhere in 
      # Extract_data.Rmd for air. For SOMEC, both aerial and ship data pass through
      # here and all aerial surveys are 2 sided whereas ship are 1.
      TransectSides = case_when(dataset == "SOMEC" &
                                  SurveyType == "Aerial" ~ 2, .default = 1)
    ) %>% 
    # Combine multiple observers on a watch to a single value
    group_by(across(-ObserverName)) %>%
    summarise(
      ObserverName = paste(unique(ObserverName), collapse = ", "),
      .groups = "drop"
    ) %>% 
    arrange(CruiseID, Sample.Label, ObserverName, Date, StartTime)
  
  # Make sure we didn't loose any watches
  if(length(setdiff(raw.dat$WatchID, watches$WatchID)) != 0) {
    missing_ids <- setdiff(raw.dat$WatchID, watches$WatchID)
    stop("Create.survey.data: lost the following WatchIDs after collapsing watches with multiple observers ", 
         paste(missing_ids, collapse = ", "))
  }
  
  # make sure all data for a watch is consistent. Find rows with duplicate watchIDs
  # and remove these watches
  dups <- watches %>%
    group_by(WatchID) %>% 
    summarise(nrows = n()) %>% 
    filter(nrows > 1) %>% 
    pull("WatchID")
  
  # Remove watches with inconsistent watch info
  if (length(dups) > 0) {
    warning(paste0(sprintf("Removing %d watches due to inconsistent watch info between rows: ", length(dups)),
                   paste(dups, collapse = ", ")), immediate. = TRUE)
    watches <- filter(watches, !(WatchID %in% dups))
  }
  
  # clip to study area
  watches <- watches %>%
    st_as_sf(coords = c("LongStart", "LatStart"), crs = st_crs(inproj)) %>% 
    st_transform(st_crs(4326)) %>% # for ms_clip below
    select(WatchID) %>% # just keep WatchID
    ms_clip(study.area %>% st_transform(st_crs(4326))) %>%   # do the clipping -
    left_join(watches, by = "WatchID") %>%  # add other cols back in
    st_transform(outproj) %>% 
    mutate(StartTime = as.character(StartTime),
           EndTime = as.character(EndTime),) 
  
  # Save as shapefile
  if (saveshp) {
    # Will whine about discarded datum and abbreviated field names until these
    # warnings are removed from new rgdal.
    layer.name <- paste0(file.prefix, "_watches.shp")
    message(sprintf(
      "Saving watches to shapefile '%s'",
      file.path(ShapeDir, layer.name)
    ))
    
    suppressWarnings(
      st_write(
        watches,
        dsn = ShapeDir,
        layer = layer.name,
        driver = "ESRI Shapefile",
        delete_layer = TRUE
      )
    )
  }

  # Create Obs: remove non-birds, convert distances to km,  windforce is
  # converted to windspeed if necessary, convert InTransect to T/F, remove ship
  # followers, remove zeros (ie. watches where Count == NA), 
  # only keep Flyswim == W or F (not L or S). Add in sample.Label,
  # rename columns and select ones of interest
  message("Creating observations...")
  
  obs <- raw.dat %>%
    mutate(
      Region.Label = 1,
      InTransect = case_when(InTransect == -1 ~ TRUE,
                             InTransect == 0 ~ FALSE,
                             TRUE ~ NA),
      Distance = Distance / 1000
    ) %>%
    filter(
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
    mutate(DistType = assign.dist.type(.),
           FlockID = case_when(all(is.na(FlockID)) ~ 1:nrow(.),
                               TRUE ~ FlockID)) %>% 
    left_join(watches[, c("WatchID", "Sample.Label")], by = "WatchID") %>%
    rename(object = FlockID,
           size = Count,
           distance = Distance) %>% 
    select(
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
    mutate(
      Windspeed = case_when(
        is.na(Windspeed) &
          !is.na(Windforce) ~ left_join(., beaufort_conversion,
                                        by = c("Windforce" = "beaufort"))$speed.kts,
        TRUE ~ Windspeed
      ),
      # Only applicable to standard ECSAS ship surveys. Not currently used for 
      # anything.
      weights = case_when(
        DistanceCode %in% c("A", "B") ~ 2,
        DistanceCode %in% c("C", "D") ~ 1,
        TRUE ~ 0
      ),
      Visibility = case_when(Visibility > 20 ~ 20,
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
  obs <- obs %>%
    st_as_sf(coords = c("LongStart", "LatStart"), crs = st_crs(inproj)) %>% 
    st_transform(st_crs(4326)) %>% # for ms_clip below
    select(object) %>% # just keep WatchID
    ms_clip(study.area %>% st_transform(st_crs(4326))) %>%   # do the clipping -
    left_join(obs, by = "object") %>%  # add other cols back in
    st_transform(outproj)  
 
  # Save as shapefile
  if (saveshp) {
    # Will whine about discarded datum and abbreviated field names until these
    # warnings are removed from new rgdal.
    layer.name <- paste0(file.prefix, "_obs.shp")
    message(sprintf(
      "Saving obs to shapefile '%s'",
      file.path(ShapeDir, layer.name)
    ))
    
    suppressWarnings(
      st_write(
        obs,
        dsn = ShapeDir,
        layer = layer.name,
        driver = "ESRI Shapefile",
        delete_layer = TRUE
      )
    )
  }
  
  ##### After all that, now just create a single distdata for use in distance
  #sampling.This doesn't actually get used (unless we do analysis of all
  #seabirds). Rather, specific datasets for a given species/group are created by
  #create.dsm.data( in order to properly generate the samples where there were 0
  #observations of a given species.
  message("Creating distdata....")
  distdata <- obs %>% 
    st_drop_geometry() %>% 
    droplevels

  # Set dataset attribute
  distdata$Dataset <- watches$Dataset <- dataset
  
  the.data <-
    list(
      distdata = distdata,
      watches = st_drop_geometry(watches)
    )
  
  message("\nDone\n")
  the.data
}

# Convert watches object to spatialinesdataframe
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
    qq.gam(b, rep = rep, level = level, type = type, rl.col = rl.col, 
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
  kchck <- k.check(b, subsample = k.sample, n.rep = k.rep)
  if (!is.null(kchck)) {
    cat("Basis dimension (k) checking results. Low p-value (k-index<1) may\n")
    cat("indicate that k is too low, especially if edf is close to k'.\n\n")
    printCoefmat(kchck, digits = 3)
  }
}

# Do observed/expected for density response or abundance response
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


# Render a full preliminary dsm analysis Rmd for one species.
# XXX Used?
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
  render(
    file.path(RDir, "Generic_full_test.Rmd"),
    params = list(species = species),
    output_file = out.file
  )
}


# Set default convert.units, cutpoints, etc values when setting up ddf model
# list.
set.def.df.spec.values <- function(ddf.def, nm) {
  
  if(length(nm) != 1) 
    stop("set.def.df.spec.values: length of ddf name is not 1.")
  
  
  # Note that the convert_units won't actually get used because we 
  # only use ds() to compute the detection function and not
  # abundance. 
  if(str_detect(nm, "Ship")) {
    # ECSAS and SOMEC Ship surveys have set cutpoints whereas aerial uses distbegin
    # and distend columns computed in 00.01_Extract_data.Rmd.
    ddf.def$convert_units = ecsas.ship.convert.units
    ddf.def$cutpoints = ecsas.ship.ctpoints
    ddf.def$distance.centers = ecsas.ship.distance.centers
  } else if (str_detect(nm, "Aerial")) {
    ddf.def$convert_units = ecsas.air.convert.units
    ddf.def$distance.centers = ecsas.air.distance.centers
  } else
    stop("set.def.df.spec.values: unrecognized ddf name.")
  ddf.def
}

# Parse a ddf name of the form "ECSAS_Ship_F_D" and return a list of the 
# 4 elements.
# dataset - ECSAS, Quebec, ???
# platform_class - ship/aerial
# behav - fly, swim
# dist_type - distances/no_distances
parse.df.name <- function(nm){
  res <- str_split_1(nm, fixed("_"))
  list(dataset = res[1], platform_class = res[2], behav = res[3], dist_type = res[4])
}


# Three functions to peel off specific aspects of a ddf spec name
# ECSAS, Quebec, etc.
get.dataset <- function(nm) {
  parse.df.name(nm)$dataset
}

#Ship or Aerial
get.platform.class <- function(nm) {
  parse.df.name(nm)$platform_class
}

#F or W
get.behav <- function(nm) {
  parse.df.name(nm)$behav
}

# D or N
get.dist_type <- function(nm) {
  parse.df.name(nm)$dist_type
}

# Process all ddf specs for a given species. Called from Generic_1_ddf_fitting.Rmd
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
    imap(
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

# Save ddf specs (multiple specs - one for each dataset, platform_class and behav)
# for a single species.
save.ddf.specs <- function(df.mod.specs, species) {
  filename <- file.path(RDataDir, paste0(species, dfModListSuffix))
  message(sprintf("Saving ddf specs for %s in '%s'", species, filename))
  save(df.mod.specs, file = filename)
}


# Do ddf fitting for a given df.spec (which will be for one
# dataset/platform_class/behav for 1 species), extracting the correct data, building
# either all candidate ddfs, or just the final one. If final.only, store the
# returned ddf in the df.mod.list structure.
#
# df.spec - the ddf specification (see my structure in analysis_settings.R)
# 
# df.spec.nm - name of the df.spec - e.g., ECSAS_Ship_W
# 
# species - char vector 
# 
# distdata - should only contain obs for given species, but we filter on species
#   anyway just in case.
# 
# Typically called once for each ddf specification in each species in df.mod.list
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
    distdata <- filter(
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
    distdata <- filter(
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
    distdata %<>% filter(!(CruiseID %in% fishingCruises))

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

# Return a dummy_ddf() for various situations
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
      df_final <- dummy_ddf(
        distdata$object %>%
          str_replace("[a-zA-Z_]+", "") %>%
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
      df_final <- dummy_ddf(
        distdata$object %>%
          str_replace("[a-zA-Z_]+", "") %>%
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

# Called to run final ddf for one specie and dataset combination. Run
# ds() analysis either for all possible models or for just the final model
# specified by df.model.
#
#
# distdata = the observation dataframe to be subsetted for the analysis
# species = char vector
# df.model = det fcn model specification used if do.final.only is true.
# dataset = chr string indicating which dataset (ECSAS, SOMEC, etc) we're analyzing. 
# do.final.only = if true then the final chosen det fcn is fitted else all 
#   candidate det fcns are run
# rerun = logical, if true any previously run candidate det fcns will be re run.
#   Previously run detection functions are those that have results in folder.
# folder = where to save the df summaries. 
# parallel = logical, should parallel processing be used to run multiple
#   candidate det fcns at the same time.
# nCores = number of (virtual) cores to use for parallel processing
# ... = other args passed on to ds().
# 
# Return list with df_final and distdata (augmented with detProb and adjSize)
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
        
        # remove unneeded cols, but be careful because some cols are needed
        # downstream by create.dsm.data) for example.
        # and remove redundant distance or distbegin/distend cols.
        distdata <-
          distdata %>% 
          dplyr::select(object, size, distbegin, distend, distance, Season, 
                        SurveyType, FlySwim, Sample.Label, WatchID, Alpha, Dataset,
                        LatStart, LongStart, 
                        all_of(all.vars(df.model$final.formula))) %>% 
          check.distdata.cols()
            
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
          pluck(1)
        
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

# Exploratory data analysis
# Produce plots distance vs various covariates
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
    walk(
      vars,
      plot.covar,
      distdata = distdata,
      species = species,
      df.spec = df.spec,
      suffix = suffix
    )  
}


# Make a plot of a covar vs distance. Note the use of get(covar) to handle non-standard
# evaluation and allows us to pass a column name in covar.
plot.covar <-
  function(covar = NULL,
           distdata = NULL,
           species = NULL,
           df.spec = NULL,
           suffix = "") # to identify dataset, platform_class, and behaviour
{
  title <- paste(species, suffix, sep = " ")
  if (covar == "size") {
    p <- ggplot(data = as.data.frame(distdata), aes(cut(distance, df.spec$cutpoints, right = FALSE), size))
    p <- p + geom_boxplot(varwidth = TRUE)
    p <- p + labs(x = "Distance Category", y = "Size", title = title)
    print(p)
  } else if (is.factor(distdata[, covar][[1]])) {
    p <- ggplot(data = as.data.frame(distdata), aes(get(covar), distance))
    p <- p + geom_violin(draw_quantiles = c(.25, .5, .75), scale = "count")
    p <- p + scale_y_continuous(labels = as.character(df.spec$distance.centers), 
                                breaks = df.spec$distance.centers)
    p <- p + labs(y = "Distance Category", x = covar, title = title)
    print(p)
  } else {
    p <- ggplot(data = as.data.frame(distdata), aes(get(covar), distance))
    p <- p + geom_point()
    p <- p + geom_smooth(method = "lm")
    p <-
      p + scale_y_continuous(
        labels = as.character(df.spec$distance.centers),
        breaks = df.spec$distance.centers
      )
    p <- p + labs(y = "Distance Category", x = covar, title = title)
    print(p)
  }
}

# The next 3 functions are used to convert a dataframe that contains 
# coords of endpoints of a line (ie transect or watch) into an sf line object
#
# Make a linestring from from a matrix of coords for 2 points
make.line <- function(xy2){
  st_linestring(matrix(xy2, nrow=2, byrow=TRUE))
}

# convert dataframe of coords to lines
make.lines <- function(df, names=c("LongStart","LatStart","LongEnd","LatEnd"), crs){
  m = as.matrix(df[,names])
  lines = apply(m, 1, make.line, simplify=FALSE)
  st_sfc(lines, crs = crs)
}

# Convert a dataframe to lines. Used in Extract_data.Rmd
sf.pts.to.lines <- function(df, names=c("LongStart","LatStart","LongEnd","LatEnd"), crs){
  geom = make.lines(df, names, crs)
  df = st_sf(df, geometry=geom)
  df
}

# Assign season based on date ignoring year.  Works across new year boundary
#
# Inputs:
#   dat - object containing data to be assigned. 
#   season.def - definition of the seasons.
#   datefield - Character string giving the name of the date field in dat.
#   
# Outputs:
#   Input dataframe with a Season field added (if needed) set to  correct season.
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
  dates <- pull(dat, datefield)
  dat$monthday <- month(dates) * 100 + day(dates)
  
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
        between(dat$monthday, season.def[[i]]["from"], season.def[[i]]["to"])
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


# Convert a dataframe to sf and save as shapefile. Called from Extract_data.Rmd
df.to.shapefile <- function(df,
                            coords = c("LongStart", "LatStart"),
                            crs = st_crs(4326),
                            out.proj,
                            dsn = ShapeDir,
                            layer)
{
  st_as_sf(df,
           coords = coords,
           crs = crs,
           remove = FALSE) %>%
    st_transform(out.proj) %>%
    st_write(
      dsn = dsn,
      layer = layer,
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )
}

# Get distances and SD between successive GPS points along a watch. Called from
# Extract_data.Rmd.
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

# get the sum of lengths between GPS positions for a given watch. Return value
# is in km. Called from Extract_data.Rmd.
get.gps.length <- function(watch, posns) {
  # get positions in this watch  
  posns <- filter(posns,
                  between(posns$datetime, watch$WatchStartTime, watch$WatchEndTime))
  
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


# make a raster from an sf object element for a given season. Called from
# dsm.pred()
make.season.raster <- function(season, obj, variable) {
  v <- filter(obj, Season == season) %>% 
    # Convert sf to SpatVector
    vect()
  
  # Create a raster template with the same extent and resolution. Assumes
  # predgridCellLength is in km and raster projection units are metres.
  r <- rast(v, resolution = predgridCellLength * 1000)
    
  # Rasterize, using an attribute field (e.g., "ID")
  ret <- rasterize(v, r, field = variable)
    
  ret
}

# copy objects from one env to another. Useful after running a background job
# whose results are returned in an environment that you save in an .rda file.
# Reloading the .rda into a clean session will give one object (the environment)
# in the global env. So use:
#
# copy.env(name_of_loaded_env, .GlobalEnv)
copy.env <- function(src, dst) {
  for(n in ls(src, all.names=TRUE)) 
    assign(n, get(n, src), dst)
}

## Create initial segment data from watches
# inproj - proj4string or EPSG of input watches
# outproj - proj4string or EPSG of output segdata
# depth and depht.g rasters are created by predgrid creation and are passed in.
#
# Note internally all extraction of environmental values from rasters uses
# segProj and all saved rasters are in segProj. So, watch data coords will need
# to be in segProj to do the extractino
# 
# Not using rxtractogon any more. Instead, netcdf files have been downloaded by hand 
# and rasters extracted from them. Those rasters are utilized here.
create.segdata <- function(the.data,
                           study.area,
                           inproj,
                           outproj,
                           scale.factors,
                           verbose = FALSE) {
  
  # create initial segdata and reproject
  if (verbose) message("Creating initial segdata from watches")
  segdata <- the.data$watches %>% 
    rename(Effort = WatchLenKm) %>%
    mutate(TransectID = as.character(TransectID),
           year = lubridate::year(Date),
           yday = lubridate::yday(Date),
           MonthYear = format(Date, "%Y-%m"),
           segment.area = Effort * TotalWidthKm) %>% 
    st_as_sf(coords = c("LongStart", "LatStart"), crs = inproj, remove = FALSE) %>% 
    st_transform(outproj) %>% 
    cbind(st_coordinates(.)) %>% 
    rename(x = X, y = Y)
  
  ###---------------------------------------------------------------------------
  #### Load rasters 
  if (verbose) message("Loading environmental rasters")
  
  # get the dates of monthly rasters needed (ie. sst, sst.g, etc) so we can read
  # the needed files into a big SpatRaster
  dates.needed <- segdata$Date %>% 
    as.character %>% 
    str_sub(end = -4) %>% 
    paste0("-16") %>% 
    unique %>% 
    str_sort
  
  ### Depth and other static rasters
  # Depth
  depth <- rast(file.path(predLayerStudyAreaDir, "depth.img"))
  
  # Depth gradient
  depth.g <- rast(file.path(predLayerStudyAreaDir, "depth.g.img"))
  
  # SST  
  files <- file.path(predLayerStudyAreaDir, "sst", paste0("sst.", dates.needed, ".img"))
  sst <- rast(files)
  
  # SST gradient
  files <- file.path(predLayerStudyAreaDir, "sst", paste0("sst.g.", dates.needed, ".img"))
  sst.g <- rast(files)
  
  ###---------------------------------------------------------------------------
  ### Extract raster values at segdata locations  

  ### Extract depth values
  if (verbose) message("Extracting depth at watch locations")
  segdata <-
    terra::extract(
      depth,
      vect(segdata),
      bind = TRUE
    ) %>%  
    st_as_sf 
  
  ### Extract depth gradient values
  if (verbose) message("Extracting depth gradient at watch locations")
  segdata <-
    terra::extract(
      depth.g,
      vect(segdata),
      bind = TRUE
    ) %>%  
    st_as_sf 
  
  ### Extract SST values
  if (verbose) message("Extracting sst at watch locations")
  needed.layers <- paste0("sst.", segdata$MonthYear, "-16")
  segdata <-
    terra::extract(
      sst,
      vect(segdata),
      layer = needed.layers,
      bind = TRUE
    ) %>% 
    st_as_sf %>% 
    rename(sst = value) %>% 
    select(-layer) # Note that layer is off by one even though the sst values is correct

  ### Extract SST gradient values
  if (verbose) message("Extracting sst gradient at watch locations")
  needed.layers <- paste0("sst.g.", segdata$MonthYear, "-16")
  segdata <-
    terra::extract(
      sst.g,
      vect(segdata),
      layer = needed.layers,
      bind = TRUE
    ) %>% 
    st_as_sf %>% 
    rename(sst.g = value) %>% 
    select(-layer) # Note that layer is off by one even though the sst values is correct
  
  ###---------------------------------------------------------------------------
  ## Add scaled versions of all preds.
  if (verbose) message("Scaling covars")
  
  segdata %<>%
    mutate(depth.sc = (depth - scale.factors$depth_mean)/scale.factors$depth_sd,
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

  st_write(
    segdata,
    dsn = ShapeDir,
    layer = "segdata.shp",
    driver = "ESRI Shapefile",
    delete_layer = TRUE
  )
    
  if (verbose) print("Done.")
  segdata
}

# NOTE: no longer used due to issues with timeouts/lags/delays/etc. 
# but kept for posterity in case I return to this approach.
#
#
# Extract remote sensing data given by dataset and parameter in a polygon given
# by xcoord, ycoord (decimal degrees) for date given by tcoord and convert to a
# raster with projection prj.
#
# note: if tcoord is a length 2 vector, rxtractogon will return a whole series
# of data from tcoord[1] to tcoord[2], however that is not the intended use
# here.
#
# NOTE: to get the right monthly raster, tcoord must be of form "year-month-15"
# (or so), b/c rxtractogon() will choose the monthly dataset with the
# closest data so specifying "year-month" or "year-month-01" may get the previous
# month, which will cause no end of tears when trying to debug....
#
#
# this is designed to be used as: map(tcoords, rxtractogon.rast, blah, blah,
# blah, ...) to extract one raster per tcoord
# #
# Params:
# tcoord - time coordinate for dataset to download. Not all datasets (e.g. ETOPO depth)
#     have a tcoord in which case it should be null.
#
# dataset - name of ERDDAP dataset to access (e.g. jplMURSST41mday)
#
# parameter - name of paramter to get from dataset (e.g. sst)
#
# xcoord, ycoord - outline of area required (it will extract the bounding box)
#
# plotit - plot downloaded data
#
# saveit - save downloaded data to RasterDir 
# 
# folder - 
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
    rxtractogon(
      rerddap::info(dataset),
      parameter =  parameter,
      xcoord = xcoord,
      ycoord = ycoord,
      tcoord = tcoord
    )
  layername <- names(dat)[1]
  dat <- pluck(dat, 1) # use 1 since it is always 1st element, but not always called same as value of parameter
  if (length(dim(dat)) > 2) # remove useless third dimension
    dat <- dat[,,1]
  
  rast <- dat %>%
    t %>%                   # rxtracto returns matrix in odd order with x and y transposed and south to north so: transpose
    .[nrow(.):1, ] %>%       # ... and reverse order of rows
    raster(
      xmn = min(xcoord),
      xmx = max(xcoord),
      ymn = min(ycoord),
      ymx = max(ycoord),
      crs = latlongproj
    )
  
  if (plotit)
    plot(rast, main = sprintf("%s from %s for %s", parameter, dataset, tcoord))
  
  if (saveit) {
    
    if (is.null(tcoord))
      filename <- file.path(folder, paste(layername, "img", sep = "."))
    else
      filename <- file.path(folder, paste(layername, tcoord, "img", sep = "."))
    
    message("Saving downloaded ERDDAP raster to ", filename)
    if (!dir.exists(folder))
      dir.create(folder, recursive = TRUE)
    writeRaster(rast, filename = filename, format = "HFA", overwrite = TRUE)
  }
  
  rast
}

# XXX No longer used??
# function stolen from http://r-sig-geo.2731867.n2.nabble.com/Run-focal-function-on-a-multi-layer-raster-td7589931.html
# to perform focal() on RasterStack
multi.focal <- function(x, w = matrix(1, nrow = 3, ncol = 3), ...) {
  
  if (is.character(x)) {
    x <- brick(x)
  }
  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    focal(x[[ind]], w = w, ...)
  }
  
  n <- seq(nlayers(x))
  list <- lapply(X = n, FUN = fun, x = x, w = w, ...)
  
  out <- stack(list)
  return(out)
}

# Reads all raster files in folder with names matching pattern and returns a raster
# stack.
# 
# Useful for debugging and to avoid downloading from ERDDAP again.
# 
# e.g. x <- recreate.sst.mnth.from.files(RasterDir, "sst.+img$")
recreate.sst.mnth.from.files <- function(folder, pattern){
  files <- list.files(folder, pattern = pattern, full.names = T)
  r <- map(files, raster) %>% 
    stack
  names(r) <- basename(files) %>% 
    str_replace(fixed(".img"), "")
  r
}

# Extract data from a 3d matrix of data retrieved from a netCDF file
# with ncvar_get() (arg. dat). Create a raster and name it according to the.date. 
# 
# Note: the data in the y-dimension need to be flipped since they are "upside down"
#    in netcdf data.
# 
# 
# dat - 3D matrix of values. Dimensions are (lon, lat, month) (may have lon and lat 
#   reversed) 
# the.date - the date (16 of the month for "sst")
# index -  index of which month to pick off 
# datname - the name of the layer (e.g. "sst") that dat represents. only used
#    to label the plot.
# x - x coords of dat
# y - y coords of dat
# inproj - projection of dat
# outproj - projection of return raster dat
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
      crs = CRS(inproj)
    ) %>%
    raster::flip(direction = "y") %>% 
    rast
  
  # is reprojection/resampling required
  if (!is.null(to)) {
    if (!is.null(outproj))
      warning("create.ncdf.rast: both 'outproj' and 'to' are provided, ignoring outproj",
              immediate. = TRUE)
    
    res <- terra::project(res, rast(to), threads = TRUE) %>% 
      terra::mask(study.area)
  } else if (!is.null(outproj))
    res <- terra::project(res, outproj, threads = TRUE) %>% 
      terra::mask(study.area)
  
  plot(res, main = paste(datname, the.date))
  res
}



# Extract data from one NetCDF file
# dataset - character string giving name of dataset to extract. (eg. "sst")
# XXX TODO: figure out how to replace this with stars::read_ncdf()
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
      ncdf4::ncvar_get(x, "time")  %>%  # returns seconds since 01/01/1970
      `/`(86400) %>%  # convert to days since 01/01/1970
      as_date(origin = lubridate::origin)
  else
    dates <- ""

  # get rasters for each date
  res <-
    map2(
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
    res <- rast(res)
    names(res) <- paste(dataset, dates, sep = ".")
  } else {
    res <- res[[1]]
    names(res) <- dataset
  }
  res
}

# Reads all raster files in folder with names matching pattern and returns a raster
# or rasterstack.Called from Extract_env_rasters.Rmd
# 
# folder - where to find the files
# pattern - filename pattern to match e.g. "^jplMURSST41mday.+nc$")
# variable - the name of the variable in the NetCDF to extract, e.g. "sst"
# inproj - projection of the input netcdf
# outproj - projection to reproject to
# to - a raster to be used to reproject/resample/match extent to. If both
#     "outproj" and "to" are provided, "outproj" is ignored with a warning.
import.netCDF <-
  function(folder, pattern, variable, inproj, outproj = NULL, to = NULL) {
    
  files <- list.files(folder, pattern = pattern, full.names = T)
  
  if (length(files) == 0) {
    stop("import.netCDF: no matching filenames")
  } else {
    # reading one or more
    res <-
      map(
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
      res <- rast(res)
    } else
      res <- res[[1]]
    
  } 
  res
}


# Function to take a SpatRaster of monthly rasters (all months), extract the monthly
# rasters of interest and take their mean. Usefule for creating, for example,
# mean jan, feb, etc monthly averages. Used by Extract_env_rasters.rmd.
rast.monthly.mean <- function(mnth, r){
  message("Getting monthly means for month ", mnth)
  r.mnths <- names(r) %>%
    str_split_fixed(fixed("."), n = Inf) %>% 
    extract(, 2) %>% 
    str_split_fixed(fixed("-"), n = Inf) %>%  
    extract(, 2) %>% 
    as.integer
  
  sel <- r[[which(r.mnths == mnth)]]
  mean(sel)
}

# Produce a dotchart of a dataframe column. Used by Generic_1.5_dsm_EDA.Rmd
do.dotchart <- function(varname, dat) {
  if (varname %in% names(dat))
    dotchart(dat[, varname],
             main = varname,
             xlab = "Values of variable",
             ylab = "Order of the data")
}

# Compute seasonal mean for dynamic variable "var" (and its scaled version) from
# monthly values stored in separate columns in dat.
# Called once for each season and dynamic variable by get.seas.mean()
# var - variable of interest (e.g., "sst")
# dat - normally predgrid containing columns of monthly values variable var
# start, end - start and end MONTHS to have values averaged across.
get.seas.mean.var <- function(var, dat, start, end) {
  var.names <- paste(var, start:end, sep = ".")
  var.names.sc <- paste0(var.names, "_sc")
  ret = list(dat %>%
               select(all_of(var.names)) %>%
               rowMeans,
             dat %>%
               select(all_of(var.names.sc)) %>%
               rowMeans)
  names(ret) <- c(var, paste0(var, "_sc"))
  ret
}


# Compute species-specific seasonal means for dynamic vars listed in dyn.vars
# mean of appropriate months for each dynamic variable (and their scaled
# verions) listed in dyn.vars. Called once for each season and species from 
# create.seasonal.predgrid()
# 
# seas - character string given season
# season.spec - species-specific season boundaries, 
# dat - normally predgrid containing columns of monthly values for various dynamic 
#     variable
# dyn.vars - char vector of names of dynamic variables (e.g. "sst", "sst.g")
get.seas.mean <- function(seas, season.spec, dat, dyn.vars) {
  start <- season.spec[[seas]]["from"] %/% 100
  end <- season.spec[[seas]]["to"] %/% 100
  dat <- filter(dat, as.character(Season) == seas) 
  new.cols <- map(dyn.vars, get.seas.mean.var, dat = dat, start = start, end = end)
  cbind(dat, new.cols)
}


# Save current prediction htmls to folder. Useful for saving and comparing to
# subsequent improved iterations.
# 
# XXX Perhaps not needed anymore now that i'm saving results htmls in git LFS?
save.prediction.htmls <- function(folder) {
  if (!dir.exists(folder))
    dir.create(folder, recursive = TRUE)
  
  names(spec.grps) %>% 
    map(function(species) {
      filename <- file.path(ResultsDir, species, paste0(species, "_3_prediction.html"))
      message(sprintf("Copying '%' in '%s'", filename, folder))
      file.copy(filename, folder)
    })
}


# Reclassify values in a terra raster, r, according to class.arg, project to
# layer, to, and save in filename. Used by Extract_env_rasters.rmd
# 
# Note: class.arg can be either a 3-col matrix (from, to, becomes), a 2-col 
#     matrix (is, becomes) or a 1-col matrix (or vector) which specifies cut points.
#     normally used to convert Inf to NA with class.arg = cbind(Inf, NA)
reclassify.project.save <- function(r, class.arg, to, filename){
  outdir <- unique(dirname(filename))
  
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE )
  
  r %>%
    classify(rcl = class.arg) %>%
    terra::project(y = to,
                   method = "bilinear",
                   threads = TRUE) %>%
    terra::mask(to) %>% 
    terra::writeRaster(filename = filename,
                       overwrite = TRUE)
}

# extract names for dynamic variables (and their gradients if used) from
# env_covars table
dynamic.env.covar.names <- function(){
  dyn_vars <- filter(env_covar_spec, var_type == "dynamic")
  grads <- dyn_vars$var_name[dyn_vars$do_gradient]
  if (length(grads) > 0)
    c(dyn_vars$var_name, paste0(grads, ".g"))
  else
    dyn_vars$var_name
}

# Create seasonal predgrid - used by Generic_3_prediction.rmd
# 
# We will need 1 copy of the basic predgrid for each season. Seasonal dynamic 
# predictors are formed as the average of monthly values within each season. Seasons
# are defined on a per-species basis.
# 
# Once this seasonal predgrid is created, we will again need to make multiple copies
# of it (just like segdata - 1 for Fly and 1 for Swim). 
# Update May 30, 2025 - there should be as many copies as there are levels of the
#  platform variable in the dsm model. This can be found in ddftype_to_platform().
#
create.seasonal.predgrid <- function(species, predgrid) {
  # Get species-specific season setting and create predgrid.all.seas with 4 seasons
  season.spec <- seasons[[species]]
  ret <-
    rbind(predgrid, predgrid, predgrid, predgrid) %>%
    mutate(Season = as.factor(rep(season.names, each = nrow(predgrid))))
  
  # Get seasonal means for dynamic variables
  match.re <- c("[0-9]$", "[0-9]_sc$")
  p.geom <- st_geometry(ret) # save geometry
  ret <- season.names %>%
    map_dfr(
      get.seas.mean,
      season.spec = season.spec,
      dat = st_drop_geometry(ret),
      dyn.vars = dynamic.env.covar.names()
    ) %>%
    # Remove monthly values from predgrid to make it smaller now that we have
    # seasonal means computed, add area, and make back into sf object.
    select(!matches(match.re)) %>%
    cbind(p.geom) %>% # add geometry back in
    st_sf
  
  # Rename scaled columns to match what was in the model specs. This is harmless
  # if they are already named correctly with .sc.
  nms <- gsub( "_sc", ".sc", names(ret), fixed = TRUE)
  names(ret) <- nms
  
  # Save the seasonal prediction grid for GIS mapping. When doing multiple
  # species, there will be a separate predgrid for each species since it's
  # possible for the season boundaries to be species specific.  
  st_write(ret, 
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
    list_rbind(names_to = "platform") %>% 
    st_sf

  
  # From multiddf paper code:
  # create an extra column to account for the variance propagation model
  # the variance propagation adds a random effect covariate named "XX"
  # which we can safely give the value 0 now the variance has been propagated.
  # pred$XX <- matrix(0, nrow=nrow(pred), ncol=3)
  ret
}

# Do prediction maps for all four seasons, for a single model, and species.
# Called from Generic_3_prediction.rmd 
# 
# dat - a dataframe with 3 sets of seasonal predicions: 1 set for flying, 
#     1 set for swimming, and 1 set for combined.
#     
# subs - which subset of predictions to plot: flying, swimming or combined
do.pred.maps <-
  function(dat,
           model,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           ...) {
    
    
  # Get data subset
  subs <- match.arg(subs)
  dat <- filter(dat, subset == subs)
  segdata <- model$data %>% 
    st_as_sf
  
  message(sprintf("%s, %s: Doing %s abundance prediction map for",
                  species, modname, subs))
  
    
  ret <- season.names %>% 
    map(do.pred.map, dat, segdata, modname, species, subs, ...)

  names(ret) <- season.names
  ret
}

# Produce a leaflet map of Density predictions from model modname for one 
# species, season, and subset (Combined, F, or W). Called from do.pred.maps()
# 
# dat - predgrid (sf polygons) with NHat and Dens for all seasons
# segdata - segdata for all seasons
# samp_n - number of polygons from dat to plot. If NA, plot all. Otherwise, draw a 
#    sample of samp_n from the rows of dat (after filtering by season)
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
    filter(Season == season) %>%
    st_transform(latlongproj) %>%
    mutate(lDens = case_when(Dens == 0 ~ 0,
                             TRUE ~ log(Dens))) %>% 
    select(Dens, geometry) %>% 
    ms_simplify()
  
  segdata <- segdata %>% 
    filter(Season == season) %>% 
    get.combined.segdata() %>% 
    st_transform(latlongproj) 
  
  # Plot only a sample of the polygons for efficiency? Typically used for
  # testing.
  if (!is.na(samp_n)) {
    index <- sample(1:nrow(dat), size = samp_n)
    dat <- dat[index,]
  }

  # Remove ridiculously large densities b/c they mess up the legend and swamp
  # everything else
  dat <- mutate(dat,
                Dens = case_when(Dens > MAX_DENS_VALUE ~ NA,
                                 TRUE ~ Dens))
  
  if ((n.na <- sum(is.na(dat$Dens))) > 0)
    message("Warning: ", n.na, " cells larger than ", MAX_DENS_VALUE, " were converted to NA")
  
  groups <- c("est abund", "Pred Dens")
  m <-
    leaflet(
      data = dat,
      options = leafletOptions(preferCanvas = TRUE)
    ) %>%
    # Options help to speed up rendering. NOTE - dont't use addProviderTiles
    # if you want to save the map and reload in a subsequent R session - it won't
    # work.
    addTiles(options = tileOptions(updateWhenZooming = FALSE,
                                 updateWhenIdle = FALSE)) %>% 
    addMapPane("density", zIndex = 410) %>%
    addMapPane("abund", zIndex = 420) %>%
    # Predicted density
    addPolygons(
      fillColor = ~ pal_pred(Dens),
      color = ~ pal_pred(Dens),
      fillOpacity = 1.0,
      opacity = 1.0,
      weight = 1,
      group = "Pred Dens",
      options = pathOptions(pane = "density")
    ) %>%
    # 0 Abund
    addCircles(
      lng = ~ LongStart,
      lat = ~ LatStart,
      stroke = FALSE,
      radius = rep(1000, times = nrow(filter(segdata, estAbund == 0))),
      color = "white",
      fillColor = "white",
      fillOpacity = 0.1,
      opacity = 0.9,
      group = "est abund",
      data = filter(segdata, estAbund == 0),
      options = pathOptions(pane = "abund")
    ) %>%
    # Est Abund
    addCircles(
      lng = ~ LongStart,
      lat = ~ LatStart,
      radius = ~ estAbund * (max.circ.radius/max(estAbund)),#  scale so largest is 100km
      color = "black",
      weight = 1,
      popup = ~ htmlEscape(paste0("Abund ", round(estAbund, 3), ", raw ", round(rawCount, 3))),
      group = "est abund",
      data = filter(segdata, estAbund != 0),
      options = pathOptions(pane = "abund")
    ) %>%
    addLegend(
      pal = pal_pred,
      opacity = 1,
      values = ~ Dens,
      # values = class_intervals$brks, # only for discrete color scale
      title = paste(species, season)
    ) %>%
    addMouseCoordinates() %>%
    addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = FALSE)) %>%
    addLayersControl(overlayGroups = groups,
                     options = layersControlOptions(collapsed = FALSE)) %>%
    hideGroup(c("est abund"))
  
  m
}

# Find the species grp that a species Alpha code belongs to. Used by
# Extract_data.Rmd.
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



# Take segdata that has copies for each platform and combine them by summing
# estimated densities, abundances, counts, etc.
#
# Used by do.pred.map() and create.species.shapefiles()
#
# segdata - assumed to be in format returned by create.dsm.data
#
get.combined.segdata <- function(segdata){
  # Note that the species-specific segdata has multiple copies for each original
  # segment - currently up to 4 per aerial and 4 per ship-based segments: one
  # for all combinations of [A]erial/[S]hip, W[ater]/F[lying], and
  # [L]ine/S[trip] transect with Sample.Label suffixes as specified by
  # ddftype_levels defined in analysis_settings.R, so we first need to chop these
  # suffixes to do the filter and then sum the estDen[sities].
  segdata <-
    select(
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
    mutate(Sample.Label = str_sub(Sample.Label, 1, nchar(Sample.Label) - 
                                    nchar(ddftype_levels[1]) - 1)) %>%
    arrange(Sample.Label)
  
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
    mutate(
      estDens = zoo::rollapply(segdata$estDens, ncopies, by = ncopies, sum),
      estAbund = zoo::rollapply(segdata$estAbund, ncopies, by = ncopies, sum),
      rawCount = zoo::rollapply(segdata$rawCount, ncopies, by = ncopies, sum)
    ) %>% 
    select(-platform) # No longer makes any sense since it will have value of first row in group
 
  res 
}

# Sum fly and water detection-corrected segment densities for a given species in
# each segment (optionally limited to only those sample labels in sample.labs)
# and save as a set of seasonal shapefiles in folder. Called from Generic_2_dsm.Rmd
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
      segdata <- filter(segdata, Sample.Label %in% sample.labs)
    
    # save as seasonal shapefiles if required
    for (seas in season.names) {
      layer <- paste(spec, seas, "segdata", sep = "_")
      
      # if shapefile doesn't exist or recreateSpecSegShapefiles is TRUE then
      # save the shapefile
      if (!file.exists(paste0(folder, "/", layer, ".shp")) ||
          recreateSpecSegShapefiles) {
        message(sprintf("Creating segdata shapefile for %s %s", spec, seas))
        
        filter(segdata, Season == seas) %>%
          st_write(
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
        mutate(object = as.character(object)) %>%
        st_as_sf(
          coords = c("LongStart", "LatStart"),
          crs = st_crs(4326),
          remove = FALSE
        ) %>%
        st_write(
          dsn = ShapeDir,
          layer = layer,
          driver = "ESRI Shapefile",
          delete_layer = TRUE
        )
    }
  }

# Return names of prediction raster file 
# for a single species/season pair. Used by 03b_Save_chosen_model_predictions.Rmd
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

# Return names of variance raster file 
# for a single species/season pair. Used by 04.703b_Save_chosen_model_variance.Rmd
# Should be combined with previous function.
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

#
# Copy final model predictions html summary to a subfolder (ie
# basename(predVersionDir)) of the species-specific "Prediction summary" folder.
# Used by 03b_Save_chosen_model_predictions.Rmd. 
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

# undebug all debugged functions. This is tricky and sometimes doesn't work.
undebug.all <- function(where=search()) {
  aa <- all_debugged(where)
  lapply(aa$env,undebug)
  ## now debug namespaces
  invisible(mapply(function(ns,fun) {
    undebug(getFromNamespace(fun,ns))
  },names(aa$ns),aa$ns))
}


# Create  2x2 table of maps - one panel for each season.
# XXX Not currently used
create.annual.map.grid <- function(maps) {
  if(length(maps) != 4)
    stop(sprintf("create_annual_map_grids: maps argument does contains %d maps - should be 4."),
         length(maps))
  
  res <-
    tagList(tags$table(
      style = "width:100%",
      tags$tr(tags$td(tagList(maps$Spring)),
              tags$td(tagList(maps$Summer))),
      tags$tr(tags$td(tagList(maps$Fall)),
              tags$td(tagList(maps$Winter)))
    ))
  
  res
}

# Save a list of maps  (typically 4 seasonal prediction maps) from model
# modname to a summary folder. Used by Generic_3_prediction.Rmd.
save.map <- function(maps, modname, species) {
  dirname <- here(ResultsDir, species, "Prediction summaries")
  if (!dir.exists(dirname))
    dir.create(dirname, recursive = TRUE)
  timestr <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- here(dirname, paste0(modname, "_", timestr, ".html"))
  message(sprintf("%s, %s: Saving map in %s.", species, modname, filename))
  list(h2(paste0(species, "_", modname, "_", timestr)),
       leafsync::sync(maps)) %>%
    tagList %>%
    save_html(file = filename)
}

# Make quick and dirty prediction maps with ggplot. Called from
# Generic_3_prediction.Rmd
do.pred.maps.ggplot <-
  function(dat,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           ...) {
    
    
    # Get data subset
    subs <- match.arg(subs)
    dat <- filter(dat, subset == subs)
    
    message(sprintf("%s, %s: Doing %s abundance prediction map for",
                    species, modname, subs))
    
    ret <- season.names %>% 
      map(do.pred.map.ggplot, dat, modname, species, subs, ...)
    
    ret <- wrap_plots(ret) + plot_annotation(title = modname)
    ret
  }

# Plot a predicton ggplot for one model with 4 seasons. Called from
# do.pred.maps.ggplot()
do.pred.map.ggplot <-
  function(season,
           dat,
           modname,
           species,
           subs = c("Combined", "F", "W"),
           samp_n = NA) {
    
    message(sprintf("\t%s",season))
    dat <- dat %>%
      filter(Season == season) %>%
      select(Dens, geometry) %>% 
      ms_simplify()
    
    
    # Plot only a sample of the polygons for efficiency? Typically used for
    # testing.
    if (!is.na(samp_n)) {
      index <- sample(1:nrow(dat), size = samp_n)
      dat <- dat[index, ]
    }
    
    # Remove ridiculously large densities b/c they mess up the legend and swamp
    # everything else
    dat <- mutate(dat,
                  Dens = case_when(Dens > MAX_DENS_VALUE ~ NA,
                                   TRUE ~ Dens))

    ret <- ggplot(dat = dat) +
      geom_sf(aes(fill = Dens), color = NA)  +
      scale_fill_continuous(
        type = "viridis"
        # breaks = class_intervals$brks,
        # labels = round(class_intervals$brks, 4)
      ) +
      ggtitle(season)
    ret
  }


# Quick and dirty leaflet map for watches. Used when doing debugging of watch
# data in Extract_data.Rmd
watch.map <- function(dat) {
  leaflet(dat) %>%
    addTiles() %>%
    addCircles(
      lng = ~ LongStart,
      lat = ~ LatStart,
      radius = 1,
      popup = ~ htmlEscape(paste0("WatchID ", WatchID))
    )
}

# Recursively create a pathname if it doesn't exist. Return pathname invisibly
# so this function can be used in a pipeline.
create.dir.if.needed <- function(pathname){
  if (!dir.exists(pathname)){
    dir.create(pathname, recursive = TRUE)
    message(sprintf("Creating needed folder %s", pathname))
  }
  return(invisible(pathname))
  
}


# Create an initial generic ddf model list for each species group, optionaly
# writing it to dfModlistLoc (default false). Called from
# Create_final_ddf_model_specs.Rmd and Generic_1_ddf_fitting.Rmd
init.df.mod.list <- function(saveit = FALSE) {
  
  # Create one copy of the basic list for each species group
  df.mod.list <- rep(list(def.ddf.list), length(spec.grps))
  names(df.mod.list) <- names(spec.grps)
  
  # Set default convert_units value based on whether survey is aerial or ship.
  # This may get updated dynamically by later processing.
  df.mod.list %<>%
    map( ~ imap(., set.def.df.spec.values))
  
  # Sometimes we want to save it (ie if it didn't already exist) but other times
  # we're just called to return a initial structure that can be modified (eg in
  # 01.02_Create_final_ddf_model_specs.Rmd)
  if (saveit)
    save(df.mod.list, file = dfModlistLoc)
  
  df.mod.list
}

# Take a dataframe that contains at least DistMeth and FlySwim and assign
# the distance type (None, Perp., or Radial) by looking it up in lkpDistMeth
# in ECSAS database.
# 
# Used by Extract_data.Rmd and create.survey.data()
assign.dist.type <- function(dat) {
  distmeth <- ECSAS.get.table(ecsas.path = ECSAS.Path, "lkpDistMeth")
  DistType <- left_join(dat, distmeth, by = c("DistMeth" = "DistMethCode")) %>%
    mutate(DistType = as.factor(
      case_when(
        FlySwim == "F" ~ DistMethFlyHow,
        FlySwim == "W" ~ DistMethWaterHow,
        TRUE ~ NA_character_
      )
    )) %>%
    pull(DistType)
  
  if (any(is.na(DistType)))
    warning("assign.dist.type: ",
            sum(is.na(DistType)),
            " rows could not be assigned a distance type",
            immediate. = TRUE)
  
  DistType
}

# Generate the cononical model name given a key, formula and adj term.
# 
# Used by do.det.fcn().
# 
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
      str_replace(as.character(form)[2], fixed(" + "), ".")
  }
  # Deal with adjustment that might be NULL with ifelse()
  nam <- paste0(key, ifelse(is.null(form), "", paste0(".", form)))
  nam
}


## A few functions to make backing up and cleaning up DSM summaries easier

# Make backup copy of current DSM summaries for species labeling them with their
# modification date/time.
backup.dsm.summary <- function(species){
  folder <- here(ResultsDir, "Backups", species, "DSM Summaries")
  if (!dir.exists(folder))
    dir.create(folder, recursive = TRUE)

  message("Backing up ", species, " DSM summaries to ", folder)
    
  # src files
  src <- list.files(path = here(ResultsDir, species, "DSM Summaries"),
                      pattern = "*.Rdata", full.names = T)
  
  if (length(src) == 0){
    message("\tNo files to backup - quitting.")
    return
  }
  
  ## create dst filenames
  
  # add modification dates to files
  mtimes <- str_replace_all(file.info(src)$mtime, " ", "_") %>% 
    str_replace_all(":", "") %>% 
    str_replace("\\..*$", "") # remove trailing milliseconds
  stopifnot(length(src) == length(mtimes))
  dst <- src %>% 
    basename() %>% 
    tools::file_path_sans_ext() %>% 
    here(folder, .) %>%  
    paste0("_", mtimes, ".", tools::file_ext(src))
  
  file.copy(src, dst, overwrite = TRUE, copy.date = TRUE) %>% 
    invisible
}

# Do DSM summaries backup for all species
backup.dsm.summaries <- function(){
  names(spec.grps) %>% 
    walk(backup.dsm.summary)
}


# Remove all current DSM summaries for all species
remove.dsm.summaries <- function(){
  names(spec.grps) %>% 
    walk(\(species){
      message("Removing DSM summaries for ", species)
      files <- list.files(path = here(ResultsDir, species, "DSM Summaries"),
                 pattern = "*.Rdata", full.names = T)
      file.remove(files)
      })
}


# Get the dates of a dynamic env covariate given a variable name. Useful to see
# what files we need to download.
# Returns a dataframe with two columns: date, filename (typcially containing
# 12 rows - 1 for each month)
get.covar.netCDF.dates <- function(var.name){
  # list files in the folder here(predLayerDir, "NetCDF", var.name)
  # get dates associated with files, convert to char and return
  # 
  files <- list.files(here(predLayerDir, "NetCDF", var.name),
                      ".*\\.nc$",
                      full.names = TRUE)
  
  suppressWarnings(
    dates <- map_dfr(
      files,
      \(filenm) {
        data.frame(
          date = stars::read_ncdf(filenm, var = "time", proxy = TRUE) %>%
            stars::st_get_dimension_values("time") %>%
            str_replace(" UTC", ""),
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

# called once with each row of env_covars to get netcdf if needed, extract
# data from it, reproject and clip, etc.
get.env.covar <- function(env_covar_spec,
                          dates.needed,
                          verbose = TRUE
) {
  
  
  var_name <- env_covar_spec$var_name
  
  message("get.env.covar: getting covariate '", var_name, "'")
  
  # Dynamic or static covar
  if (env_covar_spec$var_type == "static"){
    env_dat <- stars::read_ncdf(here(
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
    files <- filter(netcdf.dates.have, date %in% dates.needed) %>% 
      pull(filename) %>% 
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
    if(is.na(st_crs(env_dat)))
      st_crs(env_dat) <- env_covar_spec$CRS
  } else
    stop("get.env.covar: variable ",
         var_name,
         ": illegal var_type '",
         env_covar_spec$var_type, "'")
  
  if (verbose)
    message("\tProjecting and clipping to study area...", appendLF = FALSE)
  
  # Re-project, clip, etc 
  final <- st_transform(env_dat, segProj) %>% 
    `[`(study.area) %>% 
    setNames(var_name) # b/c previous processing steps loose the name
  
  if (verbose)
    message("done.")
  
  final
}



# called once with each row of env_covars to get netcdf if needed, extract
# data from it, reproject and clip, etc.
# 
do.env.covar <- function(env_covar_spec,
                         segdata.dates.needed,
                         predgrid.dates.needed) {
  var_name <- env_covar_spec$var_name
  
  # Dynamic or static covar
  if (env_covar_spec$var_type == "static"){
    env.dat <- stars::read_ncdf(here(
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
    files <- filter(netcdf.dates.have, date %in% all.dates.needed) %>% 
      pull(filename) %>% 
      unique()
    
    # Read all needed netcdf files
    res <- map(files, \(filenm) {
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
  final <- st_transform(env.dat, segProj) %>% 
    `[`(study.area)

  final
}


# Create subproject folders when a new subproject is added
create_subproject_folders <- function(subproj) {
  create.dir.if.needed(file.path(GenericRDataDir, subproj)) 
  create.dir.if.needed(file.path(GenericShapeDir, subproj)) 
  create.dir.if.needed(here("Results", subproj))
  create.dir.if.needed(file.path(predLayerDir, "Study area resolution & extent", subproj)) 
  create.dir.if.needed(file.path(GISDir, "Predictions", SubProject)) 
  create.dir.if.needed(file.path(GISDir, "Rasters", SubProject)) 
  
}


# make a raster from an sf object
make.raster <- function(sfobj, variable){
  # Convert sf to SpatVector
  v <- vect(sfobj)
  
  # Create a raster template with the same extent and resolution. Assumes
  # predgridCellLength is in km and raster projection units are metres.
  r <- rast(v, resolution = predgridCellLength * 1000)
  
  # Rasterize, using an attribute field (e.g., "ID")
  ret <- rasterize(v, r, field = variable)
  
  ret
}

# Check if bam/gam model was fitted with discrete = TRUE
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


do_oneoff_render <- function(rmdfile, species) {
  source(here::here("R/analysis settings.r"), echo = T)
  
  rmarkdown::render(
    here(RDir, rmdfile),
    params = list(species = species),
    output_file = here(
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

# run a single ds model with specification given in dataframe mod
# 
#
# This Function is used by do.ds() to run many candidate models and by do.det.fcn 
# to re-run final model if needed.
# 
# Params:
# mod - model specification with 4 cols: label, key, adj, form
# data - the ovservation data
# folder - folder to store results. If null, then just return results
# logFileConn - destination for logging messsage. Will either be "" (i.e. stdout)
#   or an open file connection. Note, in the latter case it will not work when 
#   using parallel processing - logging output is lost.
# verbose - provide extra debugging info
# min_data_size_limit - minimum number of observations needed in order to attempt
#   fitting a ddf. If nrow(data) < min_data_size_limit model fitting is not attempted
#   and a failed model dummy file is created.
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
      str_replace_all(mod$form, " ", "") != "~1") {
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
        ch <- ddf.gof(model$ddf, qq = FALSE)$chisquare$chi1
        summ <- summary(model)
        cat(paste0("Date: ", date(), "\n"), file = modOutFileConn)
        cat(paste0("Model name: ", mod$label, "\n"), file = modOutFileConn)
        cat(paste0("Key: ", summ$ds$key, "\n"), file = modOutFileConn)
        cat(paste0("Formula: ", model$ddf$ds$aux$ddfobj$scale$formula),
            "\n",
            file = modOutFileConn)
        cat(paste0("AIC: ",  round(model$ddf$criterion, 3), "\n"), file = modOutFileConn)
        #ddf.gof(model$ddf, qq=FALSE)$dsgof$CvM$p,
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


# For a given key and vector of covars, generate all possible combinations of n
# items from covars.
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

# run a series of ds models using all combinations of covariates (including none).
#
# inputs:
#   data - same as ds
#   key - key function(s) to use
#   covars - a character vector of covariate names
#   runModels - for debugging. Should models actually be run or just return the dataframe with a list of what models would be run.
#       Does not do logging.
#   parallel - use parallel processing via doSNOW library to speed things up?
#   nCores - number of cores to use in the parallel cluster
#   folder - the name of a folder where summaries of each model will be written (one per file named
#       modDat$label). This is especially useful for monitoring progress when parallel==TRUE,
#       since the progress bar does not work. If summary files already exist in folder for some models,
#       then these will not be re-run, unless rerun == TRUE.
#   logfile - If provided all logging information will be sent to this file. Otherwise it goes to stdout, which will be lost if parallel == T.
#        Note that all message directly from do.ds() (but not run.ddf.model()) are also sent to stdout.
#   rerun - should models that have already been run and output saved in folder be run again.
#     Only applies when folder == T
#   verbose - enable extra output
#   models - dataframe indicating the models to be run. Possibly created from a
#     previous invocation of do.ds with runModels == F.
#   incl.adj - If TRUE(default), include key+adjustment only models in model set.
## Works by generating a dataframe of models to be run first and then running them.
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
    modDat <- map_dfr(key, function(k, covars) {
      if (k == "unif")
        gen.form.N(0, k, NULL)
      else
        map_dfr(0:length(covars), gen.form.N, k, covars)
    }, covars = covars)
  

  # append the adjustment only models
  if (incl.adj)
    modDat <- bind_rows(adj.models, modDat)

  # setup logfile connection. If it's "" (the default) then logfileConn will be
  # "", which corresponds to stdout.
  if (logfile != "") {
    logfileConn <- file(logfile, "at") # open for append.
    cat(paste0("Logging info to ", logfile, "\n"))
  } else
    logfileConn <- logfile

  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

  cat(paste0("\n\ndoDS: ", date(), "\n"), file = logfileConn)

  if (parallel & logfile != "" & runModels)
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
        cl <- makeCluster(min(nCores, nrow(modDat)), type = "SOCK")
        registerDoSNOW(cl)
        
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
        stopCluster(cl)
        
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
          map(
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


# stole this from
# https://github.com/gsk3/taRifx/blob/master/R/Rfunctions.R#L1161. I'm using it
# to turn NA's into "" so as.numeric wont issue warning about NA's introduced by
# coercion
destring <- function(x,keep="0-9.e+-") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep = ""),"",x)) )
}

# Check a ddf. Returns original distdata augmented with det. probs and adjSize
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
    cat(sprintf("Number of detection probs < 0.15: %d\n", nrow(filter(
      distdata, detProb < .15
    ))))
    cat(sprintf("range of size: %s\n", paste(range(distdata$size), collapse = " - ")))
    cat(sprintf("range of adjusted size: %s\n", paste(round(
      range(distdata$adjSize), 2
    ), collapse = " - ")))
    
    # plot hist of det probs
    p <- ggplot(data = distdata, aes(x = detProb))
    p <- p + geom_histogram(binwidth = 0.1)
    p <- p + scale_x_continuous(breaks = seq(0, 1, .1))
    print(p)
  }
  
  # GOF testing
  message("GOF testing")
  if (!("fake_ddf" %in% class(model))){
    print(ddf.gof(model$ddf, asp = 1))
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

# XXXX Might be useful but needs to be updated to use current ddftype column
# implementation instead of "det.fcn.type" left over from subArctic analysis
# 
# find segments with observations with more than one df type
get.trouble <- function(distdata, segdata) {
  distdata %>%
    group_by(Sample.Label) %>%
    summarize(ndet = n_distinct(det.fcn.type)) %>%
    filter(ndet > 1) %>%
    left_join(segdata, by = "Sample.Label") %>%
    left_join(distdata, by = "Sample.Label") %>%
    select(Sample.Label, ndet, det.fcn.type, df.type, distance, size, TransectType, Observer, object) %>%
    as.data.frame()
}


# Extract distance sampling data from a detection function model
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

# XXXX Could be useful but needs some arg changes - behav probably needs removing??
# and args to check.det.fcn need to be updated.
# 
# check top model for given species, behav, and model det fcn type (e.g.
# "normal"). If list.only is TRUE, just print out the filename of the top model.
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

# Check all top models for given species and behav
# XXXX Could be useful but needs some arg changes - 
# folder naming and layout has changed....
check.top.models <- function(species, behav, ...){
  folder <- here("R", species, "DF Summaries", behav)

  if (!dir.exists(folder)) {
    message(sprintf("Folder %s does not exist!", folder))
    return()
  }

  dirs <- list.dirs(folder, recursive = FALSE)
  walk(dirs, check.top.model, species, behav, ...)
}


# function to summarize a dsm model
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



# call function func (either dsm.var.gam or predict) and return either the
# predictions (for predict) or a list containing the pred.var and pred
# (dsm.var.gam). Is predict even necessary if I'm getting the prediction from
# dsm.var.gam anyway?
# Need to make off.set either a single value, vector of values, or character
# name of field in dat to get off.set from
apply.dsm.var <- function(dat, this.dsm){
  res <- dsm_var_gam(this.dsm, dat, map(dat, ".my.off.set"))
  
  list(pred.var = res$pred.var, pred = unlist(res$pred))
}

# break a predgrid into chunks and apply func in parallel (if nchunks > 1)
#
# this.dsm - the dsm_final
# df - the predgrid
# nchunks - number of chunks to break df up into. This applies even if not 
#    doing parallel. Useful b/c doing dsm_var_gam on entire large predgrid
#    fails trying to allocate many GB of RAM whereas doing it in chunks (even
#    if not in parallel) will succed.
# off.set = the offset for the gam (normally predgrid cell area). Can be:
#   - single value in which case it is used for all cells
#   - a vector the same length as nrow(df), in which case it will need to be 
#     split into chucks same as df 
# parallel = do chunks in parallel?
# nodes = number of cluster nodes for parallel processing

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
    # by either map() of parLapply() below.
    dat.split <- list(df)
  }

  # split each chunk in dat.split into sublists with 1 cell per element.
  dat.split <- map(dat.split, ~ split(.x, 1:nrow(.x)))

  if (parallel) {
    
    cl <- makeCluster(nodes)
    
    # Could just execute the things we need instead.
    # clusterEvalQ(cl, source(here::here("R/analysis settings.R")))
    clusterEvalQ(cl, {
      library(dsm)
      library(purrr)
    })
    
    # Need envir arg or else it won't find data objects when being rendered. 
    # takes about 1 min
    clusterExport(
      cl,
      c(
        "apply.dsm.var",
        "dat.split",
        "this.dsm"
      ),
      envir = environment()
    )

    # Run the function
    system.time(res <- parLapply(
      cl,
      dat.split,
      apply.dsm.var,
      this.dsm = this.dsm
      )
    )
    
    stopCluster(cl)
  } else {  # non-parallel version
    # apply the function to the chunks serially with map
    print(system.time(
      res <- map(dat.split, \(chunk) {
        apply.dsm.var(chunk, this.dsm)
        gc()
      }
      )
    ))
  }
  
  res
}


# XXX Not currently used
# Given a DSM and a dataframe to predict to, produce a density estimate with
# measures of uncertainty. Previously Used dsm.var.gam to get an ABUNDANCE
# estimate and then turn that into density estimate. Now just gets density
# estimate directly bu supplying 1 as the off.set
get.dens.est <- function(dsm_final, predgrid) {

  # use dsm.var.gam to get estimated abundance params
  densEst <- summary(dsm.var.gam(dsm_final, predgrid, off.set = predgridCellArea))

  # The estimates in dsm.var.gam are summed for the entire predgrid study area
  # so we need to divide by the number of predgrid cells * cellarea.
  densEst %<>%
    map_at(c("pred.est", "se"), ~ .x/(predgridCellArea * nrow(predgrid)))

  #calculate a lognormal CI for the density est.
  cv.square <- densEst$cv^2
  asymp.ci.c.term <- exp(1.96*sqrt(log(1 + cv.square))) # stolen from print.summary.dsm.var.R. Also in buckland et al 2001.
  list(pred.est = densEst$pred.est, CV = densEst$cv, SE = densEst$se, CI = c("%5" = densEst$pred.est/asymp.ci.c.term,
                                                                             "Mean" = densEst$pred.est, "95%" = densEst$pred.est * asymp.ci.c.term))
}

# XXX not currently used
# pretty-print the output from get.dens.est
print.dens.est <- function(densEst){
  cat("Density estimate:\n\n")

  cat("Approximate asymptotic confidence interval:\n")
  print(densEst$CI)

  cat("\n")
  cat("Point estimate                 :", densEst$pred.est,"\n")
  cat("Standard error                 :", densEst$SE,"\n")
  cat("Coefficient of Variation       :", densEst$CV,"\n")
}


# Render one of the Generic_x_xxx.Rmd files passing species as a param
# and save resulting knitted output.
# 
# Eg. do.generic.render("ATPU", "Generic_2_dsm.Rmd")
#
#
# We may be called from in parallel with other ongoing renders(), this could
# cause problems overwriting intermediate files, so we create a separate folder
# for each render.
do.generic.render <- function(species, file){
  
  suffix <- str_replace(file, "^Generic", "") %>%
    str_replace("Rmd$", "html")
  out.file <- file.path(ResultsDir, species, paste0(species, suffix))
  
  # Make sure output dir exists
  if (!dir.exists(dirname(out.file)))
    dir.create(dirname(out.file), recursive = TRUE)
  render_file <- file.path(RDir, file)

  message(sprintf("Rendering generic file %s for %s to %s", file, species, out.file))
  rmarkdown::render(
    render_file,
    params = list(species = species),
    intermediates_dir = tempdir(),
    output_file = out.file
  )

  return(species)
  
  callr::r(function(file, species, out.file) {
    rmarkdown::render(
      file,
      params = list(species = species),
      intermediates_dir = tempdir(),
      output_file = out.file
    )
  },
  args = list (
    file = render_file,
    species = species,
    out.file = out.file
  ))
}


# Render the extrapolation analysis Rmd for one dataset
do.extrapolation <- function(spill, dataset, debug = FALSE){
  browser(expr = debug)

  out.file <- file.path(here(), paste(spill, dataset, "0_extrapolation.html", sep = "_"))
  message(sprintf("Rendering generic extrapolation  for %s to %s", dataset, out.file))
  render(file.path(here(), "Generic_0_extrapolation.Rmd"),
                  params = list(spill = spill, dataset = dataset),
                  output_file = out.file)
}


# fitted dsm model checking
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
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), main = "s(x.sc,y.sc) (response scale)",
            type = "response", asp = 1, plot.type = "contour")
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = 0, phi = 45,
            main = "s(x.sc,y.sc) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = 60, phi = 45,
            main = "s(x.sc,y.sc) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x.sc","y.sc"), theta = -60, phi = 45,
            main = "s(x.sc,y.sc) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
  }
  
  # Shouldn't these be seasonal?
  if (any(grepl("s(x, y", as.character(dsm_final$formula), fixed = T))) {
    vis.gam(dsm_final,  view = c("x","y"), main = "s(x, y) (response scale)",
            type = "response", asp = 1, plot.type = "contour")
    vis.gam(dsm_final,  view = c("x","y"), theta = 0, phi = 45,
            main = "s(x,y) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x","y"), theta = 60, phi = 45,
            main = "s(x,y) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
    
    vis.gam(dsm_final,  view = c("x","y"), theta = -60, phi = 45,
            main = "s(x,y) (response scale)", type = "response",
            asp = 1, ticktype = "detailed")
  }
  
  
  # Sometimes whines about S3 methods.
  message("Gratia checks")
  suppressWarnings(gratia::appraise(dsm_final))
  
  # Gam checks from MGCV
  par(mfrow = c(1,1))
  message("MGCV checks")
  try(my.gam.check(dsm_final))
  message("rqgam_check():")
  rqgam_check(dsm_final)
  
  # Remove "dsm" class to make DHARMa happy
  message("DHARMa checks")
  simmod <- dsm_final
  class(simmod) <- class(simmod)[-1] # Simulate residuals doesn't like dsm class
  sims <- simulateResiduals(fittedModel = simmod)
  plot(sims)
  testResiduals(sims)
  testZeroInflation(sims)
  # May fail if x,y locations are not unique
  res <- try(testSpatialAutocorrelation(sims, segdata$x, segdata$y))
  
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
    walk(function(termlab, sims, segdata) {
      plotResiduals(sims, segdata[, termlab], xlab = termlab)
    }, sims = sims, segdata = segdata))

  if (!brief) {
    # Concurvity
    message("Concurvity checks\nEach term with whole of rest of model")
    try(print(concurvity(dsm_final) %>% round(digits = 3)))
    message(
      "Concurvity of pairwise terms ('estimate' measure presented)\nEach row shows how terms in columns depend on the term in that row."
    )
    try(print(concurvity(dsm_final, full = FALSE)[["estimate"]] %>% round(digits = 3)))
    message("Plot is non-symmetric, showing how terms on y-axis depend on terms on the x-axis")
    try(vis_concurvity(dsm_final))
  }
  
  par(mfrow = c(1,1))
  # check observed vs expected. See Miller et al 2021 pg 11 (of 18) for 
  # reccommendation to use the "platform" variable to aggregate
  # by.
  message("Observed vs expected plot")
  try(oe.dens(dsm_final, covar = "ddftype", plotit = T))
  try(oe.dens(dsm_final, covar = "depth", plotit = T))
  # oe.dens(dsm_final, covar = "depth.g", plotit = T)
  # oe.dens(dsm_final, covar = "sst", plotit = T)
  # oe.dens(dsm_final, covar = "sst.g", plotit = T)
  # # oe.dens(dsm_final, covar = "year", plotit = T)
  # # oe.dens(dsm_final, covar = "yday", plotit = T)

  # check zero infl - not sure this is right
  (OD_dsm_final <- sum(resid(dsm_final, type = "pearson"))/dsm_final$df.res)

  # Bubble plot
  message("Doing bubbleplot")
  mydata <-
    data.frame(
      E = resid(dsm_final, type = "pearson"),
      x = segdata$x / 1000,
      y = segdata$y / 1000
    )
  coordinates(mydata) <- ~ x + y
  
  print(bubble(
    mydata,
    "E",
    col = c("black", "red"),
    main = "Residuals",
    xlab = "X-coords",
    ylab = "y-coords"
  ))
  
  if (!brief) {
    # check autocorellogram
    # create transect label as cruiseid & date & flyswim in order to avoid
    # having the duplicated segments (one for fly one for swim) together in the
    # same "Transect".
    # create segment label as %H:%M:%S
    message("Doing autocorellogram")
    dsm_final$data <- dsm_final$data %>%
      mutate(
        tr.lab = paste(
          dsm_final$data$CruiseID,
          dsm_final$data$Date,
          dsm_final$data$FlySwim,
          sep = "_"
        ),
        seg.lab = format(as_datetime(dsm_final$data$StartTime), "%H:%M:%S")
      )
    par(mfrow = c(1, 1))
    dsm_cor(
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



# Augment segdata with summed counts of observations and create zeros for
# segments where there were no observations of the species of interest.
#
# Also add fields for rawCount, density-corrected estAbund, and estDens from a
# given set of observation distdata. These are used to check if my computed
# response is equal to that from dsm, after adjusting for one-sided transects
# with convert.units. This is used in run.dsm.model() to ensure that we are
# building segdata properly, and that convert.units is doing what we think it
# is.
#
# Note that if you want to restrict the analysis to a specific species/group,
# then distdata should already have been filtered for the species of interest
# before calling this function so that this code will fill in the zeros
# properly. This happens in Generic_2_dsm.Rmd in the usual case.
augment.segdata <- function(segdata, distdata) {
  newsegdata <- distdata %>%
    group_by(Sample.Label) %>%
    summarize(estAbund = sum(adjSize), rawCount = sum(size)) %>%
    right_join(segdata, by = "Sample.Label") %>% 
    st_as_sf
  newsegdata$estAbund[is.na(newsegdata$estAbund)] <- 0
  newsegdata$rawCount[is.na(newsegdata$rawCount)] <- 0
  newsegdata$estDens <- newsegdata$estAbund/newsegdata$segment.area
  newsegdata
}


# add 'k =  xxx' modifier to smooth terms for year or yday if there are less
# than k (default 10) unique values of that variable in segdata.
#
# if the resulting k would be < 2 then just remove the term
#
# XXX This should be made generic to go through all model terms. For each term
# just convert the formula to character and use string functions to grep and
# replace the term as needed instead of trying to use update().
#
adjust.time.covars <- function(form, segdata, k = 10) {
  # check if form contains year or yday
  # check unique number of values for each
  # if either is less than default 10 then

  # check for year
  new.k <- length(unique(segdata$year))
  if (grepl("s(year", as.character(form)[3], fixed = TRUE) &&
      (new.k < k)) {
    message("Removing year from formula for lack of data")
    form <- update(form, . ~ . - s(year, bs = "ts"))

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
    form <- update(form, . ~ . - s(yday, bs = "ts"))

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

# Run a dsm model.
#
# mod.def = dataframe with elements modname, formula, and family (the latter is a list columns).
# folder = path to folder to save results in.
# ddf.obj = is the list of ddfs in the correct order
# rerun.dsms = If TRUE (default) then rerun dsm() for mod.def, otherwise
# if .Rdata file exists from previous run of this model then load results from 
# that file.
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
  # source(here::here("R/analysis settings.R"))
  # logfileConn <- file(here::here(ResultsDir,  
  #                         paste0("dsm_logfile_", Sys.getpid(), ".txt")), "at")
  # 
  # needed if being called from future_map() on a worker process. If not,
  # it messes up trying to access any sf object (segment.data) when using 
  # the s2 spherical geometry package.
  sf::sf_use_s2(FALSE) 
  
  filename <- file.path(folder, paste0(mod.def$modname, ".Rdata"))
  # Rerun dsm?
  if (rerun.dsms == FALSE && file.exists(filename)) {
    # Load saved model result from file
    message(sprintf(
      "Loading previously saved model results for %s.",
      mod.def$modname
    ))
    load(filename)
  } else {
    # Rerun the dsm
    message("Running dsm model ", mod.def$modname)
    
    # Adjust formula
    form <- adjust.time.covars(mod.def$formula[[1]], segment.data)

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
      # add it back in
      kd <- control[["keepData"]]
      if (!is.null(kd) & isTRUE(kd) & !("data" %in% names(model))){
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
    } else { # Model failed to fit. Give message and then continue on to save.
      message("run.dsm.model: dsm() failed: ", model)
    } 
    
    message(sprintf("Saving model result to %s", filename))
    save(model, file = filename, compress = FALSE)
  } # End rerun dsm
  
  model
}


# Calculate predictions for a dsm model
#
# modname = name of model to extract from mod.res and predict for
# mod.res = list of dsm model objects
# species = name of species this is for
# predgrid = prediction grid containing one row per cell (polygons) in each of 4 
#     seasons. Must contain all predictors used in the model indicated by modname 
#     (excluding offset which is calculated internally by predict.)
#
# Initially predgrid contains 400240 rows, with one copy of the cells for 
# each combination of season and flyswim (The spatial extent of the study area 
# contains 50030 cells):
# 
#
#               F     W
#     Fall   50030 50030
#     Spring 50030 50030
#     Summer 50030 50030
#     Winter 50030 50030
#
# Value: predgrid augmented with NHat and Density
# 
# At the end ret will contain:
#   - predictions for fly (200120 rows - one for each spatial cell in each season)
#   - predictions for swim (200120 rows - one for each spatial cell in each season)
#   - predictions for combined (swim+fly) (200120 rows - one for each spatial cell in each season)
#   
# 
# NOTE NOTE NOTE if you do not supply newdata arg (ie. predgrid) but rather just
# predict to the same data that was used to fit the model then the "...order of
# the results will not necessarily be the same as the segdata (segment data)
# data.frame that was supplied (it will be sorted by the Segment.Label field)."
# (From predict.dsm manual)
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
  ret <- predgrid %>% 
    mutate(subset = platform)
  
  # Predict from model for all seasons and values of platform, which
  # may or may not be a factor in the given model. So, in a no_factor model
  # we will end up with n identical copies of predictions: 1 for each value
  # of platform. In a factor model these predictions for each level of platform
  # are different.
  ret <- ret %>% 
    mutate(NHat = predict(model, newdata=ret, off.set=ret$area),
           Dens = NHat/area)
  
  # Create summed (ie Combined Nhat) across all levels of platform.
  ret <- ret %>% 
    # Peel off one copy of all seasons predgrid template. ie nrow(ret) divided
    # by the number of platforms.
    head(nrow(.) / length(unique(.$platform))) %>% 
    # Sum NHats across platforms splitting the dataframe into a list with one
    # element per platform, extracting the NHat columns and summing with reduce
    mutate(subset = "Combined",
           NHat =  split(ret, ~ platform) %>% 
             map(~ .x$NHat) %>%
             reduce(`+`),
           Dens = NHat/area,
           platform = NA) %>% 
    rbind(ret)   # tack on original rows for each platform level
  
  
  
  # if (do.plots) {
  #   # Plot combined result
  #   combined_plot <- ggplot() +
  #     geom_sf(
  #       data = filter(ret, subset == "Combined"),
  #       mapping = aes(colour = Dens, fill = Dens)
  #     ) +
  #     labs(x = "", y = "", fill = "Dens") +
  #     ggtitle("Fly+Swim Combined") +
  #     theme_minimal() +
  #     scale_colour_viridis_c(option = "E") +
  #     scale_fill_viridis_c(option = "E")
  #   print(combined_plot)
  #   
  #   # Plot individual fly and swim results
  #   ind_plot <- ggplot() +
  #     geom_sf(
  #       data = filter(ret, subset != "Combined"),
  #       mapping = aes(colour = Dens, fill = Dens)
  #     ) +
  #     labs(x = "", y = "", fill = "Dens") +
  #     facet_wrap(vars(FlySwim)) +
  #     theme_minimal() +
  #     scale_colour_viridis_c(option = "E") +
  #     scale_fill_viridis_c(option = "E")
  #   print(ind_plot)
  # }  
  
  # Save as a shapefile
  ret %>% 
    filter(subset == "Combined") %>% 
    st_write(
      dsn = ShapeDir,
      layer = paste(species, modname, "predictions", sep = "_"),
      driver = "ESRI Shapefile",
      delete_layer = TRUE
    )
  

    
  # create 4-lyr seasonal raster of combined values for plotting and
  # saving
  message("Creating and saving seasonal rasters.")
  pp_raster <- season.names %>%
    map(make.season.raster,
        obj = filter(ret, subset == "Combined"),
        variable = "Dens") %>% 
    rast
  names(pp_raster) <- season.names
  
  # produce 4 panel plot. Not necessary since this is done at end of Generic_3_prediction.rmd
  # plot(pp_raster, main = paste(season.names, paste(species, modname, sep = "_")))
  
  # Save as a raster
  writeRaster(pp_raster,
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



# render ddf 
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
  
  # Make sure output dir exists
  if (!dir.exists(dirname(out.file)))
    dir.create(dirname(out.file), recursive = TRUE)
  
  render(
    file.path(here("R"), "Generic_1_ddf_fitting.rmd"),
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

# produce a summary table of candidate det. fcns in one folder
candidate.detfcn.summary <- function(df.spec.name, species) {
  
  folder <- file.path(ResultsDir, species, paste("DF Summaries", df.spec.name, sep = "_"))
  do.ds.det.fcn.checks(folder, species = species, dsetname = df.spec.name)
  
  get.ds.res(folder) %>%
    mutate(dataset = df.spec.name, species = species) %>%
    relocate(species, dataset, Model, Key, Formula, DetProb, AIC, deltaAIC, )
}


# Check all .txt files in folder for DS summary and return a dataframe with
# model name and summary items.
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

#Load each detfcn model in folder in turn and plot it
do.ds.det.fcn.checks <- function(folder, species, dsetname) {
  fl <- list.files(path = folder, pattern = "AIC.*\\.RData", full.names = T)
  
  walk(fl, function(filename) {
    message(sprintf("Checking detection function: %s", filename))
    load(filename)
    check.det.fcn(model, 
                  species = species, 
                  str_replace(basename(folder), fixed("DF Summaries_"), ""))
  })
}


# Get various values from a file containing a dsmodel object summary.
#folder <- "C:/Users/fifieldd/Documents/Offline/R/DS Utils/Test Output"
#path <- paste(folder, "hn.size 19-Jan-2016_145149.txt", sep = "/")
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


# Takes a ddf spec and adds to $fitted.distdata the following:
# 1) ddftype as a suffix to the Sample.Labels in fitted.distdata.
# 2) ddfobj.orig as an integer version of the ddftype factor.
# 
# Used by create.dsm.data to make Sample.Labels unique in each segdata copy, and
# to assign the initial ddfobj needed by dsm(). Note the final ddfobj value must
# be sequential with no missing values and will be assigned in create.dsm.data()
# once we know which ddfs are included in the dsm (ie. ones with no obs are not
# included).
augment.distdata <- function(df.mod.spec) {
  df.mod.spec$fitted.distdata <- df.mod.spec$fitted.distdata %>%
    mutate(Sample.Label = paste(Sample.Label, df.mod.spec$ddftype, sep = "_"),
           ddfobj.orig = as.integer(df.mod.spec$ddftype))
  df.mod.spec
}  


# Create a copy of the appropriate segdata (according to surveytype - either
# aerial or ship) for a df.mod.spec. Set the Sample.Label suffix and platform
# accordingly. If I decide to compact the ddftypes down to a smaller number of
# platforms then ddftype_to_platform allows for this. 
create.segdata.copy <- function(df.mod.spec, init_segdata) {
  
  # This all only works if df.mod.spec$ddftype is a single value (which it
  # should be), but I'm paranoid.
  stopifnot(length(unique(df.mod.spec$ddftype)) == 1)
  
  segdata <- switch(
    substr(df.mod.spec$ddftype, 1, 1),
    S = filter(init_segdata, SurveyType == "Ship"),
    A = filter(init_segdata, SurveyType == "Aerial")
  ) %>%
    mutate(
      Sample.Label = paste(Sample.Label, df.mod.spec$ddftype, sep = "_"),
      ddftype_orig = df.mod.spec$ddftype,
      # Lookup the platform type for this ddftype. May or may not have a 1:1
      # correspondence and ddftype_to_platform gives flexibility.
      platform = ddftype_to_platform[as.character(df.mod.spec$ddftype)]
    )
  segdata
}

# Nov 15 2024
#   - ddfobj may have less than the full number of values (currently 4) if some
#     ddfs have no obs
#   - nonetheless, there will be segdata copies for each unique ddftype, 6 in 
#     this case as below. so need to create a factor variable in segdata to 
#     hold this info.
#     
# May 2, 2025:
#   - noted difference bewteen number of levels of ddfobj (SW, SF, SS, AW, AF, AS)
#   which applies to distdata (6 - one for each ddf that we fitted including
#   dummies), and requires 6 copies of segdata, with platform which applies to
#   the dsm model term "platform" (4 in the case of the _factor models - one for
#   each level SW, SF, AW, SF).
#     
#
# Function to create distdata, segdata, and ddfs for a given species to pass to
# dsm().
#
# Segdata needs to contain (even if there are no observations in some of these
# classes):
#
# 1. ship segdata for water birds (SW) with perp distances - ddf.obj == 1 
# 2. ship segdata for flying birds (SF) with perp distances - ddf.obj == 2 
# 3. aerial segdata for water birds (AW) with perp distances - ddf.obj == 3 
# 4. aerial segdata for flying birds (AF) with perp distances - ddf.obj == 4 
# 5. water (ship or aerial) segdata for birds with NO perp distances - ddf.obj == 5 
# 6. flying (ship or aerial) segdata for birds with NO perp distances - ddf.obj == 6
#
# Instead of hardcoding the number and nature of these, I should probably do it
# dynamically based on df.mod.specs.
#
# We need to explicitly calculate the segment.area for each segment b/c aerial
# segs are (usually) two sided and ship ones are not. Normally, this is handled
# by the convert_units= arg to dsm, but it only takes a single scalar value
# whereas we need a different one for each type of segment (aerial vs ship).
# With this setup there is no need for a convert_units= arg to dsm anymore. Note
# that under normal circumstances it would be fine to just divide the ship-based
# effort by 2 and let dsm:::make.data() do the calc of segment.area but it isn't
# that simple when one of the ddfs has no observations. In that case,
# there is currently no way to make a dummy_ddf with no obs and so you have to
# arbitrarily assign a ddf.obj to the segments in that ddftype. I just chose
# ddfobj == 1 for this. If the transect geometry differs between these two then
# make.data() cannot properly calculate the segment.area. So I do it here by
# hand and pass the segment.area to dsm(). Note that this calculation can't be
# done in Create_seg_data.Rmd b/c we don't have the ddf info at that point.
# 
# Note that perhaps the "right" solution would be for dsm to work properly with
# to dummy_ddfs that have no observations.
#
# Segment area is:
#
# For two-sided transects: Effort * (width - left) * 2 
# - for one-sided transects: Effort * (width - left)
#
# Width is taken from the meta-data for the ddf.
#
# Arguments: 
# species - char string. Needed to chose the correct season for each
#  segment b/c seasons can be species specific. 
# df.mod.specs - list of detection function specifications - 
# init_segdata - the initial segdata computed by Create_seg_data.Rmd. This is
#  subdivided into aerial and ship.
#
# Value: 
#  Returns a list of 3 elements: 
#  
#  list(distdata, segdata, ddfs) 
#  
# distdata will be augmented with ddfobj, and have appropriate suffix added to
#   Sample.Label (see segdata below)
#
# segdata will consist of multiple copies of the appropriate (aerial or ship) 
#   segdata (one for each ddftype aka observation process) augmented with 
#   variables for:
#   segment.area: 
#   ddfobj - tells dsm:::make_data() which det prob to use to modify the
#     response, this is original ddftype re-numbered sequentially from 1 to number
#     of ddfs with data
#   ddfobj.orig - the original ddftype before renumbering. 1-6 as above table 
#   Season
# Sample.Label will have  _AW, _AF, _AS, _SW, _SF, or _SS added to it.
#
# ddfs - a list containing up to n detection functions (eg, Ship_fly, ship_swim,
# aerial_fly, aerial_swim, ship_strip, aerial_strip) excluding any for which
# there were no obs.
create.dsm.data <- function(species, df.mod.specs, init_segdata) {
  
  #### Create distdata ---------------------
    
  # Create the observation data that will be passed to dsm(). Note that some of
  # these distdata components specified by df.mod.specs may have no obs and will
  # thus not be included
  distdata <- df.mod.specs %>% 
    map(augment.distdata) %>% 
    # Extract list of fitted.distdatas
    map("fitted.distdata") %>% 
    # Remove datasets with no obs - note use of base R Filter()
    Filter(function(x) nrow(x) > 0, .) %>%
    map_dfr( ~ .) %>%  # Convert list to single dataframe
    mutate(
      # Need to renumber the ddfobj values sequentially from 1 to the number
      # of ddfs, whereas they may currently have gaps in the numeric sequence.
      # For example, if there were no obs in aerial_F_D (ddfobj.orig == 4) or 
      # aerial_F_N (ddfobj.orig == 8) then the
      # ddfobj.orig would be numbered 1, 2, 3, 5, 6, 7 but should be renumbered
      # as 1:6.
      ddfobj = dense_rank(ddfobj.orig)
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
    map(create.segdata.copy, init_segdata = init_segdata) %>% 
    map_dfr(~ .) %>% # Convert to one large dataframe
    assign.season(seasons[[species]], datefield = "Date") %>% 
    select(
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
      ends_with(".sc")
    ) %>%
    # Now that platform has been created, turn it into a factor, and assign
    # ddfobj based on conversion vector created above
    mutate(platform = as.factor(platform),
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
  
  # get aerial and ship initial segdata sizes and name them "A" and "S"
  sizes <- table(init_segdata$SurveyType) %>% 
    setNames(names(.) %>% substr(1,1))
  
  # Get survey type of each ddf, and then get the number of segemnts that had a
  # ddf with no obs that have used min(newddfs) as their ddfobj. 
  # This requires that the order of ddftype_levels and def.ddf.list have the 
  # same order.
  survey_type_index <- substr(ddftype_levels, 1,1)
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
  ddfs <- map(df.mod.specs, "fitted.model") %>% 
    magrittr::extract(unique(distdata$ddfobj.orig))

  ### TODO: need to deal with segments where WhatCount was something funky
  ### like only counting gannets. In that case, if we are doing a different species
  ### then these segments will be zeros but should not be included at all. For 
  ### now I think this NOGA exception is the only case. 
  
  list(distdata = distdata, segdata = segdata, ddfs = ddfs)
}


# Get a dataframe listing all ddf models that currently have no results files
# in folder. List of models is constructed from all candidate keys, adjustment 
# terms and covars. Useful for figuring out what ddf model fitting went off into
# la-la land...
list.unfinished.ddfs <- function(folder, covars = dfCovars){
  do.ds(
    data = NULL,
    folder = folder,
    runModels = FALSE,
    covars = covars
  )
}

# Creates a failed ddf model file in the given folder for model specified by mod.
# Useful when a given ddf never finishes running and needs to be marked as failed
# by hand.
create.failed.ddf.file <- function(mod, folder) {
    file.create(file.path(folder, paste0("failed_model_", mod$label, ".txt")))
}

