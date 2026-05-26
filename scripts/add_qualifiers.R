# add_qualifiers.R
# Adds package:: qualifiers to all unqualified external function calls in
# R/DSMHelper.R, using R's own parser to correctly identify call sites
# (avoiding strings, comments, and already-qualified calls).
#
# Run once from the package root:
#   Rscript scripts/add_qualifiers.R

src_file <- "R/DSMHelper.R"
lines <- readLines(src_file)

# ---------------------------------------------------------------------------
# Parse the file and collect all function-call tokens
# ---------------------------------------------------------------------------
parsed <- parse(src_file, keep.source = TRUE)
pd     <- getParseData(parsed)

# All tokens that appear in function-call position
calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL", ]
calls <- calls[order(calls$line1, calls$col1), ]

# ---------------------------------------------------------------------------
# Identify already-qualified calls (preceded by :: or :::)
# ---------------------------------------------------------------------------
pd_sorted <- pd[order(pd$line1, pd$col1), ]

is_qualified <- logical(nrow(calls))
for (i in seq_len(nrow(calls))) {
  idx <- which(pd_sorted$id == calls$id[i])
  if (idx > 1) {
    prev_tok <- pd_sorted$token[idx - 1]
    is_qualified[i] <- prev_tok %in% c("NS_GET", "NS_GET_INT")
  }
}

# ---------------------------------------------------------------------------
# Identify internally defined functions (must NOT receive a qualifier)
# ---------------------------------------------------------------------------
fn_def_pattern <- "^\\s*([a-zA-Z][a-zA-Z0-9._]*)\\s*(<-|=)\\s*function\\s*\\("
m <- regexpr(fn_def_pattern, lines, perl = TRUE)
matched_lines <- lines[m > 0]
internal_fns <- trimws(sub("\\s*(<-|=)\\s*function.*", "",
                            regmatches(matched_lines,
                                       regexpr(fn_def_pattern, matched_lines, perl = TRUE)),
                            perl = TRUE))
message(sprintf("Found %d internally defined functions", length(internal_fns)))

# ---------------------------------------------------------------------------
# Function -> package mapping (curated; priority = order declared here).
#
# INTENTIONALLY EXCLUDED (base R generics that dispatch via S3):
#   plot, summary, print, predict, residuals, simulate, fitted, coef,
#   vcov, logLik, AIC, BIC, format, as.data.frame, as.character, as.numeric,
#   as.integer, as.logical, as.factor, nrow, ncol, length, names, dim,
#   aggregate, subset, merge, match, which, table, apply, sapply, lapply,
#   tapply, do.call, Reduce, Filter, Map, Position, Find,
#   distance, res, origin, as.raster, as.matrix, as.vector, c, t,
#   min, max, sum, mean, var, sd, range, diff, cumsum, sort, order, rev,
#   unique, duplicated, is.na, is.null, is.numeric, is.character, is.logical
# ---------------------------------------------------------------------------
fn_to_pkg <- c(
  # ---- dplyr (over stats for filter / select / lag) ----
  mutate        = "dplyr",
  select        = "dplyr",
  filter        = "dplyr",
  arrange       = "dplyr",
  group_by      = "dplyr",
  ungroup       = "dplyr",
  summarise     = "dplyr",
  summarize     = "dplyr",
  pull          = "dplyr",
  rename        = "dplyr",
  distinct      = "dplyr",
  count         = "dplyr",
  left_join     = "dplyr",
  right_join    = "dplyr",
  inner_join    = "dplyr",
  full_join     = "dplyr",
  anti_join     = "dplyr",
  semi_join     = "dplyr",
  bind_rows     = "dplyr",
  bind_cols     = "dplyr",
  slice         = "dplyr",
  slice_max     = "dplyr",
  slice_min     = "dplyr",
  slice_head    = "dplyr",
  slice_tail    = "dplyr",
  slice_sample  = "dplyr",
  case_when     = "dplyr",
  if_else       = "dplyr",
  coalesce      = "dplyr",
  lag           = "dplyr",
  lead          = "dplyr",
  transmute     = "dplyr",
  across        = "dplyr",
  any_of        = "dplyr",
  all_of        = "dplyr",
  where         = "dplyr",
  everything    = "dplyr",
  starts_with   = "dplyr",
  ends_with     = "dplyr",
  contains      = "dplyr",
  matches       = "dplyr",
  n             = "dplyr",
  n_distinct    = "dplyr",
  between       = "dplyr",
  near          = "dplyr",
  tally         = "dplyr",
  add_tally     = "dplyr",
  add_count     = "dplyr",
  relocate      = "dplyr",
  rename_with   = "dplyr",
  rowwise       = "dplyr",
  c_across      = "dplyr",
  first         = "dplyr",
  last          = "dplyr",
  nth           = "dplyr",
  row_number    = "dplyr",

  # ---- sf ----
  st_as_sf              = "sf",
  st_transform          = "sf",
  st_crs                = "sf",
  st_set_crs            = "sf",
  st_write              = "sf",
  st_read               = "sf",
  st_drop_geometry      = "sf",
  st_geometry           = "sf",
  st_set_geometry       = "sf",
  st_area               = "sf",
  st_length             = "sf",
  st_distance           = "sf",
  st_buffer             = "sf",
  st_intersection       = "sf",
  st_union              = "sf",
  st_difference         = "sf",
  st_sym_difference     = "sf",
  st_bbox               = "sf",
  st_coordinates        = "sf",
  st_linestring         = "sf",
  st_multilinestring    = "sf",
  st_point              = "sf",
  st_multipoint         = "sf",
  st_polygon            = "sf",
  st_multipolygon       = "sf",
  st_geometrycollection = "sf",
  st_sfc                = "sf",
  st_sf                 = "sf",
  st_make_valid         = "sf",
  st_simplify           = "sf",
  st_crop               = "sf",
  st_intersects         = "sf",
  st_within             = "sf",
  st_contains           = "sf",
  st_is_valid           = "sf",
  st_centroid           = "sf",
  st_nearest_points     = "sf",
  st_cast               = "sf",
  st_as_text            = "sf",
  st_as_wkt             = "sf",
  st_as_wkb             = "sf",
  st_convex_hull        = "sf",
  st_snap               = "sf",
  st_segmentize         = "sf",
  st_sample             = "sf",
  st_make_grid          = "sf",
  st_graticule          = "sf",
  st_layers             = "sf",
  read_sf               = "sf",
  write_sf              = "sf",

  # ---- ggplot2 ----
  ggplot              = "ggplot2",
  aes                 = "ggplot2",
  geom_point          = "ggplot2",
  geom_line           = "ggplot2",
  geom_bar            = "ggplot2",
  geom_col            = "ggplot2",
  geom_histogram      = "ggplot2",
  geom_boxplot        = "ggplot2",
  geom_violin         = "ggplot2",
  geom_smooth         = "ggplot2",
  geom_raster         = "ggplot2",
  geom_tile           = "ggplot2",
  geom_contour        = "ggplot2",
  geom_contour_filled = "ggplot2",
  geom_ribbon         = "ggplot2",
  geom_abline         = "ggplot2",
  geom_hline          = "ggplot2",
  geom_vline          = "ggplot2",
  geom_text           = "ggplot2",
  geom_label          = "ggplot2",
  geom_segment        = "ggplot2",
  geom_path           = "ggplot2",
  geom_polygon        = "ggplot2",
  geom_sf             = "ggplot2",
  geom_sf_text        = "ggplot2",
  geom_errorbar       = "ggplot2",
  geom_crossbar       = "ggplot2",
  geom_density        = "ggplot2",
  geom_freqpoly       = "ggplot2",
  geom_step           = "ggplot2",
  geom_area           = "ggplot2",
  stat_summary        = "ggplot2",
  aes_string          = "ggplot2",
  vars                = "ggplot2",
  scale_x_continuous  = "ggplot2",
  scale_y_continuous  = "ggplot2",
  scale_x_discrete    = "ggplot2",
  scale_y_discrete    = "ggplot2",
  scale_x_log10       = "ggplot2",
  scale_y_log10       = "ggplot2",
  scale_color_manual  = "ggplot2",
  scale_colour_manual = "ggplot2",
  scale_fill_manual   = "ggplot2",
  scale_color_gradient   = "ggplot2",
  scale_colour_gradient  = "ggplot2",
  scale_fill_gradient    = "ggplot2",
  scale_color_gradient2  = "ggplot2",
  scale_colour_gradient2 = "ggplot2",
  scale_fill_gradient2   = "ggplot2",
  scale_color_gradientn  = "ggplot2",
  scale_colour_gradientn = "ggplot2",
  scale_fill_gradientn   = "ggplot2",
  scale_colour_viridis_c = "ggplot2",
  scale_color_viridis_c  = "ggplot2",
  scale_fill_viridis_c   = "ggplot2",
  scale_colour_viridis_d = "ggplot2",
  scale_color_viridis_d  = "ggplot2",
  scale_fill_viridis_d   = "ggplot2",
  scale_size_continuous  = "ggplot2",
  scale_alpha_continuous = "ggplot2",
  labs            = "ggplot2",
  xlab            = "ggplot2",
  ylab            = "ggplot2",
  ggtitle         = "ggplot2",
  theme           = "ggplot2",
  theme_bw        = "ggplot2",
  theme_minimal   = "ggplot2",
  theme_classic   = "ggplot2",
  theme_void      = "ggplot2",
  theme_grey      = "ggplot2",
  theme_gray      = "ggplot2",
  theme_set       = "ggplot2",
  theme_get       = "ggplot2",
  element_text    = "ggplot2",
  element_line    = "ggplot2",
  element_rect    = "ggplot2",
  element_blank   = "ggplot2",
  margin          = "ggplot2",
  unit            = "ggplot2",
  facet_wrap      = "ggplot2",
  facet_grid      = "ggplot2",
  coord_fixed     = "ggplot2",
  coord_sf        = "ggplot2",
  coord_flip      = "ggplot2",
  coord_cartesian = "ggplot2",
  coord_polar     = "ggplot2",
  coord_trans     = "ggplot2",
  ggsave          = "ggplot2",
  expansion       = "ggplot2",
  guide_colorbar  = "ggplot2",
  guide_legend    = "ggplot2",
  guides          = "ggplot2",
  after_stat      = "ggplot2",
  position_dodge  = "ggplot2",
  position_jitter = "ggplot2",
  position_stack  = "ggplot2",
  position_fill   = "ggplot2",
  ggplotGrob      = "ggplot2",

  # ---- terra (only unambiguous terra-specific names) ----
  # Excluded: nrow, ncol, aggregate, distance, res, origin, as.raster,
  #           plot, values, extract (all are base/S3 generics — keep via terra methods)
  vect        = "terra",
  rast        = "terra",
  rasterize   = "terra",
  project     = "terra",
  crop        = "terra",
  mask        = "terra",
  resample    = "terra",
  focal       = "terra",
  app         = "terra",
  lapp        = "terra",
  tapp        = "terra",
  ifel        = "terra",
  classify    = "terra",
  subs        = "terra",
  ext         = "terra",
  nlyr        = "terra",
  varnames    = "terra",
  sources     = "terra",
  writeRaster = "terra",
  as.polygons = "terra",
  as.points   = "terra",
  as.lines    = "terra",
  disagg      = "terra",
  terrain     = "terra",
  patches     = "terra",
  boundaries  = "terra",
  expanse     = "terra",

  # ---- purrr ----
  map           = "purrr",
  map_dbl       = "purrr",
  map_chr       = "purrr",
  map_int       = "purrr",
  map_lgl       = "purrr",
  map_df        = "purrr",
  map_dfr       = "purrr",
  map_dfc       = "purrr",
  map2          = "purrr",
  map2_dbl      = "purrr",
  map2_chr      = "purrr",
  map2_lgl      = "purrr",
  map2_int      = "purrr",
  map2_dfr      = "purrr",
  imap          = "purrr",
  imap_dbl      = "purrr",
  imap_chr      = "purrr",
  imap_dfr      = "purrr",
  iwalk         = "purrr",
  pmap          = "purrr",
  pmap_dbl      = "purrr",
  pmap_dfr      = "purrr",
  pmap_chr      = "purrr",
  walk          = "purrr",
  walk2         = "purrr",
  reduce        = "purrr",
  reduce2       = "purrr",
  accumulate    = "purrr",
  pluck         = "purrr",
  chuck         = "purrr",
  modify        = "purrr",
  modify_if     = "purrr",
  modify_at     = "purrr",
  keep          = "purrr",
  discard       = "purrr",
  compact       = "purrr",
  flatten       = "purrr",
  flatten_dbl   = "purrr",
  flatten_chr   = "purrr",
  flatten_int   = "purrr",
  flatten_lgl   = "purrr",
  possibly      = "purrr",
  safely        = "purrr",
  quietly       = "purrr",
  set_names     = "purrr",
  list_rbind    = "purrr",
  list_cbind    = "purrr",
  has_element   = "purrr",
  detect        = "purrr",
  detect_index  = "purrr",

  # ---- tidyr ----
  pivot_longer  = "tidyr",
  pivot_wider   = "tidyr",
  gather        = "tidyr",
  spread        = "tidyr",
  nest          = "tidyr",
  unnest        = "tidyr",
  unnest_wider  = "tidyr",
  unnest_longer = "tidyr",
  separate      = "tidyr",
  separate_rows = "tidyr",
  unite         = "tidyr",
  fill          = "tidyr",
  replace_na    = "tidyr",
  drop_na       = "tidyr",
  crossing      = "tidyr",
  expand_grid   = "tidyr",
  complete      = "tidyr",
  hoist         = "tidyr",

  # ---- lubridate ----
  year          = "lubridate",
  month         = "lubridate",
  day           = "lubridate",
  hour          = "lubridate",
  minute        = "lubridate",
  second        = "lubridate",
  yday          = "lubridate",
  ymd           = "lubridate",
  mdy           = "lubridate",
  dmy           = "lubridate",
  ymd_hms       = "lubridate",
  mdy_hms       = "lubridate",
  dmy_hms       = "lubridate",
  make_date     = "lubridate",
  make_datetime = "lubridate",
  floor_date    = "lubridate",
  ceiling_date  = "lubridate",
  round_date    = "lubridate",
  as_date       = "lubridate",
  as_datetime   = "lubridate",
  isoweek       = "lubridate",
  epiweek       = "lubridate",
  quarter       = "lubridate",
  semester      = "lubridate",
  days          = "lubridate",
  weeks         = "lubridate",
  hours         = "lubridate",
  minutes       = "lubridate",
  seconds       = "lubridate",
  dseconds      = "lubridate",
  dminutes      = "lubridate",
  dhours        = "lubridate",
  ddays         = "lubridate",
  dweeks        = "lubridate",
  dyears        = "lubridate",
  parse_date_time = "lubridate",
  with_tz       = "lubridate",
  force_tz      = "lubridate",
  is.Date       = "lubridate",

  # ---- stringr ----
  str_c           = "stringr",
  str_split       = "stringr",
  str_split_1     = "stringr",
  str_split_fixed = "stringr",
  str_detect      = "stringr",
  str_replace     = "stringr",
  str_replace_all = "stringr",
  str_extract     = "stringr",
  str_extract_all = "stringr",
  str_match       = "stringr",
  str_match_all   = "stringr",
  str_pad         = "stringr",
  str_trim        = "stringr",
  str_squish      = "stringr",
  str_glue        = "stringr",
  str_length      = "stringr",
  str_sub         = "stringr",
  str_to_lower    = "stringr",
  str_to_upper    = "stringr",
  str_to_title    = "stringr",
  str_to_sentence = "stringr",
  str_count       = "stringr",
  str_starts      = "stringr",
  str_ends        = "stringr",
  str_remove      = "stringr",
  str_remove_all  = "stringr",
  str_locate      = "stringr",
  str_locate_all  = "stringr",
  str_dup         = "stringr",
  str_flatten     = "stringr",
  str_sort        = "stringr",
  str_order       = "stringr",
  str_wrap        = "stringr",
  str_trunc       = "stringr",
  str_view        = "stringr",
  fixed           = "stringr",
  regex           = "stringr",
  coll            = "stringr",

  # ---- readr ----
  write_csv   = "readr",
  read_csv    = "readr",
  write_tsv   = "readr",
  read_tsv    = "readr",
  write_rds   = "readr",
  read_rds    = "readr",
  write_delim = "readr",
  read_delim  = "readr",
  write_lines = "readr",
  read_lines  = "readr",
  read_file   = "readr",

  # ---- here ----
  here = "here",

  # ---- rmarkdown ----
  render = "rmarkdown",

  # ---- checkmate ----
  check_class             = "checkmate",
  assert_class            = "checkmate",
  expect_class            = "checkmate",
  check_numeric           = "checkmate",
  assert_numeric          = "checkmate",
  expect_numeric          = "checkmate",
  check_character         = "checkmate",
  assert_character        = "checkmate",
  expect_character        = "checkmate",
  check_logical           = "checkmate",
  assert_logical          = "checkmate",
  check_list              = "checkmate",
  assert_list             = "checkmate",
  check_data_frame        = "checkmate",
  assert_data_frame       = "checkmate",
  check_number            = "checkmate",
  assert_number           = "checkmate",
  check_int               = "checkmate",
  assert_int              = "checkmate",
  check_flag              = "checkmate",
  assert_flag             = "checkmate",
  check_string            = "checkmate",
  assert_string           = "checkmate",
  check_function          = "checkmate",
  assert_function         = "checkmate",
  check_file_exists       = "checkmate",
  assert_file_exists      = "checkmate",
  check_directory_exists  = "checkmate",
  assert_directory_exists = "checkmate",
  check_names             = "checkmate",
  assert_names            = "checkmate",
  check_count             = "checkmate",
  assert_count            = "checkmate",
  check_null              = "checkmate",
  assert_null             = "checkmate",
  check_true              = "checkmate",
  assert_true             = "checkmate",
  check_false             = "checkmate",
  assert_false            = "checkmate",
  check_subset            = "checkmate",
  assert_subset           = "checkmate",
  check_choice            = "checkmate",
  assert_choice           = "checkmate",
  test_class              = "checkmate",
  test_numeric            = "checkmate",
  test_flag               = "checkmate",
  test_string             = "checkmate",
  test_list               = "checkmate",
  test_data_frame         = "checkmate",
  test_number             = "checkmate",
  test_int                = "checkmate",
  test_character          = "checkmate",
  test_logical            = "checkmate",
  makeAssertionFunction   = "checkmate",
  makeCheckFunction       = "checkmate",
  makeTestFunction        = "checkmate",
  makeExpectationFunction = "checkmate",

  # ---- DHARMa ----
  simulateResiduals          = "DHARMa",
  testDispersion             = "DHARMa",
  testUniformity             = "DHARMa",
  testSpatialAutocorrelation = "DHARMa",
  plotResiduals              = "DHARMa",
  plotQQunif                 = "DHARMa",
  recalculateResiduals       = "DHARMa",
  testOutliers               = "DHARMa",
  createDHARMa               = "DHARMa",
  testZeroInflation          = "DHARMa",
  testQuantiles              = "DHARMa",

  # ---- mgcv ----
  gam        = "mgcv",
  bam        = "mgcv",
  gamm       = "mgcv",
  s          = "mgcv",
  te         = "mgcv",
  ti         = "mgcv",
  t2         = "mgcv",
  gam.check  = "mgcv",
  concurvity = "mgcv",
  qq.gam     = "mgcv",
  vis.gam    = "mgcv",
  k.check    = "mgcv",
  in.out     = "mgcv",
  smoothCon  = "mgcv",
  PredictMat = "mgcv",

  # ---- gratia ----
  draw             = "gratia",
  derivatives      = "gratia",
  smooth_estimates = "gratia",
  appraise         = "gratia",
  basis            = "gratia",
  variance_comp    = "gratia",

  # ---- leaflet ----
  leaflet              = "leaflet",
  addTiles             = "leaflet",
  addProviderTiles     = "leaflet",
  addPolygons          = "leaflet",
  addCircles           = "leaflet",
  addCircleMarkers     = "leaflet",
  addMarkers           = "leaflet",
  addAwesomeMarkers    = "leaflet",
  addRasterImage       = "leaflet",
  addLegend            = "leaflet",
  addLayersControl     = "leaflet",
  addMiniMap           = "leaflet",
  addScaleBar          = "leaflet",
  setView              = "leaflet",
  fitBounds            = "leaflet",
  flyTo                = "leaflet",
  clearBounds          = "leaflet",
  colorNumeric         = "leaflet",
  colorBin             = "leaflet",
  colorQuantile        = "leaflet",
  colorFactor          = "leaflet",
  providers            = "leaflet",
  hideGroup            = "leaflet",
  showGroup            = "leaflet",
  labelFormat          = "leaflet",
  highlightOptions     = "leaflet",
  popupOptions         = "leaflet",
  markerOptions        = "leaflet",
  pathOptions          = "leaflet",
  tileOptions          = "leaflet",
  layersControlOptions = "leaflet",
  leafletOutput        = "leaflet",
  renderLeaflet        = "leaflet",
  leafletProxy         = "leaflet",
  removeShape          = "leaflet",
  clearShapes          = "leaflet",

  # ---- mapview ----
  mapview        = "mapview",
  mapshot        = "mapview",
  viewExtent     = "mapview",
  viewRGB        = "mapview",
  mapviewOptions = "mapview",

  # ---- leafsync ----
  sync = "leafsync",

  # ---- rmapshaper ----
  ms_clip     = "rmapshaper",
  ms_simplify = "rmapshaper",
  ms_dissolve = "rmapshaper",
  ms_union    = "rmapshaper",
  ms_erase    = "rmapshaper",
  ms_filter   = "rmapshaper",
  ms_explode  = "rmapshaper",
  ms_lines    = "rmapshaper",

  # ---- ncdf4 ----
  nc_open   = "ncdf4",
  nc_close  = "ncdf4",
  ncvar_get = "ncdf4",
  ncvar_put = "ncdf4",
  ncatt_get = "ncdf4",
  ncatt_put = "ncdf4",
  nc_create = "ncdf4",
  ncdim_def = "ncdf4",
  ncvar_def = "ncdf4",
  nc_sync   = "ncdf4",

  # ---- tidync ----
  tidync       = "tidync",
  activate     = "tidync",
  hyper_dims   = "tidync",
  hyper_vars   = "tidync",
  hyper_tibble = "tidync",
  hyper_array  = "tidync",
  hyper_filter = "tidync",

  # ---- stars ----
  read_ncdf   = "stars",
  read_stars  = "stars",
  write_stars = "stars",
  st_as_stars = "stars",

  # ---- geosphere ----
  distHaversine         = "geosphere",
  distVincentyEllipsoid = "geosphere",
  destPoint             = "geosphere",
  bearing               = "geosphere",
  gcIntermediate        = "geosphere",
  areaPolygon           = "geosphere",
  alongTrackDistance    = "geosphere",

  # ---- sp ----
  SpatialPoints            = "sp",
  SpatialPointsDataFrame   = "sp",
  SpatialLines             = "sp",
  SpatialLinesDataFrame    = "sp",
  SpatialPolygons          = "sp",
  SpatialPolygonsDataFrame = "sp",

  # ---- raster (only functions not covered by terra above) ----
  brick           = "raster",
  stack           = "raster",
  nlayers         = "raster",
  extent          = "raster",
  projectRaster   = "raster",
  projectExtent   = "raster",
  cellStats       = "raster",
  calc            = "raster",
  overlay         = "raster",
  addLayer        = "raster",
  dropLayer       = "raster",
  getValues       = "raster",
  setValues       = "raster",
  rasterToPoints  = "raster",
  rasterToPolygons = "raster",
  pointsToRaster  = "raster",

  # ---- viridis ----
  viridis     = "viridis",
  viridis_pal = "viridis",
  magma       = "viridis",
  plasma      = "viridis",
  inferno     = "viridis",
  cividis     = "viridis",
  turbo       = "viridis",

  # ---- furrr ----
  future_map     = "furrr",
  future_map2    = "furrr",
  future_pmap    = "furrr",
  future_walk    = "furrr",
  future_map_dbl = "furrr",
  future_map_chr = "furrr",
  future_map_int = "furrr",
  future_map_lgl = "furrr",
  future_map_dfr = "furrr",
  future_imap    = "furrr",

  # ---- future ----
  plan         = "future",
  multisession = "future",
  multicore    = "future",
  sequential   = "future",

  # ---- doSNOW ----
  registerDoSNOW = "doSNOW",

  # ---- foreach ----
  foreach = "foreach",

  # ---- htmltools ----
  tagList = "htmltools",
  HTML    = "htmltools",

  # ---- knitr ----
  kable       = "knitr",
  knit_print  = "knitr",

  # ---- zoo ----
  rollapply = "zoo",
  rollmean  = "zoo",

  # ---- gstat ----
  variogram     = "gstat",
  vgm           = "gstat",
  fit.variogram = "gstat",
  krige         = "gstat",
  idw           = "gstat",

  # ---- plyr (ddply etc.; arrange deliberately excluded — dplyr wins above) ----
  ddply = "plyr",
  ldply = "plyr",
  dlply = "plyr",
  llply = "plyr",
  adply = "plyr",
  daply = "plyr",
  laply = "plyr",
  aaply = "plyr",
  alply = "plyr",
  join  = "plyr",

  # ---- geodist ----
  geodist = "geodist",

  # ---- rerddapXtracto ----
  rxtractogon = "rerddapXtracto",
  rxtracto    = "rerddapXtracto",

  # ---- ECSASconnect (Suggests) ----
  ECSAS.find.suspicious.posn = "ECSASconnect"
)

# ---------------------------------------------------------------------------
# Determine which calls to qualify
# ---------------------------------------------------------------------------
to_qualify <- calls[
  !is_qualified &
  calls$text %in% names(fn_to_pkg) &
  !calls$text %in% internal_fns,
]
message(sprintf("Found %d call sites to qualify", nrow(to_qualify)))

changes <- table(to_qualify$text)
message("\nCall counts by function:")
print(sort(changes, decreasing = TRUE))

# ---------------------------------------------------------------------------
# Apply replacements bottom-up (last line first preserves column positions)
# ---------------------------------------------------------------------------
to_qualify <- to_qualify[order(-to_qualify$line1, -to_qualify$col1), ]

for (i in seq_len(nrow(to_qualify))) {
  fn  <- to_qualify$text[i]
  pkg <- fn_to_pkg[[fn]]
  ln  <- to_qualify$line1[i]
  cl  <- to_qualify$col1[i]

  line <- lines[ln]
  lines[ln] <- paste0(
    substr(line, 1, cl - 1),
    pkg, "::",
    substr(line, cl, nchar(line))
  )
}

writeLines(lines, src_file)
message(sprintf("\nDone. %d qualifiers added to %s", nrow(to_qualify), src_file))
