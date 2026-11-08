test_that("add.base.tiles() adds the Esri ocean tiles, not OpenStreetMap (#123)", {
  m <- add.base.tiles(leaflet::leaflet())
  calls <- Filter(function(x) x$method == "addTiles", m$x$calls)
  expect_length(calls, 1)
  url <- calls[[1]]$args[[1]]
  expect_match(url, "arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base", fixed = TRUE)
  expect_false(grepl("openstreetmap", url, fixed = TRUE))
  expect_equal(calls[[1]]$args[[4]]$maxNativeZoom, 13)
})

test_that("add.base.tiles() refuses a non-leaflet object", {
  expect_error(add.base.tiles(list()))
})
