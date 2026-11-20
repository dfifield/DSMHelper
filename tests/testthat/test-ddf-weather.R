# Weather covariates for the ddf candidates (Atlantic-DSM #168):
# gen.form.N() groups of alternatives, n.candidate.formulas(),
# create.model.name() with three terms, and add.weather.covars().

test_that("a plain character vector gives every subset, as before groups", {
  f <- dplyr::bind_rows(lapply(0:2, gen.form.N, "hr", c("Season", "size")))
  expect_setequal(f$label, c("hr", "hr.Season", "hr.size", "hr.Season.size"))
  expect_equal(nrow(f), n.candidate.formulas(c("Season", "size")))
})

test_that("a group of alternatives never puts two of them in one formula", {
  cv <- list("Season", "size", c("Windspeed", "SeaStateGrp"))
  f <- dplyr::bind_rows(lapply(0:3, gen.form.N, "hr", cv))
  # 3 x 2 x 2 = 12: Season in/out, size in/out, weather none/Windspeed/SeaStateGrp
  expect_equal(nrow(f), 12L)
  expect_equal(n.candidate.formulas(cv), 12L)
  expect_false(any(grepl("Windspeed", f$label) & grepl("SeaStateGrp", f$label)))
  expect_true("hr.Season.size.Windspeed" %in% f$label)
  expect_true("hr.SeaStateGrp" %in% f$label)
  expect_equal(anyDuplicated(f$label), 0L)
  # the formula and the label agree, term for term
  expect_equal(f$form[f$label == "hr.Season.size.SeaStateGrp"], "~ Season + size + SeaStateGrp")
})

test_that("create.model.name() rebuilds a three-term label", {
  expect_equal(create.model.name("hr", ~ Season + size + Windspeed, NULL), "hr.Season.size.Windspeed")
  expect_equal(create.model.name("hn", ~ Season + size, NULL), "hn.Season.size")
  expect_equal(create.model.name("hn", ~ 1, "cos"), "hn.cos")
  expect_equal(create.model.name("hr", ~ 1, NULL), "hr")
})

wbs <- data.frame(SeaState = 0:9, wind_kt = c(6, 10, 12, 15.8, 20, 25, 30, 35, 35, 44))
bft <- data.frame(beaufort = 0:11, speed.kts = c(0, 2, 5, 8.5, 13.5, 19, 24.5, 30.5, 37, 44, 51.5, 59.5))

test_that("add.weather.covars() fills wind in order and groups sea state", {
  t0 <- as.POSIXct("2020-06-01 10:00", tz = "UTC")
  d <- data.frame(
    WatchID = paste0("w", 1:6),
    SeaState = c(1, NA, 4, NA, 2, NA),
    Windspeed = c(8, 18, NA, NA, NA, NA),
    Windforce = c(NA, NA, NA, 3, NA, NA),
    Date = as.Date("2020-06-01"), ObserverName = "A",
    ObsTime = t0 + c(0, 10, 20, 30, 40, 41) * 60)
  out <- add.weather.covars(d, wbs, 3, 15, bft)
  expect_equal(out$WindSource, c("recorded", "recorded", "seastate", "beaufort", "seastate", "neighbour"))
  expect_equal(out$Windspeed, c(8, 18, 20, 8.5, 12, 12))
  # row 2: no sea state, 18 kt >= 15 -> high; row 4: no sea state, 8.5 kt -> low;
  # row 6: neither -> copies row 5 (nearest in time, 1 min): low
  expect_equal(as.character(out$SeaStateGrp), c("0-2", "3+", "3+", "0-2", "0-2", "0-2"))
  expect_equal(levels(out$SeaStateGrp), c("0-2", "3+"))
  expect_equal(nrow(out), nrow(d))
})

test_that("add.weather.covars() stops, naming the watch, when nothing can fill a row", {
  d <- data.frame(WatchID = c("w1", "w2"), SeaState = c(NA, 2), Windspeed = c(NA, 9),
                  Date = as.Date(c("2020-06-01", "2020-06-02")), ObserverName = "A")
  expect_error(add.weather.covars(d, wbs), "w1")
})
