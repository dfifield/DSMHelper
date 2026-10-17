# The reuse test behind 02.05's family CV parts and 02.10's saved checks. Both
# used to reuse a saved result because the file existed, which is how 02.15
# once rolled up results from before the study-area buffer as current
# (Atlantic DSM issue #71).

# A result and its input, written in a chosen order. file.mtime() resolution on
# Windows/NTFS is fine-grained, but sleep past a second to be safe everywhere.
write.pair <- function(dir, input.first = TRUE) {
  input  <- file.path(dir, "suite.RData")
  result <- file.path(dir, "result.rds")
  first  <- if (input.first) input else result
  second <- if (input.first) result else input
  saveRDS(1, first)
  Sys.sleep(1.1)
  saveRDS(1, second)
  list(input = input, result = result)
}

prov <- list(mgcv = "1.9-3", seed = 20260820, thresholds = list(kindex = 0.9))

test_that("a current result is reused", {
  f <- write.pair(withr::local_tempdir())
  expect_identical(saved.result.stale.reason(f$result, f$input, prov, prov),
                   NA_character_)
})

test_that("a missing result is stale", {
  d <- withr::local_tempdir()
  expect_identical(
    saved.result.stale.reason(file.path(d, "none.rds"), file.path(d, "x"),
                              prov, prov),
    "no saved result")
})

test_that("a result saved without provenance is stale", {
  f <- write.pair(withr::local_tempdir())
  expect_match(saved.result.stale.reason(f$result, f$input, NULL, prov),
               "unrecorded")
})

test_that("differing settings are named, scalars with both values", {
  f <- write.pair(withr::local_tempdir())
  now <- prov
  now$mgcv <- "1.9-4"
  now$thresholds$kindex <- 0.8
  out <- saved.result.stale.reason(f$result, f$input, prov, now)
  expect_match(out, "mgcv 1.9-3 -> 1.9-4", fixed = TRUE)
  expect_match(out, "thresholds", fixed = TRUE)
  expect_no_match(out, "seed", fixed = TRUE)
})

test_that("a field present on only one side counts as a difference", {
  f <- write.pair(withr::local_tempdir())
  now <- c(prov, list(DHARMa = "0.5.0"))
  expect_match(saved.result.stale.reason(f$result, f$input, prov, now),
               "DHARMa unset -> 0.5.0", fixed = TRUE)
})

test_that("an input newer than the result makes it stale", {
  f <- write.pair(withr::local_tempdir(), input.first = FALSE)
  expect_match(saved.result.stale.reason(f$result, f$input, prov, prov),
               "suite.RData is .* newer")
})

test_that("the provenance records carry the installed versions", {
  cv <- get.family.cv.provenance(folds = 5, block.size.m = 1e5, nsim = 200,
                                 seed = 20260820, balance.column = "SurveyType")
  expect_identical(cv$mgcv, as.character(utils::packageVersion("mgcv")))
  expect_identical(cv$blockCV, as.character(utils::packageVersion("blockCV")))
  expect_identical(cv$balance_column, "SurveyType")

  ck <- get.dsm.check.provenance(finalists.only = TRUE, brief = FALSE,
                                 k.check.seed = 20260828,
                                 spatial.autocorr.max.n = 5000,
                                 spatial.autocorr.seed = 20260828)
  expect_identical(ck$DHARMa, as.character(utils::packageVersion("DHARMa")))
  expect_identical(ck$thresholds, default.check.thresholds())
})
