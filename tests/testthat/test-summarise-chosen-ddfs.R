# summarise.chosen.ddfs() reads the candidate summary files 01.01 writes, finds
# the chosen model by name and ranks it against the other candidates.

write.candidate <- function(folder, label, aic, p, nobs = 100) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    paste0("Model name: ", label),
    "Key: Hazard-rate",
    "Formula: ~1",
    paste0("AIC: ", aic),
    "Chi_scores: 0.1, 0.2",
    "Chisquare p: 0.5",
    paste0("Det prob: ", p),
    "SE(p): 0.01",
    paste0("CV(p): ", round(0.01 / p, 3)),
    paste0("Number of observations : ", nobs),
    "N in covered region  250 20 0.08"
  ), file.path(folder, paste0("AIC_", aic, "_", label, ".txt")))
}

spec <- function(...) {
  utils::modifyList(list(strip = FALSE, final.key = NA, final.formula = ~1,
                         final.adj = NULL, notes = NA), list(...))
}

test_that("the chosen model is found, ranked and compared with the AIC-best", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "SPP", "DF Summaries_ECSAS_Ship_W_D")
  write.candidate(f, "hr.Season", 100, 0.40)
  write.candidate(f, "hn.cos", 103.5, 0.45)

  dml <- list(SPP = list(
    ECSAS_Ship_W_D = spec(final.key = "hn", final.adj = "cos", notes = "fits A better"),
    ECSAS_Ship_F_D = spec(strip = TRUE),
    ECSAS_Ship_W_N = spec(strip = TRUE)
  ))
  res <- summarise.chosen.ddfs(dml, dir)

  # _N types are left out; the strip D type keeps a row
  expect_equal(nrow(res), 2)
  w <- res[res$ddftype == "ECSAS_Ship_W_D", ]
  expect_equal(w$status, "fitted")
  expect_equal(w$model, "hn.cos")
  expect_equal(w$behaviour, "Water")
  expect_equal(w$p, 0.45)
  expect_equal(w$chisq, 0.3)
  expect_equal(w$gof_p, 0.5)
  expect_equal(w$rank, 2L)
  expect_equal(w$delta_aic, 3.5)
  expect_equal(w$best_model, "hr.Season")
  expect_equal(w$best_p, 0.40)
  expect_equal(w$notes, "fits A better")
  expect_equal(res$status[res$ddftype == "ECSAS_Ship_F_D"], "strip")
})

test_that("a choice with no matching candidate file is reported, not dropped", {
  dir <- withr::local_tempdir()
  write.candidate(file.path(dir, "SPP", "DF Summaries_ECSAS_Ship_W_D"), "hr", 100, 0.4)
  dml <- list(SPP = list(ECSAS_Ship_W_D = spec(final.key = "hn", final.formula = ~ size),
                         ECSAS_Ship_F_D = spec()))
  res <- summarise.chosen.ddfs(dml, dir)
  expect_equal(res$status, c("not a candidate", "no choice"))
  expect_equal(res$best_model[1], "hr")
})
