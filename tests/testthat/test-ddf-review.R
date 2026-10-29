test_that("label parts invert create.model.name", {
  expect_equal(.review.label.parts("hn.cos"), list(key = "hn", formula = "~1", adj = "cos"))
  expect_equal(.review.label.parts("hr.hp")$adj, "herm")
  expect_equal(.review.label.parts("unif.sp")$adj, "poly")
  expect_equal(.review.label.parts("hr.Season.size"),
               list(key = "hr", formula = "~ Season + size", adj = NA_character_))
  expect_equal(.review.label.parts("hr"), list(key = "hr", formula = "~1", adj = NA_character_))
  for (lab in c("hn.cos", "hr.hp", "hr.Season.size", "hr.size", "hr")) {
    p <- .review.label.parts(lab)
    expect_equal(create.model.name(p$key, stats::as.formula(p$formula),
                                   if (is.na(p$adj)) NULL else p$adj), lab)
  }
})

test_that("ddf.choices.text aligns columns and round-trips through tribble()", {
  ch <- data.frame(species = c("BLKI", "Shearwaters"), ddftype = c("Ship_W_D", "Ship_F_D"),
                   key = c("hn", "strip"), formula = c("~1", "~ Season"),
                   adj = c("cos", "poly"), notes = c("say \"why\"", ""),
                   stringsAsFactors = FALSE)
  txt <- ddf.choices.text(ch)
  expect_equal(txt[1], "ddf.choices <- tribble(")
  expect_equal(txt[length(txt)], ")")
  env <- new.env()
  env$tribble <- tibble::tribble
  eval(parse(text = txt), envir = env)
  got <- as.data.frame(env$ddf.choices)
  expect_equal(got$key, c("hn", "strip"))
  expect_equal(got$formula, c("~1", "~1"))      # strip rows carry no model
  expect_true(is.na(got$adj[2]))
  expect_equal(got$notes[1], "say \"why\"")
  # every data line has the same column starts
  pos <- vapply(txt[2:4], function(l) as.integer(regexpr("~key|\"hn\"|\"strip\"", l)), 1L)
  expect_equal(length(unique(pos)), 1)
})

test_that("flag.ddf.candidates raises the checklist codes", {
  cands <- tibble::tibble(
    species = "BLKI", ddftype = "ECSAS_Ship_W_D",
    model = c("hr.cos", "hn.cos", "hn.Season", "hr.size"),
    key = c("hr", "hn", "hn", "hr"), formula = c("~1", "~1", "~ Season", "~ size"),
    adj = c("cos", "cos", NA, NA), failed = c(FALSE, FALSE, FALSE, TRUE),
    AIC = c(100, 100.5, 104, NA), delta_aic = c(0, 0.5, 4, NA),
    p = c(0.37, 0.30, 0.60, NA), gof_p = c(0.2, 0.01, NA, NA),
    chi_first_max = c(FALSE, TRUE, NA, NA), nadj = c(1L, 0L, NA, NA),
    opt_flag = c("", "AIC +5.0 at nadj 1", "", ""), n_par = c(3L, 2L, 4L, NA),
    mono = c(TRUE, FALSE, NA, NA), hr_spike = c(FALSE, FALSE, FALSE, NA),
    n_p_low = c(0L, 0L, 3L, NA), min_season_n = c(50L, 50L, 12L, NA), n = 300)
  ch <- data.frame(species = "BLKI", ddftype = "ECSAS_Ship_W_D", model = "hn.Season",
                   notes = "", stringsAsFactors = FALSE)
  fl <- flag.ddf.candidates(cands, ch)
  f <- setNames(fl$candidates$flags, fl$candidates$model)
  expect_equal(f[["hr.cos"]], "")
  expect_match(f[["hn.cos"]], "nadj0")
  expect_match(f[["hn.cos"]], "optimiser")
  expect_match(f[["hn.cos"]], "gof")
  expect_match(f[["hn.cos"]], "bin1")
  expect_match(f[["hn.cos"]], "not-monotone")
  expect_match(f[["hn.Season"]], "p-outlier")
  expect_match(f[["hn.Season"]], "params>3")
  expect_match(f[["hn.Season"]], "low-p")
  expect_match(f[["hn.Season"]], "sparse-season")
  expect_equal(f[["hr.size"]], "failed")
  t <- fl$ddftypes
  expect_equal(t$best_model, "hr.cos")
  expect_match(t$choice_flags, "no-note")      # 4 AIC behind, no note
  expect_match(t$choice_flags, "p-spread")     # 0.37 vs 0.30 within 2 AIC
  ch$model <- "hr.Season"
  expect_match(flag.ddf.candidates(cands, ch)$ddftypes$choice_flags, "not-a-candidate")
  ch$model <- NA
  expect_match(flag.ddf.candidates(cands, ch)$ddftypes$choice_flags, "no-choice")
})

test_that("collect.ddf.candidates reads what run.ddf.model saves", {
  skip_if_not_installed("Distance")
  library(Distance)  # run.ddf.model() calls ds() unqualified; the pipeline has it attached
  set.seed(114)
  x <- abs(rnorm(500, sd = 0.1)); x <- x[x <= 0.3]
  cuts <- c(0, 0.05, 0.1, 0.2, 0.3)
  bin <- cut(x, cuts, include.lowest = TRUE, labels = FALSE)
  dat <- data.frame(distbegin = cuts[bin], distend = cuts[bin + 1], object = seq_along(x),
                    size = sample(1:3, length(x), TRUE),
                    Season = sample(c("Fall", "Spring"), length(x), TRUE))
  res <- withr::local_tempdir()
  folder <- file.path(res, "BLKI", "DF Summaries_ECSAS_Ship_W_D")
  dir.create(folder, recursive = TRUE)
  mods <- data.frame(label = c("hn.cos", "hr.size"), key = c("hn", "hr"),
                     adj = c("cos", ""), form = c("~1", "~size"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(mods)))
    suppressWarnings(run.ddf.model(mods[i, ], data = dat, folder = folder, verbose = FALSE,
                                   truncation = 0.3, max_adjustments = 2))
  file.create(file.path(folder, "failed_model_hn.Season.txt"))
  dml <- list(BLKI = list(ECSAS_Ship_W_D = list(), ECSAS_Ship_W_N = list()))
  cands <- collect.ddf.candidates(res, dml)
  expect_setequal(cands$model, c("hn.cos", "hr.size", "hn.Season"))
  expect_true(cands$failed[cands$model == "hn.Season"])
  hn <- cands[cands$model == "hn.cos", ]
  expect_equal(hn$adj, "cos")
  expect_false(is.na(hn$nadj))
  expect_equal(length(hn$obs[[1]]), 4)
  expect_equal(sum(hn$obs[[1]]), hn$n)
  expect_equal(hn$bins[[1]], c("0-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.3"))
  expect_match(cands$season_n[cands$model == "hr.size"], "Fall")
  expect_equal(min(cands$delta_aic, na.rm = TRUE), 0)
})
