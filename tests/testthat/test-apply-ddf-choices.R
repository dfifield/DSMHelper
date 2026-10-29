test_that("apply.ddf.choices writes the choices and refuses the mistakes 01.02 used to make silently", {
  spec <- list(strip = FALSE, final.key = NA, final.formula = stats::as.formula("~1"),
               final.adj = NULL, final.nadj = NULL, notes = NA)
  strip.spec <- spec; strip.spec$strip <- TRUE
  dml <- list(BLKI = list(ECSAS_Ship_W_D = spec, ECSAS_Ship_F_D = spec,
                          ECSAS_Ship_W_N = strip.spec))
  ch <- tibble::tribble(
    ~species, ~ddftype,     ~key,    ~formula,          ~adj,  ~notes,
    "BLKI",   "Ship_W_D",   "hn",    "~1",              "cos", "huge SE",
    "BLKI",   "Ship_F_D",   "hr",    "~ Season + size", NA,    "")

  out <- apply.ddf.choices(dml, ch)
  w <- out$BLKI$ECSAS_Ship_W_D
  f <- out$BLKI$ECSAS_Ship_F_D
  expect_equal(w$final.key, "hn")
  expect_equal(w$final.adj, "cos")
  expect_equal(w$notes, "huge SE")
  expect_equal(all.vars(f$final.formula), c("Season", "size"))
  expect_null(f$final.adj)
  expect_true(is.na(f$notes))
  expect_identical(out$BLKI$ECSAS_Ship_W_N, strip.spec)

  # strip
  ch2 <- ch; ch2$key[2] <- "strip"; ch2$formula[2] <- "~1"
  expect_true(apply.ddf.choices(dml, ch2)$BLKI$ECSAS_Ship_F_D$strip)

  # mistakes
  bad <- ch; bad$ddftype[1] <- "Ship_W"
  expect_error(apply.ddf.choices(dml, bad), "ddftype must be")
  bad <- ch; bad$species[1] <- "BKLI"
  expect_error(apply.ddf.choices(dml, bad), "not in spec.grps")
  bad <- ch; bad$key[1] <- "hnn"
  expect_error(apply.ddf.choices(dml, bad), "is not hn, hr, unif or strip")
  bad <- ch; bad$adj[2] <- "cos"
  expect_error(apply.ddf.choices(dml, bad), "both covariates and adjustment")
  bad <- ch; bad$formula[2] <- "~ Season +"
  expect_error(apply.ddf.choices(dml, bad), "does not parse")
  expect_error(apply.ddf.choices(dml, ch[1, ]), "no choice recorded for: BLKI Ship_F_D")
  expect_error(apply.ddf.choices(dml, rbind(ch, ch[1, ])), "more than one row")

  # a ddftype this SubProject does not survey is dropped with a warning
  extra <- rbind(ch, tibble::tibble(species = "BLKI", ddftype = "Aerial_W_D",
                                    key = "strip", formula = "~1", adj = NA, notes = ""))
  expect_warning(o <- apply.ddf.choices(dml, extra), "row ignored: BLKI Aerial_W_D")
  expect_null(o$BLKI$ECSAS_Aerial_W_D)
})
