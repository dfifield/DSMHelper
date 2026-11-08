# refit.dsm.exact() and load.final.dsm() - the exact final refit of issue #117.


# tw()'s internals call ldTweedie() unqualified, so mgcv must be attached, as it
# always is in the pipeline.
suppressPackageStartupMessages(library(mgcv))
toy.dsm <- function(fam = mgcv::tw()) {
  set.seed(3)
  n <- 4000
  d <- data.frame(x = runif(n), z = runif(n), area = runif(n, 0.5, 1.5))
  d$off.set <- log(d$area)
  d$y <- rpois(n, exp(-1 + sin(2 * pi * d$x) + d$z) * d$area)
  m <- mgcv::bam(y ~ s(x) + s(z) + offset(off.set), data = d, family = fam,
                 method = "fREML", discrete = TRUE,
                 control = mgcv::gam.control(keepData = TRUE))
  m$data <- d   # dsm() keeps the segment data here
  m$ddf <- list(placeholder = "ddf")
  class(m) <- c("dsm", class(m))
  m
}

test_that("refit.dsm.exact returns an undiscretised dsm with the ddf carried over", {
  disc <- toy.dsm()
  exact <- refit.dsm.exact(disc, nthreads = 1)
  expect_s3_class(exact, "dsm")
  expect_s3_class(exact, "bam")
  expect_true(is_discrete_gam(disc))
  expect_false(is_discrete_gam(exact))
  expect_identical(exact$ddf, disc$ddf)
  # a non-discrete bam() drops $data; 03.00's check.platform.combination() needs it
  expect_identical(exact$data, disc$data)
  expect_s3_class(exact$data, "data.frame")
  expect_identical(length(coef(exact)), length(coef(disc)))
  info <- attr(exact, "exact.refit")
  expect_false(info$at_iter_limit)
  expect_length(info$warnings, 0)
  # same model, so the fits agree closely
  expect_gt(cor(fitted(disc), fitted(exact)), 0.99)
})

test_that("refit.dsm.exact refuses a model fitted without keepData", {
  disc <- toy.dsm()
  disc$data <- NULL
  expect_error(refit.dsm.exact(disc), "keepData")
})

test_that("load.final.dsm returns the model and refuses a missing, stale or mismatched file", {
  dir <- withr::local_tempdir()
  src <- file.path(dir, "dsm_tw_allpred.Rdata")
  fin <- file.path(dir, "SP.dsm_final.RData")
  model <- 1
  save(model, file = src)
  Sys.setFileTime(src, Sys.time() - 60)

  expect_error(load.final.dsm(fin, "dsm_tw_allpred", src), "render 02.50")

  dsm_final <- structure(list(a = 1), modname = "dsm_tw_allpred")
  save(dsm_final, file = fin)
  expect_identical(load.final.dsm(fin, "dsm_tw_allpred", src)$a, 1)
  expect_error(load.final.dsm(fin, "dsm_nb_allpred", src), "now names dsm_nb_allpred")

  Sys.setFileTime(src, Sys.time() + 60)
  expect_error(load.final.dsm(fin, "dsm_tw_allpred", src), "older than its candidate")
})
