# Cross-validation helpers for the sdmTMB DSM track (Atlantic-DSM issue #130,
# Phase 2, step 05.00).
#
# 05.00 chooses, per species, the fixed Matern range and the family of an
# sdmTMB DSM out of sample, on the same spatial-block folds and with the same
# deciding score (the pooled pattern score) as the mgcv family choice in 02.05,
# so the two verdicts are on one footing. These functions are the pieces of that
# step that have rules worth testing: the fishing term, the fold merge, and the
# two verdicts (range by "longest within tolerance of the best", family by
# 02.05's margin rule). Fitting and scoring stay in the step itself.

#' Add the fishing multiplier to an sdmTMB DSM formula
#'
#' Fishing-affected species (Atlantic-DSM issue #107) carry `Season:fishing`,
#' one log-multiplier per season for segments surveyed while the ship was
#' fishing; predictions set `fishing = 0`. The mgcv track adds the same term at
#' fit time rather than in `dsm.mod.specs`, so one spec serves every species;
#' this is the sdmTMB counterpart for `sdmtmb.mod.specs`.
#'
#' @param formula A model formula from `sdmtmb.mod.specs`.
#' @return The formula with `Season:fishing` added, keeping its environment.
#' @examples
#' add.sdmtmb.fishing.term(abundance.est ~ 0 + Season + platform)
#' @export
add.sdmtmb.fishing.term <- function(formula) {
  checkmate::expect_formula(formula)
  # A spec that already has a fishing term would get it twice (update() would
  # merge them, but a spec written that way means the spec and this step
  # disagree about who adds it), so refuse it.
  if ("fishing" %in% all.vars(formula))
    stop("the formula already uses `fishing`; sdmtmb.mod.specs should not include it, it is added per species")
  # update() with ". ~ . + Season:fishing" keeps the response and every
  # right-hand-side term and appends the interaction. It writes `0 +` as
  # `- 1`, which is the same model (no global intercept).
  stats::update(formula, . ~ . + Season:fishing)
}

#' Combine cross-validation folds into fewer, larger folds
#'
#' 05.00 uses 02.05's spatial-block folds merged in pairs (10 into 5), so that
#' every segment is still held out exactly once, in the same blocks as the mgcv
#' comparison, while the number of fits halves. Fold `k` of `n` goes to merged
#' fold `ceiling(k / (n / n.merged))`: with 10 into 5, folds 1-2 become 1,
#' 3-4 become 2, and so on. Which pairs are joined carries no meaning, because
#' blockCV assigns blocks to folds at random (`selection = "random"`).
#'
#' @param folds Integer fold of each row, 1..n (from [assign.blocks()]).
#' @param n.merged Number of folds wanted. Must divide the number of folds.
#' @return Integer vector of merged folds, the same length as `folds`.
#' @examples
#' combine.cv.folds(1:10, 5)
#' @export
combine.cv.folds <- function(folds, n.merged) {
  checkmate::expect_integerish(folds, lower = 1, any.missing = FALSE, min.len = 1)
  checkmate::expect_count(n.merged, positive = TRUE)
  n <- max(folds)
  # Every fold 1..n must be present: a missing one means the caller passed a
  # subset of rows, and the merge would then silently give unequal groups.
  if (!setequal(unique(folds), seq_len(n)))
    stop("folds must use every value 1..", n)
  if (n %% n.merged != 0)
    stop(sprintf("%d folds cannot be merged evenly into %d", n, n.merged))
  # Size of each merged group, e.g. 10 / 5 = 2 original folds per merged one.
  per <- n %/% n.merged
  as.integer(ceiling(folds / per))
}

#' Rank cross-validated candidates and decide whether the best is decisive
#'
#' 02.05's rule, on one pooled score per candidate (lower is better): the best
#' is decisive when the runner-up is worse by more than `tol` relative AND by
#' more than `abs.tol` absolute. The absolute arm matters for the pattern score
#' (a mean |log O/E|), where two good candidates at 0.02 and 0.04 differ by 50%
#' and mean nothing.
#'
#' @param scores Data frame with one row per candidate: `candidate` (any id)
#'   and `score`.
#' @param tol Relative margin (`dsm.options$family.cv.tol`, 0.05).
#' @param abs.tol Absolute margin (`dsm.options$family.cv.pattern.abs.tol`,
#'   0.05).
#' @return One-row tibble: `best`, `runner_up`, `best_score`, `second_score`,
#'   `rel_margin`, `decisive`. With one candidate, `decisive` is `FALSE`:
#'   nothing was compared.
#' @examples
#' rank.cv.candidates(data.frame(candidate = c("tw", "dl"), score = c(0.40, 0.50)),
#'                    tol = 0.05, abs.tol = 0.05)
#' @export
rank.cv.candidates <- function(scores, tol, abs.tol) {
  checkmate::expect_data_frame(scores, min.rows = 1)
  checkmate::expect_names(names(scores), must.include = c("candidate", "score"))
  checkmate::expect_number(tol, lower = 0, upper = 1)
  checkmate::expect_number(abs.tol, lower = 0)
  if (anyDuplicated(scores$candidate)) stop("one row per candidate expected")
  # Non-finite scores (a candidate that produced no usable predictions) sort
  # last and cannot be best.
  s <- scores[order(!is.finite(scores$score), scores$score), , drop = FALSE]
  two <- nrow(s) > 1
  out <- tibble::tibble(
    best = as.character(s$candidate[1]),
    runner_up = if (two) as.character(s$candidate[2]) else NA_character_,
    best_score = s$score[1],
    second_score = if (two) s$score[2] else NA_real_)
  # Relative margin as 02.05 defines it: the gap as a share of the runner-up.
  out$rel_margin <- (out$second_score - out$best_score) / out$second_score
  out$decisive <- two && is.finite(out$rel_margin) && out$rel_margin > tol &&
    (out$second_score - out$best_score) > abs.tol
  out
}

#' Choose a fixed Matern range from cross-validated scores
#'
#' The rule settled for #130 Phase 2: among the ranges whose pooled score is
#' within tolerance of the best (relative `tol` OR absolute `abs.tol`, i.e.
#' not decisively worse by 02.05's rule), take the LONGEST. A longer range
#' spreads isolated detections less, and trial 1b found 100-200 km
#' indistinguishable on NL ATPU, so a tie should resolve towards the smoother
#' surface rather than to whichever range happened to score 1% better.
#'
#' The chosen range is flagged `at_edge` when it is the shortest or longest
#' range tried and the ladder could be extended past it (`extend`): the best
#' value may then lie outside what was tried.
#'
#' @param scores Data frame, one row per range: `range_km`, `score` (lower is
#'   better).
#' @param tol,abs.tol As in [rank.cv.candidates()].
#' @param extend Ranges that may be added past either end of the ladder, e.g.
#'   `c(50, 450)`. A chosen range at an end is flagged only if `extend` holds
#'   a value beyond it that has not been tried.
#' @return One-row tibble: `range_km` (chosen), `best_range_km` (lowest
#'   score), `tied` (the ranges within tolerance, comma-separated, ascending),
#'   `best_score`, `chosen_score`, `decisive` (`TRUE` when no other range is
#'   within tolerance of the best, so nothing was tied), `at_edge`,
#'   `extend_to` (the untried range to add next, or `NA`).
#' @examples
#' choose.sdmtmb.range(data.frame(range_km = c(75, 150, 300), score = c(0.50, 0.40, 0.41)),
#'                     tol = 0.05, abs.tol = 0.05, extend = c(50, 450))
#' @export
choose.sdmtmb.range <- function(scores, tol, abs.tol, extend = numeric()) {
  checkmate::expect_data_frame(scores, min.rows = 1)
  checkmate::expect_names(names(scores), must.include = c("range_km", "score"))
  checkmate::expect_numeric(scores$range_km, lower = 0, any.missing = FALSE, unique = TRUE)
  checkmate::expect_number(tol, lower = 0, upper = 1)
  checkmate::expect_number(abs.tol, lower = 0)
  checkmate::expect_numeric(extend, lower = 0, any.missing = FALSE)
  s <- scores[is.finite(scores$score), , drop = FALSE]
  if (!nrow(s)) stop("no range has a finite score")
  best <- s[which.min(s$score), , drop = FALSE]
  # "Within tolerance" is the complement of 02.05's decisive rule applied to
  # each range against the best: a range is tied unless it is worse by more
  # than tol relative AND more than abs.tol absolute.
  gap <- s$score - best$score
  rel <- gap / s$score
  is_tied <- !(rel > tol & gap > abs.tol)
  chosen <- s[is_tied, , drop = FALSE]
  chosen <- chosen[which.max(chosen$range_km), , drop = FALSE]
  # At an edge of what was tried, and is there an untried value beyond it?
  lo <- min(scores$range_km)
  hi <- max(scores$range_km)
  beyond <- if (chosen$range_km == hi) extend[extend > hi] else
    if (chosen$range_km == lo) extend[extend < lo] else numeric()
  # The nearest untried value past the edge: the smallest one above the top,
  # or the largest one below the bottom.
  extend_to <- if (!length(beyond)) NA_real_ else
    if (chosen$range_km == hi) min(beyond) else max(beyond)
  tibble::tibble(range_km = chosen$range_km, best_range_km = best$range_km,
                 tied = paste(sort(s$range_km[is_tied]), collapse = ", "),
                 best_score = best$score, chosen_score = chosen$score,
                 decisive = sum(is_tied) == 1, at_edge = !is.na(extend_to),
                 extend_to = extend_to)
}
