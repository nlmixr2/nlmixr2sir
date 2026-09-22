# Rank-deficient retained vectors: refuse by default, repair on request.
#
# The default stays abort. A covariance of m vectors in p dimensions has rank
# at most m - 1, and a deficient one means some direction carries no support;
# flooring its eigenvalue does not recover the missing information, it
# fabricates variance that the next iteration then proposes along. That is why
# this is opt-in and loud rather than silent and automatic.
#
# But refusing outright is too blunt for a real case. QR model N029 -- a
# one-compartment Michaelis-Menten model on single-dose oral data -- has a
# `fit$cov` that is numerically singular (smallest eigenvalue 3.22e-10) and a
# proposal with a negative eigenvalue (-1.86e-07). Every draw then lies in an
# 8-dimensional subspace of 9, so the retained vectors are deficient before any
# resampling choice is made. Raising nResample from 200 to 500 changed nothing.
# For that fit an analyst may legitimately want the run to proceed with the
# degenerate direction pinned to a token variance, having been told.

.sirDeficientMat <- function(n = 40L, seed = 1L) {
  set.seed(seed)
  a <- rnorm(n)
  b <- rnorm(n)
  # third column is an exact combination of the first two: rank 2 of 3
  m <- cbind(p1 = a, p2 = b, p3 = a + b)
  m
}

test_that("refusing is the default", {
  expect_equal(runSIRControl()$rankDeficiency, "abort")
  expect_equal(runSIRControl(rankDeficiency = "repair")$rankDeficiency, "repair")
  expect_error(runSIRControl(rankDeficiency = "sometimes"))
})

test_that("a rank-deficient retained set is refused by default", {
  expect_error(
    nlmixr2sir:::.sirCheckProposalRank(.sirDeficientMat()),
    "rank deficient"
  )
})

test_that("repair proceeds, but says what it is doing", {
  w <- tryCatch(
    {
      nlmixr2sir:::.sirCheckProposalRank(.sirDeficientMat(),
                                         onDeficient = "repair")
      NULL
    },
    warning = function(x) conditionMessage(x)
  )
  expect_false(is.null(w))
  expect_match(w, "rank.deficient")
  # It must not be possible to read this as harmless.
  expect_match(w, "INVENTS", fixed = TRUE)
})

test_that("repair reports the deficiency rather than hiding it", {
  r <- withCallingHandlers(
    nlmixr2sir:::.sirCheckProposalRank(.sirDeficientMat(),
                                       onDeficient = "repair"),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(r$rank, 2L)
  expect_equal(r$nParams, 3L)
  expect_equal(r$deficient, 1L)
})

test_that("the refusal names the other cause, not just too few vectors", {
  # The message used to say only "repeated or collinear draws are the likely
  # cause; increase nResample". On N029 that advice was wrong and actively
  # misleading -- the proposal was already degenerate before sampling, and more
  # vectors could not help. It must name that possibility.
  err <- tryCatch(
    nlmixr2sir:::.sirCheckProposalRank(.sirDeficientMat()),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "fit$cov", fixed = TRUE)
  expect_match(err, "rankDeficiency")
})

test_that("a repaired proposal is usable where an unrepaired one is not", {
  # The substance: after repair the covariance must actually be positive
  # definite, so sampling can proceed.
  mat <- .sirDeficientMat()
  expect_error(
    sirUpdateProposal(mat, boxcox = FALSE),
    "rank deficient"
  )
  up <- withCallingHandlers(
    sirUpdateProposal(mat, boxcox = FALSE, rankDeficiency = "repair"),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_true(isTRUE(up$rankRepaired))
  expect_false(inherits(try(chol(up$covMat), silent = TRUE), "try-error"))
  # The invented direction is a token, not a real one: it must be far smaller
  # than the directions the sample does support.
  ev <- eigen(up$covMat, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), 0)
  expect_lt(min(ev) / max(ev), 1e-6)
})

test_that("an unrepaired run is unaffected by the option existing", {
  mat <- .sirDeficientMat()
  full <- cbind(mat[, 1:2], p3 = rnorm(nrow(mat)))
  a <- sirUpdateProposal(full, boxcox = FALSE)
  b <- sirUpdateProposal(full, boxcox = FALSE, rankDeficiency = "repair")
  expect_equal(a$covMat, b$covMat)
  expect_false(isTRUE(b$rankRepaired))
})
