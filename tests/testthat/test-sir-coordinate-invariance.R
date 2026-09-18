# Deliberately no skip_on_cran(): every test in this file is pure arithmetic on
# synthetic inputs -- no model fit, no fixture, no file I/O -- and the whole
# file runs in well under a second. skip_on_cran() is for tests that are slow,
# need external resources, or are fragile across platforms, and none of that
# applies here. Skipping would leave CRAN running no check at all of the
# package's mathematical core, which is the part most worth protecting.
#
# B1: rank and positive-definite repair must not depend on parameter units
#
# Pharmacometric parameters do not share a unit. A clearance, a log-scale
# THETA, a small OMEGA element and a residual SD can differ by many orders of
# magnitude. Two samples that differ only by a per-parameter unit conversion
# describe the same statistical problem and must receive the same rank and
# repair decisions.
#
# The previous implementation tested eigenvalues of the covariance in its raw
# coordinates, so the answer moved with the units. The earlier scale test only
# ever multiplied the WHOLE matrix by one scalar, which cannot detect this: a
# global rescale leaves eigenvalue *ratios* untouched, while an independent
# per-column rescale does not. These tests use a diagonal D with entries
# spanning many orders of magnitude.

.ciDiag <- function(p, seed = 7L) {
  set.seed(seed)
  diag(10^seq(-6, 6, length.out = p))
}

test_that("numerical rank is invariant to independent column rescaling", {
  set.seed(11)
  x <- cbind(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  d <- .ciDiag(3L)
  scaled <- x %*% d
  colnames(scaled) <- colnames(x)

  # Both are full rank: rescaling a column cannot destroy independence.
  expect_silent(.sirCheckProposalRank(x))
  expect_silent(.sirCheckProposalRank(scaled))
})

test_that("the reported rank itself is unchanged by column rescaling", {
  set.seed(12)
  x <- cbind(a = rnorm(60), b = rnorm(60), c = rnorm(60), d = rnorm(60))
  d <- .ciDiag(4L)
  scaled <- x %*% d
  colnames(scaled) <- colnames(x)

  expect_equal(
    .sirCheckProposalRank(x)$rank,
    .sirCheckProposalRank(scaled)$rank
  )
})

test_that("the two-column 1e-6 counterexample is accepted", {
  # Reported in the follow-up review: a full-rank sample rejected as rank one
  # purely because the second column was expressed in different units.
  set.seed(11)
  x <- cbind(a = rnorm(100), b = 1e-6 * rnorm(100))
  expect_silent(.sirCheckProposalRank(x))
  expect_equal(.sirCheckProposalRank(x)$rank, 2L)
})

test_that("genuine rank deficiency is still rejected after rescaling", {
  set.seed(13)
  a <- rnorm(50)
  b <- rnorm(50)
  # c is an exact linear combination: deficient in any units.
  x <- cbind(a = a, b = b, c = 2 * a - 3 * b)
  expect_error(.sirCheckProposalRank(x), "rank")

  d <- .ciDiag(3L)
  scaled <- x %*% d
  colnames(scaled) <- colnames(x)
  expect_error(.sirCheckProposalRank(scaled), "rank")
})

test_that("an exactly constant column is reported by name", {
  set.seed(14)
  x <- cbind(a = rnorm(40), fixedone = rep(2.5, 40L), c = rnorm(40))
  expect_error(.sirCheckProposalRank(x), "fixedone")
})

test_that("a small but independent column is not treated as constant", {
  # The distinction the raw-eigenvalue test could not make: tiny variance is
  # not the same thing as no variance.
  set.seed(15)
  x <- cbind(a = rnorm(80), tiny = 1e-9 * rnorm(80))
  expect_silent(.sirCheckProposalRank(x))
  expect_equal(.sirCheckProposalRank(x)$rank, 2L)
})

# Positive-definite repair ----------------------------------------------------

test_that("PD repair commutes with independent column rescaling", {
  set.seed(16)
  p <- 4L
  a <- matrix(rnorm(p * p), p, p)
  s <- crossprod(a)
  # Push one eigenvalue slightly negative: a repair is genuinely needed.
  e <- eigen(s, symmetric = TRUE)
  e$values[p] <- -1e-14 * max(e$values)
  s <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  s <- (s + t(s)) / 2

  d <- .ciDiag(p)
  expect_equal(
    .sirEnsurePosDef(d %*% s %*% d),
    d %*% .sirEnsurePosDef(s) %*% d,
    tolerance = 1e-8
  )
})

test_that("PD repair preserves marginal variances exactly", {
  # The defect: diag(c(1, 1e-14)) was returned as roughly diag(c(1, 1e-12)),
  # inflating the second variance a hundredfold because the first happened to
  # be 1. A repair must never change a parameter's marginal variance on the
  # strength of another parameter's units.
  m <- diag(c(1, 1e-14))
  out <- .sirEnsurePosDef(m)
  expect_equal(diag(out), diag(m), tolerance = 1e-20)

  set.seed(17)
  p <- 3L
  a <- matrix(rnorm(p * p), p, p)
  s <- crossprod(a) * 1e-8
  s[1L, 1L] <- s[1L, 1L] * 1e10
  repaired <- .sirEnsurePosDef(s)
  expect_equal(diag(repaired), diag(s), tolerance = 1e-10)
})

test_that("PD repair still returns a positive-definite matrix", {
  set.seed(18)
  p <- 5L
  a <- matrix(rnorm(p * p), p, p)
  s <- crossprod(a)
  e <- eigen(s, symmetric = TRUE)
  e$values[c(p - 1L, p)] <- c(-1e-12, -1e-10) * max(e$values)
  s <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  s <- (s + t(s)) / 2

  out <- .sirEnsurePosDef(s)
  expect_false(inherits(try(chol(out), silent = TRUE), "try-error"))
  expect_true(isSymmetric(unname(out)))
})

test_that("PD repair leaves a well-conditioned matrix untouched", {
  set.seed(19)
  p <- 3L
  a <- matrix(rnorm(p * p), p, p)
  m <- crossprod(a) + diag(p)
  expect_equal(.sirEnsurePosDef(m), m, tolerance = 1e-12)
  # Including when the parameters differ wildly in scale.
  d <- .ciDiag(p)
  scaled <- d %*% m %*% d
  expect_equal(.sirEnsurePosDef(scaled), scaled, tolerance = 1e-10)
})

test_that("PD repair tolerates a zero-variance coordinate", {
  # A Wishart fallback produces a zero proposal variance for an OMEGA element
  # estimated at zero. Such a coordinate cannot be standardized, and must not
  # acquire uncertainty it was never given.
  s <- diag(c(0.5, 0, 0.25))
  out <- .sirEnsurePosDef(s)
  expect_equal(diag(out), c(0.5, 0, 0.25), tolerance = 1e-12)
  expect_true(all(out[2L, ] == 0))
  expect_true(all(out[, 2L] == 0))
})
