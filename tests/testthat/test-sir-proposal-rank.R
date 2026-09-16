# R5: rank-deficient proposals, and a scale-equivariant PD repair ------------
#
# The empirical covariance of m retained vectors in p dimensions has rank at
# most m - 1, so a full-rank proposal needs m > p. Forcing such a matrix
# positive definite does not recover the missing information: it invents
# variation in directions the retained sample never supported, and the
# proposal then draws along them.
#
# The old repair floored every eigenvalue at the ABSOLUTE value
# sqrt(.Machine$double.eps). That is not scale-equivariant: on one and the same
# singular problem, expressed in different units, the floor is 8.7e-3 of the
# largest eigenvalue at one scale and 8.7e-15 at another. Pharmacometric
# parameters genuinely span those scales.
#
# Policy: a genuinely rank-deficient retained sample is an error (as in PsN),
# and the positive-definite repair exists only to clean up floating-point
# roundoff on a matrix that is already full rank -- relative to that matrix's
# own scale.

.r5Singular <- function(m, p, scale = 1) {
  set.seed(1)
  stats::cov(matrix(stats::rnorm(m * p), m, p) * scale)
}

test_that(".sirEnsurePosDef is scale-equivariant", {
  # The statistical property that the absolute floor violated: repairing a
  # covariance then rescaling must equal rescaling then repairing.
  set.seed(11)
  p <- 4L
  a <- matrix(stats::rnorm(p * p), p, p)
  m <- crossprod(a)
  # Nudge one eigenvalue to a tiny negative value: pure roundoff, full rank.
  e <- eigen(m, symmetric = TRUE)
  e$values[p] <- -1e-18
  m <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  m <- (m + t(m)) / 2

  for (s in c(1e-3, 1, 1e3)) {
    expect_equal(
      .sirEnsurePosDef(s^2 * m),
      s^2 * .sirEnsurePosDef(m),
      tolerance = 1e-8,
      info = paste("scale", s)
    )
  }
})

test_that(".sirEnsurePosDef leaves a well-conditioned matrix alone", {
  set.seed(12)
  p <- 3L
  a <- matrix(stats::rnorm(p * p), p, p)
  m <- crossprod(a) + diag(p)
  expect_equal(.sirEnsurePosDef(m), m, tolerance = 1e-12)

  # And at a very small scale, where the old absolute floor would have
  # swamped every eigenvalue.
  small <- m * 1e-10
  expect_equal(.sirEnsurePosDef(small), small, tolerance = 1e-20)
})

test_that(".sirEnsurePosDef does not invent rank at any scale", {
  # A singular matrix must not come back with usable variance in its null
  # space, whatever the units.
  for (s in c(1e-3, 1, 1e3)) {
    cm <- .r5Singular(3L, 4L, scale = s)
    fixed <- .sirEnsurePosDef(cm)
    ev <- eigen(fixed, symmetric = TRUE)$values
    expect_lt(min(ev) / max(ev), 1e-9, label = paste("scale", s))
  }
})

test_that(".sirCheckProposalRank rejects a rank-deficient sample", {
  # 3 vectors cannot support a 4-parameter covariance.
  mat <- matrix(stats::rnorm(3L * 4L), 3L, 4L,
                dimnames = list(NULL, c("a", "b", "c", "d")))
  expect_error(.sirCheckProposalRank(mat), "rank")

  # Duplicated vectors are deficient even when there are enough rows.
  dup <- matrix(rep(c(1, 2, 3), each = 4L), nrow = 4L, ncol = 3L,
                dimnames = list(NULL, c("a", "b", "c")))
  expect_error(.sirCheckProposalRank(dup), "rank")
})

test_that(".sirCheckProposalRank accepts a full-rank sample", {
  set.seed(13)
  mat <- matrix(stats::rnorm(20L * 3L), 20L, 3L,
                dimnames = list(NULL, c("a", "b", "c")))
  expect_silent(.sirCheckProposalRank(mat))
})

test_that("sirUpdateProposal refuses a rank-deficient retained sample", {
  mat <- matrix(stats::rnorm(3L * 5L), 3L, 5L,
                dimnames = list(NULL, letters[1:5]))
  expect_error(
    sirUpdateProposal(mat, boxcox = FALSE),
    "rank"
  )
})

test_that("runSIR rejects nResample too small for the parameter count", {
  skip_on_cran()
  fit <- theoFit()
  p <- nrow(.sirParamSpace(fit))
  dir <- withr::local_tempdir()

  # A full-rank empirical covariance needs more retained vectors than
  # parameters. Catch it before any sampling or model evaluation happens.
  expect_error(
    .sirQuiet(runSIR(
      fit,
      nSamples = 40L,
      nResample = p,
      directory = dir,
      control = runSIRControl(recover = FALSE, workers = 1L)
    )),
    "nResample"
  )
  expect_false(file.exists(file.path(dir, "raw_results.csv")))
})

test_that("a raw-results proposal with too few vectors is an error, not a warning", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  path <- sirRawResultsPath()

  # Skip most of the file so fewer vectors survive than there are parameters.
  # PsN validates raw-results rank before use; the old code warned and carried
  # on with a forced-positive-definite matrix.
  total <- nrow(utils::read.csv(path))
  expect_error(
    .sirProposalFromRawResults(
      fit,
      ps,
      rawresInput = path,
      offsetRawres = total - 2L,
      boxcox = FALSE,
      capCorrelation = 0.8
    ),
    "rank|Too few"
  )
})

test_that("nSamples too small for the parameter count fails before evaluation", {
  skip_on_cran()
  fit <- theoFit()
  np <- nrow(.sirParamSpace(fit))
  # nResample is large enough to clear the existing check, so this can only be
  # the nSamples one. A candidate must be drawn before it can be retained.
  expect_error(
    runSIR(
      fit,
      nSamples = np,
      nResample = np + 1L,
      control = runSIRControl(workers = 1L, saveFiles = FALSE)
    ),
    "nSamples"
  )
})

test_that("a resampling cap that cannot supply nResample fails before evaluation", {
  skip_on_cran()
  fit <- theoFit()
  np <- nrow(.sirParamSpace(fit))
  # capResampling = 1 retains each candidate at most once, so nSamples draws
  # can never fill more than nSamples slots.
  expect_error(
    runSIR(
      fit,
      nSamples = np + 2L,
      nResample = np + 5L,
      control = runSIRControl(workers = 1L, capResampling = 1, saveFiles = FALSE)
    ),
    "capResampling"
  )
})

test_that("too few usable candidates aborts with a feasibility message", {
  skip_on_cran()
  fit <- theoFit()
  np <- nrow(.sirParamSpace(fit))

  # The schedule is feasible as requested: 20 samples, np + 1 resamples. It is
  # infeasible as realized, because only np candidates evaluate. Before the
  # check this surfaced inside .sirCheckProposalRank() as advice to increase
  # nResample, which was not the problem.
  testthat::local_mocked_bindings(
    sirEvalOFV = function(fit, paramSamples, ...) {
      n <- nrow(paramSamples)
      c(rep(NA_real_, n - np), rep(fit$objf, np))
    }
  )

  expect_error(
    suppressWarnings(sirRunIteration(
      fit,
      mu = .sirProposalMu(fit),
      proposalCov = sirGetProposalCov(fit),
      nSamples = 20L,
      nResample = np + 1L,
      iterNum = 1L,
      directory = NULL
    )),
    "usable candidate"
  )
})
