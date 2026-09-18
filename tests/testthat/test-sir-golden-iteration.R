# Deliberately no skip_on_cran(): every test in this file is pure arithmetic on
# synthetic inputs -- no model fit, no fixture, no file I/O -- and the whole
# file runs in well under a second. skip_on_cran() is for tests that are slow,
# need external resources, or are fragile across platforms, and none of that
# applies here. Skipping would leave CRAN running no check at all of the
# package's mathematical core, which is the part most worth protecting.
#
# S3: a frozen end-to-end fixture for one SIR iteration
#
# What this is, precisely: a regression fixture that pins the whole component
# chain -- proposal density, Box-Cox Jacobian, importance ratios, normalized
# weights, resampling, proposal update, weight diagnostics -- against values
# computed once and frozen. Everything is synthetic arithmetic: no model fit, no
# NONMEM, no PsN, and no RNG in the inputs.
#
# What this is NOT: a claim of numerical identity with PsN. The review of 0.3
# was right that the "PsN oracle" label had been stretched too far. The genuine
# PsN oracle values live in `test-sir-psn-oracles.R`, and the cores this fixture
# exercises are checked there against PsN's own unit tests -- the correlated
# multivariate-normal density, the sample/resample count adjustment, and the
# RSE-to-variance conversion. This file's job is different: it catches a change
# anywhere in the chain, including in the parts where nlmixr2sir deliberately
# differs from PsN (the Box-Cox Jacobian).
#
# If a value here changes, that is either a bug or a deliberate change that
# needs its reason recorded in NEWS.md. It should never be updated casually.

.goldMu <- c(a = 1.0, b = 2.0, c = 3.0)

.goldCov <- matrix(
  c(
    0.40, 0.05, 0.02,
    0.05, 0.30, 0.01,
    0.02, 0.01, 0.20
  ),
  nrow = 3L,
  ncol = 3L,
  byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
)

# Six candidates within about 1.5 SD of the centre, so the weights are spread
# rather than degenerate: a fixture in which one candidate takes all the mass
# would pin almost nothing.
.goldSamples <- matrix(
  c(
    1.00, 2.00, 3.00,
    1.45, 1.70, 3.20,
    0.60, 2.35, 2.80,
    1.25, 2.20, 3.15,
    0.75, 1.80, 2.90,
    1.10, 2.45, 3.05
  ),
  ncol = 3L,
  byrow = TRUE,
  dimnames = list(NULL, c("a", "b", "c"))
)

.goldDofv <- c(0.00, 1.20, -0.60, 0.35, 2.10, -0.15)

# lambda 0.5 / 0 / 1 covers the three branches of the transform: a power, the
# log case, and the identity that must contribute no Jacobian at all.
.goldBoxcox <- data.frame(
  param = c("a", "b", "c"),
  lambda = c(0.5, 0, 1),
  delta = c(2, 3, 0),
  stringsAsFactors = FALSE
)

test_that("golden: relative proposal density and importance ratios", {
  w <- sirCalcWeights(.goldSamples, .goldMu, .goldCov, dOFV = .goldDofv)

  expect_equal(
    w$relPDF,
    c(
      1, 0.57308826190944, 0.569307391128196,
      0.842643358644274, 0.865886450711239, 0.710933662071744
    ),
    tolerance = 1e-10
  )
  expect_equal(
    w$importance_ratio,
    c(
      1, 0.957638940754906, 2.37105442263974,
      0.996218640018486, 0.404138150936207, 1.51615292451275
    ),
    tolerance = 1e-10
  )
  expect_equal(
    w$prob_resample,
    c(
      0.138022356187296, 0.132175582979699, 0.327258518061047,
      0.137500443973055, 0.0557800998173925, 0.209262998981509
    ),
    tolerance = 1e-10
  )
  expect_equal(sum(w$prob_resample), 1, tolerance = 1e-12)
})

test_that("golden: Box-Cox log-Jacobian, relative to the centre", {
  lj <- .sirBcLogJacobian(.goldSamples, .goldBoxcox) -
    .sirBcLogJacobian(
      matrix(.goldMu, nrow = 1L, dimnames = list(NULL, names(.goldMu))),
      .goldBoxcox
    )
  expect_equal(
    lj,
    c(
      0, -0.0080055674694921, 0.0038917733465218,
      -0.0792420669900498, 0.0843276830150699, -0.102572607652548
    ),
    tolerance = 1e-10
  )
  # The centre contributes exactly zero by construction.
  expect_equal(lj[1L], 0, tolerance = 1e-15)
})

test_that("golden: weights with the Box-Cox Jacobian applied", {
  lj <- .sirBcLogJacobian(.goldSamples, .goldBoxcox) -
    .sirBcLogJacobian(
      matrix(.goldMu, nrow = 1L, dimnames = list(NULL, names(.goldMu))),
      .goldBoxcox
    )
  w <- sirCalcWeights(
    .goldSamples, .goldMu, .goldCov,
    dOFV = .goldDofv, logJacobian = lj
  )
  expect_equal(
    w$prob_resample,
    c(
      0.134103374823921, 0.12945483596675, 0.316731351629618,
      0.144613475368461, 0.0498134372392397, 0.22528352497201
    ),
    tolerance = 1e-10
  )

  # And it genuinely differs from the no-Jacobian weights: this is the
  # deliberate divergence from PsN, so the fixture has to show it moving.
  plain <- sirCalcWeights(.goldSamples, .goldMu, .goldCov, dOFV = .goldDofv)
  expect_gt(max(abs(w$prob_resample - plain$prob_resample)), 1e-3)
})

test_that("golden: updated proposal covariance", {
  up <- sirUpdateProposal(.goldSamples, boxcox = FALSE, capCorrelation = 0.8)
  expect_equal(
    as.numeric(up$covMat),
    c(
      0.09875, -0.0285, 0.037848822086471,
      -0.0285, 0.0906666666666667, -0.0131666666666667,
      0.037848822086471, -0.0131666666666667, 0.0226666666666667
    ),
    tolerance = 1e-10
  )
  expect_equal(dimnames(up$covMat), list(c("a", "b", "c"), c("a", "b", "c")))
  # Full rank, so no positive-definite repair was needed.
  expect_false(up$posDefAdjusted)
})

test_that("golden: weight degeneracy diagnostics", {
  w <- sirCalcWeights(.goldSamples, .goldMu, .goldCov, dOFV = .goldDofv)
  d <- .sirWeightDiagnostics(w$prob_resample, nSuccessful = 6L)
  expect_equal(d$ess, 4.77492239275491, tolerance = 1e-10)
  expect_equal(d$essFraction, 0.795820398792484, tolerance = 1e-10)
  expect_equal(d$maxWeight, 0.327258518061047, tolerance = 1e-10)
  expect_equal(d$perplexity, 5.29887450655741, tolerance = 1e-10)
  expect_equal(d$nNonNegligible, 6L)
})

test_that("golden: resampling obeys its documented properties", {
  # Deliberately property-based rather than a frozen index vector. The exact
  # selection depends on R's RNG stream, which is not a stable contract across
  # R versions, so pinning it would make the fixture fail for a reason that has
  # nothing to do with this package.
  w <- sirCalcWeights(.goldSamples, .goldMu, .goldCov, dOFV = .goldDofv)

  set.seed(20260914)
  rs <- sirResample(.goldSamples, w, m = 4L, capResampling = 1)
  expect_equal(nrow(rs$samples), 4L)
  expect_true(all(rs$selectionOrder %in% seq_len(nrow(.goldSamples))))
  # cap = 1 is sampling without replacement: no candidate may repeat.
  expect_equal(anyDuplicated(rs$selectionOrder), 0L)
  expect_equal(sum(rs$resampleCounts), 4L)
  # Every returned row is one of the candidates.
  for (i in seq_len(nrow(rs$samples))) {
    expect_true(any(apply(
      .goldSamples, 1L, function(r) isTRUE(all.equal(unname(r), unname(rs$samples[i, ])))
    )))
  }

  set.seed(20260914)
  again <- sirResample(.goldSamples, w, m = 4L, capResampling = 1)
  expect_equal(rs$selectionOrder, again$selectionOrder)
})

test_that("golden: summary quantiles over the frozen candidate set", {
  probs <- .sirPercentileProbs()
  qa <- unname(stats::quantile(.goldSamples[, "a"], probs = probs, names = FALSE))
  expect_equal(
    qa,
    unname(stats::quantile(
      c(1.00, 1.45, 0.60, 1.25, 0.75, 1.10),
      probs = probs, names = FALSE
    )),
    tolerance = 1e-12
  )
  # PsN's percentile set, which the summary reports.
  expect_equal(
    .sirPercentileLabels(),
    c(2.5, 5, 10, 30, 50, 70, 90, 95, 97.5),
    tolerance = 1e-12
  )
})
