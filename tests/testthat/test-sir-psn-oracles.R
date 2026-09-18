# Deliberately no skip_on_cran(): every test in this file is pure arithmetic on
# synthetic inputs -- no model fit, no fixture, no file I/O -- and the whole
# file runs in well under a second. skip_on_cran() is for tests that are slow,
# need external resources, or are fragile across platforms, and none of that
# applies here. Skipping would leave CRAN running no check at all of the
# package's mathematical core, which is the part most worth protecting.
#
# Numeric cross-validation against PsN, using the oracle values in PsN's own
# unit tests. These exercise the pure maths, so they need neither a model fit
# nor NONMEM -- which is what makes them worth porting.
#
# Source: PsN test/unit/tool/sir.t lines 405-490, against the model
# test_files/mox_sir_block2.mod (5 THETA; a 3x3 OMEGA holding a 1x1 block and
# a 2x2 block, so four estimated elements: (1,1), (2,2), (3,2), (3,3)).

# The PsN final estimates for that model, in PsN's parameter order.
.moxValues <- c(
  theta = c(32.8872, 20.9156, 0.296626, 0.0992828, 0.3337),
  omega = c(
    `1,1` = 0.409882,
    `2,2` = 1.24558,
    `3,2` = 0.136766,
    `3,3` = 0.218255
  )
)

# A stand-in fit with mox_sir_block2's structure and estimates. Unclassed:
# nlmixr2FitCore defines its own `$` method that a list cannot satisfy.
.moxFit <- function() {
  om <- matrix(
    0,
    3L,
    3L,
    dimnames = list(
      c("eta.1", "eta.2", "eta.3"),
      c("eta.1", "eta.2", "eta.3")
    )
  )
  om[1L, 1L] <- 0.409882
  om[2L, 2L] <- 1.24558
  om[3L, 2L] <- 0.136766
  om[2L, 3L] <- 0.136766
  om[3L, 3L] <- 0.218255

  list(
    nsub = 74L,
    cov = NULL,
    theta = c(
      th1 = 32.8872,
      th2 = 20.9156,
      th3 = 0.296626,
      th4 = 0.0992828,
      th5 = 0.3337
    ),
    omega = om,
    iniDf = data.frame(
      ntheta = c(1:5, NA, NA, NA, NA),
      neta1 = c(rep(NA, 5L), 1L, 2L, 3L, 3L),
      neta2 = c(rep(NA, 5L), 1L, 2L, 2L, 3L),
      name = c(
        "th1",
        "th2",
        "th3",
        "th4",
        "th5",
        "eta.1",
        "eta.2",
        "(eta.2,eta.3)",
        "eta.3"
      ),
      lower = c(rep(-Inf, 5L), rep(-Inf, 4L)),
      est = c(
        32.8872,
        20.9156,
        0.296626,
        0.0992828,
        0.3337,
        0.409882,
        1.24558,
        0.136766,
        0.218255
      ),
      upper = rep(Inf, 9L),
      fix = rep(FALSE, 9L),
      err = rep(NA_character_, 9L),
      stringsAsFactors = FALSE
    )
  )
}

test_that("the mox fixture reproduces PsN's parameter order", {
  ps <- .sirParamSpace(.moxFit())
  # PsN's parameter_hash order is theta, then OMEGA by column then row:
  # (1,1), (2,2), (3,2), (3,3).
  expect_equal(
    ps$sirName,
    c(
      "th1",
      "th2",
      "th3",
      "th4",
      "th5",
      "eta.1",
      "eta.2",
      "eta.3:eta.2",
      "eta.3"
    )
  )
  expect_equal(unname(ps$est), unname(.moxValues))
})

# --- setup_inflation() oracles, sir.t:417-432 -------------------------------

test_that("PsN inflation oracle 1: scalar omega recycles across the block", {
  ps <- .sirParamSpace(.moxFit())
  v <- .sirInflationVector(ps, thetaInflation = 1, omegaInflation = 2)
  expect_equal(unname(v), c(1, 1, 1, 1, 1, 2, 2, 2, 2))
})

test_that("PsN inflation oracle 2: per-THETA vector with scalar omega", {
  ps <- .sirParamSpace(.moxFit())
  v <- .sirInflationVector(
    ps,
    thetaInflation = c(1, 2, 3, 4, 5),
    omegaInflation = 3
  )
  expect_equal(unname(v), c(1, 2, 3, 4, 5, 3, 3, 3, 3))
})

test_that("PsN inflation oracle 3: off-diagonal is sqrt(i) * sqrt(j)", {
  ps <- .sirParamSpace(.moxFit())
  v <- .sirInflationVector(
    ps,
    thetaInflation = c(1, 2, 3, 4, 5),
    omegaInflation = c(6, 16, 25)
  )
  # PsN: [1,2,3,4,5,6,16,20,25]. The off-diagonal (3,2) is sqrt(16)*sqrt(25).
  expect_equal(unname(v), c(1, 2, 3, 4, 5, 6, 16, 20, 25))
  expect_equal(v[["eta.3:eta.2"]], sqrt(16) * sqrt(25))
})

test_that("PsN inflation oracle 4: all ones is a no-op", {
  ps <- .sirParamSpace(.moxFit())
  v <- .sirInflationVector(
    ps,
    thetaInflation = rep(1, 5L),
    omegaInflation = 1
  )
  # PsN signals "no inflation" by returning an empty vector; this returns all
  # ones, which is the same thing applied.
  expect_true(all(v == 1))
})

test_that("PsN illegal-inflation cases are rejected", {
  ps <- .sirParamSpace(.moxFit())
  # sigma inflation given for a model with no residual-error parameters
  expect_error(.sirInflationVector(ps, sigmaInflation = 2), "no estimated")
  # four values for five THETAs
  expect_error(
    .sirInflationVector(ps, thetaInflation = c(1, 2, 3, 4)),
    "one value per"
  )
  # four values for three OMEGA diagonals
  expect_error(
    .sirInflationVector(ps, omegaInflation = c(1, 2, 3, 4)),
    "one value per"
  )
  # a negative factor
  expect_error(.sirInflationVector(ps, thetaInflation = -11), "thetaInflation")
})

# --- RSE-to-variance oracles, sir.t:459-486 ---------------------------------

test_that("PsN get_offdiagonal_variance type 1 oracle", {
  # cmp_float(get_offdiagonal_variance(type = 1, covariance = 0.03,
  #   rse_i = 20, rse_j = 40, var_i = 3, var_j = 0.5),
  #   (0.03**2 + 1.5) / (25 + 2.5**2 + 1))
  n <- (100 / 20)^2 + (100 / 40)^2 + 1
  expected <- (0.03^2 + 3 * 0.5) / n
  expect_equal(expected, (0.03^2 + 1.5) / (25 + 2.5^2 + 1))

  # The same formula as implemented here, via a two-eta fixture.
  fit <- .moxFit()
  ps <- .sirParamSpace(fit)
  v <- .sirRseVariance(ps, rseTheta = 20, rseOmega = c(20, 20, 40))
  off <- v[["eta.3:eta.2"]]
  nn <- (100 / 40)^2 + (100 / 20)^2 + 1
  expect_equal(
    off,
    (0.136766^2 + 1.24558 * 0.218255) / nn,
    tolerance = 1e-12
  )
})

test_that("PsN setup_variancevec diagonals are (rse * estimate / 100)^2", {
  fit <- .moxFit()
  ps <- .sirParamSpace(fit)
  v <- .sirRseVariance(ps, rseTheta = 20, rseOmega = 10)

  # PsN oracle 1's diagonal entries, with rse_theta = 20 and rse_omega = 10.
  expect_equal(v[["th1"]], (32.8872 * 0.2)^2, tolerance = 1e-12)
  expect_equal(v[["th2"]], (20.9156 * 0.2)^2, tolerance = 1e-12)
  expect_equal(v[["th3"]], (0.296626 * 0.2)^2, tolerance = 1e-12)
  expect_equal(v[["th4"]], (0.0992828 * 0.2)^2, tolerance = 1e-12)
  expect_equal(v[["th5"]], (0.3337 * 0.2)^2, tolerance = 1e-12)
  expect_equal(v[["eta.1"]], (0.409882 * 0.1)^2, tolerance = 1e-12)
  expect_equal(v[["eta.2"]], (1.24558 * 0.1)^2, tolerance = 1e-12)
  expect_equal(v[["eta.3"]], (0.218255 * 0.1)^2, tolerance = 1e-12)
})

test_that("PsN setup_variancevec oracle 2: per-parameter RSE vectors", {
  fit <- .moxFit()
  ps <- .sirParamSpace(fit)
  v <- .sirRseVariance(
    ps,
    rseTheta = c(20, 30, 40, 30, 10),
    rseOmega = c(15, 30, 20),
    rseSigma = NULL
  )
  expect_equal(v[["th1"]], (32.8872 * 0.20)^2, tolerance = 1e-12)
  expect_equal(v[["th2"]], (20.9156 * 0.30)^2, tolerance = 1e-12)
  expect_equal(v[["th3"]], (0.296626 * 0.40)^2, tolerance = 1e-12)
  expect_equal(v[["th4"]], (0.0992828 * 0.30)^2, tolerance = 1e-12)
  expect_equal(v[["th5"]], (0.3337 * 0.10)^2, tolerance = 1e-12)
  expect_equal(v[["eta.1"]], (0.409882 * 0.15)^2, tolerance = 1e-12)
  expect_equal(v[["eta.2"]], (1.24558 * 0.30)^2, tolerance = 1e-12)
  expect_equal(v[["eta.3"]], (0.218255 * 0.20)^2, tolerance = 1e-12)
})

test_that("PsN illegal setup_variancevec cases are rejected", {
  ps <- .sirParamSpace(.moxFit())
  # a zero RSE
  expect_error(
    .sirRseVariance(ps, rseTheta = c(20, 30, 40, 30, 10), rseOmega = 0),
    "rseOmega"
  )
  # four omega values for three diagonals
  expect_error(
    .sirRseVariance(
      ps,
      rseTheta = c(20, 30, 40, 30, 10),
      rseOmega = c(1, 2, 3, 4)
    ),
    "one value per"
  )
  # four theta values for five THETAs
  expect_error(
    .sirRseVariance(
      ps,
      rseTheta = c(20, 40, 30, 10),
      rseOmega = c(1, 2, 3)
    ),
    "one value per"
  )
})

# PsN's multivariate-normal density oracle ------------------------------------
#
# Source: PsN test/unit/tool/sir.t, the mvnpdf_cholesky block (approx. lines
# 238-262). These are the values PsN itself checks against Matlab's mvnpdf.
#
# This is the oracle that matters most: sirCalcWeights() divides the likelihood
# ratio by exactly this density, so an error here misprices every candidate.
# It is deliberately *correlated* -- a diagonal covariance cannot distinguish a
# Cholesky solve from its transpose, which is why the existing relPDF tests
# passed while the quadratic form was wrong.

.psnMvnSigma <- matrix(
  c(
    3.0, 0.1, 0.2,
    0.1, 8.0, 0.3,
    0.2, 0.3, 2.0
  ),
  nrow = 3L,
  ncol = 3L,
  byrow = TRUE
)

test_that("PsN oracle: relative MVN density for a correlated covariance", {
  mu <- c(1, 2, 3)
  x <- matrix(c(0, 0, 0), nrow = 1L)

  res <- sirCalcWeights(x, mu = mu, covMat = .psnMvnSigma, dOFV = 0)

  expect_equal(res$relPDF[1L], 0.08373785511747776, tolerance = 1e-12)
})

test_that("relative MVN density equals the explicit Mahalanobis form", {
  mu <- c(1, 2, 3)
  pts <- matrix(
    c(
      0, 0, 0,
      1, 2, 3,
      -2, 5, 1.5,
      4.25, -1, 0.75
    ),
    ncol = 3L,
    byrow = TRUE
  )

  res <- sirCalcWeights(
    pts,
    mu = mu,
    covMat = .psnMvnSigma,
    dOFV = rep(0, nrow(pts))
  )

  inv <- solve(.psnMvnSigma)
  expected <- apply(pts, 1L, function(z) {
    d <- z - mu
    exp(-0.5 * as.numeric(t(d) %*% inv %*% d))
  })

  expect_equal(res$relPDF, expected, tolerance = 1e-12)
})

test_that("relative MVN density matches mvtnorm::dmvnorm up to the constant", {
  skip_if_not_installed("mvtnorm")
  set.seed(20260914)
  p <- 4L
  a <- matrix(stats::rnorm(p * p), p, p)
  sigma <- crossprod(a) + diag(p) # positive definite by construction
  mu <- stats::rnorm(p)
  pts <- matrix(stats::rnorm(6L * p), ncol = p)

  res <- sirCalcWeights(pts, mu = mu, covMat = sigma, dOFV = rep(0, nrow(pts)))

  # dmvnorm carries the normalising constant; the density at mu is that
  # constant, so the ratio is exactly relPDF.
  expected <- exp(
    mvtnorm::dmvnorm(pts, mean = mu, sigma = sigma, log = TRUE) -
      mvtnorm::dmvnorm(t(mu), mean = mu, sigma = sigma, log = TRUE)
  )

  expect_equal(res$relPDF, as.numeric(expected), tolerance = 1e-10)
})

test_that("a diagonal covariance cannot detect the transposed solve", {
  # Pins the reason the bug survived: with a diagonal covariance the Cholesky
  # factor is its own transpose, so both solves agree. This test documents the
  # blind spot rather than guarding behaviour.
  mu <- c(1, 2, 3)
  x <- matrix(c(0, 0, 0), nrow = 1L)
  d <- diag(c(3.0, 8.0, 2.0))

  res <- sirCalcWeights(x, mu = mu, covMat = d, dOFV = 0)

  expect_equal(
    res$relPDF[1L],
    exp(-0.5 * sum((x[1L, ] - mu)^2 / c(3.0, 8.0, 2.0))),
    tolerance = 1e-12
  )
})
