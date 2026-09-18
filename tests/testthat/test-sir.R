# Step 1: sirGetProposalCov ----------------------------------------------------

test_that("sirGetProposalCov returns fit$cov values under SIR names", {
  skip_on_cran()
  result <- sirGetProposalCov(theoFit())
  ps <- .sirParamSpace(theoFit())
  expect_equal(unname(result), unname(theoFit()$cov), tolerance = 1e-10)
  expect_equal(rownames(result), ps$sirName)
  expect_equal(colnames(result), ps$sirName)
})

test_that("sirGetProposalCov inflates each kind with its own factor", {
  skip_on_cran()
  orig <- theoFit()$cov
  ps <- .sirParamSpace(theoFit())
  result <- sirGetProposalCov(
    theoFit(),
    thetaInflation = 2,
    sigmaInflation = 3,
    omegaInflation = 4,
    capCorrelation = 1
  )
  expected <- c(theta = 2, sigma = 3, omegaDiag = 4)[ps$kind]
  expect_equal(
    unname(diag(result)),
    unname(diag(orig) * expected),
    tolerance = 1e-10
  )
  # Correlations must be preserved
  expect_equal(
    unname(cov2cor(result)),
    unname(cov2cor(orig)),
    tolerance = 1e-10
  )
})

test_that("sirGetProposalCov applies sigmaInflation to residual error only", {
  skip_on_cran()
  orig <- theoFit()$cov
  ps <- .sirParamSpace(theoFit())
  result <- sirGetProposalCov(theoFit(), sigmaInflation = 9, capCorrelation = 1)
  isSigma <- ps$kind == "sigma"
  expect_equal(
    unname(diag(result)[isSigma]),
    unname(diag(orig)[isSigma] * 9),
    tolerance = 1e-10
  )
  expect_equal(
    unname(diag(result)[!isSigma]),
    unname(diag(orig)[!isSigma]),
    tolerance = 1e-10
  )
})


test_that("sirGetProposalCov caps correlations at capCorrelation", {
  skip_on_cran()
  # Inflate theta 10x to push correlations toward 1, then cap at 0.3
  result <- sirGetProposalCov(
    theoFit(),
    thetaInflation = 10,
    capCorrelation = 0.3
  )
  corr <- cov2cor(result)
  off <- corr[lower.tri(corr)]
  expect_true(all(off >= -0.3 - 1e-10 & off <= 0.3 + 1e-10))
})

test_that("sirGetProposalCov with capCorrelation = 1 leaves correlations unchanged", {
  skip_on_cran()
  orig <- theoFit()$cov
  result <- sirGetProposalCov(theoFit(), capCorrelation = 1)
  expect_equal(
    unname(cov2cor(result)),
    unname(cov2cor(orig)),
    tolerance = 1e-10
  )
})

test_that("sirGetProposalCov result is symmetric", {
  skip_on_cran()
  result <- sirGetProposalCov(
    theoFit(),
    thetaInflation = 2,
    capCorrelation = 0.8
  )
  expect_equal(result, t(result), tolerance = 1e-14)
})

test_that("sirGetProposalCov renames fit$cov rows to SIR parameter names", {
  skip_on_cran()
  result <- sirGetProposalCov(blockFit())
  ps <- .sirParamSpace(blockFit())
  expect_equal(rownames(result), ps$sirName)
  expect_equal(colnames(result), ps$sirName)
  expect_true("eta.cl:eta.ka" %in% rownames(result))
})

test_that("sirGetProposalCov errors when fit has no covariance matrix", {
  skip_on_cran()
  expect_error(sirGetProposalCov(theoFitNoCov()), "covariance matrix")
})

test_that("sirGetProposalCov errors on non-fit input", {
  skip_on_cran()
  expect_error(sirGetProposalCov(list()), class = "error")
})

# Step 2: sirSampleTheta -------------------------------------------------------

# Fixed 3-parameter setup reused across Step 2 tests
.mu3 <- c(tka = 0.45, tcl = 0.98, tv = 3.47)
.cov3 <- matrix(
  c(0.065, -0.004, 0.004, -0.004, 0.006, -0.002, 0.004, -0.002, 0.003),
  nrow = 3L,
  dimnames = list(names(.mu3), names(.mu3))
)

test_that("sirSampleTheta returns n rows with no bounds", {
  skip_on_cran()
  set.seed(1)
  res <- sirSampleTheta(.mu3, .cov3, n = 200L)
  expect_equal(nrow(res$samples), 200L)
  expect_equal(ncol(res$samples), 3L)
  expect_equal(res$nRejected, 0L)
})

test_that("sirSampleTheta column names match mu", {
  skip_on_cran()
  set.seed(1)
  res <- sirSampleTheta(.mu3, .cov3, n = 50L)
  expect_equal(colnames(res$samples), names(.mu3))
})

test_that("sirSampleTheta respects lower and upper bounds", {
  skip_on_cran()
  set.seed(42)
  lo <- .mu3 - 0.5
  hi <- .mu3 + 0.5
  res <- sirSampleTheta(.mu3, .cov3, n = 200L, lower = lo, upper = hi)
  expect_true(all(res$samples >= rep(lo, each = nrow(res$samples))))
  expect_true(all(res$samples <= rep(hi, each = nrow(res$samples))))
})

test_that("sirSampleTheta nRejected increases with tight bounds", {
  skip_on_cran()
  set.seed(7)
  # Wide bounds: expect few rejections
  res_wide <- sirSampleTheta(
    .mu3,
    .cov3,
    n = 200L,
    lower = .mu3 - 10,
    upper = .mu3 + 10
  )
  # Tight bounds: expect more rejections
  res_tight <- suppressWarnings(
    sirSampleTheta(
      .mu3,
      .cov3,
      n = 200L,
      lower = .mu3 - 0.01,
      upper = .mu3 + 0.01
    )
  )
  expect_gt(res_tight$nRejected, res_wide$nRejected)
})

test_that("sirSampleTheta is reproducible with set.seed", {
  skip_on_cran()
  set.seed(99)
  r1 <- sirSampleTheta(.mu3, .cov3, n = 50L)
  set.seed(99)
  r2 <- sirSampleTheta(.mu3, .cov3, n = 50L)
  expect_equal(r1$samples, r2$samples)
})

test_that("sirSampleTheta warns and returns fewer rows when bounds exclude all draws", {
  skip_on_cran()
  # Bounds set to a tiny region far from mu — all draws will be rejected
  lo <- .mu3 + 1e6
  hi <- .mu3 + 1e6 + 1e-9
  expect_warning(
    res <- sirSampleTheta(.mu3, .cov3, n = 10L, lower = lo, upper = hi),
    regexp = "within bounds"
  )
  expect_equal(nrow(res$samples), 0L)
  expect_equal(res$nRejected, 100L) # 10 * 10 max attempts, all rejected
})

test_that("sirSampleTheta scalar bounds are recycled to length p", {
  skip_on_cran()
  set.seed(3)
  res <- sirSampleTheta(.mu3, .cov3, n = 100L, lower = -100, upper = 100)
  expect_equal(nrow(res$samples), 100L)
  expect_true(all(res$samples >= -100 & res$samples <= 100))
})

# Step 3: sirSampleOmegaSigma --------------------------------------------------

# 2×2 OMEGA — lower triangle has 3 elements: [1,1], [2,1], [2,2]
.omega2 <- matrix(
  c(0.5, 0.1, 0.1, 0.3),
  nrow = 2L,
  dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
)
# Diagonal omegaCovMat built from plausible SEs
.omega2_cov <- diag(c(0.05^2, 0.02^2, 0.04^2))

test_that("sirSampleOmegaSigma returns n matrices", {
  skip_on_cran()
  set.seed(1)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 50L)
  expect_length(res$samples, 50L)
})

test_that("sirSampleOmegaSigma all samples are positive definite", {
  skip_on_cran()
  set.seed(2)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 100L)
  pd_ok <- vapply(
    res$samples,
    function(m) {
      tryCatch(
        {
          chol(m)
          TRUE
        },
        error = function(e) FALSE
      )
    },
    logical(1L)
  )
  expect_true(all(pd_ok))
})

test_that("sirSampleOmegaSigma all samples are symmetric", {
  skip_on_cran()
  set.seed(3)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 50L)
  sym_ok <- vapply(
    res$samples,
    function(m) isTRUE(all.equal(m, t(m))),
    logical(1L)
  )
  expect_true(all(sym_ok))
})

test_that("sirSampleOmegaSigma preserves dimnames from omegaEst", {
  skip_on_cran()
  set.seed(4)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 10L)
  dn_ok <- vapply(
    res$samples,
    function(m) {
      identical(dimnames(m), dimnames(.omega2))
    },
    logical(1L)
  )
  expect_true(all(dn_ok))
})

test_that("sirSampleOmegaSigma nRejected is non-negative integer", {
  skip_on_cran()
  set.seed(5)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 50L)
  expect_true(is.integer(res$nRejected) || is.numeric(res$nRejected))
  expect_gte(res$nRejected, 0L)
})

test_that("sirSampleOmegaSigma is reproducible with set.seed", {
  skip_on_cran()
  set.seed(77)
  r1 <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 30L)
  set.seed(77)
  r2 <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 30L)
  expect_equal(r1$samples, r2$samples)
})

test_that("sirSampleOmegaSigma warns and returns fewer matrices when budget exhausted", {
  skip_on_cran()
  # Mean is a non-PD configuration (large off-diagonal, tiny diagonals);
  # near-zero variance pins draws close to the mean → ~0% PD rate.
  bad_omega <- matrix(
    c(1e-4, 1.0, 1.0, 1e-4),
    2L,
    2L,
    dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
  )
  bad_cov <- diag(c(1e-10, 1e-10, 1e-10))
  expect_warning(
    res <- sirSampleOmegaSigma(bad_omega, bad_cov, n = 50L),
    regexp = "positive definite"
  )
  expect_lt(length(res$samples), 50L)
})

test_that("sirSampleOmegaSigma works with 1x1 OMEGA (single eta)", {
  skip_on_cran()
  omega1 <- matrix(0.4, 1L, 1L, dimnames = list("eta.ka", "eta.ka"))
  cov1 <- matrix(0.01, 1L, 1L)
  set.seed(9)
  res <- sirSampleOmegaSigma(omega1, cov1, n = 30L)
  expect_length(res$samples, 30L)
  pd_ok <- vapply(res$samples, function(m) m[1L, 1L] > 0, logical(1L))
  expect_true(all(pd_ok))
})

# Step 4: sirEvalOFV -----------------------------------------------------------

test_that("sirEvalOFV returns OFV close to fit$objf at true estimates", {
  skip_on_cran()
  mu <- .sirProposalMu(theoFit())
  param_mat <- matrix(mu, nrow = 1L, dimnames = list(NULL, names(mu)))

  ofvs <- sirEvalOFV(theoFit(), param_mat)

  expect_length(ofvs, 1L)
  expect_false(is.na(ofvs[1L]))
  expect_lt(abs(ofvs[1L] - theoFit()$objf), 1)
})

test_that("sirEvalOFV returns NA for a clearly invalid parameter vector", {
  skip_on_cran()
  mu <- .sirProposalMu(theoFit())
  bad <- mu
  bad["tka"] <- 1e10 # absurdly large; model should fail or produce NA
  param_mat <- matrix(bad, nrow = 1L, dimnames = list(NULL, names(bad)))

  ofvs <- suppressWarnings(sirEvalOFV(theoFit(), param_mat))

  expect_length(ofvs, 1L)
  expect_true(is.na(ofvs[1L]) || is.finite(ofvs[1L])) # NA or a number (not NaN/Inf)
})

test_that("sirEvalOFV returns a vector with one entry per row", {
  skip_on_cran()
  mu <- .sirProposalMu(theoFit())
  # Three rows: true estimates, slight perturbations
  param_mat <- rbind(
    mu,
    mu + 0.01,
    mu - 0.01
  )
  dimnames(param_mat) <- list(NULL, names(mu))

  ofvs <- sirEvalOFV(theoFit(), param_mat)

  expect_length(ofvs, 3L)
  expect_true(all(!is.na(ofvs)))
})

test_that("sirEvalOFV returns the same OFVs with future workers", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  # Workers are fresh R processes that library() nlmixr2sir. Under
  # devtools::load_all() there is no installed copy for them to load -- or an
  # older one, which would test stale code -- so this only means something
  # against an installed package (R CMD check). An installed package has
  # Meta/package.rds; a source tree does not.
  skip_if_not(
    file.exists(file.path(
      getNamespaceInfo("nlmixr2sir", "path"), "Meta", "package.rds"
    )),
    "future workers need an installed nlmixr2sir"
  )
  plan_before <- future::plan()
  on.exit(future::plan(plan_before), add = TRUE)

  mu <- .sirProposalMu(theoFit())
  param_mat <- rbind(
    mu,
    mu + 0.01
  )
  dimnames(param_mat) <- list(NULL, names(mu))

  ofv_seq <- suppressMessages(sirEvalOFV(theoFit(), param_mat, workers = 1L))
  ofv_par <- suppressMessages(sirEvalOFV(theoFit(), param_mat, workers = 2L))

  expect_equal(ofv_par, ofv_seq, tolerance = 1e-8)
})

# rxode2::ini() evaluates the OMEGA line as lotri({...}) in the caller's
# environment. lotri is an Imports of rxode2 and nlmixr2est, so it is never
# attached; without nlmixr2sir importing it, every OFV evaluation returned NA
# in a session where the user had not also attached nlmixr2.
test_that("lotri is imported, so OFV evaluation works without library(nlmixr2)", {
  skip_on_cran()
  expect_true(
    exists(
      "lotri",
      envir = parent.env(asNamespace("nlmixr2sir")),
      inherits = FALSE
    )
  )
})

test_that("sirEvalOFV keeps the underlying error when evaluation fails", {
  skip_on_cran()
  nms <- .sirParamSpace(theoFit())$sirName
  bad <- matrix(
    NA_real_,
    nrow = 1L,
    ncol = length(nms),
    dimnames = list(NULL, nms)
  )
  res <- suppressMessages(sirEvalOFV(theoFit(), bad, workers = 1L))
  expect_true(is.na(res))
  expect_type(attr(res, "evalErrors"), "character")
})

test_that("sirEvalOFV handles diagonal omega columns for multi-ETA models", {
  skip_on_cran()
  fit <- threeEtaFit()
  proposal <- .sirInitialProposal(
    fit,
    .sirProposalMu(fit),
    sirGetProposalCov(fit),
    capCorrelation = 0.8
  )
  param_mat <- matrix(
    proposal$mu,
    nrow = 1L,
    dimnames = list(NULL, names(proposal$mu))
  )

  eta_cols <- c("eta.ka", "eta.cl", "eta.v")
  expect_setequal(intersect(eta_cols, colnames(param_mat)), eta_cols)

  ofvs <- sirEvalOFV(fit, param_mat, workers = 1L)

  expect_length(ofvs, 1L)
  expect_false(is.na(ofvs[1L]))
  expect_lt(abs(ofvs[1L] - fit$objf), 1)
})

test_that("sirEvalOFV OFV increases away from the estimates", {
  skip_on_cran()
  mu <- .sirProposalMu(theoFit())
  big_offset <- mu + 1 # large shift in all thetas
  param_mat <- rbind(mu, big_offset)
  dimnames(param_mat) <- list(NULL, names(mu))

  ofvs <- suppressWarnings(sirEvalOFV(theoFit(), param_mat))

  expect_true(!is.na(ofvs[1L]))
  # OFV at true values should be lower (better fit)
  if (!is.na(ofvs[2L])) expect_lt(ofvs[1L], ofvs[2L])
})

test_that("sirEvalOFV errors when no columns match fit parameters", {
  skip_on_cran()
  bad_mat <- matrix(1, nrow = 1L, dimnames = list(NULL, "not_a_param"))
  expect_error(sirEvalOFV(theoFit(), bad_mat), "match")
})

# Step 5: sirCalcWeights -------------------------------------------------------

# Fixed 3-parameter setup reused across Step 5 tests (same as Step 2)
.wt_mu <- c(tka = 0.45, tcl = 0.98, tv = 3.47)
.wt_cov <- matrix(
  c(0.065, -0.004, 0.004, -0.004, 0.006, -0.002, 0.004, -0.002, 0.003),
  nrow = 3L,
  dimnames = list(names(.wt_mu), names(.wt_mu))
)

test_that("sirCalcWeights: relPDF = 1 at the proposal mean", {
  skip_on_cran()
  samp <- matrix(.wt_mu, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = 0)
  expect_equal(res$relPDF[1L], 1, tolerance = 1e-12)
})

test_that("sirCalcWeights: prob_resample sums to 1", {
  skip_on_cran()
  set.seed(1)
  samp <- mvtnorm::rmvnorm(50L, mean = .wt_mu, sigma = .wt_cov)
  dofv <- rnorm(50L, mean = 0, sd = 2)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = dofv)
  expect_equal(sum(res$prob_resample), 1, tolerance = 1e-12)
})

test_that("sirCalcWeights: negative dOFV (better fit) gives IR > 1 at mu", {
  skip_on_cran()
  # At mu: relPDF = 1, so IR = exp(-0.5 * dOFV); dOFV < 0 → IR > 1
  samp <- matrix(.wt_mu, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = -2)
  expect_gt(res$importance_ratio[1L], 1)
})

test_that("sirCalcWeights: NA dOFV gives prob_resample = 0", {
  skip_on_cran()
  set.seed(2)
  samp <- mvtnorm::rmvnorm(5L, mean = .wt_mu, sigma = .wt_cov)
  dofv <- c(1, NA, 2, NA, 0.5)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = dofv)
  expect_equal(res$prob_resample[c(2L, 4L)], c(0, 0))
  expect_true(is.na(res$importance_ratio[2L]))
  expect_equal(sum(res$prob_resample), 1, tolerance = 1e-12)
})

test_that("sirCalcWeights: output has correct structure", {
  skip_on_cran()
  set.seed(3)
  samp <- mvtnorm::rmvnorm(10L, mean = .wt_mu, sigma = .wt_cov)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = rep(0, 10L))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10L)
  expect_named(
    res,
    c(
      "sample_id",
      "dOFV",
      "likelihood_ratio",
      "relPDF",
      "importance_ratio",
      "prob_resample"
    )
  )
})

test_that("sirCalcWeights: relPDF decreases away from mu", {
  skip_on_cran()
  far <- matrix(.wt_mu + 5, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  near <- matrix(.wt_mu + 0.01, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  samp <- rbind(far, near)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = c(0, 0))
  expect_lt(res$relPDF[1L], res$relPDF[2L])
})

test_that("sirCalcWeights: positive dOFV (worse fit) at mu gives IR < 1", {
  skip_on_cran()
  samp <- matrix(.wt_mu, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = 4)
  expect_lt(res$importance_ratio[1L], 1)
})

test_that("sirCalcWeights: errors on non-PD covMat", {
  skip_on_cran()
  bad_cov <- matrix(c(1, 2, 2, 1), 2L) # not PD
  samp <- matrix(c(1, 2), nrow = 1L)
  expect_error(
    sirCalcWeights(samp, c(1, 2), bad_cov, dOFV = 0),
    "positive definite"
  )
})

test_that("sirCalcWeights: log-space normalization handles extreme weights", {
  skip_on_cran()
  cov_mat <- diag(2)
  mu <- c(a = 0, b = 0)
  samp <- matrix(
    c(100, 100, 0, 0),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, names(mu))
  )
  res <- sirCalcWeights(samp, mu, cov_mat, dOFV = c(-1000, 1000))
  expect_equal(sum(res$prob_resample), 1, tolerance = 1e-12)
  expect_true(all(is.finite(res$prob_resample)))
  expect_gt(res$prob_resample[1L], 0.999)
})

# Step 6: sirResample ----------------------------------------------------------

# Helper: build a weights data frame with uniform probs for n samples
.uniform_weights <- function(n) {
  data.frame(prob_resample = rep(1 / n, n))
}

# Helper: build a weights data frame with all weight on row 1
.spike_weights <- function(n) {
  p <- c(1, rep(0, n - 1L))
  data.frame(prob_resample = p)
}

test_that("sirResample returns m rows", {
  skip_on_cran()
  set.seed(1)
  samp <- matrix(rnorm(30L), nrow = 10L)
  res <- sirResample(samp, .uniform_weights(10L), m = 5L)
  expect_equal(nrow(res$samples), 5L)
  expect_equal(ncol(res$samples), 3L)
})

test_that("sirResample resampleCounts sums to m", {
  skip_on_cran()
  set.seed(2)
  samp <- matrix(rnorm(30L), nrow = 10L)
  res <- sirResample(samp, .uniform_weights(10L), m = 6L)
  expect_equal(sum(res$resampleCounts), 6L)
  expect_length(res$resampleCounts, 10L)
})

test_that("sirResample with capped replacement can repeatedly select row 1", {
  skip_on_cran()
  samp <- matrix(1:20, nrow = 5L)
  res <- sirResample(samp, .spike_weights(5L), m = 5L, capResampling = 5L)
  # Every resampled row should equal samp[1, ]
  expect_true(all(sweep(res$samples, 2L, samp[1L, ], "==") == 1L))
  expect_equal(res$resampleCounts[1L], 5L)
  expect_true(all(res$resampleCounts[-1L] == 0L))
})

test_that("sirResample cap is respected: no sample exceeds capResampling", {
  skip_on_cran()
  set.seed(3)
  samp <- matrix(rnorm(50L), nrow = 10L)
  cap <- 3L
  res <- sirResample(samp, .uniform_weights(10L), m = 20L, capResampling = cap)
  expect_true(all(res$resampleCounts <= cap))
  expect_equal(sum(res$resampleCounts), 20L)
  expect_equal(nrow(res$samples), 20L)
})

test_that("sirResample cap = 1 samples without replacement", {
  skip_on_cran()
  set.seed(4)
  samp <- matrix(rnorm(30L), nrow = 10L)
  res <- sirResample(samp, .uniform_weights(10L), m = 10L, capResampling = 1)
  expect_equal(nrow(res$samples), 10L)
  expect_true(all(res$resampleCounts <= 1L))
  expect_equal(sum(res$resampleCounts), 10L)
})

test_that("sirResample with cap forces spread when one sample dominates", {
  skip_on_cran()
  set.seed(5)
  samp <- matrix(rnorm(50L), nrow = 10L)
  # Uniform weights + cap = 3: row 1 cannot be selected more than 3 times
  res <- sirResample(samp, .uniform_weights(10L), m = 20L, capResampling = 3L)
  expect_lte(res$resampleCounts[1L], 3L)
  expect_equal(sum(res$resampleCounts), 20L)
  expect_equal(nrow(res$samples), 20L)
})

test_that("sirResample selected rows all come from original samples", {
  skip_on_cran()
  set.seed(6)
  samp <- matrix(seq_len(20L), nrow = 5L)
  res <- sirResample(samp, .uniform_weights(5L), m = 5L)
  # Every resampled row must exactly match one of the original rows
  for (i in seq_len(nrow(res$samples))) {
    row_match <- apply(samp, 1L, function(r) identical(r, res$samples[i, ]))
    expect_true(any(row_match))
  }
})

test_that("sirResample returns sample-order metadata", {
  skip_on_cran()
  set.seed(7)
  samp <- matrix(seq_len(20L), nrow = 5L)
  res <- sirResample(samp, .uniform_weights(5L), m = 5L)
  expect_length(res$sampleOrder, 5L)
  expect_length(res$selectionOrder, 5L)
  expect_equal(sum(res$resampleCounts), 5L)
})

# Step 7: sirBoxCox / sirBoxCoxInverse -----------------------------------------

test_that("sirBoxCox with lambda = 1 is a linear (affine) transformation of x", {
  skip_on_cran()
  x <- c(1, 2, 4, 8, 16)
  res <- sirBoxCox(x, lambda = 1)
  # transformed = (x + delta)^1 - 1 = x + delta - 1, linear in x
  expect_equal(res$lambda, 1)
  expect_equal(cor(res$transformed, x), 1, tolerance = 1e-12)
})

test_that("sirBoxCox with lambda = 0 returns log(x + delta)", {
  skip_on_cran()
  x <- c(1, 2, 4, 8, 16)
  delta <- abs(min(x)) + 1e-6
  res <- sirBoxCox(x, lambda = 0)
  expect_equal(res$transformed, log(x + delta), tolerance = 1e-12)
})

test_that("sirBoxCoxInverse recovers x after sirBoxCox (lambda estimated)", {
  skip_on_cran()
  set.seed(1)
  x <- exp(rnorm(50L)) # log-normal: should yield lambda ≈ 0
  res <- sirBoxCox(x)
  x_back <- sirBoxCoxInverse(res$transformed, res$lambda, res$delta)
  expect_equal(x_back, x, tolerance = 1e-8)
})

test_that("sirBoxCoxInverse recovers x for fixed lambda = 1", {
  skip_on_cran()
  x <- c(0.5, 1, 2, 5, 10)
  res <- sirBoxCox(x, lambda = 1)
  x_back <- sirBoxCoxInverse(res$transformed, res$lambda, res$delta)
  expect_equal(x_back, x, tolerance = 1e-10)
})

test_that("sirBoxCoxInverse recovers x for fixed lambda = 0 (log case)", {
  skip_on_cran()
  x <- c(0.1, 0.5, 1, 5, 20)
  res <- sirBoxCox(x, lambda = 0)
  x_back <- sirBoxCoxInverse(res$transformed, res$lambda, res$delta)
  expect_equal(x_back, x, tolerance = 1e-10)
})

test_that("sirBoxCox estimated lambda improves normality vs untransformed", {
  skip_on_cran()
  set.seed(2)
  x <- rexp(100L, rate = 2) # clearly right-skewed
  res <- sirBoxCox(x)
  # Shapiro-Wilk p should be higher after transformation
  p_raw <- shapiro.test(x)$p.value
  p_bc <- shapiro.test(res$transformed)$p.value
  expect_gt(p_bc, p_raw)
})

test_that("sirBoxCox returns correct list structure", {
  skip_on_cran()
  x <- rnorm(20L) + 5
  res <- sirBoxCox(x)
  expect_named(res, c("transformed", "lambda", "delta"))
  expect_length(res$transformed, 20L)
  expect_length(res$lambda, 1L)
  expect_length(res$delta, 1L)
})

test_that("sirBoxCox delta defaults to |min(x)| + 1e-6", {
  skip_on_cran()
  x <- c(-3, 0, 1, 5)
  res <- sirBoxCox(x, lambda = 1)
  expected_delta <- abs(min(x)) + 1e-6
  expect_equal(res$delta, expected_delta, tolerance = 1e-12)
})

test_that("sirBoxCox supplied delta is respected", {
  skip_on_cran()
  x <- c(1, 2, 3)
  res <- sirBoxCox(x, lambda = 1, delta = 0.5)
  expect_equal(res$delta, 0.5)
})

test_that("sirBoxCoxInverse errors when back-transform base is non-positive", {
  skip_on_cran()
  # lambda = 2, x_transformed = -1 → base = 2*(-1) + 1 = -1 ≤ 0
  expect_error(sirBoxCoxInverse(-1, lambda = 2, delta = 0), "positive")
})

# Step 8: sirUpdateProposal ----------------------------------------------------

test_that("sirUpdateProposal covMat is symmetric", {
  skip_on_cran()
  set.seed(1)
  mat <- matrix(
    rexp(100L),
    nrow = 20L,
    dimnames = list(NULL, c("a", "b", "c", "d", "e"))
  )
  res <- sirUpdateProposal(mat)
  expect_equal(res$covMat, t(res$covMat), tolerance = 1e-12)
})

test_that("sirUpdateProposal covMat is positive semi-definite", {
  skip_on_cran()
  set.seed(2)
  mat <- matrix(rexp(60L), nrow = 20L, dimnames = list(NULL, c("a", "b", "c")))
  res <- sirUpdateProposal(mat)
  eigs <- eigen(res$covMat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("sirUpdateProposal preserves column names", {
  skip_on_cran()
  set.seed(3)
  nms <- c("tka", "tcl", "tv")
  mat <- matrix(rnorm(60L) + 5, nrow = 20L, dimnames = list(NULL, nms))
  res <- sirUpdateProposal(mat)
  expect_equal(rownames(res$covMat), nms)
  expect_equal(colnames(res$covMat), nms)
})

test_that("sirUpdateProposal boxcoxParams has correct structure when boxcox = TRUE", {
  skip_on_cran()
  set.seed(4)
  nms <- c("p1", "p2")
  mat <- matrix(rexp(40L), nrow = 20L, dimnames = list(NULL, nms))
  res <- sirUpdateProposal(mat, boxcox = TRUE)
  expect_s3_class(res$boxcoxParams, "data.frame")
  expect_named(res$boxcoxParams, c("param", "lambda", "delta"))
  expect_equal(nrow(res$boxcoxParams), 2L)
  expect_equal(res$boxcoxParams$param, nms)
})

test_that("sirUpdateProposal boxcoxParams is NULL when boxcox = FALSE", {
  skip_on_cran()
  set.seed(5)
  mat <- matrix(rnorm(40L), nrow = 20L)
  res <- sirUpdateProposal(mat, boxcox = FALSE)
  expect_null(res$boxcoxParams)
})

test_that("sirUpdateProposal boxcox = FALSE matches cov() directly", {
  skip_on_cran()
  set.seed(6)
  mat <- matrix(rnorm(60L), nrow = 20L, dimnames = list(NULL, c("a", "b", "c")))
  res <- sirUpdateProposal(mat, boxcox = FALSE)
  expected <- cov(mat)
  expect_equal(res$covMat, expected, tolerance = 1e-12)
})

test_that("sirUpdateProposal lambdas are retrievable and finite", {
  skip_on_cran()
  set.seed(7)
  mat <- matrix(
    rexp(80L),
    nrow = 20L,
    dimnames = list(NULL, c("a", "b", "c", "d"))
  )
  res <- sirUpdateProposal(mat, boxcox = TRUE)
  expect_true(all(is.finite(res$boxcoxParams$lambda)))
  expect_true(all(res$boxcoxParams$lambda >= -3 & res$boxcoxParams$lambda <= 3))
})

# Step 9: sirRunIteration ------------------------------------------------------
# These tests require OFV evaluation so use tiny nSamples (8) for speed.

test_that("sirRunIteration returns correct list structure", {
  skip_on_cran()
  expect_named(
    iter1(),
    c(
      "resampledMat",
      "newMu",
      "newCov",
      "proposalRepair",
      "referenceOfv",
      "newReferenceOfv",
      "iterSummary",
      "boxcoxState",
      "rawResults"
    )
  )
})

test_that("sirRunIteration resampledMat has nResample rows", {
  skip_on_cran()
  expect_equal(nrow(iter1()$resampledMat), 8L)
})

test_that("sirRunIteration iterSummary has expected columns and values", {
  skip_on_cran()
  s <- iter1()$iterSummary
  expect_named(
    s,
    c(
      "iter",
      "nSamples",
      "nAttempted",
      "nDrawAttempts",
      "nCollected",
      "nSuccessful",
      "nFailed",
      "nResample",
      "nResampled",
      "ess",
      "essFraction",
      "maxWeight",
      "perplexity",
      "nNonNegligible",
      "posDefAdjusted",
      "proposalRepaired",
      "proposalRepairMagnitude",
      "minDOFV",
      "meanDOFV",
      "nNegativeDOFV",
      "thetaRejected",
      "omegaRejected",
      "sigmaRejected",
      "inverseRejected"
    )
  )
  expect_equal(s$iter, 1L)
  expect_equal(s$nSamples, 16L)
  expect_equal(s$nAttempted, 16L)
  expect_true(s$nCollected <= 16L && s$nCollected >= 1L)
  expect_gte(s$nFailed, 0L)
  expect_true(is.numeric(s$minDOFV))
  expect_gte(s$sigmaRejected, 0L)
})

test_that("sirRunIteration newMu is a named numeric vector matching proposal params", {
  skip_on_cran()
  mu_names <- colnames(iter1()$resampledMat)
  expect_named(iter1()$newMu, mu_names)
  expect_true(is.numeric(iter1()$newMu))
})

test_that("sirRunIteration newCov is symmetric and positive semi-definite", {
  skip_on_cran()
  m <- iter1()$newCov
  expect_equal(m, t(m), tolerance = 1e-12)
  eigs <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("sirRunIteration boxcoxState has correct structure when boxcox = TRUE", {
  skip_on_cran()
  bc <- iter1()$boxcoxState
  expect_s3_class(bc, "data.frame")
  expect_named(bc, c("param", "lambda", "delta"))
  expect_equal(bc$param, colnames(iter1()$resampledMat))
})

test_that("sirRunIteration rawResults has dOFV column", {
  skip_on_cran()
  expect_true("dOFV" %in% names(iter1()$rawResults))
  expect_equal(nrow(iter1()$rawResults), iter1()$iterSummary$nCollected + 1L)
  expect_equal(iter1()$rawResults$sample_id[[1L]], 0L)
  expect_equal(iter1()$rawResults$dOFV[[1L]], 0)
})

test_that("sirRunIteration rawResults has PsN-like SIR metadata", {
  skip_on_cran()
  expect_true(all(
    c(
      "sample_id",
      "likelihood_ratio",
      "probability_resample",
      "resamples",
      "sample_order"
    ) %in%
      names(iter1()$rawResults)
  ))
  expect_equal(
    sum(iter1()$rawResults$resamples),
    iter1()$iterSummary$nResampled
  )
})

test_that(".sirBuildRawResults expands capped raw-result rows", {
  skip_on_cran()
  param_mat <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )
  weights <- data.frame(
    likelihood_ratio = c(1, 2, 3),
    relPDF = 1,
    importance_ratio = c(1, 2, 3),
    prob_resample = c(0.2, 0.3, 0.5)
  )
  resampled <- list(
    resampleCounts = c(2L, 0L, 1L),
    sampleOrder = c("2;4", "", "1")
  )
  raw <- nlmixr2sir:::.sirBuildRawResults(
    paramMat = param_mat,
    weights = weights,
    dOFV = c(0.1, 0.2, 0.3),
    resampled = resampled,
    mu = c(a = 0, b = 0),
    capResampling = 2
  )
  expect_equal(nrow(raw), 7L)
  expect_equal(raw$sample_id[[1L]], 0L)
  expect_equal(raw$resamples[raw$sample_id == 1L], c(1L, 1L))
  expect_equal(raw$sample_order[raw$sample_id == 1L], c("2", "4"))
  expect_equal(raw$resamples[raw$sample_id == 2L], c(0L, 0L))
})

test_that("sirRunIteration recentering: newMu shifts when a better sample exists", {
  skip_on_cran()
  # Use a mu perturbed away from estimates so samples near true values get dOFV < 0
  mu_perturbed <- .sirProposalMu(theoFit()) + 0.5
  prop_cov <- sirGetProposalCov(theoFit())
  set.seed(1)
  res <- .sirQuiet(
    sirRunIteration(
      theoFit(),
      mu = mu_perturbed,
      proposalCov = prop_cov,
      nSamples = 16L,
      nResample = 8L,
      iterNum = 1L,
      recenter = TRUE,
      boxcox = FALSE,
      directory = NULL
    )
  )
  # newMu should differ from mu_perturbed if any dOFV < 0 was found
  min_dofv <- res$iterSummary$minDOFV
  if (!is.na(min_dofv) && min_dofv < 0) {
    expect_false(isTRUE(all.equal(res$newMu, mu_perturbed, tolerance = 1e-6)))
  } else {
    expect_equal(res$newMu[names(mu_perturbed)], mu_perturbed, tolerance = 1e-6)
  }
})

test_that("sirRunIteration keeps raw results in memory when directory is provided", {
  skip_on_cran()
  mu <- .sirProposalMu(theoFit())
  prop_cov <- sirGetProposalCov(theoFit())
  tmp_dir <- tempfile("sir_test_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  set.seed(9)
  res <- .sirQuiet(
    sirRunIteration(
      theoFit(),
      mu = mu,
      proposalCov = prop_cov,
      nSamples = 16L,
      nResample = 8L,
      iterNum = 2L,
      directory = tmp_dir
    )
  )
  expect_s3_class(res$rawResults, "data.frame")
  expect_true("dOFV" %in% names(res$rawResults))
  expect_false(file.exists(file.path(
    tmp_dir,
    "raw_results_sir_iteration2.csv"
  )))
})

test_that("sirRunIteration chained: iter 2 accepts boxcoxState from iter 1", {
  skip_on_cran()
  set.seed(77)
  res2 <- .sirQuiet(
    sirRunIteration(
      theoFit(),
      mu = iter1()$newMu,
      proposalCov = iter1()$newCov,
      nSamples = 16L,
      nResample = 8L,
      iterNum = 2L,
      boxcoxState = iter1()$boxcoxState
    )
  )
  expect_named(
    res2,
    c(
      "resampledMat",
      "newMu",
      "newCov",
      "proposalRepair",
      "referenceOfv",
      "newReferenceOfv",
      "iterSummary",
      "boxcoxState",
      "rawResults"
    )
  )
  expect_equal(res2$iterSummary$iter, 2L)
})

test_that("sirRunIteration new proposal includes omega and sigma columns", {
  skip_on_cran()
  expect_true("eta.ka" %in% rownames(iter1()$newCov))
  expect_true("add.sd" %in% rownames(iter1()$newCov))
  expect_true("eta.ka" %in% iter1()$boxcoxState$param)
  expect_true("add.sd" %in% iter1()$boxcoxState$param)
})

test_that("sir initial proposal applies omega and sigma inflation", {
  skip_on_cran()
  mu <- .sirProposalMu(theoFit())
  prop_cov <- sirGetProposalCov(theoFit())
  base <- nlmixr2sir:::.sirInitialProposal(
    theoFit(),
    mu,
    prop_cov,
    capCorrelation = 1
  )
  inflated <- nlmixr2sir:::.sirInitialProposal(
    theoFit(),
    mu,
    prop_cov,
    omegaInflation = 4,
    sigmaInflation = 9,
    capCorrelation = 1
  )
  expect_equal(
    inflated$covMat["eta.ka", "eta.ka"],
    4 * base$covMat["eta.ka", "eta.ka"],
    tolerance = 1e-8
  )
  expect_equal(
    inflated$covMat["add.sd", "add.sd"],
    9 * base$covMat["add.sd", "add.sd"],
    tolerance = 1e-8
  )
})

# Gap fix: .sirParamSpace fixed-parameter handling, .sirFallbackSe,
# .sirReconstructOmega ---------------------------------------------------------

# Minimal stand-in for a fit, so fixed-parameter handling can be tested
# without a model that takes a minute to converge.
.fakeFit <- function(fix) {
  # A plain list, deliberately unclassed: nlmixr2FitCore defines its own `$`
  # method that a list stand-in cannot satisfy.
  list(
    nsub = 20L,
    cov = NULL,
    theta = c(tka = 0.4),
    omega = matrix(
      c(0.4, 0.05, 0.05, 0.2),
      2L,
      2L,
      dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
    ),
    iniDf = data.frame(
      ntheta = c(1L, NA, NA, NA),
      neta1 = c(NA, 1L, 2L, 2L),
      neta2 = c(NA, 1L, 1L, 2L),
      name = c("tka", "eta.ka", "(eta.ka,eta.cl)", "eta.cl"),
      lower = c(-Inf, -Inf, -Inf, -Inf),
      est = c(0.4, 0.4, 0.05, 0.2),
      upper = c(Inf, Inf, Inf, Inf),
      fix = c(FALSE, fix),
      err = c(NA_character_, NA, NA, NA),
      stringsAsFactors = FALSE
    )
  )
}

test_that(".sirParamSpace excludes fixed omega elements", {
  skip_on_cran()
  ps <- .sirParamSpace(.fakeFit(c(FALSE, TRUE, TRUE)))
  om <- ps[ps$kind %in% c("omegaDiag", "omegaOffdiag"), ]
  expect_equal(om$sirName, "eta.ka")
  expect_equal(om$neta1, 1L)
  expect_equal(om$neta2, 1L)
})

test_that(".sirParamSpace keeps only THETA when every omega is fixed", {
  skip_on_cran()
  ps <- .sirParamSpace(.fakeFit(c(TRUE, TRUE, TRUE)))
  expect_equal(ps$sirName, "tka")
  expect_equal(ps$kind, "theta")
})

test_that(".sirParamSpace marks covName NA for a fit with no cov at all", {
  skip_on_cran()
  ps <- .sirParamSpace(.fakeFit(c(FALSE, FALSE, FALSE)))
  expect_true(all(is.na(ps$covName)))
  expect_equal(nrow(ps), 4L)
})

test_that(".sirFallbackSe follows the Wishart formula for omega", {
  skip_on_cran()
  fit <- .fakeFit(c(FALSE, FALSE, FALSE))
  ps <- .sirParamSpace(fit)
  se <- .sirFallbackSe(fit, ps)
  df <- 20L - 1L
  om <- fit$omega
  expect_equal(se[["eta.ka"]], sqrt(2 * om[1, 1]^2 / df), tolerance = 1e-12)
  expect_equal(se[["eta.cl"]], sqrt(2 * om[2, 2]^2 / df), tolerance = 1e-12)
  expect_equal(
    se[["eta.cl:eta.ka"]],
    sqrt((om[1, 1] * om[2, 2] + om[2, 1]^2) / df),
    tolerance = 1e-12
  )
})

test_that(".sirFallbackSe honours omegaDf", {
  skip_on_cran()
  fit <- .fakeFit(c(FALSE, FALSE, FALSE))
  se <- .sirFallbackSe(fit, omegaDf = 100)
  expect_equal(
    se[["eta.ka"]],
    sqrt(2 * fit$omega[1, 1]^2 / 100),
    tolerance = 1e-12
  )
})

test_that(".sirFallbackSe uses sigmaFallbackRse for residual error", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  se <- .sirFallbackSe(fit, ps, sigmaFallbackRse = 50)
  expect_gt(se[["add.sd"]], 0)
})

test_that(".sirReconstructOmega sets the diagonal from a named vector", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  mat <- .sirReconstructOmega(ps, c("eta.ka" = 0.99), theoFit()$omega)
  expect_equal(mat["eta.ka", "eta.ka"], 0.99, tolerance = 1e-12)
})

test_that(".sirReconstructOmega leaves unspecified elements unchanged", {
  skip_on_cran()
  fit <- .fakeFit(c(FALSE, FALSE, FALSE))
  ps <- .sirParamSpace(fit)
  mat <- .sirReconstructOmega(
    ps,
    c("eta.ka" = 0.8, "eta.cl:eta.ka" = 0.05),
    fit$omega
  )
  expect_equal(mat[1L, 1L], 0.8, tolerance = 1e-12)
  expect_equal(mat[2L, 1L], 0.05, tolerance = 1e-12)
  expect_equal(mat[1L, 2L], 0.05, tolerance = 1e-12)
  expect_equal(mat[2L, 2L], 0.2, tolerance = 1e-12) # unchanged
})

test_that("sirRunIteration rawResults contains sigma column (add.sd)", {
  skip_on_cran()
  expect_true("add.sd" %in% names(iter1()$rawResults))
})

test_that("sirRunIteration rawResults contains omega column (eta.ka)", {
  skip_on_cran()
  expect_true("eta.ka" %in% names(iter1()$rawResults))
})

test_that("sirRunIteration resampledMat contains sigma column (add.sd)", {
  skip_on_cran()
  expect_true("add.sd" %in% colnames(iter1()$resampledMat))
})

# Step 11: sirSummary ----------------------------------------------------------

test_that("sirSummary returns final SIR summary statistics", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  expect_s3_class(s, "data.frame")
  # PsN's percentile set, from prediction intervals 0, 40, 80, 90, 95.
  expect_named(
    s,
    c(
      "param",
      "estimate",
      "mean",
      "sd",
      "rse",
      "rse_sd_scale",
      "p2.5",
      "p5",
      "p10",
      "p30",
      "p50",
      "p70",
      "p90",
      "p95",
      "p97.5"
    )
  )
  expect_equal(s$param, colnames(iter1()$resampledMat))
  expect_equal(nrow(s), ncol(iter1()$resampledMat))
})

test_that("sirSummary computes empirical SD, RSE, and percentiles", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  first_param <- s$param[[1L]]
  x <- iter1()$resampledMat[, first_param]
  expect_equal(
    s$sd[s$param == first_param],
    stats::sd(x),
    tolerance = 1e-12
  )
  expect_equal(
    s$p50[s$param == first_param],
    unname(stats::quantile(x, probs = 0.5)),
    tolerance = 1e-12
  )
  expect_equal(
    s$rse[s$param == first_param],
    stats::sd(x) / abs(s$estimate[s$param == first_param]) * 100,
    tolerance = 1e-12
  )
})

test_that("sirSummary attaches empirical covariance and correlation matrices", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  expect_equal(
    attr(s, "covMatrix"),
    stats::cov(iter1()$resampledMat),
    tolerance = 1e-12
  )
  expect_equal(
    attr(s, "corMatrix"),
    stats::cor(iter1()$resampledMat),
    tolerance = 1e-12
  )
})

# Step 10: runSIR --------------------------------------------------------------

# The sample-count adjustment is now checked against PsN's own oracle values
# in test-sir-psn-parity.R. The test that stood here asserted 112 for a 10%
# loss and no adjustment at exactly 0.95 turnout -- neither of which is what
# PsN does -- so it was superseded rather than corrected.

test_that("runSIR runs end-to-end and writes Step 10 artifacts", {
  skip_on_cran()
  tmp_dir <- tempfile("sir_run_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  set.seed(20260420)
  res <- .sirQuiet(
    runSIR(
      theoFit(),
      nSamples = c(16L, 16L),
      nResample = c(8L, 8L),
      directory = tmp_dir,
      control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L, boxcox = TRUE)
    )
  )

  expect_s3_class(res, "nlmixr2SIR")
  expect_s3_class(res, "data.frame")
  expect_true(all(
    c("param", "estimate", "sd", "rse", "p2.5", "p97.5") %in%
      names(res)
  ))
  expect_true(file.exists(file.path(tmp_dir, "raw_results.csv")))
  expect_true(file.exists(file.path(tmp_dir, "raw_results.rds")))
  expect_true(file.exists(file.path(tmp_dir, "raw_results_header.json")))
  expect_true(file.exists(file.path(tmp_dir, "summary_iterations.csv")))
  expect_true(file.exists(file.path(tmp_dir, "sample_rejection_summary.txt")))
  expect_true(file.exists(file.path(tmp_dir, "sir_results.csv")))
  expect_true(file.exists(file.path(tmp_dir, "sir_state.rds")))
  expect_true(file.exists(file.path(tmp_dir, "sir_seed.rds")))

  iter_summary <- attr(res, "iterationSummary")
  expect_equal(nrow(iter_summary), 2L)
  expect_equal(iter_summary$nSamples, c(16L, 16L))
  expect_equal(iter_summary$nResample, c(8L, 8L))

  iterations <- attr(res, "iterations")
  expect_null(iterations[[2L]]$boxcoxState)

  raw <- nlmixr2utils::readRawResults(tmp_dir)
  resampled <- attr(res, "resampledMat")
  expect_equal(raw$sample[[1L]], 0L)
  expect_equal(raw$role[[1L]], "reference")
  expect_equal(raw$source, rep("sir", nrow(raw)))
  expect_equal(sum(raw$role == "sample"), nrow(resampled))

  parsed <- nlmixr2utils::parseRawResultsParams(
    raw,
    theoFit()
  )
  expect_equal(length(parsed), nrow(resampled))
  expect_equal(names(parsed), paste0("sample_", seq_len(nrow(resampled))))
})

test_that("runSIRControl rejects invalid workers before running", {
  skip_on_cran()
  expect_error(runSIRControl(objfStencil = FALSE, workers = 0L), "workers")
})

test_that("runSIR rejects settings passed outside the control object", {
  skip_on_cran()
  expect_error(
    runSIR(theoFit(), nSamples = 1L, nResample = 1L, workers = 1L),
    "runSIRControl"
  )
})

test_that("runSIR requires a control object of the right class", {
  skip_on_cran()
  expect_error(
    runSIR(theoFit(), nSamples = 1L, nResample = 1L, control = list()),
    "runSIRControl"
  )
})

# Step 12: S3 print and plot methods ------------------------------------------

test_that("print.nlmixr2SIR returns object invisibly and prints tables", {
  skip_on_cran()
  printed <- utils::capture.output(ret <- print(sirObj()))
  expect_identical(ret, sirObj())
  expect_match(paste(printed, collapse = "\n"), "param")
  expect_match(paste(printed, collapse = "\n"), "estimate")
  expect_match(paste(printed, collapse = "\n"), "nAttempted")
})

test_that("plot.nlmixr2SIR returns parameter distribution plot", {
  skip_on_cran()
  p <- plot(sirObj())
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Parameter value")
})

test_that("plot.nlmixr2SIR returns dOFV diagnostic plot", {
  skip_on_cran()
  p <- plot(sirObj(), type = "dofv", bins = 10L)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "dOFV")
})

test_that("plot.nlmixr2SIR returns resampling diagnostic plot", {
  skip_on_cran()
  p <- plot(sirObj(), type = "resampling")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "Probability resample")
})

test_that("sirSummary reports rse_sd_scale only for OMEGA diagonals", {
  skip_on_cran()
  fit <- blockFit()
  ps <- .sirParamSpace(fit)
  # blockFit has 7 estimated parameters, so 8 retained vectors would give a
  # covariance of rank exactly 7 -- no margin for a repeated draw.
  it <- .sirQuiet(sirRunIteration(
    fit,
    mu = .sirProposalMu(fit),
    proposalCov = sirGetProposalCov(fit),
    nSamples = 24L,
    nResample = 12L,
    iterNum = 1L,
    recenter = FALSE,
    boxcox = FALSE,
    directory = NULL
  ))
  s <- sirSummary(it$resampledMat, fit)
  kind <- ps$kind[match(s$param, ps$sirName)]

  # A variance has an SD-scale counterpart; RSE(sqrt(v)) ~= RSE(v) / 2.
  diag_rows <- which(kind == "omegaDiag")
  expect_gt(length(diag_rows), 0L)
  expect_equal(
    s$rse_sd_scale[diag_rows],
    s$rse[diag_rows] / 2,
    tolerance = 1e-12
  )

  # An OMEGA off-diagonal is a covariance: it may be negative or zero, and has
  # no square root, so there is no SD scale to convert to.
  off_rows <- which(kind == "omegaOffdiag")
  expect_gt(length(off_rows), 0L)
  expect_true(all(is.na(s$rse_sd_scale[off_rows])))

  # THETAs and residual-error SDs are already on their reported scale.
  other_rows <- which(!kind %in% c("omegaDiag", "omegaOffdiag"))
  expect_true(all(is.na(s$rse_sd_scale[other_rows])))
})
