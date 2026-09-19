# nlmixr2sir

`nlmixr2sir` provides Sampling Importance Resampling (SIR) for `nlmixr2`
population PK/PD models.

SIR is a simulation-based approach for parameter uncertainty estimation. Rather
than relying only on the asymptotic covariance matrix around the maximum
likelihood estimate, it repeatedly:

1. Samples parameter vectors from a proposal distribution centered on the MLE.
2. Evaluates each sampled vector by fixing parameters and recomputing the OFV.
3. Weights each sample by the ratio of target likelihood to proposal density.
4. Resamples according to those importance weights.
5. Updates the proposal from the empirical covariance of the resampled vectors.

The result is a set of empirical draws from the parameter uncertainty
distribution that can be summarized with nonparametric intervals, covariance
matrices, and diagnostic plots.

## The details

`runSIR()` implements the SIR workflow for `nlmixr2` Where practical, it follows
the same process as Perl-speaks-NONMEM.

`nlmixr2sir` builds the initial proposal from `fit$cov` by default, with
optional inflation and correlation capping, and can take it instead from
relative standard errors, a supplied covariance matrix, or the parameter
vectors in a raw-results file -- see *Requirements and Practical Notes*. THETA,
OMEGA, and residual-error parameters are sampled, and parameter-space
constraints such as bounds and positive-definite OMEGA matrices are enforced.

Since `nlmixr2est` 7, `foceiControl(covFull = TRUE)` is the default and `fit$cov`
covers THETA, residual error, and OMEGA together. `nlmixr2sir` uses that matrix
directly (`omegaFallback = "cov"`, the default), which means the initial
proposal carries the **correlations between THETA and OMEGA** rather than
treating the two as independent blocks.

When `fit$cov` is present but does not carry OMEGA -- a fit run with
`covFull = FALSE`, or a partial covariance -- `nlmixr2sir` fills the missing
block automatically with a Wishart-style approximation, and
`omegaFallback = "wishart"` forces that route even when OMEGA is available.

This automatic fallback completes an **incomplete** covariance. It cannot
replace an **absent** one. A fit with `covMethod = ""`, or whose covariance step
failed, has `fit$cov == NULL` and carries no THETA uncertainty at all; the
Wishart approximation derives OMEGA uncertainty from the OMEGA estimates and
the subject count, but nothing in the fit supplies THETA. `runSIR()` therefore
stops and asks for one of `rseTheta`, `covmatInput`, or `rawresInput` -- see
[Requirements and Practical Notes](#requirements-and-practical-notes) below. Inventing THETA uncertainty
from a default assumed RSE would fabricate exactly the quantity SIR is there to
measure. The fallback takes the free lower-triangular OMEGA elements and, with
`omegaDf = nSubjects - 1` by default, approximates diagonal SEs as
`sqrt(2 * omega^2 / df)` and off-diagonal SEs as
`sqrt((omega[i, i] * omega[j, j] + omega[i, j]^2) / df)`. That route gives a
block-diagonal proposal, so it discards the THETA-OMEGA correlations the
default route keeps. The route actually taken is reported in the run log.

Only the free lower-triangular elements are sampled directly. Each proposed
vector is reconstructed into an OMEGA matrix, and non-positive-definite draws
are discarded. After each SIR iteration, the next proposal covariance is
updated from the empirical covariance of the retained samples.

Sampled vectors are re-evaluated by fixing the population parameters and
recomputing the objective function against the data, without any estimation
being performed. Importance ratios are computed, weighted resampling is
performed, and the proposal for the next iteration is updated.

The SIR tool supports recentering, Box-Cox proposal updates, recovery from
saved state, extending a finished run with further iterations, iteration
summaries, and diagnostic plots.

The package is designed to work alongside `nlmixr2utils`, which provides the
shared worker-plan helpers and core infrastructure.

## Installation

The package is not on CRAN. Install it from GitHub together with
`nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "nlmixr2/nlmixr2utils",
  "nlmixr2/nlmixr2sir"
))
```

Using `remotes`:

```r
remotes::install_github("nlmixr2/nlmixr2utils")
remotes::install_github("nlmixr2/nlmixr2sir")
```

## Basic Use

```r
library(nlmixr2)
library(nlmixr2sir)

one_cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.00
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

fit <- nlmixr2(
  one_cmt,
  data = nlmixr2data::theo_sd,
  est = "focei",
  control = list(print = 0L),
  table = list(npde = TRUE, cwres = TRUE)
)

sir <- runSIR(
  fit,
  nSamples = c(1000L, 1000L, 1000L, 2000L, 2000L),
  nResample = c(200L, 400L, 500L, 1000L, 1000L)
)

print(sir)
plot(sir, type = "parameters")
```

Everything that tunes *how* the run behaves lives in `runSIRControl()`:

```r
sir <- runSIR(
  fit,
  nSamples = c(1000L, 1000L, 1000L, 2000L, 2000L),
  nResample = c(200L, 400L, 500L, 1000L, 1000L),
  control = runSIRControl(
    thetaInflation = 2,
    workers = 4,
    rxThreads = 2
  )
)
```

## Diagnostics

```r
plot(sir, type = "convergence")   # dOFV vs reference chi-square, per iteration
plot(sir, type = "intervals")     # proposal vs SIR interval, per parameter
plot(sir, type = "rsecor")        # RSE / correlation, with CI asymmetry
plot(sir, type = "parameters")    # resampled parameter distributions
```

`type = "convergence"` is the most informative diagnostic. For each iteration
it draws the empirical dOFV quantile curve for the proposal and for the
retained SIR distribution against a reference chi-square on the number of
estimated parameters. Convergence reads as the SIR curve settling onto the
reference. If the first iteration's proposal falls below the reference for more
than a quarter of the quantiles, **this plot** warns: the proposal is too
narrow, and resampling cannot recover from that -- restart with inflation. The
check runs when the convergence plot is drawn, not during `runSIR()`, so a run
you never plot will not raise it.

The chi-square reference rests on regular likelihood asymptotics, so read it as
evidence rather than a certificate. It can mislead for variance components on a
boundary, weakly identified or multimodal parameters, non-smooth likelihoods,
or a small subject count -- and percentile intervals from a likelihood-weighted
sample do not automatically have nominal frequentist coverage. The conditions
are set out in the [technical reference](docs/sir-technical-reference.md).

Each iteration also reports importance-weight diagnostics -- effective sample
size (ESS), its fraction of the usable samples, the largest single weight, and
perplexity. `runSIR()` warns on three separate conditions, because they mean
different things and have different remedies:

- **ESS below 100.** The retained distribution rests on about that many
  effectively independent points, whatever `nResample` says, and percentile
  intervals drawn from so few are dominated by resampling noise. ESS grows
  roughly in proportion to `nSamples`, so more samples genuinely help here.
- **Efficiency (ESS/n) below 10%.** Most draws carry negligible weight. This
  ratio is a property of the *proposal*, not of the sample size: asymptotically
  it converges to a constant fixed by the proposal-target mismatch, so drawing
  more samples raises ESS but leaves the ratio where it is. Widen the proposal
  with the inflation controls instead. Note too that this ratio is
  *optimistically biased* when estimated from few samples, so a small run
  flatters its own proposal and the figure can fall as `nSamples` rises simply
  because the estimate is becoming honest.
- **One candidate carrying more than half the weight.**

`type = "rsecor"` annotates each parameter's RSE with the confidence-interval
asymmetry ratio `(high - median) / (median - low)`. A symmetric
normal-approximation covariance reports one standard error per parameter and
cannot express that asymmetry, which is a large part of why SIR is run at all.

## Parity with PsN

| PsN option | `nlmixr2sir` | Status |
|---|---|---|
| `-samples` | `nSamples` | supported |
| `-resamples` | `nResample` | supported |
| covariance matrix from the fit | default | supported |
| `-rse_theta` / `-rse_omega` / `-rse_sigma` | `rseTheta` / `rseOmega` / `rseSigma` | supported |
| `-covmat_input=<file>` / `=identity` | `covmatInput` | supported |
| `-rawres_input` | `rawresInput` | supported |
| `-offset_rawres` | `offsetRawres` | supported |
| `-in_filter` | `inFilter` | supported |
| `-theta_inflation` etc., scalar or vector | `thetaInflation` etc. | supported |
| `-inflate_only_diagonal` semantics | always applied | supported |
| `-recenter` | `recenter` | supported |
| `-boxcox` | `boxcox` | partial — deliberate divergence, see below |
| `-cap_resampling` | `capResampling` | supported |
| `-cap_correlation` | `capCorrelation` | supported |
| `-add_iterations` | `addIterations` | supported |
| sample / resample count adjustment | automatic | supported, oracle-tested |
| dOFV vs chi-square plot | `plot(type = "convergence")` | supported |
| CI-by-iteration plot | `plot(type = "intervals")` | supported |
| RSE / correlation plot | `plot(type = "rsecor")` | supported |
| `empirical_statistics()` output | `sirSummary()` | supported |
| `<model>_sir.cov` | `<fitName>_sir.cov` | supported |
| estimation methods accepted | the deterministic ladder (`fo`/`foi`/`foce`/`focei`/`focep`/`laplace`/`agq` + `m…`/`i…`) and the importance-sampling family (`imp`/`impmap`/`qrpem`) | matches PsN on IMP/IMPMAP; still a deliberate difference on the rest, where PsN warns and continues; see below |
| draw-attempt budget | `10 x nSamples` | deliberate difference — PsN uses `2000 x nSamples`; see below |
| OMEGA/SIGMA block adjustment after prolonged rejection | — | deliberate difference — not implemented; see below |
| `-auto_rawres` | — | not implemented |
| `-print_iter` | — | not implemented |
| `-fast_posdef_checks` | — | not implemented |
| `rplots_level = 2` extras | — | not implemented |
| `-mceta`, `-copy_data`, `-problems_per_file`, `-nm_version` | — | not applicable (NONMEM execution) |

Numeric parity for the sample/resample adjustment, the inflation vector and
the RSE-to-variance conversion is checked against oracle values taken from
PsN's own unit tests.

PsN is used here as a comparator and a source of numerical oracles, not as the
specification. `nlmixr2sir` defines and tests its own statistical contract, and
the differences listed below are the known deliberate ones rather than an
exhaustive catalogue of every divergence.

### Deliberate differences from PsN

**Box-Cox includes the change-of-variables Jacobian.** With `boxcox = TRUE`,
candidates are drawn on a transformed scale and mapped back, so the density
induced on the original parameter scale is
`q_x(x) = q_y(T(x)) * |det J_T(x)|`. nlmixr2sir divides the likelihood by
`q_x`; PsN divides by `q_y`, omitting the Jacobian. Omitting it retains a
sample from `L(x)|det J_T(x)|` rather than `L(x)`, so the answer depends on
which smooth parameterization the model happens to be written in — exactly for
the skewed and weakly identified parameters Box-Cox is meant to help with.

The size of the effect is not subtle. Importance-sampling a known Gamma(3, 1)
target through a Box-Cox proposal recovers a mean of 3.00 and a second moment
of 12.00 with the Jacobian (truth 3 and 12), against 2.25 and 7.31 without it.
`tests/testthat/test-sir-boxcox-jacobian.R` holds that simulation plus the
analytic one-dimensional checks.

So SIR here targets the normalized likelihood **on nlmixr2's own parameter
scale**. The Jacobian is what makes Box-Cox an internal *proposal*
transformation: the retained distribution does not move when the internal
transform changes, because the induced proposal density is divided out
correctly. That is the estimand; it is a different one from PsN's, and it is
why `-boxcox` is marked partial above.

It does **not** make the result invariant to re-expressing the *model*. Flat
normalized likelihood is a choice of base measure: rewrite the model in another
parameterization and define the estimand against a flat measure there, and the
two disagree by the Jacobian of that reparameterization. Two models that are
reparameterizations of each other can give different SIR intervals. That is a
property of the estimand, not a defect.

**Only deterministic objectives are accepted, and the refusal is an error.**
`runSIR()` accepts the conditional-estimation ladder — `fo`, `foi`, `foce`,
`focei`, `focep`, `laplace`, `agq`, and their `m…`/`i…` mu-referencing
variants — plus the importance-sampling family `imp`, `impmap` and `qrpem`.
Candidates on the ladder are scored by re-evaluating the fit's own method at
fixed population parameters, carrying the fit's own control settings, so the
candidate surface and the `fit$objf` reference are the same function.

A stochastic or non-FOCEi-family fit is refused. On a SAEM fit of `theo_sd` the
stored objective is 208.512 (Gaussian quadrature) against 205.820 from a FOCEi
re-evaluation at the same estimates — a 2.69-unit gap that is a different
function, not noise. It does not cancel: under `recenter = TRUE` the centre
scores dOFV ≈ −2.69 and the run "finds" a better optimum made entirely of the
offset. The preflight rejects such a fit before any directory is created.

`imp`/`impmap`/`qrpem` are accepted, but on a different basis from the ladder,
and it is worth being precise about why. Their candidates are **not** scored by
re-running importance sampling. nlmixr2est recomputes the objective of every
imp-family fit as a nested FOCEi evaluation at the converged estimates, because
the in-C++ finalize leaves the eta-Hessian without its data term — so `fit$objf`
on such a fit has never been the importance-sampling objective. (That one lives
in `fit$env$impObj`, and `runSIR()` never reads it.) SIR therefore scores these
candidates directly as FOCEi, matching the calculation that produced the
reference. Measured on `theo_sd` at nlmixr2est 7.1.0, the two routes agree
bit-for-bit — 193.6046289649 with one random effect, 116.8319956005 with three
— and scoring as FOCEi is 6–9× faster per candidate, because it skips an E-step
whose result is discarded.

`runSIR()` warns once per run when it takes this route, naming both the method
you fitted with and the method it will score with, and the run fingerprint
records `evalMethod` (what the objectives were produced with) alongside
`estMethod` (what the fit was run with).

Verified against nlmixr2est 7.1.0. No minimum version is declared for it,
because SIR never uses the expectation-only evaluation path that 7.1.0 added —
what it relies on is the objective recompute, which is older. If a version ever
published the raw importance-sampling objective as `fit$objf` instead, the
per-run preflight would catch it: the FOCEi re-evaluation at the centre would
miss by about 19.96 units on a one-eta model, thousands of times the tolerance,
and the run would be refused before any sampling.

`npag`/`npb` do not have a normal `Omega` to propose from, and the variational
methods optimise a bound rather than the marginal likelihood. `saem` is refused
as above.

PsN's equivalent check lives in `set_maxeval_zero()`, which sets `MAXEVAL=0`
for classical methods and `EONLY=1` for `IMP`/`IMPMAP`, but for anything else —
`SAEM` included — only sets an internal failure flag and prints a warning. That
flag is discarded by its caller, so the run proceeds with evaluation models
that still carry the original method. PsN admits the IMP family where this
package cannot, and is more permissive at the edge; this package refuses rather
than warns. Adding a method here means validating an evaluator that reproduces
its objective, not adding a string to a list.

**The draw-attempt budget is `10 x nSamples`**, where PsN uses
`2000 x nSamples`. `runSIR()` is called from an interactive R session, where a
proposal bad enough to reject 99.95% of draws is better reported quickly than
ground through. On exhaustion the run warns with the attempted and successful
counts and continues on however many samples it collected.

**OMEGA/SIGMA blocks are not adjusted after prolonged rejection**, as PsN's are.
A model whose OMEGA block sits near the positive-definite boundary will
therefore reject more draws here than under PsN, and may exhaust the budget
where PsN would have continued. Widen the proposal with the inflation controls.

**`rse` is a percentage** where PsN reports a fraction.

**`rse_sd_scale` halves the RSE of OMEGA diagonals only.** PsN halves
everything that is not a NONMEM THETA, which catches `$SIGMA` because NONMEM
parameterises residual error as a variance, whereas nlmixr2 parameterises it on
the standard-deviation scale already. Off-diagonals are `NA`: a covariance can
be negative or zero and has no standard-deviation counterpart.

## Requirements and Practical Notes

`runSIR()` does not require a successful covariance step, but it does require
a proposal. When `fit$cov` is unavailable the run stops and names these three
routes; pick one and supply it:

```r
# from relative standard errors
runSIR(fit, control = runSIRControl(rseTheta = 30))

# from a diagonal proposal, widened by inflation
runSIR(fit, control = runSIRControl(covmatInput = "identity",
                                    thetaInflation = 0.05))

# seeded from the parameter vectors in a raw-results file
runSIR(fit, control = runSIRControl(rawresInput = "raw_results.csv"))
```

### Run directories, recovery, and running without files

A directory created by `runSIR()` carries a `sir_manifest.dcf` file recording
what produced it -- package version, fit, schedule, parameters, estimation
method, data rows, dependency versions, and seed. That manifest is also what
marks the directory as this package's: `runSIR()` refuses to clear a non-empty
directory that does not have one, so pointing `directory` at something valuable
cannot destroy it.

`recover = TRUE` (the default) resumes a run from its saved state, but only
after checking that the state belongs to the run being asked for. The fingerprint
covers the model, data, parameter set, estimates, objective, estimation method,
schedule, and statistical controls; a mismatch stops the run and names the
fields that changed rather than handing back another run's result. Worker and
thread settings are deliberately excluded -- they do not change the answer.

To run without touching the filesystem at all:

```r
runSIR(fit, control = runSIRControl(saveFiles = FALSE))
```

Nothing is written and no directory is created. The result is returned as usual
and `setCov()` registration still works, but recovery, `addIterations`, and
per-iteration seed reproduction all need the saved state; use `set.seed()`
before the call to make such a run reproducible.

For practical use:

* Use larger production schedules than toy examples; the default PsN-style
  schedule is usually a good starting point.
* Review the convergence diagnostic to make sure the proposal is not too
  narrow or too wide.
* Use `workers` and `rxThreads` to parallelize OFV evaluation when runs are
  large enough to justify it. Whenever `workers > 1`, `workers * rxThreads`
  must not exceed the machine's core count.
* `nlmixr2est::setCov(fit, "sir")` switches the fit's reported uncertainty to
  the SIR result after a run.

## Acknowledgments

The SIR methodology used here is based primarily on the method of [Dosne 
et al](https://link.springer.com/article/10.1007/s10928-016-9487-8). The technical implementation is based heavily on the 
[PsN tool](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/sir_userguide.pdf).

## References 

* Dosne, A.-G., Bergstrand, M., Harling, K., & Karlsson, M.O. (2016).
  Improving the estimation of parameter uncertainty distributions in nonlinear
  mixed effects models using sampling importance resampling. *Journal of
  Pharmacokinetics and Pharmacodynamics*, 43(6), 583-596.

For a fuller worked example, see the package vignette:
`vignette("runSIR", package = "nlmixr2sir")`.

For a specification of the implementation -- the proposal sources, the
importance-ratio and resampling maths, the Box-Cox proposal update, the
diagnostics, and the deliberate differences from PsN -- see
[`docs/sir-technical-reference.md`](docs/sir-technical-reference.md).

## Credit where it's due

`nlmixr2sir` is based on the [PsN implementation](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/sir_userguide.pdf) written by Lars Lindbom, 
Niclas Jonsson, Pontus Pihlgren, Mats Karlsson, Andrew Hooker, Kajsa Harling, 
Rikard Nordgren and Svetlana Freiberga.