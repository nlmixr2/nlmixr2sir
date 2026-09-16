# `nlmixr2sir` technical reference

This document describes the sampling importance resampling (SIR) procedure
implemented by `nlmixr2sir`. Source
references name functions rather than line numbers so that they remain useful
as the package evolves.

## Scope and statistical target

A SIR run produces an empirical sample from the parameter uncertainty
distribution of a fitted nonlinear mixed-effects model, without assuming that
distribution is multivariate normal. Given a fit with parameter estimate
$\widehat\psi$ and objective function value $\mathrm{OFV}(\widehat\psi)$, each
iteration:

1. draws $M$ candidate vectors $\psi_1,\ldots,\psi_M$ from a proposal
   distribution $g$;
2. evaluates each candidate's objective function with the population
   parameters fixed, giving
   $\Delta\mathrm{OFV}_i = \mathrm{OFV}(\psi_i) - \mathrm{OFV}(\widehat\psi)$;
3. forms importance ratios comparing the target density to $g$;
4. resamples $m < M$ vectors without replacement in proportion to those
   ratios; and
5. rebuilds $g$ from the empirical covariance of the resampled vectors,
   feeding the next iteration.

The retained vectors after the final iteration are the product. They
support nonparametric intervals, an empirical covariance, and asymmetric intervals.

The method is based on that of [Dosne, Bergstrand, Harling & Karlsson
(2016)](https://doi.org/10.1007/s10928-016-9487-8), which introduced SIR for
this purpose and characterised its behaviour against the covariance step,
bootstrap, and log-likelihood profiling. The iterative, self-tuning form
implemented here (repeated rounds with the proposal rebuilt each time, sample
counts compensated for evaluation failures) follows [Dosne, Bergstrand &
Karlsson (2017)](https://doi.org/10.1007/s10928-017-9542-0). See
[Literature](#literature).

The output is not a Bayesian posterior (there is no prior,
and the target is the likelihood surface around the maximum likelihood
estimate), it is not a bootstrap (the data are never resampled), and it does
not correct a misspecified model. A SIR interval is a statement about
parameter uncertainty under the fitted model, conditional on that model being
the right one.

Saying "the likelihood surface" is not quite enough to pin the target down.
Once a likelihood is normalized and treated as a density over parameters, a
base measure has been chosen, and a measure that is flat in $\psi$ is not flat
in a transform of $\psi$. **The target here is the normalized likelihood with
respect to a flat measure on nlmixr2's own parameter scale** — the scale the
`ini` block is written in. That choice is what the Box-Cox Jacobian in
[Proposal update and Box-Cox](#proposal-update-and-box-cox) enforces: the
Jacobian makes Box-Cox an internal *proposal* transformation while preserving
the normalized-likelihood target on nlmixr2's original parameter scale.

It is worth being precise about what that does and does not buy, because the
two are easy to conflate. If $x$ is the nlmixr parameter and $z = T(x)$ is only
a proposal coordinate, the induced proposal density is

$$
q_x(x) = q_z(T(x))\,\left|T'(x)\right|,
$$

and dividing by that density targets $L(x)\,dx$. The answer therefore does not
depend on the internal transform -- Box-Cox changes the efficiency of the
proposal, not the estimand.

Re-expressing the **model** is a different matter. If the model is rewritten in
$y = h(x)$ and the estimand is again defined against a flat measure $dy$, then
mapping that target back to $x$ introduces $\left|h'(x)\right|$. Flat
normalized likelihood is a choice of base measure, not a parameterization-free
object, so it is **not** invariant to redefining the model's parameter scale.
Two models that are reparameterizations of each other can give different SIR
intervals, and that is a property of the estimand rather than a defect in the
implementation.

## Public interface

| Function | Purpose |
| --- | --- |
| [`runSIR()`](../R/sir-run.R) | Run the iterative SIR workflow and write its artifacts |
| [`runSIRControl()`](../R/sir-control.R) | Construct and validate the run settings |
| [`sirSummary()`](../R/sir-results.R) | Empirical summary of a resampled matrix |
| `print()`, `plot()` | S3 methods on the returned `nlmixr2SIR` object |

`runSIR()` takes the fit, the sampling schedule, and where output goes;
everything that tunes *how* the run behaves lives in `runSIRControl()`.

```r
sir <- runSIR(
  fit,
  nSamples  = c(1000, 1000, 1000, 2000, 2000),
  nResample = c(200, 400, 500, 1000, 1000),
  control   = runSIRControl(thetaInflation = 2, workers = 4, rxThreads = 2)
)
```

## The parameter vector

[`.sirParamSpace()`](../R/sir-paramspace.R) is the single source of truth: one
row per estimated, non-fixed parameter, carrying every name that parameter is
known by.

| kind | SIR name | `fit$cov` rowname | raw-results column |
| --- | --- | --- | --- |
| `theta` | `tka` | `tka` | `tka` |
| `sigma` | `add.sd` | `add.sd` | `add.sd` |
| `omegaDiag` | `eta.ka` | `om.eta.ka` | `omega(eta.ka,eta.ka)` |
| `omegaOffdiag` | `eta.cl:eta.ka` | `cov.eta.cl.eta.ka` | `omega(eta.cl,eta.ka)` |

Rows are ordered THETA, then residual error, then the OMEGA lower triangle by
column and then row. That ordering reproduces `rownames(fit$cov)` exactly, so
the fitted covariance can be consumed without reordering, and coincides
with PsN's `parameter_hash` order, which facilitates positional comparison against
PsN output if needed.

OMEGA off-diagonals are ordered
`(neta1, neta2)`, not alphabetically, and eta names may contain periods. This means
`cov.eta.cl.eta.ka` cannot be split back unambiguously. Every name is derived
from the eta index pair, and is never parsed. `kind == "sigma"` is decided
structurally from `iniDf$err`, not from absence from `fit$cov`.

## The initial proposal

Four sources are available, dispatched in this order by
[`.sirResolveInitialProposal()`](../R/sir-proposal-input.R). They are mutually
exclusive and validated at control construction, rather than resolved silently
by precedence.

### From the fitted covariance

This is the default. Since `nlmixr2est` 7, `foceiControl(covFull = TRUE)` is the
default and `fit$cov` spans THETA, residual error, and OMEGA jointly, so

$$
  g_1 = N\!\left(\widehat\psi,\ \widehat V\right),
  \qquad \widehat V = \texttt{fit\$cov},
$$

including the estimated THETA-OMEGA correlations. 

### Wishart-style fallback

When `fit$cov` is present but does not include OMEGA (`covFull = FALSE`, or a
partial covariance), or when `omegaFallback = "wishart"` forces the route even
though OMEGA is available, OMEGA uncertainty is approximated with
$df = n_{\text{sub}} - 1$ by default:

$$
  \operatorname{Var}(\widehat\Omega_{jj}) = \frac{2\widehat\Omega_{jj}^2}{df},
  \qquad
  \operatorname{Var}(\widehat\Omega_{jk})
    = \frac{\widehat\Omega_{jj}\widehat\Omega_{kk}
            + \widehat\Omega_{jk}^2}{df}.
$$

This route gives a block-diagonal proposal. The route taken is reported
in the run log and on the returned object.

**It completes an incomplete covariance; it does not replace an absent one.**
A fit with `covMethod = ""`, or whose covariance step failed, has
`fit$cov == NULL`. The approximation above needs only the OMEGA estimates and
the subject count, so OMEGA uncertainty would still be available -- but nothing
in such a fit supplies THETA uncertainty, and normalising a likelihood over a
THETA whose scale was assumed rather than estimated would fabricate the very
quantity being reported. `runSIR()` therefore stops, naming
[`rseTheta`](#from-relative-standard-errors), `covmatInput` and `rawresInput`
as the routes that do carry the missing information.

### From relative standard errors

`rseTheta`, `rseOmega` and `rseSigma` can be used to build a diagonal proposal with no
covariance step at all. Diagonal
variances are $(\mathrm{rse}\cdot\widehat\psi_j/100)^2$. OMEGA off-diagonals
use

$$
  N = \left(\frac{100}{\mathrm{rse}_j}\right)^2
    + \left(\frac{100}{\mathrm{rse}_k}\right)^2 + 1,
  \qquad
  \operatorname{Var}(\widehat\Omega_{jk})
    = \frac{\widehat\Omega_{jk}^2
            + \widehat\Omega_{jj}\widehat\Omega_{kk}}{N}.
$$ 

Each argument is used for the whole class or one value per estimated
element of it. A scalar `rseTheta` fills in an unset `rseOmega`
and `rseSigma`, a vector `rseTheta` does not, and setting `rseOmega` without
`rseTheta` is an error. Inflation cannot be combined with this route; the RSE
already states the width.

### Supplied directly, or from previous parameter vectors

`covmatInput` accepts a matrix, a NONMEM-style `.cov` file, or `"identity"`;
the last together with inflation is the cheap "any diagonal proposal" route.
`rawresInput` seeds the first proposal from the parameter vectors in a
canonical raw-results file, taking their empirical mean
and covariance, with `offsetRawres` and `inFilter` narrowing which rows are
used. Any canonical raw-results file works, including one written by
`nlmixr2boot`.

## Inflation, correlation capping, and positive-definiteness

Inflation multiplies each parameter's proposal *variance*, preserving
correlations, and is applied to the initial proposal only - from the second
iteration the proposal is the previous empirical covariance and re-inflating
it each round would be unhelpful. In [`.sirInflationVector()`](../R/sir-proposal.R), each argument is a scalar or one value per
*diagonal* element of its class, and an OMEGA off-diagonal is never given a
factor directly but derives

$$
  c_{jk} = \sqrt{c_{jj}}\,\sqrt{c_{kk}},
$$

which leaves the correlation unchanged when the two factors are equal. The
rescaling is implemented as
$\Sigma_{jk} \mapsto \Sigma_{jk}\sqrt{c_j c_k}$, algebraically identical to
rescaling standard deviations but requiring no `cov2cor()`, which would fail
on a variance of exactly zero, as the Wishart fallback produces for an OMEGA
element estimated at zero.

[`.sirCapCovCorrelation()`](../R/sir-utils.R) then clamps every off-diagonal
correlation to $\pm$ `capCorrelation` (default 0.8) while holding the
standard deviations fixed, and [`.sirEnsurePosDef()`](../R/sir-utils.R)
symmetrises and conditions the matrix. Capping a correlation can itself destroy
positive-definiteness, so the order matters.

The repair works in **standardized coordinates**, not raw ones. It divides out
the marginal standard deviations, floors the eigenvalues of the resulting
*correlation* matrix at `relTol * max(eigenvalue)` (with `relTol = 1e-12`),
renormalizes to a unit diagonal, and maps back with the original standard
deviations. Marginal variances therefore come back exactly as they went in, and
only the correlation structure is conditioned.

That matters because pharmacometric parameters do not share a unit. An
eigenvalue test in raw coordinates measures the spread of the units as much as
the spread of the information, so the same statistical problem gets a different
verdict depending on how a parameter happens to be expressed: an absolute floor
of $\sqrt{\varepsilon}$ once turned `diag(c(1, 1e-14))` into roughly
`diag(c(1, 1e-12))`, inflating one parameter's variance a hundredfold because a
*different* parameter happened to have variance one. A zero-variance coordinate
carries no uncertainty, cannot be standardized, and is held at zero rather than
being given uncertainty by the repair.

The repair exists to clean roundoff on an already-full-rank matrix. Genuine
rank deficiency is not repaired; it aborts.

## Sampling and rejection

[`.sirSampleFullProposal()`](../R/sir-proposal.R) draws from the multivariate
normal proposal in batches until $M$ valid vectors are collected or a budget
of `maxAttemptFactor * M` draws is exhausted. A draw is rejected if:

- back-transformation from the Box-Cox scale produces any non-finite value
  (`inverseRejected`). This is tested with `!is.finite()`, not `is.na()`:
  inversion can overflow to `Inf` without raising an error, and an infinite
  value would otherwise pass the bounds test too, since `Inf > Inf` is `FALSE`.
  A non-finite draw is an inverse failure whatever else is true of it, so it is
  classified here, before any parameter-specific check;
- any THETA or residual-error parameter falls outside its `iniDf` bounds
  (`thetaRejected`, `sigmaRejected`); or
- the reconstructed OMEGA matrix is not positive-definite, tested by Cholesky
  (`omegaRejected`).

OMEGA elements are constrained by the positive-definiteness test rather than
by element-wise bounds, which keeps the four rejection counts separable — they
are reported per iteration and written to `sample_rejection_summary.txt`. A
run that rejects heavily in one category is diagnosable; a single pooled count
would not be.

### The draw-attempt budget

`maxAttemptFactor` is `10`, so the budget is `10 * M` draws for `M` requested
samples. PsN uses `2000 * M`. This is a deliberate choice rather than an
oversight: `runSIR()` is called interactively from an R session, where a
proposal bad enough to reject 99.95% of draws is better reported quickly than
ground through. On exhaustion the run warns with the attempted and successful
counts and proceeds on however many samples it did collect, so a marginal case
still produces a result -- it just says so.

PsN additionally adjusts OMEGA and SIGMA blocks after prolonged rejection.
`nlmixr2sir` does not. The consequence is concrete: a model whose OMEGA block
sits near the positive-definite boundary will reject more draws here than under
PsN, and may exhaust the budget where PsN would have continued. Widening the
proposal with the inflation controls is the remedy.

### Requested, attempted, collected, successful, usable, retained

Six counts appear in the iteration summary and they are deliberately distinct:

| Count | Column | Meaning |
| --- | --- | --- |
| Requested | `nSamples` | What the schedule asked for |
| Attempted | `nAttempted` | After the turnout adjustment of the previous iteration |
| Draw attempts | `nDrawAttempts` | Raw draws made, including rejected ones |
| Collected | `nCollected` | Draws that survived rejection |
| Successful | `nSuccessful` | Collected draws with a usable objective value |
| Retained | `nResampled` | Distinct vectors actually kept by the resampler |

A seventh, *usable*, sits between successful and retained: a finite dOFV does
not guarantee a finite importance ratio, so only candidates with non-zero
resampling probability can actually be drawn. Conflating any of these makes a
failure report point at the wrong cause -- which is why an infeasible retained
set now aborts naming how many candidates could be *scored*, rather than
letting the rank check report the final count and advise raising `nResample`.

## Objective function evaluation

[`sirEvalOFV()`](../R/sir-eval.R) sets each sampled vector into the model with
`rxode2::ini()` and re-evaluates **the fit's own estimation method** with
`maxOuterIterations = 0`, so no estimation occurs — the population parameters
are fixed at the proposed values and the inner problem is solved.

Both halves of that matter. The method is taken from the fit rather than
hard-coded: an `fo` fit evaluated as `focei` scored 103.870 against its own
127.982, *below* FOCEi's own minimum, because the etas were being estimated
rather than held at zero. And the **whole** control object comes from the fit,
with only the evaluation fields overridden
(`maxOuterIterations`, `calcTables`, `covMethod`, `compress`, `print`).

That second point replaced an earlier design that rebuilt the control from a
hand-picked list of likelihood-relevant fields. `foceiControl()` has 150
arguments, so such a list fails *open*: a setting nobody listed is silently
dropped and the candidate is scored on a different surface. `agqLow`/`agqHi`
were lost exactly that way — an AGQ fit with `agqLow = -100` was re-evaluated
with the default `-Inf`, agreeing at the centre to 1e-08 but differing by
about 6490 OFV units at `tka = -20`, which enters the weight as
$\exp(-\Delta\mathrm{OFV}/2)$. Copying the control and overriding only the
evaluation fields inverts the failure mode: an unrecognised setting is
preserved rather than lost.

Note what this means for the preflight. The centre check and the stencil
cannot catch a defect of that shape — the stencil perturbs by a thousandth of
each estimate, and integration bounds only bite far from the mode. The
protection against it is structural (carry everything) plus the regression
tests that evaluate off-centre candidates, not the preflight. Evaluation is parallelised across `workers`, each worker using
`rxThreads` rxode2 threads; whenever `workers > 1`, `workers * rxThreads` must
not exceed the core count, since each worker is a separate process with its
own thread pool.

Failed evaluations return `NA` rather than aborting the run, but the
underlying error messages are retained and surfaced if *every* evaluation
fails. A configuration fault fails all samples identically, and reporting only
"all evaluations failed" hides the cause.

### Only validated estimation methods are accepted

`.sirSupportedEstimationMethods` holds the deterministic conditional-estimation
ladder — `fo`, `foi`, `foce`, `focei`, `focep`, `laplace`, `agq`, and the
`m…`/`i…` mu-referencing variants of each — and
[`.sirCheckObjective()`](../R/sir-preflight.R) rejects anything else before a
directory is created or a single candidate is drawn.

All of these run on the same FOCEi engine and differ only in settings the
evaluator now carries, so one evaluator reproduces each one's own objective.
Verified on `theo_sd`, re-evaluating each fit at its own estimates:

| Method | Stored | Re-evaluated | Absolute difference |
| --- | --- | --- | --- |
| `focei`/`foce`/`focep`/`laplace` | 116.8042 | 116.8042 | ~4e-06 |
| `agq` | 118.4833 | 118.4833 | 2.8e-07 |
| `mfocei`/`ifocei` | 116.8569 | 116.8569 | 1.9e-06 |
| `fo`/`foi` | 127.9822 | 127.9822 | ~1e-13 |

The restriction that remains is not conservatism for its own sake. A fit whose
objective came from a different likelihood approximation would have its
candidates scored on one surface and its reference `fit$objf` taken from
another. Measured on a SAEM fit of `theo_sd`:

| Quantity | OFV |
| --- | --- |
| stored `fit$objf` (Gaussian quadrature) | 208.512 |
| FOCEi reevaluation at the same estimates | 205.820 |
| difference | **2.69** |

That is a different function, not numerical noise, and it is not a constant
that cancels: with `recenter = TRUE` the centre itself scores dOFV around
$-2.69$, so the run would immediately "find" a better optimum manufactured
entirely out of the offset.

Adding a method to the allowlist means validating an evaluator that reproduces
*its* objective, not adding a string.

A second reason to be careful here is stochasticity. A deterministic evaluator
returns the same OFV for the same vector, so a dOFV difference is signal. An
MCMC or Monte-Carlo E-step does not, and the resulting noise propagates
straight into $\exp(-\Delta\mathrm{OFV}/2)$ and hence into the weights, with
nothing in the SIR machinery to account for it. Supporting such a method
properly would require deciding what the noise does to the importance weights,
which is a design question rather than a configuration one.

## Importance weights

For proposal $g = N(\mu,\Sigma)$ with Cholesky factor $L$,
[`sirCalcWeights()`](../R/sir-weights.R) computes a *relative* proposal
density, normalised so that $\psi_i = \mu$ gives exactly 1:

$$
  \log \mathrm{relPDF}_i
    = -\tfrac12 \left\| L^{-\mathsf T}(\psi_i - \mu) \right\|^2 .
$$

The likelihood ratio relative to the fit is
$\log \mathrm{LR}_i = -\tfrac12 \Delta\mathrm{OFV}_i$, so the importance ratio
and resampling probability are

$$
  \log \mathrm{IR}_i = \log \mathrm{LR}_i - \log \mathrm{relPDF}_i,
  \qquad
  p_i = \frac{\exp(\log \mathrm{IR}_i - \max_k \log \mathrm{IR}_k)}
             {\sum_j \exp(\log \mathrm{IR}_j - \max_k \log \mathrm{IR}_k)} .
$$

Everything is carried on the log scale and the maximum is subtracted before
exponentiating, because raw importance ratios overflow readily. Normalising
constants common to all samples cancel in $p_i$ and are never formed.

## Resampling

[`sirResample()`](../R/sir-weights.R) draws $m$ vectors **without
replacement** with probability proportional to $p_i$. Sampling without
replacement is what distinguishes SIR here from a naive importance sample: it
prevents a single high-weight vector from dominating the retained set, at the
cost of requiring $m < M$.

`capResampling` relaxes this. A value $c > 1$ expands each candidate into $c$
slots before drawing, so a vector may be selected up to $c$ times — limited
replacement, with the cap bounding how far any one vector can dominate.

## Proposal update and Box-Cox

[`sirUpdateProposal()`](../R/sir-boxcox.R) rebuilds the proposal from the
empirical covariance of the retained vectors. With `boxcox = TRUE` each
column is first transformed by

$$
  x^{(\lambda)} = \begin{cases}
    \dfrac{(x+\delta)^{\lambda} - 1}{\lambda}, & \lambda \neq 0 \\[2ex]
    \log(x+\delta), & \lambda = 0
  \end{cases}
$$

with $\delta = |\min x| + 10^{-6}$ guaranteeing positivity, and $\lambda$
chosen on $[-3,3]$ to maximise the correlation between the sorted transformed
values and normal scores — a normality-of-fit criterion rather than a
profile likelihood. The covariance is then taken on the transformed scale, so
the next iteration proposes in a space where the parameters are closer to
normal, and draws are back-transformed before evaluation.

### The change-of-variables Jacobian

Drawing on the transformed scale means the normal density of the draws is
$q_y$, not the density induced on the original parameter scale. Those differ by
the Jacobian of the transform:

$$
  q_x(x) = q_y\bigl(T(x)\bigr)\,\bigl|\det J_T(x)\bigr|.
$$

The likelihood in the numerator of the importance ratio is a function of $x$,
so the weight must divide by $q_x$:

$$
  w(x) \;\propto\; \frac{L(x)}{q_y(T(x))\,\bigl|\det J_T(x)\bigr|}.
$$

Box-Cox is applied one coordinate at a time, so $J_T$ is diagonal and

$$
  \log\bigl|\det J_T(x)\bigr|
    = \sum_j (\lambda_j - 1)\,\log(x_j + \delta_j),
$$

which is what [`.sirBcLogJacobian()`](../R/sir-boxcox.R) computes. It is passed
to [`sirCalcWeights()`](../R/sir-weights.R) relative to the proposal centre, so
`relPDF` remains 1 at the centre and keeps its meaning: the density the weight
divides by, relative to that centre.

**This is a deliberate divergence from PsN**, which evaluates the transformed
normal density without a Jacobian (`lib/tool/sir.pm`). Omitting the term
retains a sample from $L(x)\,|\det J_T(x)|$, so the answer depends on the
parameterization the model happens to be written in. Importance-sampling a
known Gamma(3, 1) target through a Box-Cox proposal recovers mean 3.00 and
second moment 12.00 with the Jacobian, against 2.25 and 7.31 without it, where
the truth is 3 and 12.

The Box-Cox transform is deliberately **not** applied on the final iteration:
the last proposal is built on the original scale so the delivered vectors and
their covariance need no back-transformation.

If `recenter = TRUE` and any sample has $\Delta\mathrm{OFV} < 0$ — meaning a
proposed vector fits better than the reported maximum likelihood estimate, so
the fit was not fully converged — the proposal is recentred on that vector.
With `recenter = FALSE` the same condition warns instead.

## Sample count adjustment

Evaluation failures shrink the usable sample, so both counts are compensated
per iteration ([`R/sir-iterate.R`](../R/sir-iterate.R)). With
turnout $t$ = successful / requested:

$$
  M' = \begin{cases}
    \mathrm{round}(M / t_{\text{prev}}), & t_{\text{prev}} \le 0.95 \\
    M, & \text{otherwise}
  \end{cases}
  \qquad
  m' = \begin{cases}
    \mathrm{round}(m \cdot t), & |t - 1| \ge 0.05 \\
    m, & \text{otherwise}
  \end{cases}
$$

The attempted count compensates for loss only; the resample count scales on
gain *or* loss. Turnout for the resample adjustment is measured against the
originally requested sample count, not the compensated attempted count.

`round()` here is round-half-away-from-zero, implemented as
[`.sirRound()`](../R/sir-utils.R). R's own `round()` is round-half-to-even and
disagrees on exact halves — `round(20.5)` is 21 in nlmixr2sir and 20 in R. 

## Diagnostics

### Convergence: dOFV against a reference chi-square

The most informative SIR diagnostic. For a
quantile grid $q$ stopping short of 1 so the reference stays finite,
[`.sirDofvCurves()`](../R/sir-convergence.R) draws three curves per iteration:

| curve | definition |
| --- | --- |
| reference | $\chi^2_{p}$ quantiles, $p$ = number of estimated parameters |
| proposal | empirical $\Delta\mathrm{OFV}$ quantiles over all evaluated samples |
| SIR | empirical $\Delta\mathrm{OFV}$ quantiles over the resampled subset |


#### What the reference curve does and does not establish

The chi-square reference is a consequence of regular likelihood asymptotics:
that $2(\ell(\widehat\psi) - \ell(\psi))$ is approximately $\chi^2_p$ near a
well-identified interior maximum, with enough subjects for the approximation
to hold. Read the curve as evidence to interpret, not as a certificate. It can
mislead when:

- a variance component sits **on or near a boundary** (an OMEGA element
  estimated at or close to zero), where the asymptotic distribution is not
  chi-square;
- a parameter is **weakly identified**, so the likelihood is flat in some
  direction and the quadratic approximation never applies;
- the likelihood is **multimodal**, where a single centre describes none of
  the modes;
- the likelihood is **non-smooth** in the parameters, for instance through
  hard bounds or discrete model switches; or
- the **subject count is small** relative to the number of parameters.

Agreement with the reference is therefore evidence that the importance sample
has settled, not proof that the interval is correct. In the same vein, the
percentile intervals in [`sirSummary()`](../R/sir-results.R) are quantiles of
a likelihood-weighted retained sample. They are not guaranteed to have
nominal frequentist coverage, and they inherit every one of the conditions
above. Where coverage matters and these conditions are in doubt, a simulation
study on the model at hand is the only way to establish it.

Under the asymptotic theory the SIR curve should approach the reference.
Convergence reads as that curve settling onto it across iterations, with a
resampling-noise band on the last two iterations obtained by repeating the
weighted resampling and taking the 2.5th and 97.5th percentiles of the
resulting curves.

If the first iteration's proposal curve falls *below* the reference for more
than a quarter of the quantiles, the proposal is too narrow. The vectors SIR
would need were never drawn, and resampling cannot manufacture them.

The warning is emitted by [`plot()`](../R/sir-methods.R) when
`type = "convergence"` is drawn, **not** by `runSIR()`. A run that is never
plotted will not raise it, so treat the convergence plot as part of checking a
run rather than as optional decoration. This check is based on PsN's.

### Intervals by iteration, and CI asymmetry

`plot(type = "intervals")` shows the proposal and SIR interval per parameter
per iteration; uncertainty that is still moving between the last two
iterations means the run has not settled.

`plot(type = "rsecor")` draws RSE% on the diagonal and correlations off it,
annotating the diagonal with

$$
  \text{asymmetry} = \frac{P_{\text{high}} - P_{\text{med}}}
                          {P_{\text{med}} - P_{\text{low}}},
$$

banded at 0.5, 1, 1.25 and 2. A symmetric normal-approximation covariance
reports one standard error per parameter and cannot express that a parameter's
upper interval half is twice its lower. Showing it is a large part of why SIR
is run at all.

## Summaries and artifacts

[`sirSummary()`](../R/sir-results.R) reports `estimate`, `mean`, `sd`, `rse`,
`rse_sd_scale`, at percentiles 2.5, 5, 10, 30, 50, 70, 90, 95,
97.5, derived from prediction intervals 0, 40, 80, 90 and 95. Empirical
covariance, correlation, and standard-deviation/correlation matrices are
attached as attributes and written to disk.

While summaries are structured to be similar to PsN's in order to facilitate comparisons, two differences must be highlighted. First, `rse` is a percentage where PsN
reports a fraction; the returned object records this in an `rseUnits`
attribute rather than leaving it implicit. Second, `rse_sd_scale` halves the RSE
of OMEGA **diagonals** only, whereas PsN halves everything that is not a NONMEM
THETA. NONMEM parameterises residual error
as a variance, whereas nlmixr2 parameterises it on the standard-deviation
scale, so halving `add.sd` would rescale a quantity that needs no
rescaling. OMEGA off-diagonals are reported as `NA`: the delta-method relation
$\operatorname{RSE}(\sqrt v) \approx \tfrac12 \operatorname{RSE}(v)$ needs a
positive variance, and a covariance can be negative or zero and has no
standard-deviation counterpart. Use the empirical correlations in the
`sdCorMatrix` attribute for off-diagonal uncertainty.

A run writes `sir_results.csv`, `summary_iterations.csv`, `<fitName>_sir.cov` and `.sdcorr`, a canonical
`raw_results.*` set, `sample_rejection_summary.txt`, and `sir_state.rds`.
`runSIR()` also registers the empirical covariance so that
`nlmixr2est::setCov(fit, "sir")` switches the fit's reported uncertainty to
the SIR result, skipping registration with a message if the parameters do not
match `fit$cov` or the covariance is not positive-definite.

## Persistence, resume and extension

State is written after every iteration, under a versioned schema. Alongside
it goes a **run fingerprint**: the model text, a digest of the dataset and its
row count, the estimated parameter set, the parameter *schema* (kinds and
bounds, not just names), a digest of the **resolved initial proposal** -- its
covariance and mean after parsing and validation -- the estimates, the
objective and estimation method, the cumulative sample/resample schedule, a
state and algorithm version, and a digest of the statistical controls. Worker
and thread settings are excluded, because they do not change the answer and so
must not invalidate a saved run.

Digesting the *resolved* proposal rather than the control object is what covers
every input route at once. A control digest only ever captured a path string,
so replacing the file at that path left the fingerprint unchanged, and a
changed covariance could be reused as though it were the original.

Comparison is field by field, so a mismatch names what moved rather than
reporting an opaque hash difference. On the recovery path it **fails closed**:
a field that cannot be digested blocks reuse rather than being skipped, because
recovery is exactly the situation where inability to establish identity must
not be read as permission. A fresh run is more forgiving, since an unusual fit
should not be blocked because one field would not serialize.

Package versions are recorded and warned about but **not** enforced. A
dependency bump does not by itself invalidate a result, and a blanket
requirement that every version match would be unnecessarily strict; but a
version change can alter proposal construction, Box-Cox estimation or
random-number behaviour, so it is never silent.

`recover = TRUE` resumes from the last completed iteration, returning the
stored result unchanged if the schedule was already finished -- but only after
the saved fingerprint matches the current one. A mismatch stops the run and
names the fields that moved. Without that check, pointing a different fit at
an existing directory returns a stale result labelled as the new run's, which
is a provenance failure rather than a caching one.

`addIterations = TRUE` appends further iterations to a completed run, carrying
the existing iterations over rather than recomputing them. It exempts the
schedule field from the identity comparison, which it deliberately changes, and
nothing else.

The schedule it *stores* is **cumulative**: the prior schedule plus the
extension, covering every iteration the result contains. Storing the extension
alone was a provenance failure rather than a cosmetic one -- a two-iteration run
extended by one saved three completed iterations beside an identity describing
a single iteration, so a later plain recovery presenting that one-iteration
schedule matched and was handed the three-iteration result.

Directories created by `runSIR()` carry a `sir_manifest.dcf` manifest. It is
human-readable provenance, and it is also the ownership marker.

Ownership is established **before anything is written** -- before the seed file
and before the manifest itself -- and on every path, not only when overwriting.
That ordering is the point: writing the marker must not be what creates the
ownership it later checks for. Under the default `recover = TRUE` an existing
non-empty directory comes back in resume mode, and guarding only the overwrite
path once let such a directory acquire a manifest and so become eligible for
recursive deletion by the next run.

The marker is validated by content, not by filename: the manifest is parsed and
its package, prefix and state version checked. A foreign or malformed
`sir_manifest.dcf` does not establish ownership. A non-empty directory without
a valid manifest is refused outright rather than claimed, and failure to write
a manifest is fatal, because the marker participates in the deletion policy.

### Covariance repair provenance

Two distinct repairs are recorded, because they happen at different points and
mean different things:

- `initialProposalRepair` on the result, and `proposalRepaired` per iteration,
  describe the covariance that iteration actually *drew from*. Iteration 1's is
  the run's initial-proposal record, and a resumed run reads it back rather
  than losing it, since its iteration 1 is not re-run.
- `posDefAdjusted` per iteration describes the repair applied to the *empirical
  update* that iteration produced for the next one.

Each record carries whether a repair was applied, the method, the threshold,
and the magnitude -- the largest absolute change to any entry. The flag says a
repair happened; the magnitude says whether it mattered.

Seeding is managed per iteration through `nlmixr2utils::withRunSeed()`, so a
resumed run reproduces the stream it would have had.

`runSIRControl(saveFiles = FALSE)` turns all of this off: no directory, no
files, no state. The result is returned as usual, but recovery,
`addIterations`, and per-iteration seeding are unavailable, and a single
`set.seed()` before the call is what makes the run reproducible.

## Interpretation checklist and limitations

- **Read the convergence plot first.** If the iteration-1 proposal sits below
  the reference chi-square, the result is not usable; restart with inflation.
- **Check the sample-to-resample ratio.** Roughly 5:1 is the working default.
  Heavy rejection in any one category points at bounds, positive-definiteness,
  or a proposal on the wrong scale.
- **Check that the last two iterations agree.** If intervals are still moving,
  add iterations or increase sample counts.
- SIR characterises uncertainty under the fitted model. It cannot diagnose
  structural misspecification, and a tight SIR interval around a wrong model
  is still wrong.
- The target is the likelihood surface near $\widehat\psi$. With a
  poorly-identified parameter the surface may be flat or multimodal, and SIR
  will report that flatness faithfully rather than resolving it.
- Default schedules are expensive: the objective is evaluated once per sample,
  which at the default 7,000 samples is minutes to hours depending on the
  model. `workers` and `rxThreads` are the practical levers.

## Difference from PsN

The algorithm, the sample-count adjustment rules, the inflation semantics, the
RSE-to-variance conversion and the principal diagnostics follow those implemented by PsN, and the
numeric cores are checked against reference values extracted from PsN's own unit
tests (`test/unit/tool/sir.t`). The implementations diverge in execution:
PsN generates and runs NONMEM control streams, whereas `nlmixr2sir` calls
`rxode2::ini()` and `nlmixr2est::nlmixr2()` directly and keys everything by
canonical named parameter schemas rather than by position.

Four PsN options are not implemented: `-auto_rawres`, `-print_iter`,
`-fast_posdef_checks`, and the `rplots_level = 2` extras (bin-exhaustion
diagnostics and inverse-Wishart degrees-of-freedom estimation). The
NONMEM-execution options — `-mceta`, `-copy_data`, `-problems_per_file`,
`-nm_version` and similar — have no analogue. A per-option parity matrix is
kept in the [README](../README.md).

**PsN is a comparator, not a specification.** It is an independently developed
implementation of the same method, which makes it valuable in two narrower
roles: as a source of exact numerical oracles for primitives both packages
share, and as a source of workflow ideas and edge cases. It is not the
normative definition of what `nlmixr2sir` should do. Where the two differ,
the question is whether *this* package is mathematically coherent and does what
it documents -- not whether it matches PsN. `nlmixr2sir` defines and tests its
own statistical contract.

The differences below are the ones that are known and deliberate. That is not a
claim to have enumerated every divergence exhaustively; two independent
implementations of a stochastic method will differ in ways neither author has
catalogued. The deliberate differences are:

- **The Box-Cox change-of-variables Jacobian**, described under
  [The change-of-variables Jacobian](#the-change-of-variables-jacobian). This
  is an *algorithmic* divergence, not a reporting one: it changes the retained
  distribution whenever `boxcox = TRUE`, which is the default. PsN omits the
  term; including it is what makes the target the original-scale normalized
  likelihood rather than a parameterization-dependent tilt of it.
- **Unsupported estimation methods are refused rather than warned about.**
  PsN builds its evaluation models through `set_maxeval_zero()`
  (`lib/model.pm`), which handles three cases: a classical method (`FO`,
  `FOCE`, `FOCEI`, `Laplace`, or no `METHOD` at all) gets `MAXEVAL=0`; `IMP`
  and `IMPMAP` get `EONLY=1`; and anything else -- `SAEM` included -- sets an
  internal `$success = 0` and prints

  > `METHOD in last $EST was not classical nor IMP/IMPMAP. Cannot set`
  > `MAXEVAL=0 or EONLY=1.`

  That return value is discarded by the caller
  (`create_maxeval_zero_models_array()`), and `sir.pm` has no `METHOD` guard of
  its own, so the run continues and builds evaluation models that still carry
  the original method. The PsN source records the open question in a comment:
  *"if other method return error no success / should this be changed to replace
  last est with something classical?"*

  `nlmixr2sir` aborts instead, in the preflight, before any directory is
  created. See [Objective function evaluation](#objective-function-evaluation)
  for the measured size of the problem.

  The comparison is not exactly like for like -- PsN's classical path covers
  several NONMEM methods, and PsN additionally admits the IMP family, which
  `nlmixr2sir` cannot evaluate at fixed parameters at all. The difference is in
  what happens at the edge: PsN warns and proceeds, `nlmixr2sir` refuses.

- **Stochastic candidate evaluation is not supported at all.** PsN does admit
  one stochastic evaluator, `IMP`/`IMPMAP` under `EONLY=1`, whose Monte-Carlo
  E-step returns an OFV *estimate*. Nothing in PsN's SIR quantifies or
  compensates for that sampling noise; the dOFVs inherit it and so do the
  weights. `nlmixr2sir` has no equivalent path, which sidesteps the question
  rather than answering it.

- The two reporting differences described under
  [Summaries and artifacts](#summaries-and-artifacts) — `rse` as a percentage,
  and the OMEGA-only `rse_sd_scale` rule.
- The off-diagonal RSE rule following PsN's code rather than its documentation,
  described under
  [From relative standard errors](#from-relative-standard-errors).

This comparison is based on PsN's
[`tool::sir`](https://github.com/UUPharmacometrics/PsN/blob/master/lib/tool/sir.pm),
its diagnostic template
[`sir_default.R`](https://github.com/UUPharmacometrics/PsN/blob/master/R-scripts/sir_default.R),
and the
[SIR user guide](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/sir_userguide.pdf).

## Literature

1. Dosne A-G, Bergstrand M, Harling K, Karlsson MO. Improving the estimation
   of parameter uncertainty distributions in nonlinear mixed effects models
   using sampling importance resampling. *Journal of Pharmacokinetics and
   Pharmacodynamics*. 2016;43:583-596.
   [doi:10.1007/s10928-016-9487-8](https://doi.org/10.1007/s10928-016-9487-8).
   This is the original SIR paper and the primary reference for the method
   implemented here: the importance-ratio construction, the use of $\Delta$OFV
   against a reference chi-square as the convergence criterion, and the
   comparison against the covariance step, bootstrap, and log-likelihood
   profiling.

2. Dosne A-G, Bergstrand M, Karlsson MO. An automated sampling importance
   resampling procedure for estimating parameter uncertainty. *Journal of
   Pharmacokinetics and Pharmacodynamics*. 2017;44:509-520.
   [doi:10.1007/s10928-017-9542-0](https://doi.org/10.1007/s10928-017-9542-0).
   This develops the iterative, self-tuning procedure that this package
   implements — multiple rounds with the proposal rebuilt from each round's
   resampled vectors, and sample counts adjusted for evaluation failures.

3. Rubin DB. Using the SIR algorithm to simulate posterior distributions. In:
   *Bayesian Statistics 3*. Oxford University Press; 1988:395-402. The
   origin of sampling importance resampling as a general technique, on which
   the pharmacometric application above builds.

4. Box GEP, Cox DR. An analysis of transformations. *Journal of the Royal
   Statistical Society, Series B*. 1964;26:211-252.
   [doi:10.1111/j.2517-6161.1964.tb00553.x](https://doi.org/10.1111/j.2517-6161.1964.tb00553.x).
   The transformation used to make the inter-iteration proposal closer to
   normal.

5. Lindbom L, Ribbing J, Jonsson EN. Perl-speaks-NONMEM (PsN) — a Perl module
   for NONMEM related programming. *Computer Methods and Programs in
   Biomedicine*. 2004;75:85-94.
   [doi:10.1016/j.cmpb.2003.11.003](https://doi.org/10.1016/j.cmpb.2003.11.003).
   The reference implementation this package is checked against.

## Implementation index

- Orchestration, persistence, and return construction:
  [`runSIR()`](../R/sir-run.R)
- Control validation: [`runSIRControl()`](../R/sir-control.R)
- Parameter vector and naming bridge:
  [`.sirParamSpace()`, `.sirProposalMu()`](../R/sir-paramspace.R)
- Proposal construction, fallback uncertainty, inflation, and sampling:
  [`sirGetProposalCov()`, `.sirFallbackSe()`, `.sirInitialProposal()`,
  `.sirInflationVector()`, `.sirSampleFullProposal()`](../R/sir-proposal.R)
- Alternative proposal sources:
  [`.sirRseVariance()`, `.sirProposalFromCovmatInput()`,
  `.sirProposalFromRawResults()`,
  `.sirResolveInitialProposal()`](../R/sir-proposal-input.R)
- Objective function evaluation: [`sirEvalOFV()`](../R/sir-eval.R)
- Weights, resampling, and per-iteration raw results:
  [`sirCalcWeights()`, `sirResample()`,
  `.sirBuildRawResults()`](../R/sir-weights.R)
- Box-Cox and proposal update:
  [`sirBoxCox()`, `sirUpdateProposal()`](../R/sir-boxcox.R)
- One iteration, and the PsN sample-count adjustments:
  [`sirRunIteration()`, `.sirAdjustedAttemptedSamples()`,
  `.sirAdjustedResamples()`](../R/sir-iterate.R)
- Convergence diagnostic:
  [`.sirDofvCurves()`, `.sirDofvNoise()`,
  `.sirProposalTooNarrow()`](../R/sir-convergence.R)
- Interval and RSE/correlation diagnostics:
  [`.sirIterationIntervals()`, `.sirRseCorData()`](../R/sir-diagnostics.R)
- Summaries and on-disk artifacts:
  [`sirSummary()`, `.sirWriteIterationSummary()`,
  `.sirWriteCovMatrices()`](../R/sir-results.R)
- `setCov()` registration:
  [`.sirCovAsFitCov()`, `.sirRegisterCov()`](../R/sir-setcov.R)
- Matrix helpers: [`.sirCapCovCorrelation()`, `.sirEnsurePosDef()`,
  `.sirRound()`](../R/sir-utils.R)
- Print and plot methods: [`R/sir-methods.R`](../R/sir-methods.R)
