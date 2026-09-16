# nlmixr2sir (development version)

## Diagnostics and provenance

* **Every iteration now reports importance-weight degeneracy.** Effective
  sample size (Kish, `1 / sum(p^2)`), its fraction of the usable samples, the
  largest single weight, perplexity, and the count of non-negligible weights
  are added to the iteration summary. `runSIR()` warns when the effective
  sample size falls below 10% of the usable samples, or one candidate carries
  more than half the weight. Log-scale normalization stops the weights
  overflowing but says nothing about whether they are informative: a run resting
  on two candidates normalizes perfectly well and previously reported nothing
  unusual.

* **The convergence-plot noise band now uses the run's own resampler.** It
  called a local `sample.int()` that was never given the run's `capResampling`
  and treated any cap above one as unlimited replacement, whereas
  `sirResample()` expands each candidate into a finite number of slots. The
  band therefore described a different algorithm from the one that produced the
  retained sample. Replicates now go through `sirResample()` with the run's
  actual cap.

* **Results carry their full provenance.** The effective control object,
  schedule, proposal source, reference-OFV history, run fingerprint, and
  whether any covariance needed positive-definite repair are attached to the
  returned object and persisted in the run state, alongside the existing seed
  and `sir_manifest.dcf` manifest.

* A frozen end-to-end fixture (`test-sir-golden-iteration.R`) pins the whole
  component chain -- proposal density, Box-Cox Jacobian, importance ratios,
  normalized weights, resampling properties, proposal update, and weight
  diagnostics -- against values computed once from synthetic inputs. It is a
  regression fixture, not a PsN comparison; the genuine PsN oracle values remain
  in `test-sir-psn-oracles.R`.

* Terminology: "SIR posterior" is now "retained SIR distribution" throughout the
  code and documentation, matching the technical reference's statement that the
  output is not Bayesian. The technical reference and README now set out, next
  to the convergence diagnostic itself, the conditions under which the
  chi-square reference can mislead -- boundary variance components, weak
  identification, multimodality, non-smooth likelihoods, small subject counts --
  and note that percentile intervals from a likelihood-weighted sample do not
  automatically have nominal frequentist coverage.

## Run identity, feasibility and robustness

* **Added iterations now keep a cumulative schedule.** `addIterations = TRUE`
  stored only the extension schedule, so a two-iteration run extended by one
  saved three completed iterations beside an identity describing a single
  iteration. A later plain `recover` request presenting that one-iteration
  schedule matched and was handed the three-iteration result. The stored
  schedule, the result's `schedule` attribute and the saved fingerprint now all
  describe every iteration the result contains.

* **Infeasible schedules are rejected before any model is evaluated.**
  `runSIR()` checked only `nResample <= nParameters`. Two further
  impossibilities are knowable in advance: `nSamples <= nParameters` (a
  candidate must be drawn before it can be retained), and
  `nResample > nSamples * capResampling` (under limited replacement each
  candidate fills at most `capResampling` slots). The first previously surfaced
  only after a full round of model evaluations, advising that `nResample` be
  raised when `nResample` was not the problem; the second was silently clamped.

* **An infeasible retained set now aborts naming its cause.** Failed
  evaluations and the turnout and cap clamps can leave fewer usable candidates
  than the proposal update needs, on a schedule that was feasible as requested.
  This reached the rank check, which saw only the final count and advised
  raising `nResample`. It now reports how many candidates could actually be
  scored.

* **Non-finite inverse Box-Cox results are classified as inverse failures.**
  The filter used `is.na()`, and inversion can overflow to `Inf` without
  raising an error. Since `Inf > Inf` is `FALSE`, such a value also passed the
  bounds test, so it was charged to whichever check happened to reject it next
  -- or, absent one, entered the sample and reached objective evaluation.

* **The Box-Cox shift now protects the centre it will transform.** The shift
  was chosen from the retained sample alone. Under `recenter = TRUE` the next
  centre is the best candidate, which need not be among the retained rows and
  can lie below their minimum, making `mu + delta <= 0` and aborting a run that
  was proceeding normally.

* **Initial proposal covariance repair is recorded.** Only the repair of later
  empirical updates was stored. The repair of the covariance that iteration 1
  draws from is now kept on the result and in the run state, with its method,
  threshold and magnitude, and survives a resume. The iteration summary
  distinguishes `proposalRepaired` (the covariance that iteration drew from)
  from `posDefAdjusted` (the empirical update it produced for the next one).

* The convergence noise band is built from normalized resampling probabilities
  rather than raw importance ratios. `importance_ratio = exp(log_ir)` can
  overflow to `Inf` on a strongly favoured candidate, and filtering on it
  discarded exactly the candidate carrying the weight.

## Documentation

* **Corrected an incorrect claim about parameterization invariance.** The
  README and technical reference stated that including the Box-Cox
  change-of-variables Jacobian makes the retained distribution invariant to
  re-expressing the *model* in another smooth parameterization. It does not.
  The Jacobian makes Box-Cox an internal *proposal* transformation, so the
  result does not depend on that transform; but flat normalized likelihood is a
  choice of base measure, and two models that are reparameterizations of each
  other can give different SIR intervals. Both documents now say so.

* The technical reference described a `sqrt(.Machine$double.eps)` eigenvalue
  floor that the code no longer uses, and attributed the narrow-proposal
  warning to `runSIR()` when it is emitted by `plot(type = "convergence")`.
  Both corrected.

* **`runSIR()` now accepts the deterministic estimation ladder, not only
  `focei`**: `fo`, `foi`, `foce`, `focei`, `focep`, `laplace`, `agq`, and the
  mu-referencing `m...`/`i...` variants of each. Candidates are scored by
  re-evaluating the fit's own method at fixed population parameters, so the
  candidate surface and the `fit$objf` reference are the same function. Each
  method was verified by re-evaluating a real fit at its own estimates rather
  than by being listed.

* **The evaluator now carries the fit's whole control object**, overriding only
  the evaluation fields (`maxOuterIterations`, `calcTables`, `covMethod`,
  `compress`, `print`). It previously rebuilt the control from a hand-picked
  list of likelihood-relevant settings, which failed open: `foceiControl()` has
  150 arguments, and one that was not on the list was silently dropped.
  `agqLow` and `agqHi` were lost that way, so an AGQ fit with non-default
  integration bounds agreed with its own objective at the centre to 1e-08 but
  differed by about 6490 OFV units at an off-centre candidate -- a difference
  that enters the importance weight as `exp(-dOFV/2)`. Neither the centre check
  nor the stencil can detect this: the stencil perturbs by a thousandth of each
  estimate, and integration bounds only bite far from the mode.

* **The restriction that remains is documented as a deliberate
  difference from PsN**, with the measurement behind it: a SAEM fit of
  `theo_sd` stores an objective of 208.512 against 205.820 from FOCEi
  re-evaluation at the same estimates, and that 2.69-unit gap does not cancel
  under `recenter = TRUE`. PsN's `set_maxeval_zero()` handles classical methods
  and `IMP`/`IMPMAP` and only warns for anything else -- and discards the
  failure flag, so its SIR proceeds. `runSIR()` aborts in the preflight
  instead. The absence of any stochastic-evaluator path, and what its sampling
  noise would do to the importance weights, is documented alongside it.

* PsN is now described as a comparator and a source of numerical oracles rather
  than as the specification, and the claim to enumerate every deliberate
  difference is replaced. The draw-attempt budget (`10 * nSamples` against
  PsN's `2000 * nSamples`) and the absent OMEGA/SIGMA block adjustment are
  documented as deliberate choices.

## Correctness fixes

* **Importance weights were wrong for correlated proposals.** `sirCalcWeights()`
  solved against R's Cholesky factor instead of its transpose, giving the wrong
  Mahalanobis distance and so the wrong proposal density for every candidate.
  On PsN's own `mvnpdf_cholesky` oracle the relative density was
  `0.0772293088557175` where PsN and the explicit Mahalanobis form both give
  `0.0837378551174778`. Because `fit$cov` is the default proposal source and
  population covariances are essentially always correlated, this affected the
  main path: **retained samples, intervals, RSEs, and every diagnostic derived
  from them are wrong in runs made with earlier versions and should be rerun.**
  The existing tests could not catch it because all of them used a diagonal
  covariance, for which a Cholesky factor equals its own transpose; PsN's
  density oracle and a `mvtnorm` cross-check are now in the suite.

* **Box-Cox proposals now include the change-of-variables Jacobian**, and this
  is a deliberate divergence from PsN. Candidates are drawn on a transformed
  scale and mapped back, so the density induced on the original parameter scale
  is `q_x(x) = q_y(T(x)) * |det J_T(x)|`; the weight now divides by `q_x`.
  Omitting the term, as PsN does, retains a sample from `L(x)|det J_T(x)|`, so
  the answer depends on which smooth parameterization the model is written in --
  worst for the skewed and weakly identified parameters Box-Cox exists to help
  with. Importance-sampling a known Gamma(3, 1) target through a Box-Cox
  proposal recovers a mean of 3.00 and second moment 12.00 with the Jacobian
  (truth 3 and 12), against 2.25 and 7.31 without it. SIR now targets the
  normalized likelihood on nlmixr2's own parameter scale, and the retained
  distribution is invariant to re-expressing the model in another smooth
  parameterization. `boxcox = TRUE` is the default, so **Box-Cox runs made with
  earlier versions should be rerun**, and `-boxcox` is now marked *partial* in
  the README parity matrix.

* **Recovery now checks that the saved run is the run being asked for.**
  `recover = TRUE` is the default, and a completed result used to be returned
  solely because a state file existed in the directory and recorded enough
  iterations. Nothing verified it belonged to the fit in hand, so pointing a
  different model, dataset, schedule, or set of statistical controls at the same
  directory returned a stale result labelled as the new run's. State now carries
  a run fingerprint -- model, data, parameter set, estimates, objective,
  estimation method, schedule, statistical controls, and a state-format version
  -- and a mismatch aborts naming the fields that changed. `addIterations`
  exempts the schedule, which it deliberately changes. Parallelism settings are
  excluded: they do not change the answer, so they must not invalidate a run.

* **`runSIR()` will no longer delete a directory it does not recognise.** An
  explicitly supplied directory in overwrite mode was removed with
  `unlink(recursive = TRUE)` with no check of what it contained. Directories
  created by `runSIR()` now carry a `sir_manifest.dcf` file, which records what
  produced them and is also what permits them to be cleared. A non-empty
  directory without one is refused.

* **`runSIRControl(saveFiles = FALSE)` runs entirely in memory.** No directory
  is created and nothing is written; the result is returned as usual and
  `setCov()` registration still happens. Recovery, `addIterations`, and
  per-iteration seed reproduction all need the saved state, so they are
  unavailable -- seed such a run with `set.seed()` beforehand to reproduce it.
  `addIterations = TRUE` with `saveFiles = FALSE` is rejected by
  `runSIRControl()`.

* **Rank-deficient proposals are now an error rather than silently
  regularized.** The empirical covariance of `m` retained vectors in `p`
  dimensions has rank at most `m - 1`, so a full-rank proposal needs
  `m > p`. Forcing such a matrix positive definite does not recover the missing
  information -- it fabricates variance in directions the retained sample never
  supported, and the next iteration then proposes along them. `runSIR()` now
  rejects `nResample <= <number of estimated parameters>` before any model is
  evaluated, and the proposal update rejects a retained sample whose numerical
  rank is deficient (repeated or collinear draws). PsN stops here too.

* **The positive-definite repair is now scale-relative.** Eigenvalues were
  floored at the *absolute* value `sqrt(.Machine$double.eps)`, which is not
  scale-equivariant: on one and the same singular problem expressed in
  different units, that floor was 8.7e-3 of the largest eigenvalue at one scale
  and 8.7e-15 at another -- dominant in one parameterization, negligible in
  another. Pharmacometric parameters genuinely span those scales. The floor is
  now relative to the matrix's own largest eigenvalue, and exists only to clean
  up floating-point roundoff on a covariance that is already full rank.

* A raw-results proposal with fewer vectors than parameters is now an error.
  It previously warned and continued with a forced positive-definite matrix.

* When fewer samples have non-zero resampling probability than the requested
  resample count, the count is reduced with a warning naming the cause, instead
  of aborting inside the resampler with a message about probabilities. A finite
  dOFV does not guarantee a finite importance ratio, so the turnout adjustment
  alone could leave the count too high.

* **The documented automatic covariance fallback now matches what the code
  does.** The README and technical reference said a failed covariance step or
  `covMethod = ""` fell back automatically to the Wishart-style OMEGA
  approximation. It does not, and cannot: that approximation needs only the
  OMEGA estimates and the subject count, but such a fit carries no THETA
  uncertainty at all, and assuming a THETA scale would fabricate the quantity
  SIR reports. The fallback completes an *incomplete* `fit$cov` (`covFull =
  FALSE`, or a partial covariance) and that path is now covered end to end by
  tests. An *absent* `fit$cov` stops the run, and the error now names
  `rseTheta`, `covmatInput` and `rawresInput` instead of advising a re-run with
  a covariance step -- unhelpful guidance for the models SIR exists to serve.

* `runSIR()` now verifies before sampling that it can reproduce the fit's own
  objective at the fit's own estimates, and aborts naming both values if not.
  Importance sampling assumes one fixed target; candidates are scored by a
  fresh FOCEi evaluation, which is not guaranteed to be the surface that
  produced `fit$objf`. The tolerance is `runSIRControl(objfTolerance =)`.

* `runSIR()` rejects fits whose objective it cannot reproduce, and says so before doing any work.
  A SAEM fit's objective comes from Gaussian quadrature: on `theo_sd` it is
  208.512 against 205.820 from a FOCEi reevaluation, a 2.69 unit gap. That is a
  different likelihood, not numerical noise, and with `recenter = TRUE` the run
  would have "found" a better optimum from the offset alone.

* Recentring now moves the dOFV reference as well as the proposal centre, and
  the reference is persisted in the run state so recovery and `addIterations`
  resume from it. Previously every dOFV stayed measured against the original
  `fit$objf`, leaving later iterations, the negative-dOFV counts, and the
  chi-square convergence diagnostic pinned to an optimum the run had already
  superseded. Within a single iteration nothing changes, because a constant
  dOFV shift cancels in the normalised weights.

* Empirical proposal summaries no longer include the synthetic centre row that
  `runSIR()` writes to the raw results. Raw results gain a `role` column
  (`"reference"` / `"sample"`), and the dOFV curves, interval, and RSE
  diagnostics all exclude the reference row, as PsN does. The centre is not a
  draw from the proposal, and counting it biased quantiles, intervals,
  covariances, and RSEs — by up to 2.6 percentage points of RSE on an
  eight-sample iteration.

* `sirSummary()` reports `rse_sd_scale` only for OMEGA diagonals, and `NA` for
  off-diagonals. The delta-method relation `RSE(sqrt(v)) ~= RSE(v)/2` needs a
  positive variance; an off-diagonal is a covariance, which can be negative or
  zero and has no standard-deviation counterpart. Use the empirical
  correlations in the `sdCorMatrix` attribute instead.

# nlmixr2sir 0.3

* SIR now derives its parameter vector from a single internal description of the fit instead of re-deriving THETA, sigma, and OMEGA names independently in each function. This fixes `runSIR()` aborting with `Assertion on 'mu' failed: Contains missing values` against nlmixr2est 7, where `fit$cov` reports OMEGA alongside THETA and the old proposal mean came back `NA` for every OMEGA element.

* OFV evaluation no longer fails in a session where `nlmixr2sir` was attached but `nlmixr2` was not. `rxode2::ini()` evaluates the OMEGA line it is given as a `lotri({...})` call in the caller's environment, and `lotri` is an Imports rather than a Depends of both rxode2 and nlmixr2est, so it was never visible; `nlmixr2sir` now imports it. Every sample previously returned `NA` and `runSIR()` aborted with `All SIR OFV evaluations failed`.

* That abort now reports the first underlying evaluation error, instead of discarding it. A configuration problem fails every sample identically, which the old message hid.

* The package now requires nlmixr2est >= 7.0.0, nlmixr2utils >= 0.3, and rxode2 >= 5.0.0.

* `runSIR()` now takes its run settings through `control = runSIRControl()` rather than as flat arguments, cutting its signature from 21 arguments to 6. Passing a setting directly to `runSIR()` is an error naming `runSIRControl()`.

* `runSIR()` registers its empirical covariance with the fit, so `nlmixr2est::setCov(fit, "sir")` switches the fit's reported uncertainty to the SIR result. Registration is skipped, with a message, if the parameters do not match `fit$cov` exactly or the covariance is not positive definite.

* `runSIRControl()` objects deparse back into reproducible source through `rxode2::rxUiDeparse()`, emitting only the arguments that differ from the defaults.

* `plot()` gains three diagnostics. `type = "convergence"` is the dOFV-versus-chi-square plot: per iteration, the empirical dOFV quantile curve for the proposal and for the SIR posterior against a reference chi-square on the number of estimated parameters, with a resampling-noise band on the last two iterations. Convergence reads as the SIR curve settling onto the reference, and a proposal that falls below the reference for more than a quarter of the quantiles now warns and recommends inflation. `type = "intervals"` compares the proposal and SIR interval per parameter per iteration. `type = "rsecor"` draws the RSE/correlation matrix with the diagonal annotated by the confidence-interval asymmetry ratio, which a symmetric normal approximation cannot show.

* `sirSummary()` now matches PsN's output: it adds `mean` alongside the median, reports PsN's percentile set (2.5, 5, 10, 30, 50, 70, 90, 95, 97.5, from prediction intervals 0/40/80/90/95) in place of the previous set, and adds `rse_sd_scale`. Note `p25` and `p75` are no longer reported. `rse` remains a percentage where PsN reports a fraction, and the returned object now records that in an `rseUnits` attribute.

* `runSIR()` writes `<fitName>_sir.cov` and `<fitName>_sir.sdcorr`, and `summary_iterations.csv` now leads with PsN's column names so either file can be read by PsN-literate tooling.

* `runSIR()` no longer requires a successful covariance step. `runSIRControl()` gains three alternative proposal sources, following PsN: `rseTheta`/`rseOmega`/`rseSigma` build a diagonal proposal from relative standard errors, `covmatInput` takes a covariance matrix, a NONMEM-style `.cov` file, or `"identity"`, and `rawresInput` seeds the first proposal from the parameter vectors in a canonical raw-results file (with `offsetRawres` and `inFilter`). These exist for the models whose covariance step fails, which is where SIR is most wanted.

* `runSIRControl()` inflation arguments accept a vector as well as a scalar: one value per estimated THETA, OMEGA diagonal, or residual-error parameter. OMEGA off-diagonals derive `sqrt(infl_i) * sqrt(infl_j)` from the two diagonals they connect, as PsN does.

* The sample and resample count adjustments now match PsN exactly, and are checked against the oracle values in PsN's own unit tests. The attempted-sample count previously used `ceiling()` where PsN uses round-half-away-from-zero, and the resample count used `floor()` and a strict `>` where PsN uses rounding and `>=`; both also clamped in ways PsN does not.

* `runSIR()` gains `rxThreads`, controlling rxode2 OpenMP threads per worker during OFV evaluation. `nlmixr2utils` 0.3 requires it whenever `workers > 1`, so parallel SIR runs previously aborted on most multicore machines.

* `runSIR()` now applies `sigmaInflation` to residual-error parameters. They were classified as THETA, which made the argument silently unreachable.

* `runSIR()` now takes OMEGA uncertainty from `fit$cov` by default, via the new `omegaFallback = "cov"`. The initial proposal therefore carries the correlations between THETA and OMEGA, which the previous block-diagonal construction discarded even when they were available. `omegaFallback = "wishart"` forces the old approximation, and it is still used automatically when `fit$cov` does not carry OMEGA.

# nlmixr2sir 0.2

* `runSIR()` now writes a final canonical shared `raw_results.*` artifact, uses the common `nlmixr2utils` run-state and seeding helpers, and no longer emits per-iteration raw-results CSV files.

# nlmixr2sir 0.1

* Initial package split from `nlmixr2extra`, providing `runSIR()`,
  `sirSummary()`, S3 print/plot methods, tests, and the SIR vignette as a
  standalone package depending on `nlmixr2utils`.
