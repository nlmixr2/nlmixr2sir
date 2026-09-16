# Follow-up critical review of the latest updates

**Date:** 2026-09-16  
**Scope:** Commits after `origin/main` reviewed on 2026-09-15, including the
expanded deterministic estimation-method support, the single-quantile dOFV
noise fix, the PsN comparison updates, and the IMP-family rationale.

## Summary

**Verdict: Request changes.**

The IMP-family rationale and segfault issue are sound, and the full test suite
passes, but the expanded AGQ support can silently calculate incorrect
importance weights.

## Critical issues (blocking)

### 1. AGQ evaluation drops likelihood bounds

[`R/sir-eval.R`](../R/sir-eval.R) omits `agqLow` and `agqHi` from
`.sirLikelihoodControlFields`, while
[`R/sir-preflight.R`](../R/sir-preflight.R) includes `agq` among the supported
estimation methods. Consequently, `.sirEvalControl()` reconstructs an AGQ
control with the default bounds, `-Inf` and `Inf`, rather than the bounds used
by the fitted model.

This is a demonstrated objective-surface mismatch, not merely missing
plumbing. With a real fit using
`agqControl(nAGQ = 2, agqLow = -100, agqHi = Inf)`:

| Check | Result |
| --- | ---: |
| Stored fit OFV | 188.378 |
| Preflight absolute difference | 1.216185e-08 |
| Preflight stencil minimum dOFV | -4.118186e-05 |
| Current evaluator OFV at `tka = -20` | 3439.366 |
| OFV with the fitted AGQ bounds preserved | 2317.800 |
| Difference | **1121.566** |

The center and stencil preflight therefore pass at their default tolerances,
while an off-center candidate is scored on a materially different surface.
That corrupts its dOFV and hence its importance weight.

Required remediation:

1. Carry `agqLow` and `agqHi` into the fixed-parameter evaluator.
2. Audit every method-specific likelihood control for all newly admitted
   methods rather than assuming the shared FOCEi fields are exhaustive.
3. Add a regression test using non-default AGQ bounds and an off-center
   candidate.
4. Remove `agq` from the allowlist until this is complete if the broader
   control audit cannot be completed immediately.

## Required changes

### 2. User-facing documentation still describes FOCEi-only behavior

The implementation now accepts the deterministic conditional-estimation
ladder, but several passages still claim that only `focei` is accepted or that
every candidate is evaluated through a hard-coded FOCEi call:

- [`README.md`](../README.md), including the parity table and the “Only
  `focei` fits are accepted” section.
- [`docs/sir-technical-reference.md`](sir-technical-reference.md), including
  the objective-evaluation description, supported-method section, and PsN
  comparison.
- [`NEWS.md`](../NEWS.md), which contains multiple FOCEi-only entries and no
  correct release note for the expanded support.
- The unsupported-method diagnostic in
  [`R/sir-preflight.R`](../R/sir-preflight.R), which still says that candidates
  are evaluated with `est = "focei"`.

These are direct contradictions of the implementation and must be corrected
before release.

### 3. The allowlist test does not validate the allowlist

The test named “the deterministic ladder is accepted, not just focei” in
[`tests/testthat/test-sir-objective-preflight.R`](../tests/testthat/test-sir-objective-preflight.R)
mostly checks that strings are present in `.sirSupportedEstimationMethods`.
That is a tautological test of the constant, not evidence that the evaluator
reproduces those methods' objectives.

Real objective evaluation currently covers FO and FOCEi. The other admitted
methods and their method-specific controls are not protected by durable
regression tests, which is why the AGQ defect escaped.

Add tests for each distinct likelihood configuration, at minimum:

- FO and FOI;
- FOCE, FOCEI, and FOCEP;
- Laplace;
- AGQ, including non-default `agqLow` and `agqHi`;
- the distinct `muModel` variants; and
- off-center candidate values, not only reproduction at the fitted center.

## Verification performed

- Full `devtools::test()` suite: passed.
- Targeted convergence and objective-preflight tests: passed.
- `git diff --check`: clean.
- No repository files were modified during the review itself.

A green suite does not reduce the severity of the AGQ finding: the failing
configuration lies outside the present test coverage and was reproduced
directly.

## IMP-family issue draft

The rationale in [`ISSUE-nlmixr2est-segfault.md`](../ISSUE-nlmixr2est-segfault.md)
is technically coherent:

- `nIter = 0` reproducibly crashes the IMP-family methods rather than issuing
  a controlled error;
- positive `nIter` values perform full EM iterations and therefore move the
  population parameters;
- no exposed expectation-only control equivalent to NONMEM/PsN `EONLY=1` was
  identified; and
- the disagreement between `fit$objf` and `fit$env$impObj` is correctly noted
  as an additional design question.

No changes to that issue draft are required by this review.

