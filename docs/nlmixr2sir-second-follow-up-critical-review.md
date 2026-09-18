# Second follow-up critical review

**Date:** 2026-09-17

**Reviewed range:** `85518d1..1b72a80`

This review covers the evaluator-control remediation, tracking of the previous
follow-up review, and the redesign of the importance-weight degeneracy warning.

## Summary

**Verdict: Request changes.**

The previous blocking AGQ defect is fixed correctly. The evaluator preserves
the fitted control, the off-center AGQ regression coverage is substantive, and
the complete test suite passes.

No blocking correctness issues remain, but the new ESS guidance contains
misleading advice and one test does not exercise its claimed branch.

## Critical issues (blocking)

None.

## Required changes

### 1. Low efficiency does not imply that the proposal should be widened

[`R/sir-weight-diagnostics.R`](../R/sir-weight-diagnostics.R) always recommends
widening when the ESS fraction falls below its threshold. ESS/$n$ measures the
magnitude of proposal-target mismatch, not its direction.

For a standard-normal target, the following exact asymptotic efficiencies are
obtained from centered normal proposals:

| Proposal SD | Asymptotic ESS fraction |
| ---: | ---: |
| 10 | 14.11% |
| 20 | 7.07% |
| 40 | 3.53% |

At SD 20 the warning fires, but widening to SD 40 makes efficiency worse.
The warning should recommend adjusting the proposal's scale or center using
the convergence diagnostics. It should recommend widening specifically only
when those diagnostics show that the proposal is too narrow.

### 2. The large-but-inefficient test does not enter its advertised branch

The test named “absolute ESS is what a large but inefficient run is judged on”
in
[`tests/testthat/test-sir-weight-diagnostics.R`](../tests/testthat/test-sir-weight-diagnostics.R)
constructs an ESS of approximately 402 and an ESS fraction of approximately
20.1%. The fraction is above the 10% warning threshold, so the conditional
accepts the no-warning path despite the comment claiming that the efficiency
warning fires.

Replace it with weights guaranteed to satisfy both premises:

```r
ess > .sirEssWarn
essFraction < .sirEssFractionWarn
```

Assert those premises and the exact warning content without an `if` branch.

### 3. Tracked documentation contradicts the repaired implementation

- [`docs/nlmixr2sir-follow-up-critical-review.md`](nlmixr2sir-follow-up-critical-review.md)
  is not pinned to the reviewed commit and describes the now-fixed AGQ defect
  as currently blocking. Pin it to `85518d1` and add a resolution notice
  referencing `f6bc538`.
- [`NEWS.md`](../NEWS.md) still says candidates are scored with a fresh FOCEi
  evaluation, although the evaluator now uses the fit's own method.
- [`tests/testthat/setup.R`](../tests/testthat/setup.R) retains the obsolete
  explanation that larger fixtures necessarily warn harder. That explanation
  described the old fraction-only warning, not the new absolute-ESS criterion.

## Strong suggestions

### Fail closed when evaluator control reconstruction is impossible

[`R/sir-eval.R`](../R/sir-eval.R) still has two fallback paths:

- failure to build the method-specific reference control falls back to the raw,
  explicitly unnormalized override values; and
- a fit without a usable control falls back to method defaults.

Either path can recreate the objective-surface mismatch that the AGQ fix was
designed to eliminate. Abort instead when the fitted control or normalized
override values cannot be recovered.

### Describe the fixed ESS cutoffs as heuristics

The thresholds of 100 effective samples and 200 input samples are useful
diagnostic conventions, not universal reliability boundaries. Reliability
depends on the weight distribution and its tail behavior. The documentation
should label the cutoffs accordingly and avoid implying a general statistical
guarantee.

## Verification

- Targeted objective-preflight, AGQ, and weight-diagnostic tests passed.
- The complete `devtools::test()` suite passed.
- The prior blocking AGQ reproduction is now protected by an off-center test
  with non-default integration bounds.
- `git diff --check 85518d1..1b72a80` reports trailing whitespace and an extra
  end-of-file blank line in the first follow-up review document.

