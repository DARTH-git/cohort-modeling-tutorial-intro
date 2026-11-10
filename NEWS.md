# cohort-modeling-tutorial-intro 0.3.2 (2026-11-10)

## Bug fixes

* `check_sum_of_transition_array()` (`R/Functions.R`) compared the sum of each
  cycle's transition probabilities in a transition *array* with the number of
  states using exact equality, so rounding error could make it reject a valid
  array. It now uses a numerical tolerance (`all.equal()`), as in darthtools.
  This tutorial only checks transition matrices, so its results are not
  affected; the same bug stopped the PSA of the state-residence model in
  [cohort-modeling-tutorial-timedep](https://github.com/DARTH-git/cohort-modeling-tutorial-timedep)
  (fixed there in v0.3.0).

## Documentation

* Documented an erratum in the published article's parameter table (Table 2),
  which lists the cost of treatment B (`c_trtB`) as $12,000 instead of the
  $13,000 stated in the text and used in the code, the PSA and the published
  results. The code was already correct. See `manuscript/README.md` and the
  README (closes #1). Reported by @conallaque.

## Tests

* Added `tests/testthat/test-psa-params.R`, which checks that every PSA
  distribution is centred on its base-case value (within 2%), so the base
  case and the PSA cannot drift apart.
* Added `tests/testthat/test-transition-checks.R` for the transition matrix
  and array checks.

# cohort-modeling-tutorial-intro 0.3.1 (2026-11-10)

## Documentation

* `darthtools` v0.4.1 fixes the same `gen_wcc()` indexing error
  ([DARTH-git/darthtools#16](https://github.com/DARTH-git/darthtools/pull/16)).
  The README and code comments now point to it. The tutorial keeps its own
  copy of `gen_wcc()` in `R/Functions.R`, so it still runs without
  `darthtools`.
* Corrected the description of the old weights in NEWS.md and in the
  `gen_wcc()` documentation: they summed to n - 2/3 for an even number of
  cycles and to n - 1/3 for an odd number (such as the tutorial's 75), not
  always to n - 2/3.
* Added ORCIDs for Fernando Alarid-Escudero and Alan Yang to `CITATION.cff`
  and `DESCRIPTION`.

# cohort-modeling-tutorial-intro 0.3.0 (2026-11-10)

This release fixes a bug in the Simpson's 1/3rd within-cycle correction (WCC).
**It changes the numerical results reported in the published article**, but
not its conclusions.

## Bug fixes

* `gen_wcc(method = "Simpson1/3")` (`R/Functions.R`) assigned the 4/3 and 2/3
  weights by position in the 1-based R vector instead of by time
  t = 0, 1, ..., n. It returned `1/3, 2/3, 4/3, 2/3, ..., 1/3` instead of
  `1/3, 4/3, 2/3, 4/3, ..., 1/3`, so the weights summed to n - 2/3 (even n)
  or n - 1/3 (odd n) instead of n. For example, `gen_wcc(4)` integrated f(t) = 1 over [0, 4] as 3.33
  instead of 4. Weights are now indexed by time (closes #2).

  The same error is in the WCC description, the equation and the
  `within-cycle-vector` code chunk of the published manuscript. The files in
  `manuscript/` are kept unchanged as the archival source of the article (see
  `manuscript/README.md`).

  Effect on the base case of the time-independent Sick-Sicker model
  (`analysis/cSTM_time_indep.R`, 75 annual cycles, 3% discount):

  | Strategy         | QALYs v0.2.1 | QALYs v0.3.0 | Costs v0.2.1 | Costs v0.3.0 | ICER v0.2.1 | ICER v0.3.0 |
  |------------------|-------------:|-------------:|-------------:|-------------:|------------:|------------:|
  | Standard of care |        20.71 |        21.03 |     $151,580 |     $152,282 |          -- |          -- |
  | Strategy B       |        22.18 |        22.50 |     $259,100 |     $260,222 |     $72,988 |     $73,279 |
  | Strategy AB      |        23.14 |        23.46 |     $378,875 |     $380,378 |    $125,764 |    $125,302 |
  | Strategy A       |        21.50 |        21.83 |     $284,805 |     $285,889 |   Dominated |   Dominated |

  QALYs increase by about 0.32 per strategy and ICERs change by about 0.4%.
  The efficient frontier and the dominance of Strategy A are unchanged.

* The PSA model wrapper (`calculate_ce_out()` in
  `R/Functions_cSTM_time_indep.R`) now uses the corrected local `gen_wcc()`.
  It used to call `darthtools::gen_wcc()`, which has the same bug in
  darthtools <= 0.4.0. The deterministic and probabilistic analyses now use
  the same WCC, and the repository no longer needs `darthtools`.

* `calculate_ce_out()` had the discount rates swapped: costs were discounted
  with `d_e` and QALYs with `d_c`. This had no effect on the results because
  both rates are 3%, but it would have mattered for differential discounting.

* In `analysis/cSTM_time_indep.R`, the state costs under Strategy A
  (`v_c_strA`) were not scaled by `cycle_length`, unlike every other reward
  vector. This had no effect with annual cycles but would have mattered with
  monthly cycles.

## Changes in behaviour

* `gen_wcc(method = "Simpson1/3")` now supports an odd number of cycles, like
  the tutorial's 75 cycles. Simpson's 1/3rd rule needs an even number of
  intervals, so for odd `n_cycles` the composite 1/3 rule is applied over the
  first n - 3 cycles and Simpson's 3/8 rule over the last three. Both cases
  integrate polynomials up to degree three exactly, and the weights always
  sum to `n_cycles`.
* `gen_wcc()` validates its input. It now errors when `n_cycles` is not a
  single positive integer (it used to truncate, e.g., 7.5 to 7) and when
  `n_cycles = 1` is used with Simpson's rule.

## Infrastructure and reproducibility

* Added unit tests (`tests/testthat/`) for `gen_wcc()`: exact weights, exact
  integration of polynomials of degree 0-3 for odd and even n, and input
  validation. Added a regression test for the base-case CEA. Run them with
  `testthat::test_dir("tests/testthat")`.
* Added a GitHub Actions workflow (`.github/workflows/check.yaml`). It runs the
  tests and the full tutorial analysis on Ubuntu and macOS for every push and
  pull request.
* Added `renv` (`renv.lock`) to pin package versions. Run `renv::restore()` to
  reproduce the computational environment.
* Added a `DESCRIPTION` file that records the project version and
  dependencies (the repository is not an installable package).
* Added `CITATION.cff` with all authors and the preferred citation (the MDM
  article). This enables GitHub's "Cite this repository" button and complete
  Zenodo metadata.
* Added this `NEWS.md` changelog.
* Replaced deprecated `ggplot2` arguments (`size` -> `linewidth` in
  `geom_line()`, `..density..` -> `after_stat(density)`).
* Removed a stray, unnamed copy of `rate_to_prob()` from `R/Functions.R`.

## Documentation

* Rewrote the `gen_wcc()` documentation to describe the weights by time
  ("entries at odd times t = 1, 3, ... get 4/3") and to document the odd-n
  extension.
* Updated and expanded the README: a corrections notice, quick start with
  `renv`, repository structure, testing, how to report issues, and fixed
  typos and citation details (including the Zenodo DOIs).
* Added `manuscript/README.md` explaining how the archival manuscript sources
  differ from the corrected code.

## Acknowledgements

This bug was found by two attentive users of the tutorial, and we are very
grateful to both:

* **Lauren N. Steimle** (H. Milton Stewart School of Industrial and Systems
  Engineering, Georgia Institute of Technology) found the indexing error while
  preparing lecture notes on within-cycle corrections for her graduate course.
  She reported it by email in January 2025, with a proposed correction.
* **[@ambtnhlr](https://github.com/ambtnhlr)** independently reported the same
  error in #2. The report included a reproducible example showing that the
  weights did not integrate a constant correctly, a time-indexed fix, a
  suggested rewording of the documentation, and the observation that 75
  cycles is an odd number of intervals for Simpson's 1/3rd rule.

Their careful reports made this release possible.

# cohort-modeling-tutorial-intro 0.2.1 (2022-03-31)

* Second Zenodo release (DOI: 10.5281/zenodo.6400119), accompanying the
  revised manuscript accepted in *Medical Decision Making*.
* Code, manuscript text and references were updated in response to peer
  review. Rewards are scaled by `cycle_length` throughout.
* Added `rate_to_prob()`, `gen_wcc()` and the transition-matrix checks from
  `darthtools` to `R/Functions.R`, so the tutorial runs without `darthtools`.
* Added the matrix-algebra appendix.

# cohort-modeling-tutorial-intro 0.1.1 (2021-08-19)

* First Zenodo release (DOI: 10.5281/zenodo.5223093).

# cohort-modeling-tutorial-intro 0.1.0 (2021-08-19)

* First public release, accompanying the submitted manuscript.
