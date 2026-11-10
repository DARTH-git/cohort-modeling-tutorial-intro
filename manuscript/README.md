# Manuscript sources (archival)

This folder holds the R Markdown source, LaTeX output, figures and appendices
of the published article:

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P,
> Jalal H. An Introductory Tutorial on Cohort State-Transition Models in R Using
> a Cost-Effectiveness Analysis Example. *Medical Decision Making*.
> 2023;43(1):3-20. <https://doi.org/10.1177/0272989X221103163>

The files are kept **unchanged** so that they match the published article.

## Known erratum: within-cycle correction weights

The "Within-cycle correction" section of `cSTM_Tutorial_Intro.Rmd` (and the
derived `.R`, `.tex` and `.pdf` files) assigns Simpson's 1/3rd weights by
position in the 1-based R vector instead of by time t = 0, 1, ..., n:

* The equation shows **wcc = [1/3, 2/3, 4/3, ..., 1/3]**. The correct
  vector is **wcc = [1/3, 4/3, 2/3, 4/3, ..., 4/3, 1/3]**: 4/3 at odd times
  t = 1, 3, ..., and 2/3 at even times t = 2, 4, ....
* The `within-cycle-vector` code chunk applies `is_even()`/`is_odd()` to
  `v_cycles <- seq(1, n_cycles + 1)`, which are positions, not times. The fix
  is to use `v_cycles <- 0:n_cycles`.
* The tutorial uses 75 cycles. That is an odd number of intervals, so
  Simpson's 1/3rd rule alone does not apply. The corrected code uses Simpson's
  3/8 rule for the last three cycles.

As a result, the QALYs and costs reported in the article are slightly
underestimated: by about 0.32 QALYs per strategy, with ICERs off by about
0.4%. The conclusions are unchanged.

The corrected implementation is `gen_wcc()` in [`../R/Functions.R`](../R/Functions.R),
from release v0.3.0 on. It is used by
[`../analysis/cSTM_time_indep.R`](../analysis/cSTM_time_indep.R). See
[`../NEWS.md`](../NEWS.md) for details and a before/after comparison of
results.
