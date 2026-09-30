# clusterIV 0.2.0

## Weak-instrument-robust inference

* New `cjar()` implements the cluster-jackknife Anderson–Rubin test of
  Ligtenberg (2025), including analytic, grid-free confidence-set inversion,
  shifted/scaled chi-square and normal calibrations, and complete reporting of
  bounded, disjoint, and unbounded sets.
* New `cjscore()` implements Ligtenberg's cluster-jackknife score test and its
  analytically inverted confidence set. The CJAR and CJS procedures support
  both the plain and leave-two-cluster-out cross-fit variance estimators.
* New `iv_infer()` returns the CJIVE estimate, CJAR confidence set, and CJS
  point-null test in one panel while sharing preprocessing and
  leave-cluster-out calculations.
* New base-graphics methods for `cjar`, `cjscore`, and `iv_infer` display
  p-value curves and the reported confidence sets. The sets are computed by
  polynomial inversion, not inferred from the plotting grid.

## Estimation and fixed effects

* `cjive()`, `cjar()`, `cjscore()`, `iv_infer()`, and `iv_compare()` now accept
  high-dimensional fixed effects through `fixed_effects`. A base-R,
  matrix-free joint projection absorbs multiple factor dimensions without
  constructing a dummy matrix and verifies weighted orthogonality before
  estimation.
* Formula methods accept either `y ~ x | z | fe` or the controls-inside-formula
  layout `y ~ exog | fe | x ~ z`. Formula sections deliberately support a
  restricted additive grammar of bare variable names.
* Formula methods now support `subset` and `na.action`. A single complete-case
  filter is applied to every model component. The fitted-object interfaces
  (`cjive()`, `cjar()`, `cjscore()`, `iv_infer()`) store and report the number
  of omitted rows; `iv_compare()` returns a plain data frame and does not.
* `inference = "t"` uses a `t(G - 1)` reference convention for CJIVE Wald
  inference. It does not change the coefficient or standard error and is not
  presented as an exact finite-sample law.
* The dense Frisch–Waugh–Lovell path remains the default throughout the
  package. `method = "leaveout_mean"` remains an explicit special case for a
  pure grouping-instrument design and is never selected automatically.

## Diagnostics and output

* Dense-path fits from `cjive()`, `cjar()`, `cjscore()`, and `iv_infer()`
  report the clustered Montiel Olea–Pflueger effective first-stage statistic
  (`F_eff`) and its simplified-TSLS critical value. CJAR and CJS also report
  their polynomial tail diagnostics (`F_CJ` and `F_CJS`) alongside the
  existing maximum within-cluster leverage diagnostic. `iv_compare()`
  intentionally reports estimates only, without these diagnostics.
* New `summary()` methods for CJAR and CJS report the test, confidence-set
  topology, first-stage diagnostics, and design advisories together.
  `summary.cjive()` includes the same strength block.
* `nobs()` methods are available for fitted objects. Supported `tidy()` and
  `glance()` methods are registered through `generics` when that suggested
  package is installed; disjoint and unbounded sets remain available in a
  list column rather than being forced into a single interval.
* `iv_compare()` is now an S3 generic with vector and formula methods. Its
  historical `"JIVE"` row label is retained for compatibility; the
  observation-level calculation after FWL is IJIVE.

## Documentation

* Two precomputed vignettes work through the package on real, named
  datasets: `vignette("queens-workflow")` (Dube and Harish 2020; the full
  workflow, the published specification set, and the LaTeX export path) and
  `vignette("miami-bail")` (Frandsen, Leslie and McIntyre 2025; the
  judge-leniency design at scale). The data cannot be redistributed, so the
  vignettes ship with their outputs baked in and state where to obtain the
  files, verified by hash. `knitr` and `rmarkdown` enter `Suggests` for
  vignette building only.
* The README gained a Validation section summarising the external Stata
  cross-checks and the frozen dense oracles, and five FAQ entries (which function to call first,
  CJAR versus CJS, interpreting `F_CJ`/`F_CJS`, robust sets versus Wald
  intervals, and LaTeX export).

## Performance

* The cross-fit variance kernel dispatches each leave-two-cluster-out solve
  on the stacked size of the two left-out clusters, mirroring the CJIVE
  leave-out kernel: many-instrument, small-cluster pairs are solved on the
  small side through a Woodbury identity, and the leave-out Gram matrices
  are cached lazily instead of always held. Together this reduces the
  per-pair time and the memory footprint by orders of magnitude in
  many-instrument, small-cluster designs; coefficients still pass the dense
  eq. (7) oracle at the 1e-10 gate, the exact spectral guard and its
  fail-closed classification are unchanged, and the large-side path is
  unchanged bit for bit.
* `cjar()` and `cjscore()` with `variance = "crossfit"` now compute only the
  variance polynomial they consume. `iv_infer()` still shares one pass when
  both jackknife tests are requested.
* Whitening and the first-stage hat-diagonal accumulation allocate far less
  memory, which matters at high-dimensional-fixed-effects scale.

## Input and numerical safety

* Input validation is shared across the public functions. Outcomes,
  endogenous regressors, and weights must be numeric; weights must be finite
  and positive; absorbed variables, rank-deficient instrument designs, and
  numerically zero IV denominators produce explicit errors.
* Factor and character grouping instruments use reference coding whenever an
  intercept or fixed effects are partialled out, and full-level coding only
  when neither is present. Data-frame factor controls follow the same rule,
  avoiding a redundant dummy when fixed effects already span the intercept.
  Dense paths assess support from the encoded, residualized instrument matrix;
  the stricter outside-cluster support rule is applied only to the explicit
  `leaveout_mean` shortcut where that rule is mathematically required.
* Dense controls that are numerically absorbed by the joint fixed-effect span
  are removed before the second FWL step, rather than allowing pivoted QR to
  promote projection noise into a regressor. Surviving controls are normalized
  before a dimension-aware rank decision, making the result invariant to their
  units. The reported `k_controls` is now the non-redundant reference-coded
  nuisance-column upper bound rather than intercept plus every FE level.
* Numerically scaled inference preserves coefficient, standard-error, and
  test-statistic invariance under common rescaling and common rescaling of
  precision weights. Very large finite null values in CJAR/CJS and their plot
  methods are evaluated without forming overflowing powers of the null.
* Cross-fit variance estimation now rejects a leave-one- or leave-two-cluster-
  out instrument Gram whose spectral reciprocal condition number is at or
  below `sqrt(.Machine$double.eps)`, instead of returning a result that can
  depend materially on the instrument basis.
* CJAR inversion now handles two-cluster designs and variance-boundary
  tangencies consistently with pointwise test acceptance. JIVE leverage equal
  to one and unsupported plain-test covariance ranks stop or return a typed
  unavailable result rather than propagating non-finite output.

## Important changes

* `confint()` now defaults to the confidence level stored on the fitted object
  instead of silently reverting to 95%.
* Fast fixed effects make high-control designs practical, but the reported
  plug-in CJIVE sandwich can over-reject when the number of instruments and
  controls is large relative to the sample. The help pages now state this
  limitation and cite Kolesar, Min, Wang and Zhang (2026).
* The package title and description now cover both clustered IV estimation and
  weak-instrument-robust inference. The CRAN housekeeping corrections noted
  after version 0.1.0 are included.

# clusterIV 0.1.0

* Initial CRAN release with `cjive()` and `iv_compare()`.
