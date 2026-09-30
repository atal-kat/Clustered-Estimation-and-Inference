# clusterIV

Clustered instrumental variables estimation and weak-instrument-robust
inference for R.

`clusterIV` estimates CJIVE and implements cluster-jackknife Anderson–Rubin
and score inference for one endogenous regressor under one-way clustering.
CJIVE removes own-cluster first-stage terms that remain in 2SLS and
observation-level jackknife IV when observations are dependent within
clusters.

[Setup](#setup) | [Syntax](#syntax) | [Description](#description) |
[Options](#options) | [Simulated example](#simulated-example) |
[Worked examples](#worked-examples) | [Validation](#validation) |
[FAQ](#faq) | [References](#references)

## Setup

Install the released version from CRAN:

```r
install.packages("clusterIV")
```

The weak-instrument-robust inference functions documented below
(`cjar()`, `cjscore()`, `iv_infer()`), fixed-effect absorption, and the
formula interfaces require clusterIV 0.2.0 or later; version 0.1.0
contains `cjive()` and `iv_compare()` only. Check the installed version
with `packageVersion("clusterIV")`.

## Syntax

The vector interfaces expose all computational options directly:

```r
cjive(y, x, z, cluster,
      controls = NULL, weights = NULL, level = 0.95, intercept = TRUE,
      method = c("auto", "dense", "leaveout_mean"),
      fixed_effects = NULL, inference = c("asymptotic", "t"))

cjar(y, x, z, cluster,
     controls = NULL, fixed_effects = NULL, weights = NULL,
     intercept = TRUE, level = 0.95, beta0 = 0,
     calibration = c("chisq", "normal"),
     variance = c("plain", "crossfit"))

cjscore(y, x, z, cluster,
        controls = NULL, fixed_effects = NULL, weights = NULL,
        intercept = TRUE, level = 0.95, beta0 = 0,
        variance = c("plain", "crossfit"))

iv_infer(y, x, z, cluster,
         controls = NULL, fixed_effects = NULL, weights = NULL,
         intercept = TRUE, level = 0.95, beta0 = 0,
         calibration = c("chisq", "normal"),
         inference = c("asymptotic", "t"),
         tests = c("cjar", "cjscore"),
         variance = c("plain", "crossfit"))

iv_compare(y, x, z, cluster,
           controls = NULL, weights = NULL, level = 0.95,
           intercept = TRUE, fixed_effects = NULL,
           inference = c("asymptotic", "t"))
```

The same functions accept either restricted formula layout:

```r
y ~ x | z                  # controls supplied through controls =
y ~ x | z | fe1 + fe2      # absorbed fixed effects
y ~ exog1 + exog2 | x ~ z  # controls inside the formula
y ~ exog | fe1 + fe2 | x ~ z
```

Formula methods additionally accept `data`, `subset`, and
`na.action = stats::na.omit`. Formula sections are additive lists of bare
variable names. Transformed terms, multiway clustering, and multiple
endogenous regressors are not supported.

## Description

The cluster jackknife originates with Ligtenberg (2023, the first version of
the cited preprint). Frandsen, Leslie and McIntyre (2025) use it to construct
CJIVE. On the default dense path, the package residualises the outcome,
endogenous regressor, and instruments by Frisch–Waugh–Lovell before applying
the cluster jackknife. Consequently, `cjive()` and the CJIVE row of
`iv_compare()` agree for the same design. The explicitly requested
`method = "leaveout_mean"` shortcut is the documented exception.

`iv_infer()` is the main inference workflow. Its default panel contains:

- a CJIVE coefficient and cluster-robust Wald interval;
- a CJAR confidence set, valid asymptotically with weak and many instruments
  under the conditions in Ligtenberg (2025); and
- a CJS test and confidence set for the specified point null `beta0`.

The selected procedures share preprocessing, the instrument
factorisation, leave-cluster-out work, and cluster sums.

The leave-cluster-out first stage works in whitened instrument coordinates and
chooses the smaller of the cluster-size and instrument-dimension solves. It
uses triangular solves rather than an explicit inverse. High-dimensional
fixed effects are absorbed by a matrix-free joint projection and checked for
weighted orthogonality before estimation; no dummy matrix is formed.

## Options

### Core inputs

- `y` and `x` are numeric vectors. Exactly one endogenous regressor is
  supported.
- `z` is a numeric vector or matrix of excluded instruments, or a factor or
  character grouping instrument for a judge/examiner design.
- `cluster` identifies one clustering dimension. A formula must name one bare
  variable.
- `weights` supplies finite, strictly positive precision weights.

### Covariates and fixed effects

- `controls` supplies ordinary covariates. An intercept is included unless
  `intercept = FALSE`.
- `fixed_effects` supplies one or more factors to absorb without dummy
  expansion. In formulas, fixed effects can instead occupy the dedicated
  formula section.
- `method = "auto"` and `method = "dense"` use the package-wide dense FWL
  convention. `method = "leaveout_mean"` is an explicit special case for a
  pure grouping-instrument design; it is never selected automatically.

### Inference and diagnostics

- `level` sets the reported confidence level; `beta0` sets the null tested by
  CJAR and CJS.
- `variance = "plain"` uses Ligtenberg's plain variance estimator.
  `variance = "crossfit"` uses the leave-cluster-out cross-fit construction
  in Section 5.4 of Ligtenberg (2025). The cross-fit estimator solves one
  leave-out system per cluster pair — `G(G-1)/2` solves, inherent to its
  definition — with each solve dispatched on the smaller of the stacked
  cluster size and the instrument dimension.
- `calibration = "chisq"` uses the paper's shifted/scaled chi-square
  approximation for CJAR; `"normal"` uses its asymptotic normal reference.
  Neither is an exact finite-sample law.
- `inference = "t"` replaces the CJIVE normal critical value with a
  `t(G - 1)` reference convention. It does not change the coefficient or
  standard error and is not an exact finite-sample law under arbitrary
  within-cluster dependence.
- `tests` selects any subset of `"cjar"` and `"cjscore"`. Use `NULL` for a
  CJIVE-only panel.

Dense-path fits from `cjive()`, `cjar()`, `cjscore()`, and `iv_infer()`
report the maximum within-cluster leverage `maxlev` and the clustered
Montiel Olea–Pflueger effective first-stage statistic `F_eff`. CJAR and CJS
additionally report their polynomial tail diagnostics `F_CJ` and `F_CJS`.
`iv_compare()` reports estimates only, without these diagnostics. The
diagnostics describe different aspects of the design and are not automatic
model-selection rules.

## Simulated example

```r
library(clusterIV)

set.seed(42)
G <- 40L
ng <- 8L
n <- G * ng

dat <- data.frame(
  cluster = rep(seq_len(G), each = ng),             # dependence cluster
  judge = factor(sample(1:8, n, replace = TRUE)),   # grouping instrument
  w = rnorm(n)                                      # included covariate
)
cluster_shock <- rnorm(G)[dat$cluster]
dat$x <- 0.45 * as.numeric(dat$judge) +             # endogenous regressor
  0.30 * dat$w + cluster_shock + rnorm(n)
dat$y <- 1.20 * dat$x + 0.50 * dat$w +              # outcome
  cluster_shock + rnorm(n)

fit <- iv_infer(
  y ~ w | x ~ judge,      # outcome/control | endogenous ~ instrument
  data = dat,             # analysis data
  cluster = cluster,      # one-way clustering
  beta0 = 1,              # point null for CJAR and CJS
  level = 0.95,           # confidence level
  variance = "plain"      # Ligtenberg's plain variance estimator
)

fit
coef(fit$cjive)            # CJIVE point estimate
confint(fit$cjar)          # weak-ID-robust CJAR confidence set
fit$cjscore$p.value        # CJS p-value for beta = 1
plot(fit)                  # p-value curves and reported sets
```

For an estimator comparison using the same residualised data and
cluster-robust sandwich convention:

```r
iv_compare(y ~ w | x ~ judge, data = dat, cluster = cluster)
```

This returns OLS, 2SLS, observation-level IJIVE (retaining the historical row
label `"JIVE"`), and CJIVE. It matches the four-estimator layout of FLM's
Table 1; it is not an empirical replication of that table.

## Worked examples

Two precomputed vignettes run the package on real, named datasets:

- `vignette("queens-workflow")` — the full workflow on Dube and Harish
  (2020): 3,586 polity-years, 176 reign clusters, instruments that are not
  strong, dense controls plus three absorbed fixed-effect dimensions. It
  reproduces the published 2SLS anchors, walks through every diagnostic,
  runs the published specification set, and shows the `tidy()`-to-LaTeX
  export path.
- `vignette("miami-bail")` — the judge-leniency design of Frandsen, Leslie
  and McIntyre (2025): 91,421 defendants, 146 judges, over one hundred
  instrument columns, clustered by courtroom shift. It shows the workflow
  at scale, including the instrument-hygiene step every judge-dummy design
  needs.

Neither dataset can be redistributed inside the package, so the vignettes
ship with their outputs baked in; each states where to obtain the data and
verifies it by hash.

## Validation

The implementation is validated against external references and frozen
brute-force oracles rather than against itself. In brief: CJIVE agrees
with FLM's own released Stata implementation run on the deposited
Miami-Dade sample
(pointwise-identical constructed instrument; coefficient to the ado's
single-precision limit); every CJAR/CJS statistic, variance polynomial,
and confidence-set endpoint is gated at 1e-10 against literal dense
transcriptions of the papers' formulas that were written and frozen before
the fast paths existed; and at singleton clusters the tests reduce exactly
to the Mikusheva–Sun and Matsushita–Otsu statistics.

## FAQ

### Which function do I call first?

`iv_infer()`. One call returns the CJIVE estimate, the CJAR confidence
set, and the CJS test from one shared preprocessing pass, and its printed
panel carries every diagnostic discussed below. `cjive()`, `cjar()`, and
`cjscore()` are the same computations as focused standalone calls;
`iv_compare()` adds the OLS/2SLS/IJIVE comparison row layout.

### Why use CJIVE rather than 2SLS or observation-level IJIVE?

With clustered observations and many instruments, leaving out only the focal
observation does not remove first-stage terms contributed by other
observations in its cluster. CJIVE leaves out the entire cluster.

### When should `z` be a factor rather than a numeric matrix?

Use a factor or character vector for a grouping instrument such as judge or
examiner identity. With an intercept or absorbed fixed effects, the package
uses reference coding; it uses one column per level only when neither is
present. Supply a numeric matrix when the excluded instruments are already
constructed columns.

### Which test do I report, CJAR or CJS — and what if they disagree?

Report the CJAR confidence set as the headline weak-instrument-robust
result, and the CJS p-value when a specific point null is the question.
The two tests are both valid under weak and many instruments but have
different power profiles: the score test concentrates power near the null,
the AR test retains power against distant alternatives. Their confidence
sets therefore need not coincide, and neither disagreement between them
nor disagreement with the Wald interval is an error. When the CJAR and
Wald intervals disagree materially, the Wald interval is the one whose
validity is in question — check `F_eff` against its critical value.

### Why does my CJAR set differ from the 2SLS or CJIVE Wald interval?

The Wald interval assumes the estimator is approximately normal, which
requires strong instruments; the CJAR set inverts a test that is valid
without that assumption. With strong instruments the two nearly agree;
with weak or many instruments the robust set is typically wider, shifted,
or unbounded — that is information, not a bug. An unbounded set says the
data cannot rule out arbitrarily large coefficients at this level.

### What do `F_CJ` and `F_CJS` mean, and what is a "good" value?

Each compares the tail behaviour of its test's inversion polynomials with
the critical value: `F_CJ > crit` holds exactly when the CJAR confidence
set is bounded, and likewise for `F_CJS`. They are boundedness criteria
in the spirit of a first-stage F, printed next to their own critical
values — compare each against the printed critical value, not against the
conventional rule of 10, and do not interchange them with `F_eff`: the
effective F and `F_CJ` measure different objects and can diverge as the
instrument count grows.

### What does `maxlev` near one mean?

`maxlev` is the largest spectral norm of a within-cluster projection block.
A value near one means that deleting one cluster leaves a nearly singular
instrument Gram matrix. Treat it as a conditioning warning and inspect the
instrument design and cluster sizes.

### What do an unbounded or empty confidence set mean?

An unbounded set is a valid weak-identification outcome: the selected test
does not exclude arbitrarily large coefficient values. An empty set means no
coefficient value is accepted at the selected level. It is not, by itself,
proof of invalid instruments or a misspecified model; inspect the design,
assumptions, and numerical diagnostics.

### When should I use `method = "leaveout_mean"`?

Only when you deliberately want the printed leave-out group-mean formula for
a pure grouping-instrument design with intercept-only controls. The default
dense FWL path is the package-wide convention. The two forms differ through
the intercept direction by order `n_g / n` (about `1 / G` in a balanced
design); the group-mean form is never chosen automatically.

### What if I absorb many controls or fixed effects?

Fast absorption makes a many-controls design easy to fit, but it does not make
the usual plug-in sandwich reliable in every such design. In one Table 1
simulation of Kolesar, Min, Wang and Zhang (2026), with `n = 600`, CJIVE's
nominal 5% test rejects 7.9% with 50 instruments and 50 controls, and 53.4%
with 150 of each. Their L2CO/L3CO corrections are not implemented here.
Ligtenberg's Section 5.3 result does cover cluster-specific controls, including
cluster fixed effects, because they introduce dependence only within clusters;
it should not be read as a general many-controls result.

### How do I get the results into a LaTeX table?

`tidy()` and `glance()` methods are registered through the optional
`generics` package and return plain data frames — one row per inferential
procedure, with the complete confidence set retained in a `conf.set`
list-column and its topology in explicit `shape`/`n.components` columns.
They feed directly into `knitr::kable(format = "latex")` or
`modelsummary`; the queens vignette ends with a worked export. A disjoint
or unbounded set is never silently flattened into two numbers.

### What happens when every cluster is a singleton?

CJIVE reduces to improved JIVE after the package's FWL residualisation. Plain
CJAR has the Mikusheva–Sun jackknife-AR numerator/statistic relationship, and
plain CJS has the Matsushita–Otsu algebraic identity when there are no supplied
controls or the data are already partialled. The cross-fit option uses
Ligtenberg's leave-two-cluster-out variance and should not be given those
plain-variance labels.

### How has the implementation been checked?

CJIVE is tested against dense leave-cluster-out calculations, a literal R
translation of FLM's Stata/Mata algorithm, and FLM's own released Stata
implementation run on the deposited Miami-Dade data. CJAR, CJS, and
cross-fit variance are tested against direct implementations of their
defining formulas, frozen before the fast paths were written. These checks establish numerical agreement with the stated algorithms; the
cited papers supply the statistical assumptions and asymptotic results.

## References

Ackerberg, D. A. and Devereux, P. J. (2009). Improved JIVE estimators for
overidentified linear models with and without heteroskedasticity. *Review of
Economics and Statistics*, 91(2), 351–362.

Frandsen, B., Leslie, E. and McIntyre, S. (2025). Cluster Jackknife
Instrumental Variables Estimation. *Review of Economics and Statistics*.
doi:10.1162/rest.a.263. See the estimator and empirical comparison.

Kolesar, M., Min, P., Wang, W. and Zhang, Y. (2026). Cluster-Robust Inference
for Quadratic Forms. arXiv:2602.13537. See Table 1 and Sections 3–4.

Ligtenberg, J. W. (2025). Inference in clustered IV models with many and weak
instruments. arXiv:2306.08559v3. See Sections 3–5 for CJAR, CJS, controls, and
cross-fit variance.

Montiel Olea, J. L. and Pflueger, C. (2013). A robust test for weak
instruments. *Journal of Business & Economic Statistics*, 31(3), 358–369.
