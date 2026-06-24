# cjive

**Cluster-jackknife instrumental variables estimation for R.**

`cjive` computes the cluster-jackknife IV estimator (CJIVE) of Frandsen, Leslie
and McIntyre (2025) for a single endogenous regressor in a just-identified
design — the judge/examiner and shift-share settings where the instrument is
many and the errors are clustered. Each observation's first-stage value is
fitted from a regression that leaves out the observation's *entire cluster*,
which annihilates the within-cluster dependence that otherwise reintroduces the
many-instrument bias of two-stage least squares. The leave-cluster-out fits are
computed by an exact Woodbury block update — one Cholesky of the instrument Gram
matrix plus a small solve per cluster, never `G` refactorisations — so the
estimator runs comfortably on samples in the hundreds of thousands. Base R only
(`Imports: stats`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("atal-kat/clustered-estimation-and-inference")
```

## Usage

Formula interface, `y ~ x | z` (the bar separates the endogenous regressor from
the instruments); a factor on the instrument side is a judge design:

```r
library(cjive)
fit <- cjive(wage ~ incarcerated | judge_id, data = cases, cluster = ~courtroom)
summary(fit)
```

Matrix/vector interface, with a numeric instrument matrix:

```r
cjive(y, x, z, cluster = g, controls = X, weights = w)
```

`controls` are partialled out by Frisch-Waugh-Lovell (fixed effects allowed),
`weights` are optional precision weights, and `confint()`, `coef()`, `vcov()`
work as usual.

### Comparing estimators

`iv_compare()` reproduces the shape of FLM's Table 1: OLS, 2SLS, JIVE and CJIVE
on the *same* cluster-robust IV sandwich standard error, only the constructed
instrument differing between rows.

```r
iv_compare(y, x, z, cluster = g)
#>   estimator coefficient     se statistic  p.value conf.low conf.high
#> 1       OLS         ...    ...       ...      ...      ...       ...
#> 2      2SLS         ...    ...       ...      ...      ...       ...
#> 3      JIVE         ...    ...       ...      ...      ...       ...
#> 4     CJIVE         ...    ...       ...      ...      ...       ...
```

## The covariate/intercept convention

There is **one** rule: the dense Frisch-Waugh-Lovell route is the default
everywhere. Covariates (and the intercept) are partialled out globally, then the
leave-cluster-out fit runs on the residuals. A grouping-factor `z` is expanded
to a dummy design (one reference level dropped, the intercept supplying the
rest) and run through the same path, so `cjive()` and `iv_compare()` return the
identical CJIVE for any design.

FLM's printed closed form for the *pure* judge design — the leave-cluster-out
group mean of `x` — is available on explicit request via
`method = "leaveout_mean"` (grouping-factor `z`, intercept-only controls). It
differs from the default by an `O(1/n_g)` intercept term and is never selected
automatically.

## Reference

Frandsen, B., Leslie, E. & McIntyre, S. (2025). *Cluster Jackknife Instrumental
Variables Estimation.* Review of Economics and Statistics.
