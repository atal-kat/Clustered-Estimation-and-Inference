#' Cluster-jackknife instrumental variables estimation (CJIVE)
#'
#' Computes the cluster-jackknife IV estimator of Frandsen, Leslie and McIntyre
#' (2025) for one endogenous regressor and one or more excluded instrument
#' columns, with
#' cluster-robust inference.  The first-stage value for each observation is
#' fitted from a regression that leaves out the observation's entire cluster,
#' which removes the many-instrument bias that survives clustering.
#'
#' @param y Outcome (numeric vector).  Factor and character outcomes are
#'   rejected.  The generic dispatches on this argument; a formula as the
#'   first argument selects the formula method (see \code{formula}).
#' @param x Single endogenous regressor (numeric vector).  Factor and
#'   character regressors are rejected; for the print methods, a fitted
#'   \code{"cjive"} object.
#' @param z Instruments: a numeric vector/matrix, or a factor/character grouping
#'   vector for a judge design.  A grouping instrument uses reference coding
#'   whenever an intercept or any fixed effect is partialled out.  It uses one
#'   column for every level only when \code{intercept = FALSE} and no fixed
#'   effects are supplied.
#' @param cluster Cluster identifiers (length n).  In the formula methods all
#'   four forms work: a bare column name (\code{cluster = cl}), a one-sided
#'   formula (\code{cluster = ~cl}), a column name as a string
#'   (\code{cluster = "cl"}), or a full vector (\code{cluster = df$cl});
#'   a bare name is looked up in \code{data} first, then in the caller.  A
#'   cluster formula must contain exactly one bare variable; multiway cluster
#'   formulas and arithmetic combinations of cluster codes are not supported.
#' @param formula A model formula in either of two restricted layouts.  The
#'   legacy layout is \code{y ~ x | z} (endogenous regressor, then instruments),
#'   optionally \code{y ~ x | z | fe} with fixed effects after a second bar;
#'   controls enter via the \code{controls} argument.  The controls-inside-
#'   formula layout is \code{y ~ exog | endo ~ inst} or
#'   \code{y ~ exog | fe | endo ~ inst}: exogenous controls first, then
#'   (optionally) fixed effects, then the IV part \code{endo ~ inst}, which
#'   must come last.  This is not a full model-formula grammar: every section
#'   accepts only bare variable names joined by \code{+}; transformed terms,
#'   arithmetic and subtraction are rejected, apart from \code{0} or
#'   \code{-1} as intercept-removal markers in the endogenous/exogenous
#'   section.  There must be exactly one outcome and one endogenous regressor.
#'   Removing the intercept also switches a grouping instrument to full-level
#'   coding unless fixed effects are absorbed (their span already contains the
#'   intercept).  Supplying controls or fixed effects both in the formula and via
#'   \code{controls =} / \code{fixed_effects =} is an error.
#' @param data A data frame in which to evaluate the formula.
#' @param controls Optional covariates (FLM's \eqn{X}): a matrix or data frame,
#'   or a one-sided formula in the formula methods.  May be rank deficient (fixed
#'   effects are allowed).  An intercept is added unless \code{intercept = FALSE}.
#' @param fixed_effects Optional high-dimensional fixed effects to absorb: a
#'   factor, or a list/data frame of factors (one per dimension).  In the
#'   formula methods they can also enter through either formula layout
#'   or as a one-sided formula (\code{fixed_effects = ~ fe1 + fe2}).  Absorbed
#'   by a matrix-free, \code{rowsum}-based joint projection before the dense
#'   partialling out; numerically identical (to tolerance) to entering the
#'   corresponding dummies via \code{controls}, but without forming any dummy
#'   matrix.  Absorption spans the intercept within groups, and a result is
#'   returned only after normalized refinement and per-dimension weighted
#'   orthogonality checks pass.
#' @param weights Optional numeric precision weights.  Values must be finite
#'   and strictly positive; factor and character weights are rejected.  In
#'   the formula methods the same four forms as \code{cluster} work, including
#'   a bare data-frame column name.
#' @param subset Optional expression selecting the rows to use, evaluated in
#'   \code{data} (formula methods only), as in \code{\link[stats]{lm}}.
#' @param na.action How to treat missing values (formula methods only):
#'   \code{\link[stats]{na.omit}} (the default) drops every row with a missing
#'   value in any model component -- outcome, endogenous regressor,
#'   instruments, controls, fixed effects, cluster or weights -- in a single
#'   complete-cases filter,
#'   and records the count as \code{n_dropped} (reported by \code{print} when
#'   non-zero); \code{na.fail} stops on any missing value; \code{na.pass}
#'   leaves them in, so the input validation reports them.  The default
#'   methods (vector interface) do not handle missing values at all: they
#'   require clean input and stop otherwise.
#' @param level Confidence level for the reported interval: one finite numeric
#'   value strictly between 0 and 1.
#' @param intercept One non-missing logical value; partial out an intercept
#'   (default \code{TRUE}).
#' @param inference Critical values and p-values: \code{"asymptotic"} (standard
#'   normal, the default) or \code{"t"} (Student t with \eqn{G-1} degrees of
#'   freedom, a small-cluster reference convention in common use; it is not
#'   an exact finite-sample t law under arbitrary within-cluster
#'   dependence).  The coefficient and
#'   standard error are identical under both; only the interval and p-value
#'   change.
#' @param method One of \code{"auto"}, \code{"dense"}, \code{"leaveout_mean"}.
#'   \code{"auto"} and \code{"dense"} both use the dense Frisch-Waugh-Lovell
#'   block-jackknife and are the default.  \code{"leaveout_mean"} evaluates FLM's
#'   printed leave-cluster-out group-mean form and is available only for a
#'   grouping-factor \code{z} with intercept-only controls; it differs from the
#'   default by an intercept term of order \eqn{n_g/n} (about \eqn{1/G} in
#'   balanced designs) and is never selected automatically.
#' @param ... Must be empty.  Unknown, misspelled and unnamed arguments are
#'   errors.
#'
#' @return An object of class \code{"cjive"}: a list with \code{coefficient},
#'   \code{se}, \code{statistic}, \code{p.value}, \code{conf.low},
#'   \code{conf.high}, \code{term} (the endogenous-regressor label),
#'   \code{level}, \code{inference}, the diagnostics \code{n}, \code{G},
#'   \code{k} (number of instrument columns after expansion) and \code{p}
#'   (the retained v0.1.0 compatibility alias for \code{k}), \code{path}
#'   (\code{"dense"} or
#'   \code{"leaveout_mean"}), \code{maxlev} (the maximum within-cluster
#'   leverage \eqn{\max_g \lambda_{\max}(H_g)}, a conditioning diagnostic;
#'   \code{NA} on the mean path), \code{F_eff}, \code{K_eff} and
#'   \code{F_eff_crit} (the clustered Montiel Olea-Pflueger effective
#'   first-stage F, its possibly fractional effective degrees of freedom,
#'   and the simplified-TSLS critical value at \eqn{\tau = 10\%},
#'   \eqn{\alpha = 5\%}; if the estimated denominator is exactly zero with
#'   nonzero signal, \code{F_eff = Inf} while \code{K_eff} and
#'   \code{F_eff_crit} are \code{NA}; all three are \code{NA} on the
#'   \code{"leaveout_mean"} path; the details and the
#'   \sQuote{Where is the first-stage F?} FAQ live in \code{\link{cjar}}),
#'   \code{ng_max} (largest cluster size), \code{fe_dims} and
#'   \code{fe_levels} (number
#'   of absorbed fixed-effect dimensions and their total level count; 0 when
#'   none), \code{k_controls} (a raw nuisance-column upper bound: dense control
#'   columns plus the explicit intercept when there are no fixed effects, or
#'   one intercept plus \eqn{L_j-1} columns per absorbed dimension; not the
#'   effective rank-adjusted nuisance dimension, because nesting, duplicate
#'   dimensions and dense-control collinearity can reduce rank),
#'   \code{n_dropped} (rows removed by
#'   \code{na.action}; 0 on the vector interface), and the \code{call}.
#'
#' @details
#' The estimator is the covariance ratio
#' \eqn{\hat\delta = \widehat{Cov}(Y,\hat p)/\widehat{Cov}(D,\hat p)} with the
#' cluster-jackknife constructed instrument \eqn{\hat p}.  Covariates are handled
#' by Frisch-Waugh-Lovell: \eqn{Y}, \eqn{D} and each instrument are residualised
#' on the covariates (with an intercept) once, up front, then the estimator runs
#' on the residuals.  This dense route is the single convention everywhere, so
#' \code{cjive()} and \code{\link{iv_compare}} return the identical CJIVE for any
#' design.  The weighted, residualised inputs must be finite.  After dense FWL
#' or fixed-effect absorption, \code{y}, \code{x} and every instrument column
#' must retain usable variation, and the residualised instrument matrix must
#' have full numerical column rank.  These checks are scale- and
#' dimension-aware; absorbed or collinear columns are reported rather than
#' silently dropped.  Point estimation also stops before division when its IV
#' denominator is numerically zero.  The leave-cluster-out fits are computed
#' in whitened coordinates
#' (\eqn{\tilde Z = Z R^{-1}}, with \eqn{R} the Cholesky factor of
#' \eqn{Z'Z}, obtained by one triangular solve): per cluster, a
#' \eqn{k \times k} solve when \eqn{n_g > k}, while the Woodbury block update
#' on the \eqn{n_g \times n_g} projection block survives as the small-cluster
#' branch (\eqn{n_g \le k}, including singleton clusters) -- it is not
#' discarded.  The total cost is \eqn{O(nk^2 + \sum_g \min(n_g, k)^3)} time
#' and one \eqn{n \times k} whitened copy of the instruments; the
#' \code{maxlev} diagnostic is a free by-product of the same per-cluster
#' pass, since the \eqn{k \times k} and \eqn{n_g \times n_g} blocks share
#' their non-zero eigenvalues.  The fits agree numerically with the dense
#' leave-cluster-out definition, collapsing when every cluster is a singleton
#' to
#' the observation-level improved JIVE (IJIVE) of Ackerberg and Devereux
#' (2009) -- the leave-one-out first-stage fit on the covariate-partialled
#' instruments, not the original JIVE of Angrist, Imbens and Krueger (1999),
#' which would jackknife the covariates alongside the instruments.
#' The cluster jackknife originates in the 2023 first version of Ligtenberg
#' (2025); Frandsen, Leslie and McIntyre (2025) use it to construct CJIVE.
#' For a pure judge design, Ligtenberg and Woutersen (2024, Section 2 and
#' Appendix B) show the connection, up to weighting, between CJIVE on judge
#' dummies and 2SLS with a leave-cluster-out mean instrument.  The explicit
#' \code{method = "leaveout_mean"} path exposes that group-mean form; the
#' default remains the package-wide dense FWL convention described above.
#'
#' \strong{Formula layouts.}  Both supported layouts describe the same
#' designs.  The legacy layout separates the parts with bars only
#' (\code{y ~ x | z | fe}, controls via \code{controls =}); the second layout
#' (\code{y ~ exog | fe | endo ~ inst}) carries controls inside the formula.
#' Both use the restricted bare-name, additive syntax described under
#' \code{formula}; they do not implement a general R or \code{feols} formula
#' language.  Missing values are handled only on the formula path (see
#' \code{na.action}); the vector interface requires clean input.
#'
#' Fixed effects passed via \code{fixed_effects} are absorbed by a matrix-free
#' joint projection.  Each factor projector uses C-level \code{rowsum} group
#' sums.  With multiple dimensions, a conjugate-gradient solve applies
#' self-adjoint combinations of those projectors.  The result is returned only
#' after freshly normalized refinement finds no material projection correction,
#' every dimension's weighted group score and correction are below tolerance,
#' and an additional complete sweep is stable.  A single dimension is one
#' exact projection.  If the weighted FE incidence geometry or an observed
#' Krylov direction is too ill-conditioned to support the requested numerical
#' tolerance, absorption stops with a conditioning diagnosis rather than
#' returning the backward certificate as a forward-accuracy claim.  The result
#' is numerically identical, to the stated solver tolerance, to entering the
#' corresponding dummy variables through \code{controls}, but no dummy matrix
#' is ever formed.
#'
#' A dense control whose post-absorption norm is at most the fixed-effect
#' numerical tolerance relative to its original weighted norm is treated as
#' absorbed before the remaining FWL regression.  Surviving control columns
#' are normalized before pivoted QR, so numerical-rank decisions do not depend
#' on their units.  Joint collinearity among surviving controls is allowed.
#'
#' \strong{Many controls and the reported standard errors.}  Kolesar, Min, Wang
#' and Zhang (2026) show that CJIVE with the plug-in cluster-robust sandwich --
#' exactly what this package reports -- can over-reject when the instrument
#' and control dimensions are large relative to \eqn{n}.  In one of their
#' Table 1 designs with \eqn{n = 600}, rejection at a nominal 5 percent is 7.9
#' percent with 50 instruments and 50 controls and 53.4 percent with 150 of
#' each; they also report over-rejection under heterogeneous treatment effects.
#' The mechanism is the global partialling-out step, which
#' uses each observation's own cluster, so many controls reintroduce a bias the
#' first-stage cluster jackknife does not touch.  Since \code{fixed_effects}
#' makes it easy to absorb thousands of levels -- that is, to enter exactly
#' this regime -- the ratio \code{k_controls/n} is reported by
#' \code{print.cjive} as a transparency diagnostic.  No variance correction is
#' applied.
#'
#' @references
#' Ackerberg, D. A. and Devereux, P. J. (2009). Improved JIVE estimators for
#' overidentified linear models with and without heteroskedasticity.
#' \emph{Review of Economics and Statistics}, 91(2), 351--362.  The improved
#' JIVE that the singleton-cluster limit of CJIVE reduces to.
#'
#' Angrist, J. D., Imbens, G. W. and Krueger, A. B. (1999). Jackknife
#' instrumental variables estimation. \emph{Journal of Applied Econometrics},
#' 14(1), 57--67.  The original JIVE.
#'
#' Frandsen, B., Leslie, E. and McIntyre, S. (2025). Cluster Jackknife
#' Instrumental Variables Estimation. \emph{Review of Economics and Statistics}.
#'
#' Ligtenberg, J. W. (2025). Inference in clustered IV models with many and
#' weak instruments. arXiv:2306.08559v3.  The 2023 first version introduced
#' the cluster jackknife used by CJIVE.
#'
#' Ligtenberg, J. W. and Woutersen, T. (2024). Multidimensional clustering in
#' judge designs. arXiv:2406.09473.  Section 2 and Appendix B establish the
#' judge-dummy CJIVE/leave-out-mean connection up to weighting.
#'
#' Gaure, S. (2013). OLS with multiple high dimensional category variables.
#' \emph{Computational Statistics and Data Analysis}, 66, 8--18.
#'
#' Guimaraes, P. and Portugal, P. (2010). A simple feasible procedure to fit
#' models with high-dimensional fixed effects. \emph{Stata Journal}, 10(4),
#' 628--649.
#'
#' Halperin, I. (1962). The product of projection operators. \emph{Acta
#' Scientiarum Mathematicarum (Szeged)}, 23, 96--99.
#'
#' Kolesar, M., Min, P., Wang, W. and Zhang, Y. (2026). Cluster-Robust
#' Inference for Quadratic Forms. arXiv:2602.13537.
#'
#' @seealso \code{\link{iv_infer}} for the recommended one-call workflow
#'   (CJIVE estimate + CJAR confidence set + CJS point-null test with
#'   shared preprocessing); \code{\link{cjar}} and \code{\link{cjscore}} for the standalone
#'   weak-instrument-robust tests; \code{\link{iv_compare}} for estimator
#'   comparison.  Two precomputed worked examples on real data:
#'   \code{vignette("queens-workflow", package = "clusterIV")} and
#'   \code{vignette("miami-bail", package = "clusterIV")}.
#'
#' @examples
#' set.seed(1)
#' G  <- 40; ng <- 6; n <- G * ng
#' cl <- rep(seq_len(G), each = ng)
#' j  <- factor(rep(rep(1:4, length.out = ng), G))   # judge identity
#' u  <- rnorm(G)[cl]
#' x  <- as.numeric(j) + u + rnorm(n)
#' y  <- 1.5 * x + u + rnorm(n)
#' w1 <- rnorm(n); w2 <- rnorm(n)
#' fit <- cjive(y, x, j, cluster = cl)
#' print(fit)
#'
#' ## formula interface with controls supplied through the argument
#' dat <- data.frame(y = y, x = x, j = j, cl = cl, w1 = w1, w2 = w2)
#' cjive(y ~ x | j, data = dat, cluster = cl, controls = ~ w1 + w2)
#'
#' ## formula interface, controls inside the formula
#' cjive(y ~ w1 + w2 | x ~ j, data = dat, cluster = cl)
#'
#' @export
cjive <- function(y, ...) {
  .check_dispatch_dots(
    y,
    match.call(expand.dots = FALSE)$...,
    c("x", "z", "cluster", "controls", "fixed_effects", "weights",
      "intercept", "level", "method", "inference", "data", "subset",
      "na.action"),
    "cjive", parent.frame()
  )
  UseMethod("cjive")
}

# Assemble a "cjive" object from the prepared design and the inference list;
# shared by cjive.default and iv_infer() so the two produce identical fields.
# `eff` is the effective-F triple from .eff_f() (or .eff_f_na() on the
# leaveout_mean path, which has no whitened first stage).
.cjive_build <- function(d, inf, level, inference, path, maxlev, eff, call,
                         term = "x") {
  structure(c(inf[c("coefficient", "se", "statistic", "p.value",
                    "conf.low", "conf.high")],
              list(term = term, level = level, inference = inference,
                   k = d$k, p = d$k, path = path, call = call),
              .fit_common(d, eff, maxlev)),
            class = "cjive")
}

#' @rdname cjive
#' @export
cjive.default <- function(y, x, z, cluster, controls = NULL, weights = NULL,
                          level = 0.95, intercept = TRUE,
                          method = c("auto", "dense", "leaveout_mean"),
                          fixed_effects = NULL,
                          inference = c("asymptotic", "t"), ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "cjive")
  cl <- match.call()
  .check_flag(intercept, "intercept")
  .check_level(level)
  method <- match.arg(method)
  inference <- match.arg(inference)
  term <- .term_label(substitute(x))

  d <- .prep_data(y, x, z, cluster, controls, weights, intercept,
                  fixed_effects = fixed_effects)

  if (method == "leaveout_mean") {
    if (!d$grouping)
      stop("method = \"leaveout_mean\" requires a grouping-factor `z`.", call. = FALSE)
    if (!is.null(controls) || !is.null(fixed_effects))
      stop("method = \"leaveout_mean\" requires intercept-only controls.", call. = FALSE)

    phat <- .leaveout_mean(x, d$group, d$cluster, d$weights)
    # Centre y, x and p-hat (residualise on the intercept) to match the dense route.
    po <- .partial_out(as.numeric(y), as.numeric(x), phat,
                       controls = NULL, weights = d$weights,
                       intercept = intercept)
    inf <- .iv_inference(as.numeric(po$Z), po$x, po$y, d$cluster, level,
                         inference = inference)
    path <- "leaveout_mean"
    maxlev <- NA_real_
    eff <- .eff_f_na()
  } else {
    fs <- .first_stage(d$x, d$Z, R = d$R)
    lo <- .leaveout_fit(d$x, fs$Ztil, fs$t, d$groups, ncol(d$Z))
    inf <- .iv_inference(lo$phat, d$x, d$y, d$cluster, level,
                         inference = inference)
    path <- "dense"
    maxlev <- lo$maxlev
    eff <- .eff_f(fs$Ztil, fs$t, fs$e, d$cluster)
  }

  .cjive_build(d, inf, level, inference, path, maxlev, eff, cl, term)
}

#' @rdname cjive
#' @export
cjive.formula <- function(formula, data, cluster, controls = NULL,
                          weights = NULL, level = 0.95, intercept = TRUE,
                          method = c("auto", "dense", "leaveout_mean"),
                          fixed_effects = NULL,
                          inference = c("asymptotic", "t"),
                          subset, na.action = stats::na.omit, ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "cjive")
  .check_flag(intercept, "intercept")
  .check_level(level)
  .formula_method(formula, data, missing(data), substitute(cluster),
                  substitute(weights),
                  if (missing(subset)) NULL else substitute(subset),
                  controls, fixed_effects, na.action, intercept,
                  parent.frame(), match.call(),
                  fit = function(prep, intercept)
                    cjive.default(prep$y, prep$x, prep$z,
                                  cluster = prep$cluster,
                                  controls = prep$controls,
                                  fixed_effects = prep$fixed_effects,
                                  weights = prep$weights,
                                  intercept = intercept, level = level,
                                  method = method, inference = inference))
}
