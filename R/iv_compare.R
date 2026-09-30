#' Compare IV estimators on a common cluster-robust standard error
#'
#' Returns OLS, 2SLS, JIVE and CJIVE for the same design, each reported with
#' the same cluster-robust IV sandwich standard-error convention, computed
#' from that row's single constructed instrument;
#' only the constructed instrument differs between rows.  The four-row layout
#' follows Table 1 in Frandsen, Leslie and McIntyre (2025); the function does
#' not replicate that table's empirical estimates.
#'
#' @inheritParams cjive
#' @param y Outcome (numeric vector); factor and character outcomes are
#'   rejected.  A formula as the first argument selects the formula method.
#' @param x Single endogenous regressor (numeric vector); factor and character
#'   regressors are rejected.
#' @param z Instruments (numeric matrix/vector or a grouping factor).  A factor
#'   uses reference coding with an intercept or absorbed fixed effects, and all
#'   levels only when neither is present.
#' @param formula One of the two restricted formula layouts documented in
#'   \code{\link{cjive}}.  Sections accept only additive lists of bare variable
#'   names, with exactly one outcome and one endogenous regressor; \code{0} or
#'   \code{-1} removes the intercept, while other arithmetic and subtraction
#'   are rejected.
#' @param na.action How to treat missing values (formula method only):
#'   \code{\link[stats]{na.omit}} (the default) drops every row with a missing
#'   value in any model component -- outcome, endogenous regressor,
#'   instruments, controls, fixed effects, cluster or weights -- in a single
#'   complete-cases filter; \code{na.fail} stops on any missing value;
#'   \code{na.pass} leaves them in, so the input validation reports them.
#'   Unlike the fitted-object interfaces, \code{iv_compare} returns a plain
#'   data frame and does not store or report the number of dropped rows.  The
#'   default method (vector interface) does not handle missing values at all:
#'   it requires clean input and stops otherwise.
#' @param ... Must be empty.  Unknown, misspelled and unnamed arguments are
#'   errors.
#'
#' @return A data frame with one row per estimator (in the order OLS, 2SLS,
#'   JIVE, CJIVE) and columns \code{estimator}, \code{coefficient}, \code{se},
#'   \code{statistic}, \code{p.value}, \code{conf.low}, \code{conf.high}.
#'
#' @details The constructed instruments are: OLS, the residualised \eqn{x}
#'   itself; 2SLS, the full-sample fit \eqn{Z\hat\pi}; JIVE, the leave-one-out
#'   fit \eqn{(\hat x - h x)/(1 - h)}; CJIVE, the leave-cluster-out block fit.
#'   The CJIVE row equals \code{cjive(..., method = "dense")} on the same design.
#'   Weighting and residualisation precede validation: transformed values must
#'   be finite, \code{y}, \code{x} and each instrument must retain usable
#'   variation, and the residualised instrument matrix must have full numerical
#'   column rank.  Absorbed or collinear instruments are not silently dropped;
#'   a numerically zero estimator denominator is rejected before division.
#'
#'   Because covariates (and the intercept) are partialled out by
#'   Frisch-Waugh-Lovell before the jackknife -- the package-wide convention --
#'   the leverage \eqn{h} and the fit \eqn{\hat x} in the \code{"JIVE"} row are
#'   computed on the \emph{residualised} instruments.  This is the improved
#'   JIVE (IJIVE) of Ackerberg and Devereux (2009), not the original JIVE of
#'   Angrist, Imbens and Krueger (1999), which jackknifes the covariates
#'   alongside the instruments; residualising first removes the
#'   covariate-count bias that the original estimator carries, and the two
#'   coincide only when there are no covariates to partial out.  The row keeps
#'   the historical label \code{"JIVE"} for backwards compatibility with
#'   clusterIV 0.1.0; FLM label the corresponding estimator IJIVE.
#'
#' @seealso \code{\link{cjive}} for the CJIVE fit with diagnostics;
#'   \code{\link{iv_infer}} for the recommended inference workflow.  Two
#'   precomputed worked examples on real data:
#'   \code{vignette("queens-workflow", package = "clusterIV")} and
#'   \code{vignette("miami-bail", package = "clusterIV")}.
#'
#' @references
#' Ackerberg, D. A. and Devereux, P. J. (2009). Improved JIVE estimators for
#' overidentified linear models with and without heteroskedasticity.
#' \emph{Review of Economics and Statistics}, 91(2), 351--362.  The
#' partial-out-then-jackknife construction of the \code{"JIVE"} row.
#'
#' Angrist, J. D., Imbens, G. W. and Krueger, A. B. (1999). Jackknife
#' instrumental variables estimation. \emph{Journal of Applied Econometrics},
#' 14(1), 57--67.  The original JIVE that IJIVE improves upon.
#'
#' Frandsen, B., Leslie, E. and McIntyre, S. (2025). Cluster Jackknife
#' Instrumental Variables Estimation. \emph{Review of Economics and Statistics}.
#'
#' @examples
#' set.seed(2)
#' G <- 50; ng <- 5; n <- G * ng
#' cl <- rep(seq_len(G), each = ng)
#' z  <- matrix(rnorm(n * 3), n, 3)
#' u  <- rnorm(G)[cl]
#' x  <- z %*% c(1, -1, 0.5) + u + rnorm(n)
#' y  <- 2 * x + u + rnorm(n)
#' iv_compare(y, x, z, cluster = cl)
#'
#' @export
iv_compare <- function(y, ...) {
  .check_dispatch_dots(
    y,
    match.call(expand.dots = FALSE)$...,
    c("x", "z", "cluster", "controls", "fixed_effects", "weights",
      "intercept", "level", "inference", "data", "subset", "na.action"),
    "iv_compare", parent.frame()
  )
  UseMethod("iv_compare")
}

#' @rdname iv_compare
#' @export
iv_compare.default <- function(y, x, z, cluster, controls = NULL,
                               weights = NULL, level = 0.95,
                               intercept = TRUE, fixed_effects = NULL,
                               inference = c("asymptotic", "t"), ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "iv_compare")
  .check_flag(intercept, "intercept")
  .check_level(level)
  inference <- match.arg(inference)
  d <- .prep_data(y, x, z, cluster, controls, weights, intercept,
                  fixed_effects = fixed_effects)
  fs <- .first_stage(d$x, d$Z, R = d$R)

  # Each estimator differs only in the constructed instrument.
  phats <- list(
    OLS   = d$x,
    `2SLS` = .phat_2sls(fs),
    JIVE  = .phat_jive(fs, d$x),
    CJIVE = .phat_cjive(fs, d$x, d$groups)
  )

  rows <- lapply(names(phats), function(nm) {
    inf <- .iv_inference(phats[[nm]], d$x, d$y, d$cluster, level,
                         inference = inference)
    data.frame(estimator = nm,
               coefficient = inf$coefficient, se = inf$se,
               statistic = inf$statistic, p.value = inf$p.value,
               conf.low = inf$conf.low, conf.high = inf$conf.high,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname iv_compare
#' @export
iv_compare.formula <- function(formula, data, cluster, controls = NULL,
                               weights = NULL, level = 0.95,
                               intercept = TRUE, fixed_effects = NULL,
                               inference = c("asymptotic", "t"),
                               subset, na.action = stats::na.omit, ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "iv_compare")
  .check_flag(intercept, "intercept")
  .check_level(level)
  .formula_method(formula, data, missing(data), substitute(cluster),
                  substitute(weights),
                  if (missing(subset)) NULL else substitute(subset),
                  controls, fixed_effects, na.action, intercept,
                  parent.frame(), match.call(),
                  fit = function(prep, intercept)
                    iv_compare.default(prep$y, prep$x, prep$z,
                                       cluster = prep$cluster,
                                       controls = prep$controls,
                                       fixed_effects = prep$fixed_effects,
                                       weights = prep$weights,
                                       intercept = intercept, level = level,
                                       inference = inference),
                  stamp = FALSE)
}
