#' Cluster jackknife IV estimation and weak-instrument-robust inference in one call
#'
#' Runs the recommended workflow of this package -- the CJIVE point estimate
#' (\code{\link{cjive}}), the CJAR confidence set (\code{\link{cjar}}), and
#' the CJS test of a point null (\code{\link{cjscore}}) -- with the
#' expensive preprocessing shared rather than repeated.
#' The three help pages describe this workflow and previously
#' required three calls; \code{iv_infer()} performs it with one data
#' preparation, one Cholesky factorisation of the instrument Gram matrix, one
#' whitened leave-cluster-out pass (fitted values and the \code{maxlev}
#' leverage diagnostic from the same loop) and one set of per-cluster
#' instrument sums shared by every selected test kernel.  Every computational
#' result agrees with the standalone calls on the same design.
#'
#' @inheritParams cjive
#' @param formula A model formula in either restricted layout documented in
#'   \code{\link{cjive}}: \code{y ~ x | z} (optionally
#'   \code{y ~ x | z | fe}) or \code{y ~ exog | fe | endo ~ inst}.  Formula
#'   sections are additive lists of bare variable names, with exactly one
#'   outcome and one endogenous regressor; arithmetic, transformed terms and
#'   subtraction are rejected except for the \code{0}/\code{-1} intercept
#'   markers.
#' @param level One finite numeric confidence level strictly between 0 and 1
#'   for the CJIVE interval and CJAR/CJS confidence sets (default 0.95).
#' @param beta0 One finite numeric null value of the structural coefficient at
#'   which the CJAR and CJS statistics and p-values are reported (default 0).
#' @param calibration CJAR critical-value calibration, \code{"chisq"} (the
#'   default) or \code{"normal"}; see \code{\link{cjar}}.
#' @param inference CJIVE critical values and p-values, \code{"asymptotic"}
#'   (the default) or \code{"t"}; see \code{\link{cjive}}.
#' @param tests Character vector selecting the test rows of the panel (the
#'   CJIVE estimate row is always present): any subset of \code{"cjar"} and
#'   \code{"cjscore"} (the cluster jackknife tests).  When omitted, the
#'   recommended \code{"cjar"} and \code{"cjscore"} pair is selected.
#'   \code{NULL} and
#'   \code{character(0)} request a CJIVE-only panel.  Values are uniquely
#'   resolved, deduplicated and stored in the canonical order shown above;
#'   unknown or ambiguous partial values are errors.  Component fields always
#'   exist and are \code{NULL} when deselected.
#' @param variance Variance estimator for selected CJAR/CJS components,
#'   \code{"plain"} (default) or \code{"crossfit"}; see \code{\link{cjar}}.
#'   The cross-fit calculation is shared when both components are selected and
#'   skipped when neither is selected.
#' @param ... Must be empty.  Unknown, misspelled and unnamed arguments are
#'   errors.
#'
#' @return An object of class \code{"iv_infer"}: a list with
#'   \describe{
#'     \item{\code{cjive}, \code{cjar}, \code{cjscore}}{objects of the three
#'       existing classes whose computational fields match the standalone calls
#'       \code{cjive()}, \code{cjar()} and \code{cjscore()} on the same
#'       design; only \code{call} records the enclosing panel invocation.
#'       Thus every existing method (\code{coef}, \code{confint},
#'       \code{print}, \code{nobs}, ...) keeps working on the components.
#'       The CJIVE component always uses the dense path.  \code{cjar} and
#'       \code{cjscore} are \code{NULL} when deselected via \code{tests}.}
#'     \item{\code{n}, \code{G}, \code{k}}{observations, clusters, and number
#'       of instrument columns after partialling.}
#'     \item{\code{maxlev}}{the shared within-cluster leverage diagnostic
#'       (see \code{\link{cjar}}).}
#'     \item{\code{F_eff}, \code{K_eff}, \code{F_eff_crit},
#'       \code{ng_max}}{the Montiel Olea-Pflueger effective first-stage F
#'       and its simplified-TSLS critical value (identical on every
#'       component), and the largest cluster size; see \code{\link{cjar}}.
#'       If the estimated denominator is exactly zero with nonzero signal,
#'       \code{F_eff = Inf} and both \code{K_eff} and \code{F_eff_crit} are
#'       \code{NA}.}
#'     \item{\code{fe_dims}, \code{fe_levels}, \code{k_controls},
#'       \code{n_dropped}}{as in \code{\link{cjive}}.}
#'     \item{\code{tests}}{the resolved, canonical \code{tests} selector.}
#'     \item{\code{term}, \code{beta0}, \code{level}, \code{calibration},
#'       \code{inference}, \code{variance}}{the endogenous-regressor label
#'       and panel-wide reporting choices.}
#'     \item{\code{call}}{the matched call.}
#'   }
#'
#' @details
#' Inputs follow the strict shared contract in \code{\link{cjive}}.  In
#' particular, the weighted, residualised data must be finite; \code{y},
#' \code{x} and every instrument must retain usable variation after dense FWL
#' or fixed-effect absorption; and the residualised instrument matrix must
#' have full numerical column rank.  Absorbed or collinear instrument columns
#' are reported rather than silently dropped.  The CJIVE component also stops
#' before division when its IV denominator is numerically zero.
#'
#' The shared work is preprocessing, not approximation: one data
#' preparation, one instrument-Gram Cholesky factorisation, one whitened
#' leave-cluster-out pass and one set of per-cluster instrument sums feed
#' the CJIVE first stage,
#' the CJAR quartic coefficients and the CJS coefficients; each selected
#' test then performs its own
#' coefficient and inversion work.  Component results agree fieldwise with
#' the standalone functions apart from \code{call}.  The three
#' advisories the selected jackknife tests issue (fewer than 20 clusters, a
#' dominating cluster, within-cluster leverage near 1) fire once per
#' \code{iv_infer()} call, not once per component; a CJIVE-only panel emits
#' none of those test advisories.
#'
#' Reporting conventions follow the component documentation: report the CJIVE
#' coefficient for magnitude, the CJAR set for inference (it remains valid
#' under weak and many instruments, unlike the CJIVE t-interval), and use the
#' CJS p-value for the specific point null \code{beta0}, where its
#' complementary power properties are relevant.
#'
#' @seealso \code{\link{cjive}}, \code{\link{cjar}}, \code{\link{cjscore}}
#'   for the components and their full documentation (statistics, confidence
#'   set shapes, controls and fixed effects, caveats);
#'   \code{\link{iv_compare}} for estimator comparison.  Two precomputed
#'   worked examples on real data:
#'   \code{vignette("queens-workflow", package = "clusterIV")} and
#'   \code{vignette("miami-bail", package = "clusterIV")}.
#'
#' @examples
#' ## Simulated judge design, as in ?cjar
#' set.seed(42)
#' G  <- 30; ng <- 8; n <- G * ng
#' cl <- rep(seq_len(G), each = ng)            # cluster identifier
#' judge <- factor(rep(1:6, length.out = n))   # judge identity (grouping z)
#' u  <- rnorm(G)[cl]                          # cluster-level shock
#' x  <- 0.6 * as.numeric(judge) + u + rnorm(n)
#' y  <- 0.5 * x + u + rnorm(n)
#'
#' ## One call: CJIVE estimate + CJAR set + CJS point-null test
#' fit <- iv_infer(y, x, judge, cluster = cl)
#' fit
#'
#' ## The components are ordinary fitted objects
#' coef(fit$cjive)
#' confint(fit$cjar)
#'
#' ## Formula interface, controls inside the formula
#' dat <- data.frame(y = y, x = x, judge = judge, cl = cl)
#' iv_infer(y ~ 1 | x ~ judge, data = dat, cluster = cl)
#'
#' @export
iv_infer <- function(y, ...) {
  .check_dispatch_dots(
    y,
    match.call(expand.dots = FALSE)$...,
    c("x", "z", "cluster", "controls", "fixed_effects", "weights",
      "intercept", "level", "beta0", "calibration", "inference", "tests",
      "variance", "data", "subset", "na.action"),
    "iv_infer", parent.frame()
  )
  UseMethod("iv_infer")
}

#' @rdname iv_infer
#' @export
iv_infer.default <- function(y, x, z, cluster, controls = NULL,
                             fixed_effects = NULL, weights = NULL,
                             intercept = TRUE, level = 0.95, beta0 = 0,
                             calibration = c("chisq", "normal"),
                             inference = c("asymptotic", "t"),
                             tests = c("cjar", "cjscore"),
                             variance = c("plain", "crossfit"),
                             ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "iv_infer")
  cl <- match.call()
  .check_flag(intercept, "intercept")
  calibration <- match.arg(calibration)
  inference <- match.arg(inference)
  tests <- .resolve_iv_tests(tests, omitted = missing(tests))
  variance <- match.arg(variance)
  .check_beta0_level(beta0, level)
  term <- .term_label(substitute(x))

  d <- .prep_data(y, x, z, cluster, controls, weights, intercept,
                  fixed_effects = fixed_effects)

  # One Cholesky factorisation of Z'Z, shared by all three components.
  R <- d$R

  # One whitened leave-cluster-out pass: the CJIVE constructed instrument and
  # the maxlev leverage diagnostic come from the same per-cluster loop.
  fs <- .first_stage(d$x, d$Z, R = R)
  lo <- .leaveout_fit(d$x, fs$Ztil, fs$t, d$groups, d$k)

  # One set of whitened per-cluster instrument sums, shared by every selected
  # test.  With tests = NULL/character(0), the fit is deliberately CJIVE-only.
  csums <- if (length(tests)) {
    .cluster_sums(d$y, d$x, d$Z, d$cluster, R = R)
  } else NULL

  # The Montiel Olea-Pflueger effective F, from the same first-stage pass.
  eff <- .eff_f(fs$Ztil, fs$t, fs$e, d$cluster)

  jackknife_tests <- intersect(tests, c("cjar", "cjscore"))
  if (length(jackknife_tests)) {
    label <- paste(ifelse(jackknife_tests == "cjar", "CJAR", "CJS"),
                   collapse = "/")
    .advise(d$G, d$groups, d$n, lo$maxlev, label)
  }

  inf <- .iv_inference(lo$phat, d$x, d$y, d$cluster, level,
                       inference = inference)

  # Cross-fitting changes only the shared variance polynomial and is therefore
  # computed once even when both jackknife tests are requested; with a single
  # jackknife test selected, only that test's polynomial is computed.
  cf <- if (variance == "crossfit" && length(jackknife_tests)) {
    .crossfit_var(d$y, d$x, d$Z, d$cluster, R = R, Ztil = fs$Ztil,
                  which = if (length(jackknife_tests) == 2L) "both"
                          else if (jackknife_tests == "cjar") "ar"
                          else "score")
  } else NULL

  cjive_fit <- .cjive_build(d, inf, level, inference, "dense", lo$maxlev,
                            eff, cl, term)
  cjar_fit <- if ("cjar" %in% tests) {
    .cjar_build(d, csums, lo$maxlev, eff, beta0, level, calibration, cl,
                variance = variance, cf = cf, term = term)
  } else NULL
  cjscore_fit <- if ("cjscore" %in% tests) {
    .cjs_build(d, csums, lo$maxlev, eff, beta0, level, cl,
               variance = variance, cf = cf, term = term)
  } else NULL

  structure(list(
    cjive = cjive_fit, cjar = cjar_fit, cjscore = cjscore_fit,
    tests = tests, term = term, beta0 = beta0, level = level,
    calibration = calibration, inference = inference,
    variance = variance,
    n = d$n, G = d$G, k = d$k, maxlev = lo$maxlev,
    F_eff = eff$F_eff, K_eff = eff$K_eff, F_eff_crit = eff$F_eff_crit,
    ng_max = max(lengths(d$groups)),
    fe_dims = d$fe_dims, fe_levels = d$fe_levels,
    k_controls = d$k_controls, n_dropped = 0L,
    call = cl), class = "iv_infer")
}

#' @rdname iv_infer
#' @export
iv_infer.formula <- function(formula, data, cluster, controls = NULL,
                             fixed_effects = NULL, weights = NULL,
                             subset, na.action = stats::na.omit,
                             intercept = TRUE, level = 0.95, beta0 = 0,
                             calibration = c("chisq", "normal"),
                             inference = c("asymptotic", "t"),
                             tests = c("cjar", "cjscore"),
                             variance = c("plain", "crossfit"),
                             ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "iv_infer")
  .check_flag(intercept, "intercept")
  .check_beta0_level(beta0, level)
  .formula_method(formula, data, missing(data), substitute(cluster),
                  substitute(weights),
                  if (missing(subset)) NULL else substitute(subset),
                  controls, fixed_effects, na.action, intercept,
                  parent.frame(), match.call(),
                  fit = function(prep, intercept)
                    iv_infer.default(prep$y, prep$x, prep$z,
                                     cluster = prep$cluster,
                                     controls = prep$controls,
                                     fixed_effects = prep$fixed_effects,
                                     weights = prep$weights,
                                     intercept = intercept, level = level,
                                     beta0 = beta0,
                                     calibration = calibration,
                                     inference = inference, tests = tests,
                                     variance = variance),
                  components = c("cjive", "cjar", "cjscore"))
}

#' @describeIn iv_infer Print the panel: CJIVE point estimate and SE, the
#'   CJAR confidence set with its shape, the CJS p-value at \code{beta0},
#'   \code{F_CJ} against its critical value, and \code{maxlev}.  Rows
#'   deselected via \code{tests} are omitted.
#' @param digits Number of significant digits to print.
#' @export
print.iv_infer <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  fmt <- function(v) format(v, digits = digits)
  present <- c(
    "CJIVE",
    if (!is.null(x$cjar)) "CJAR",
    if (!is.null(x$cjscore)) "CJS"
  )
  if (length(present) == 1L) {
    cat("Cluster jackknife IV estimation (CJIVE only)\n")
  } else {
    cat("Cluster IV inference panel (", paste(present, collapse = " + "),
        ")\n", sep = "")
  }
  .call_line(x)
  est <- x$cjive
  ar <- x$cjar
  sc <- x$cjscore
  cat(sprintf(paste0("\n  CJIVE/Wald (H0: beta = 0): coefficient = %s   ",
                     "cluster-robust SE = %s   %s = %s   p = %s\n"),
              fmt(est$coefficient), fmt(est$se),
              if (identical(est$inference, "t")) "t" else "z",
              fmt(est$statistic), format.pval(est$p.value, digits = digits)))
  cat(sprintf("    %g%% Wald interval = [%s, %s]\n", 100 * est$level,
              fmt(est$conf.low), fmt(est$conf.high)))
  if (!is.null(ar)) {
    cat(sprintf("  CJAR (H0: beta = %s): T = %s   one-sided p = %s\n",
                fmt(ar$beta0), fmt(ar$statistic),
                format.pval(ar$p.value, digits = digits)))
    if (is.na(ar$statistic)) {
      cat("    note: Vhat(beta0) <= 0; the statistic is undefined and beta0 is accepted (conservative convention)\n")
    }
    cat(sprintf("    %g%% confidence set = %s%s\n", 100 * ar$level,
                .format_conf_set(ar$conf_set, digits),
                .topology_note(ar$conf_set)))
  }
  if (!is.null(sc)) {
    cat(sprintf("  CJS (H0: beta = %s): LM = %s   p = %s\n",
                fmt(sc$beta0), fmt(sc$statistic),
                format.pval(sc$p.value, digits = digits)))
    if (is.na(sc$statistic)) {
      cat(sprintf("    note: Vhat(beta0) = %s <= 0; the LM statistic is undefined and beta0 is accepted (conservative convention)\n",
                  fmt(sc$variance)))
    }
    cat(sprintf("    %g%% confidence set = %s%s\n", 100 * sc$level,
                .format_conf_set(sc$conf_set, digits),
                .topology_note(sc$conf_set)))
  }
  if (!is.null(ar)) {
    cat(sprintf("  F_CJ = %s   critical value = %s\n",
                fmt(ar$F_CJ), fmt(ar$crit)))
  }
  if (!is.null(sc)) {
    cat(sprintf("  F_CJS^2 = %s   critical value = %s\n",
                fmt(sc$F_CJS^2), fmt(sc$crit)))
  }
  if (identical(x$variance, "crossfit") &&
      (!is.null(ar) || !is.null(sc))) {
    cat("  variance estimator for selected jackknife tests: cross-fit ",
        "(Ligtenberg 2025, Section 5.4)\n", sep = "")
  }
  .eff_f_line(x, digits)
  cat(sprintf("  n = %d   G = %d clusters   k = %d instruments\n",
              x$n, x$G, x$k))
  .maxlev_line(x, digits)
  .fe_line(x, digits)
  .n_dropped_line(x)
  invisible(x)
}

#' @describeIn iv_infer Number of observations used.
#' @param object A fitted \code{"iv_infer"} object.
#' @export
nobs.iv_infer <- function(object, ...) object$n
