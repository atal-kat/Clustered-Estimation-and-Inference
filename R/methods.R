#' @describeIn cjive Print a concise summary of the fit.
#' @param object A fitted \code{"cjive"} object.
#' @param digits Number of significant digits to print.
#' @export
print.cjive <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Cluster-jackknife IV (CJIVE)\n")
  .call_line(x)
  cat(sprintf("\n  coefficient = %s   cluster-robust SE = %s\n",
              format(x$coefficient, digits = digits),
              format(x$se, digits = digits)))
  cat(sprintf("  %s = %s   p = %s   %g%% CI = [%s, %s]\n",
              if (identical(x$inference, "t")) "t" else "z",
              format(x$statistic, digits = digits),
              format.pval(x$p.value, digits = digits),
              100 * x$level,
              format(x$conf.low, digits = digits),
              format(x$conf.high, digits = digits)))
  k <- if (!is.null(x$k)) x$k else x$p
  cat(sprintf("  n = %d   G = %d clusters   k = %d instruments   path = %s\n",
              x$n, x$G, k, x$path))
  .n_dropped_line(x)
  if (!is.null(x$inference) && x$inference != "asymptotic") {
    cat(sprintf("  inference: t(G-1) critical values\n"))
  }
  .fe_line(x, digits)
  .maxlev_line(x, digits)
  .eff_f_line(x, digits)
  invisible(x)
}

#' @describeIn cjive Number of observations used in the fit.
#' @export
nobs.cjive <- function(object, ...) object$n

#' @describeIn cjive Build a summary object.
#' @export
summary.cjive <- function(object, ...) {
  structure(object, class = c("summary.cjive", "cjive"))
}

#' @describeIn cjive Print method for the summary object, including the
#'   instrument-strength block (the Montiel Olea-Pflueger effective F;
#'   \code{F_CJ} requires a \code{\link{cjar}} fit -- \code{\link{iv_infer}}
#'   reports both in one call).
#' @export
print.summary.cjive <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  base <- x
  base$F_eff <- NA_real_
  class(base) <- "cjive"
  print(base, digits = digits, ...)
  cat("\nCoefficients:\n")
  ct <- cbind(x$coefficient, x$se, x$statistic, x$p.value)
  dimnames(ct) <- list(.object_term(x), c("Estimate", "Std. Error",
    if (identical(x$inference, "t")) c("t value", "Pr(>|t|)")
    else c("z value", "Pr(>|z|)")))
  stats::printCoefmat(ct, digits = digits, has.Pvalue = TRUE)
  if (!is.null(x$F_eff) && !is.na(x$F_eff)) {
    cat("\n")
    .strength_block(x,
      "  F_CJ requires a cjar() fit on the same design; iv_infer() reports the\n  estimate, both strength statistics and the robust confidence set in one call.\n",
      digits)
  }
  invisible(x)
}

#' @describeIn cjive Extract the point estimate.
#' @export
coef.cjive <- function(object, ...) {
  stats::setNames(object$coefficient, .object_term(object))
}

#' @describeIn cjive Extract the cluster-robust (co)variance.
#' @export
vcov.cjive <- function(object, ...) {
  term <- .object_term(object)
  v <- matrix(object$se^2, 1L, 1L, dimnames = list(term, term))
  v
}

#' @describeIn cjive Confidence interval for the coefficient.  \code{level}
#'   defaults to the level the object was fitted at.
#' @param parm Ignored (a single coefficient is estimated).
#' @export
confint.cjive <- function(object, parm, level = object$level, ...) {
  .check_level(level)
  a <- 1 - (1 - level) / 2
  zc <- if (identical(object$inference, "t")) stats::qt(a, df = object$G - 1)
        else stats::qnorm(a)
  ci <- object$coefficient + c(-1, 1) * zc * object$se
  out <- matrix(ci, 1L, 2L, dimnames = list(
    .object_term(object),
    paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), " %")))
  out
}

# A stored term was added after v0.1.0; old serialized objects remain usable.
.object_term <- function(x) {
  if (is.character(x$term) && length(x$term) == 1L && !is.na(x$term) &&
      nzchar(x$term)) x$term else "x"
}

# Topology is read from the complete endpoint matrix, not inferred from a
# strength diagnostic or the backwards-compatible tail-classification label.
.set_reading <- function(conf_set) {
  tp <- .set_topology(conf_set)
  if (tp$n_components == 0L) {
    return("no beta is accepted at this confidence level; inspect the specification, assumptions, and design diagnostics")
  }
  if (tp$n_components == 1L && tp$unbounded_left && tp$unbounded_right) {
    return("no beta is rejected: the confidence set is the whole real line")
  }
  if (tp$unbounded_left && tp$unbounded_right) {
    return(sprintf("unbounded in both tails with %d accepted component(s); all rejected gaps are shown above",
                   tp$n_components))
  }
  if (tp$unbounded_left || tp$unbounded_right) {
    return(sprintf("unbounded to the %s with %d accepted component(s); all rejected gaps are shown above",
                   if (tp$unbounded_left) "left" else "right",
                   tp$n_components))
  }
  if (tp$n_components == 1L) {
    "one bounded accepted component; every beta outside it is rejected"
  } else {
    sprintf("%d bounded accepted components; every beta outside them is rejected",
            tp$n_components)
  }
}

.topology_note <- function(conf_set) {
  tp <- .set_topology(conf_set)
  if (tp$n_components == 0L) {
    "  <- empty: no beta is accepted at this level"
  } else if (tp$unbounded_left || tp$unbounded_right) {
    paste0("  <- unbounded ",
           if (tp$unbounded_left && tp$unbounded_right) "in both tails"
           else if (tp$unbounded_left) "to the left" else "to the right",
           "; ", tp$n_components, " accepted component(s)")
  } else if (tp$n_components > 1L) {
    paste0("  <- ", tp$n_components, " bounded accepted components")
  } else ""
}

# The advisory conditions, re-evaluated from the stored fields (same
# thresholds as .advise), as printable lines for the test summaries.
.advisory_lines <- function(x) {
  out <- character(0)
  if (!is.null(x$G) && x$G < 20L)
    out <- c(out, sprintf("only %d clusters: the calibration is asymptotic in the number of clusters", x$G))
  if (!is.null(x$ng_max) && x$ng_max / x$n > 0.2)
    out <- c(out, sprintf("the largest cluster holds %.0f%% of the sample: the asymptotics require no dominating cluster", 100 * x$ng_max / x$n))
  if (!is.null(x$maxlev) && is.finite(x$maxlev) && x$maxlev > 0.99)
    out <- c(out, sprintf("max within-cluster leverage = %.4f: a cluster nearly spans the instrument space (Assumption A2(ii))", x$maxlev))
  out
}

# The shared instrument-strength block of the summaries: the cluster
# jackknife statistic and the Montiel
# Olea-Pflueger effective F (the familiar fixed-k bias diagnostic), side by
# side so divergence is visible. `cj` carries the jackknife line already
# formatted by the caller (it differs between AR and score).
.strength_block <- function(x, cj_line, digits) {
  cat("Instrument strength:\n")
  cat(cj_line)
  if (!is.null(x$F_eff) && !is.na(x$F_eff)) {
    cat(sprintf("  F_eff = %s   vs critical value %s   (Montiel Olea-Pflueger effective F, simplified TSLS, tau = 10%%, alpha = 5%%, K_eff = %s)\n",
                format(x$F_eff, digits = digits),
                format(x$F_eff_crit, digits = digits),
                format(x$K_eff, digits = digits)))
    cat("  (Confidence-set topology is read from the endpoint matrix printed\n")
    cat("   above; no reporting branch infers it from F_eff.)\n")
  }
  invisible(NULL)
}

# One line for the Montiel Olea-Pflueger effective F, printed BESIDE the
# cluster jackknife strength statistic (never instead of it, so divergence
# between the two is visible; see the ?cjar FAQ). Silent when the statistic
# is unavailable (the leaveout_mean path). Nothing in the package branches
# on F_eff.
.eff_f_line <- function(x, digits) {
  if (is.null(x$F_eff) || is.na(x$F_eff)) return(invisible(NULL))
  cat(sprintf("  effective F (Montiel Olea-Pflueger) = %s   critical value (tau = 10%%, alpha = 5%%) = %s   K_eff = %s\n",
              format(x$F_eff, digits = digits),
              format(x$F_eff_crit, digits = digits),
              format(x$K_eff, digits = digits)))
  invisible(NULL)
}

# The design-diagnostic lines shared verbatim by every fitted-object
# printer.  Byte-identical output is part of the contract: test-output.R
# captures printed panels.
.call_line <- function(x) {
  if (is.null(x$call)) return(invisible(NULL))
  cat("Call: ", paste(deparse(x$call), collapse = " "), "\n", sep = "")
  invisible(NULL)
}
.maxlev_line <- function(x, digits) {
  if (!is.finite(x$maxlev)) return(invisible(NULL))
  cat(sprintf("  max within-cluster leverage = %s%s\n",
              format(x$maxlev, digits = digits),
              if (x$maxlev > 0.99) "  <- near 1: conditioning frontier" else ""))
  invisible(NULL)
}
.fe_line <- function(x, digits) {
  if (is.null(x$fe_dims) || x$fe_dims == 0L) return(invisible(NULL))
  cat(sprintf("  absorbed fixed effects: %d dimension(s), %d levels   raw nuisance count/n = %s\n",
              x$fe_dims, x$fe_levels,
              format(x$k_controls / x$n, digits = digits)))
  invisible(NULL)
}
.n_dropped_line <- function(x) {
  if (is.null(x$n_dropped) || x$n_dropped == 0L) return(invisible(NULL))
  cat(sprintf("  %d observation(s) dropped due to missing values (na.action)\n",
              x$n_dropped))
  invisible(NULL)
}

# Honest headers for the test printers. At G == n every cluster is a
# singleton, Pddot is the projection matrix with a zeroed diagonal, and the
# cluster jackknife tests ARE the independent-data jackknife tests: with the
# plain variance the AR side is the Mikusheva-Sun (2022) statistic with the
# plain (naive) variance -- not their cross-fit test proper -- and the score
# side is exactly the Matsushita-Otsu (2024) jackknife LM test (statistic,
# variance and chi-square(1) calibration; see the FAQ in ?cjar and
# ?cjscore). With variance = "crossfit" neither named identity holds: the
# variance is Ligtenberg's (2025, Section 5.4) leave-out construction, which
# is not Mikusheva-Sun's Phihat_2 and not Matsushita-Otsu's estimator, so
# the label names the construction instead. The header is the only field
# that changes.
.cjar_header <- function(x) {
  if (x$G == x$n) {
    if (isTRUE(x$variance_estimator == "crossfit"))
      "Jackknife Anderson-Rubin test (independent data; cross-fit variance, Ligtenberg 2025 construction)"
    else
      "Jackknife Anderson-Rubin test (independent data; Mikusheva-Sun statistic, plain variance)"
  } else "Cluster jackknife Anderson-Rubin test (CJAR)"
}
.cjs_header <- function(x) {
  if (x$G == x$n) {
    if (isTRUE(x$variance_estimator == "crossfit"))
      "Jackknife score test (independent data; cross-fit variance, Ligtenberg 2025 construction)"
    else
      "Jackknife score test (independent data; Matsushita-Otsu 2024)"
  } else "Cluster jackknife score test (CJS)"
}

# Render an m x 2 confidence-set matrix as "[a, b] U [c, Inf)" (or "{}"),
# shared by print.cjar, print.cjscore and print.iv_infer.
.format_conf_set <- function(cs, digits) {
  fmt <- function(v) format(v, digits = digits)
  if (nrow(cs) == 0L) return("{}")
  paste(apply(cs, 1L, function(r) {
    paste0(if (is.finite(r[1L])) paste0("[", fmt(r[1L])) else "(-Inf",
           ", ",
           if (is.finite(r[2L])) paste0(fmt(r[2L]), "]") else "Inf)")
  }), collapse = " U ")
}

#' @describeIn cjar Print the statistic, the confidence set with its shape,
#'   and the F_CJ / leverage diagnostics.
#' @param digits Number of significant digits to print.
#' @export
print.cjar <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  fmt <- function(v) format(v, digits = digits)
  cat(.cjar_header(x), "\n", sep = "")
  .call_line(x)
  cat(sprintf("\n  H0: beta = %s   T = %s   one-sided p = %s\n",
              fmt(x$beta0), fmt(x$statistic),
              format.pval(x$p.value, digits = digits)))
  if (is.na(x$statistic)) {
    cat("  note: Vhat(beta0) <= 0; the statistic is undefined and beta0 is accepted (conservative convention, see ?cjar)\n")
  }

  set_str <- .format_conf_set(x$conf_set, digits)
  note <- .topology_note(x$conf_set)
  cat(sprintf("  %g%% confidence set = %s%s\n", 100 * x$level, set_str, note))
  cat(sprintf("  n = %d   G = %d clusters   k = %d instruments   calibration = %s\n",
              x$n, x$G, x$k, x$calibration))
  if (isTRUE(x$variance_estimator == "crossfit")) {
    cat("  variance estimator: cross-fit (Ligtenberg 2025, Section 5.4)\n")
  }
  .n_dropped_line(x)
  cat(sprintf("  F_CJ = %s   critical value = %s\n",
              fmt(x$F_CJ), fmt(x$crit)))
  .eff_f_line(x, digits)
  .maxlev_line(x, digits)
  .fe_line(x, digits)
  invisible(x)
}

#' @describeIn cjar Number of observations used in the test.
#' @param object A fitted \code{"cjar"} object.
#' @export
nobs.cjar <- function(object, ...) object$n

#' @describeIn cjar Confidence set for the structural coefficient.
#'   \code{level} defaults to the level the object was fitted at, where the
#'   stored set is returned; at any other level the inversion is recomputed
#'   from the stored polynomial coefficients (they are level-free), without
#'   refitting.
#' @param parm Ignored (the set concerns the single structural coefficient).
#' @export
confint.cjar <- function(object, parm, level = object$level, ...) {
  .check_level(level)
  if (identical(level, object$level)) return(object$conf_set)
  crit <- .cjar_crit(level, object$k, object$calibration)
  .cjar_invert(object$coef_num, object$coef_var, object$k, crit,
               nonpos_var = if (isTRUE(object$variance_estimator == "crossfit"))
                 "accept" else "zero")$conf_set
}

#' @describeIn cjar Build a summary object: the full diagnostic picture in
#'   one block -- statistic, confidence set with a plain-language reading of
#'   its shape, the instrument-strength block (\code{F_CJ} and the Montiel
#'   Olea-Pflueger effective F side by side), design diagnostics, and the
#'   advisories that apply.
#' @export
summary.cjar <- function(object, ...) {
  structure(object, class = c("summary.cjar", "cjar"))
}

#' @describeIn cjar Print method for the summary object.
#' @export
print.summary.cjar <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  fmt <- function(v) format(v, digits = digits)
  cat(.cjar_header(x), " -- summary\n", sep = "")
  .call_line(x)
  cat(sprintf("\n  H0: beta = %s   T = %s   one-sided p = %s\n",
              fmt(x$beta0), fmt(x$statistic),
              format.pval(x$p.value, digits = digits)))
  if (is.na(x$statistic)) {
    cat("  note: Vhat(beta0) <= 0; the statistic is undefined and beta0 is accepted (conservative convention)\n")
  }
  cat("\n")
  cat(sprintf("Confidence set (%g%%, %s calibration):\n", 100 * x$level,
              x$calibration))
  cat(sprintf("  %s   shape: %s\n", .format_conf_set(x$conf_set, digits),
              x$shape))
  cat(sprintf("  reading: %s\n", .set_reading(x$conf_set)))
  if (isTRUE(x$variance_estimator == "crossfit")) {
    cat("  variance estimator: cross-fit (Ligtenberg 2025, Section 5.4)\n")
  }
  cat("\n")
  .strength_block(x,
    sprintf("  F_CJ  = %s   vs critical value %s\n",
            fmt(x$F_CJ), fmt(x$crit)),
    digits)
  cat("\nDesign:\n")
  cat(sprintf("  n = %d   G = %d clusters   k = %d instruments\n",
              x$n, x$G, x$k))
  cat(sprintf("  max within-cluster leverage = %s   raw nuisance count/n = %s\n",
              fmt(x$maxlev), fmt(x$k_controls / x$n)))
  .n_dropped_line(x)
  adv <- .advisory_lines(x)
  cat("\nAdvisories:\n")
  if (length(adv)) cat(paste0("  - ", adv, "\n"), sep = "") else cat("  none\n")
  invisible(x)
}

#' @describeIn cjscore Print the statistic, the confidence set with its
#'   shape, and the F_CJS / leverage diagnostics.
#' @param digits Number of significant digits to print.
#' @export
print.cjscore <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  fmt <- function(v) format(v, digits = digits)
  cat(.cjs_header(x), "\n", sep = "")
  .call_line(x)
  cat(sprintf("\n  H0: beta = %s   LM = %s   p = %s\n",
              fmt(x$beta0), fmt(x$statistic),
              format.pval(x$p.value, digits = digits)))
  if (is.na(x$statistic)) {
    cat(sprintf("  note: Vhat(beta0) = %s <= 0; the LM statistic is undefined and beta0 is accepted (conservative convention, see ?cjscore)\n",
                fmt(x$variance)))
  }

  set_str <- .format_conf_set(x$conf_set, digits)
  note <- .topology_note(x$conf_set)
  cat(sprintf("  %g%% confidence set = %s%s\n", 100 * x$level, set_str, note))
  cat("  (score sets can contain spurious regions; report cjar() sets, use cjscore() for point nulls)\n")
  cat(sprintf("  n = %d   G = %d clusters   k = %d instruments\n",
              x$n, x$G, x$k))
  if (isTRUE(x$variance_estimator == "crossfit")) {
    cat("  variance estimator: cross-fit (Ligtenberg 2025, Section 5.4)\n")
  }
  .n_dropped_line(x)
  cat(sprintf("  F_CJS^2 = %s (first-stage LM statistic)   critical value (chisq_1) = %s\n",
              fmt(x$F_CJS^2), fmt(x$crit)))
  .eff_f_line(x, digits)
  .maxlev_line(x, digits)
  .fe_line(x, digits)
  invisible(x)
}

#' @describeIn cjscore Number of observations used in the test.
#' @param object A fitted \code{"cjscore"} object.
#' @export
nobs.cjscore <- function(object, ...) object$n

#' @describeIn cjscore Confidence set for the structural coefficient.
#'   \code{level} defaults to the level the object was fitted at, where the
#'   stored set is returned; at any other level the inversion is recomputed
#'   from the stored coefficients (they are level-free), without refitting.
#'   Mind the spurious-region caveat in \code{?cjscore}.
#' @param parm Ignored (the set concerns the single structural coefficient).
#' @export
confint.cjscore <- function(object, parm, level = object$level, ...) {
  .check_level(level)
  if (identical(level, object$level)) return(object$conf_set)
  .cjs_invert(object$coef_score, object$coef_var,
              stats::qchisq(level, df = 1))$conf_set
}

#' @describeIn cjscore Build a summary object: the full diagnostic picture
#'   in one block, as in \code{\link{summary.cjar}}.
#' @export
summary.cjscore <- function(object, ...) {
  structure(object, class = c("summary.cjscore", "cjscore"))
}

#' @describeIn cjscore Print method for the summary object.
#' @export
print.summary.cjscore <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  fmt <- function(v) format(v, digits = digits)
  cat(.cjs_header(x), " -- summary\n", sep = "")
  .call_line(x)
  cat(sprintf("\n  H0: beta = %s   LM = %s   p = %s\n",
              fmt(x$beta0), fmt(x$statistic),
              format.pval(x$p.value, digits = digits)))
  if (is.na(x$statistic)) {
    cat(sprintf("  note: Vhat(beta0) = %s <= 0; the LM statistic is undefined and beta0 is accepted (conservative convention)\n",
                fmt(x$variance)))
  }
  cat("\n")
  cat(sprintf("Confidence set (%g%%, chisq_1 calibration):\n", 100 * x$level))
  cat(sprintf("  %s   shape: %s\n", .format_conf_set(x$conf_set, digits),
              x$shape))
  cat(sprintf("  reading: %s\n", .set_reading(x$conf_set)))
  if (isTRUE(x$variance_estimator == "crossfit")) {
    cat("  variance estimator: cross-fit (Ligtenberg 2025, Section 5.4)\n")
  }
  cat("  (score sets can contain spurious regions; report cjar() sets, use cjscore() for point nulls)\n\n")
  .strength_block(x,
    sprintf("  F_CJS^2 = %s   vs critical value %s   (first-stage LM)\n",
            fmt(x$F_CJS^2), fmt(x$crit)),
    digits)
  cat("\nDesign:\n")
  cat(sprintf("  n = %d   G = %d clusters   k = %d instruments\n",
              x$n, x$G, x$k))
  cat(sprintf("  max within-cluster leverage = %s   raw nuisance count/n = %s\n",
              fmt(x$maxlev), fmt(x$k_controls / x$n)))
  .n_dropped_line(x)
  adv <- .advisory_lines(x)
  cat("\nAdvisories:\n")
  if (length(adv)) cat(paste0("  - ", adv, "\n"), sep = "") else cat("  none\n")
  invisible(x)
}
