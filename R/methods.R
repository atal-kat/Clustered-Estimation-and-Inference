#' @describeIn cjive Print a concise summary of the fit.
#' @param object,x A fitted \code{"cjive"} object.
#' @param digits Number of significant digits to print.
#' @export
print.cjive <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Cluster-jackknife IV (CJIVE)\n")
  if (!is.null(x$call)) {
    cat("Call: ", paste(deparse(x$call), collapse = " "), "\n", sep = "")
  }
  cat(sprintf("\n  coefficient = %s   cluster-robust SE = %s\n",
              format(x$coefficient, digits = digits),
              format(x$se, digits = digits)))
  cat(sprintf("  z = %s   p = %s   %g%% CI = [%s, %s]\n",
              format(x$statistic, digits = digits),
              format.pval(x$p.value, digits = digits),
              100 * x$level,
              format(x$conf.low, digits = digits),
              format(x$conf.high, digits = digits)))
  cat(sprintf("  n = %d   G = %d clusters   p = %d instruments   path = %s\n",
              x$n, x$G, x$p, x$path))
  if (is.finite(x$maxlev)) {
    cat(sprintf("  max within-cluster leverage = %s%s\n",
                format(x$maxlev, digits = digits),
                if (x$maxlev > 0.99) "  <- near 1: conditioning frontier" else ""))
  }
  invisible(x)
}

#' @describeIn cjive Build a summary object.
#' @export
summary.cjive <- function(object, ...) {
  structure(object, class = c("summary.cjive", "cjive"))
}

#' @describeIn cjive Print method for the summary object.
#' @export
print.summary.cjive <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  class(x) <- "cjive"
  print(x, digits = digits, ...)
  cat("\nCoefficients:\n")
  ct <- cbind(Estimate = x$coefficient, `Std. Error` = x$se,
              `z value` = x$statistic, `Pr(>|z|)` = x$p.value)
  rownames(ct) <- "x"
  stats::printCoefmat(ct, digits = digits, has.Pvalue = TRUE)
  invisible(x)
}

#' @describeIn cjive Extract the point estimate.
#' @export
coef.cjive <- function(object, ...) c(x = object$coefficient)

#' @describeIn cjive Extract the cluster-robust (co)variance.
#' @export
vcov.cjive <- function(object, ...) {
  v <- matrix(object$se^2, 1L, 1L, dimnames = list("x", "x"))
  v
}

#' @describeIn cjive Confidence interval for the coefficient.
#' @param parm Ignored (a single coefficient is estimated).
#' @export
confint.cjive <- function(object, parm, level = 0.95, ...) {
  zc <- stats::qnorm(1 - (1 - level) / 2)
  ci <- object$coefficient + c(-1, 1) * zc * object$se
  out <- matrix(ci, 1L, 2L, dimnames = list(
    "x", paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), " %")))
  out
}
