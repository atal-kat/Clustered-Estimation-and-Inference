# p-value curve plots for the test objects. Base graphics only; evaluating
# the statistic at any beta is O(1) from the stored polynomial coefficients,
# so a 512-point curve costs nothing. The shaded confidence set is taken
# from the object's conf_set -- never re-derived from the plotted grid -- so
# the picture cannot silently disagree with the reported set.

# p(beta) for a cjar object, vectorised, under the object's calibration and
# the same conventions as the fit. Plain variance: Vhat clamped at 0,
# T := 0 where Vhat = 0. Cross-fit variance: p = 1 where Vhat <= 0 (the
# statistic is undefined and beta is accepted -- the conservative
# convention), so the curve cannot disagree with the reported set.
.pval_curve_cjar <- function(x, b) {
  .cjar_point_eval(
    x$coef_num, x$coef_var, x$k, b, x$calibration,
    nonpos_var = if (isTRUE(x$variance_estimator == "crossfit"))
      "accept" else "zero"
  )$p.value
}

# p(beta) for a cjscore object, vectorised, under the conservative
# non-positive-variance convention (p = 1 where Vhat <= 0). The paper's
# 1/sqrt(n) and 1/n scalings cancel in the ratio.
.pval_curve_cjs <- function(x, b) {
  .cjs_point_eval(x$coef_score, x$coef_var, b)$p.value
}

# Positive coefficient ratio raised to a fractional power.  Keep the ordinary
# arithmetic path unchanged, but when numerator / denominator overflows or
# underflows, evaluate the root on the binary-log scale and map it back with an
# exact power-of-two multiplier.  This matters after a change of outcome and
# regressor units: all stored polynomial coefficients can be finite even when
# their raw ratio is not.
.plot_root_ratio <- function(numerator, denominator, power) {
  if (!is.finite(numerator) || !is.finite(denominator) ||
      numerator <= 0 || denominator <= 0) {
    return(1)
  }

  ratio <- numerator / denominator
  if (is.finite(ratio) && ratio > 0) {
    out <- ratio^power
    if (is.finite(out) && out > 0) return(out)
  }

  log_scale <- power * (log2(numerator) - log2(denominator))
  if (!is.finite(log_scale)) return(1)
  exponent <- floor(log_scale)
  out <- .mul_pow2(2^(log_scale - exponent), exponent)
  if (is.finite(out) && out > 0) out else 1
}

# Default plot range, deterministic from the stored objects (documented in
# ?plot.cjar): the finite endpoints of conf_set padded by half their width
# when at least two distinct ones exist; otherwise centred on `centre` (the
# argmin of the statistic) with 4 times the natural beta scale `scale` on
# each side, widened to include any single finite endpoint.
.plot_range <- function(conf_set, centre, scale) {
  fe <- conf_set[is.finite(conf_set)]
  if (length(fe) >= 2L && diff(range(fe)) > 0) {
    w <- diff(range(fe))
    c(min(fe) - w / 2, max(fe) + w / 2)
  } else {
    lo <- centre - 4 * scale
    hi <- centre + 4 * scale
    if (length(fe)) {
      lo <- min(lo, min(fe) - scale)
      hi <- max(hi, max(fe) + scale)
    }
    c(lo, hi)
  }
}

# Clip the confidence set to the plotted range: each row intersected with
# [from, to], empty intersections dropped. This is the single source of the
# shaded regions, so plotting cannot redefine the reported set.
.clip_regions <- function(conf_set, from, to) {
  if (nrow(conf_set) == 0L) {
    return(matrix(numeric(0), 0L, 2L,
                  dimnames = list(NULL, c("lower", "upper"))))
  }
  lo <- pmax(conf_set[, 1L], from)
  hi <- pmin(conf_set[, 2L], to)
  keep <- lo <= hi
  matrix(c(lo[keep], hi[keep]), ncol = 2L,
         dimnames = list(NULL, c("lower", "upper")))
}

# Validate the public plotting contract before base graphics sees it.  `xlim`
# is deliberately represented by from/to, so accepting both would create two
# competing range specifications.
.validate_plot_spec <- function(from, to, n, ylim, dots) {
  if (!is.numeric(from) || length(from) != 1L || !is.finite(from) ||
      !is.numeric(to) || length(to) != 1L || !is.finite(to) || from >= to) {
    stop("`from` and `to` must be finite numeric scalars with from < to.",
         call. = FALSE)
  }
  if (!is.numeric(n) || is.factor(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a finite integer-like scalar greater than or equal to 2.",
         call. = FALSE)
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || any(!is.finite(ylim)) ||
      ylim[1L] >= ylim[2L]) {
    stop("`ylim` must contain two finite increasing numeric values.",
         call. = FALSE)
  }
  if (length(dots)) {
    nms <- names(dots)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("all additional plotting arguments in `...` must be named.",
           call. = FALSE)
    }
    reserved <- intersect(nms, c("x", "y", "xlim", "type"))
    if (length(reserved)) {
      stop("plot argument(s) ", paste(reserved, collapse = ", "),
           " conflict with the `from`/`to` p-value-curve contract.",
           call. = FALSE)
    }
  }
  as.integer(n)
}

# Shared drawing core: shaded set, p-value curve(s), alpha line, beta0 mark.
# Styles are public arguments rather than hard-coded graphics calls.
.plot_pcurve <- function(x, curves, from, to, n, main,
                         xlab = expression(beta), ylab = "p-value",
                         ylim = c(0, 1), col = "black", lty = 1L, lwd = 1.5,
                         shade.col = "grey88", shade.border = NA,
                         alpha.col = "black", alpha.lty = 3L, ...) {
  dots <- list(...)
  n <- .validate_plot_spec(from, to, n, ylim, dots)
  bs <- seq(from, to, length.out = n)
  alpha <- 1 - x$level

  # Publication styling, base graphics only: horizontal tick labels, short
  # outward ticks, an unboxed (L-shaped) frame, a non-bold right-sized title,
  # and enough right margin for the alpha label. Restored on exit so the
  # caller's graphics state is untouched. Each choice is a default the user
  # can still override through `...`.
  op <- graphics::par(c("las", "tcl", "mgp", "font.main", "cex.main",
                        "cex.lab", "cex.axis", "mar"))
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(las = 1, tcl = -0.3, mgp = c(2.3, 0.55, 0),
                font.main = 1L, cex.main = 1.05, cex.lab = 0.95,
                cex.axis = 0.85, mar = c(4.1, 4.1, 2.9, 3.4))
  if (is.null(dots[["bty"]])) dots[["bty"]] <- "n"

  do.call(graphics::plot,
          c(list(x = NA_real_, y = NA_real_, xlim = c(from, to), ylim = ylim,
                 xlab = xlab, ylab = ylab, main = main, type = "n"), dots))
  regs <- .clip_regions(x$conf_set, from, to)
  if (nrow(regs)) {
    for (i in seq_len(nrow(regs))) {
      if (regs[i, 1L] < regs[i, 2L]) {
        graphics::rect(regs[i, 1L], ylim[1L], regs[i, 2L], ylim[2L],
                       col = shade.col, border = shade.border)
      } else {
        # an isolated accepted point: a degenerate [b, b] row
        graphics::segments(regs[i, 1L], ylim[1L], regs[i, 1L], ylim[2L],
                           col = shade.col, lwd = max(1, lwd[1L]))
      }
    }
  }
  col <- rep_len(col, length(curves))
  lty <- rep_len(lty, length(curves))
  lwd <- rep_len(lwd, length(curves))
  # size line first, so the curves read on top of it
  graphics::abline(h = alpha, lty = alpha.lty, col = alpha.col, lwd = 1)
  graphics::mtext(bquote(alpha == .(alpha)), side = 4, at = alpha,
                  las = 1, cex = 0.8, line = 0.4, col = alpha.col)
  if (!is.null(x$beta0)) {
    graphics::abline(v = x$beta0, lty = 3L, col = "grey55")
  }
  for (i in seq_along(curves)) {
    graphics::lines(bs, curves[[i]](bs), lty = lty[i], col = col[i],
                    lwd = lwd[i])
  }
  if (!is.null(x$beta0)) {
    graphics::points(x$beta0, curves[[1L]](x$beta0), pch = 19, cex = 0.8,
                     col = col[1L])
  }
  if (length(curves) > 1L) {
    graphics::legend("topright", legend = names(curves), lty = lty, col = col,
                     lwd = lwd, bty = "n", cex = 0.85, seg.len = 2.4,
                     inset = c(0, 0.02))
  }
  invisible(regs)
}

#' Plot the p-value curve of a cluster jackknife test
#'
#' Draws the p-value \eqn{p(\beta)} of the CJAR (or CJS) test as a function
#' of the hypothesised coefficient, with the reported confidence set shaded,
#' a horizontal line at \eqn{\alpha = 1 - \mathrm{level}}, and \code{beta0}
#' marked.  The curve is evaluated from the polynomial coefficients stored
#' on the object, so no refitting occurs; the shaded regions are taken from
#' \code{x$conf_set} -- never re-derived from the plotted grid -- so the
#' picture cannot silently disagree with the reported set.  An unbounded set
#' stops looking like a failure here: it is a curve that never dips below
#' the threshold in the tails (Dufour 1997).
#'
#' @param x A fitted \code{"cjar"}, \code{"cjscore"} or \code{"iv_infer"}
#'   object.
#' @param from,to When supplied, finite numeric scalar plot limits satisfying
#'   \code{from < to}.  Each omitted limit is taken from a range determined by
#'   the stored coefficients: the finite endpoints of the confidence set are
#'   padded by half their width when at least two distinct ones exist.
#'   Otherwise the range is centred on the CJAR numerator-quadratic vertex or
#'   the CJS score root, with four times the coefficient-based scale on each side
#'   (\eqn{(w_0/w_4)^{1/4}} or \eqn{\sqrt{v_0/v_2}} when the required
#'   coefficients are positive, and 1 otherwise), then widened to include any
#'   single finite endpoint.
#' @param n Finite integer-like number of grid points, at least 2 (default 512).
#' @param main Plot title.  For an \code{iv_infer} object, \code{NULL} chooses
#'   a title naming the test components actually present.
#' @param xlab,ylab Axis labels.
#' @param ylim Two finite increasing y-axis limits.
#' @param col,lty,lwd Curve colour, line type, and line width; vectors are
#'   recycled when the panel contains two curves.
#' @param shade.col,shade.border Fill and border for accepted regions.
#' @param alpha.col,alpha.lty Colour and line type of the test-size line.
#' @param ... Additional named arguments passed to \code{\link[graphics]{plot}}.
#'   Every argument must be named.  Use \code{from}/\code{to}, not
#'   \code{xlim}; \code{x}, \code{y}, \code{xlim}, and \code{type} conflict
#'   with the p-value-curve contract and are rejected.
#'
#' @return Invisibly, \code{x}.
#'
#' @details \code{plot.iv_infer} plots every available jackknife test curve.
#'   With CJAR present its confidence set supplies the shading; with CJS alone,
#'   the CJS set supplies it.  The title and legend identify the plotted tests.
#'   Positive-width accepted intervals are shaded and isolated accepted points
#'   (degenerate \code{[b, b]} rows) are drawn as vertical marks.  The CJIVE
#'   point estimate and the visible part of its Wald interval are marked beneath
#'   the curves when they overlap the requested horizontal range.
#'
#' @examples
#' set.seed(42)
#' G <- 30; ng <- 8; n <- G * ng
#' cl <- rep(seq_len(G), each = ng)
#' judge <- factor(rep(1:6, length.out = n))
#' u <- rnorm(G)[cl]
#' x <- 0.6 * as.numeric(judge) + u + rnorm(n)
#' y <- 0.5 * x + u + rnorm(n)
#' plot(cjar(y, x, judge, cluster = cl))
#' plot(iv_infer(y, x, judge, cluster = cl))
#'
#' @exportS3Method plot cjar
plot.cjar <- function(x, from = NULL, to = NULL, n = 512L,
                      main = "CJAR p-value curve", xlab = expression(beta),
                      ylab = "p-value", ylim = c(0, 1), col = "black",
                      lty = 1L, lwd = 1.5, shade.col = "grey88",
                      shade.border = NA, alpha.col = "black",
                      alpha.lty = 3L, ...) {
  nc <- x$coef_num; wc <- x$coef_var
  centre <- if (nc[3L] > 0) nc[2L] / (2 * nc[3L]) else x$beta0
  scale <- .plot_root_ratio(wc[1L], wc[5L], 0.25)
  rng <- .plot_range(x$conf_set, centre, scale)
  if (is.null(from)) from <- rng[1L]
  if (is.null(to)) to <- rng[2L]
  .plot_pcurve(x, list(CJAR = function(b) .pval_curve_cjar(x, b)),
               from, to, n, main, xlab = xlab, ylab = ylab, ylim = ylim,
               col = col, lty = lty, lwd = lwd, shade.col = shade.col,
               shade.border = shade.border, alpha.col = alpha.col,
               alpha.lty = alpha.lty, ...)
  invisible(x)
}

#' @rdname plot.cjar
#' @exportS3Method plot cjscore
plot.cjscore <- function(x, from = NULL, to = NULL, n = 512L,
                         main = "CJS p-value curve", xlab = expression(beta),
                         ylab = "p-value", ylim = c(0, 1), col = "black",
                         lty = 1L, lwd = 1.5, shade.col = "grey88",
                         shade.border = NA, alpha.col = "black",
                         alpha.lty = 3L, ...) {
  sc <- x$coef_score; vc <- x$coef_var
  centre <- if (sc[2L] != 0) sc[1L] / sc[2L] else x$beta0
  scale <- .plot_root_ratio(vc[1L], vc[3L], 0.5)
  rng <- .plot_range(x$conf_set, centre, scale)
  if (is.null(from)) from <- rng[1L]
  if (is.null(to)) to <- rng[2L]
  .plot_pcurve(x, list(CJS = function(b) .pval_curve_cjs(x, b)),
               from, to, n, main, xlab = xlab, ylab = ylab, ylim = ylim,
               col = col, lty = lty, lwd = lwd, shade.col = shade.col,
               shade.border = shade.border, alpha.col = alpha.col,
               alpha.lty = alpha.lty, ...)
  invisible(x)
}

#' @rdname plot.cjar
#' @exportS3Method plot iv_infer
plot.iv_infer <- function(x, from = NULL, to = NULL, n = 512L,
                          main = NULL, xlab = expression(beta),
                          ylab = "p-value", ylim = c(-0.1, 1),
                          col = c("black", "steelblue4"), lty = c(1L, 2L),
                          lwd = 1.5, shade.col = "grey88", shade.border = NA,
                          alpha.col = "black", alpha.lty = 3L, ...) {
  ar <- x$cjar
  sc <- x$cjscore
  if (is.null(ar) && is.null(sc)) {
    stop("plot.iv_infer() requires a selected CJAR or CJS component.",
         call. = FALSE)
  }
  primary <- if (!is.null(ar)) ar else sc
  est <- x$cjive
  if (!is.null(ar)) {
    nc <- ar$coef_num; wc <- ar$coef_var
    centre <- if (nc[3L] > 0) nc[2L] / (2 * nc[3L]) else ar$beta0
    scale <- .plot_root_ratio(wc[1L], wc[5L], 0.25)
  } else {
    sco <- sc$coef_score; vc <- sc$coef_var
    centre <- if (sco[2L] != 0) sco[1L] / sco[2L] else sc$beta0
    scale <- .plot_root_ratio(vc[1L], vc[3L], 0.5)
  }
  rng <- .plot_range(primary$conf_set, centre, scale)
  # widen to keep the CJIVE Wald interval in the frame
  rng[1L] <- min(rng[1L], est$conf.low)
  rng[2L] <- max(rng[2L], est$conf.high)
  if (is.null(from)) from <- rng[1L]
  if (is.null(to)) to <- rng[2L]
  curves <- list()
  if (!is.null(ar)) curves$CJAR <- function(b) .pval_curve_cjar(ar, b)
  if (!is.null(sc)) curves$CJS <- function(b) .pval_curve_cjs(sc, b)
  if (is.null(main)) {
    main <- paste(paste(names(curves), collapse = " and "),
                  if (length(curves) > 1L) "p-value curves"
                  else "p-value curve")
  }
  .plot_pcurve(primary, curves, from, to, n, main, xlab = xlab, ylab = ylab,
               ylim = ylim, col = col, lty = lty, lwd = lwd,
               shade.col = shade.col, shade.border = shade.border,
               alpha.col = alpha.col, alpha.lty = alpha.lty, ...)
  # CJIVE point estimate and its Wald interval, beneath the curves
  ymark <- ylim[1L] + 0.04 * diff(ylim)
  ci_lo <- max(est$conf.low, from)
  ci_hi <- min(est$conf.high, to)
  if (ci_lo <= ci_hi) {
    graphics::segments(ci_lo, ymark, ci_hi, ymark, lwd = 2)
  }
  if (est$coefficient >= from && est$coefficient <= to) {
    graphics::points(est$coefficient, ymark, pch = 18, cex = 1.3)
    graphics::mtext("CJIVE estimate with Wald interval", side = 1, line = -1.1,
                    at = est$coefficient, cex = 0.7)
  }
  invisible(x)
}
