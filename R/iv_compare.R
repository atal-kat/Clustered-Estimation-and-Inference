#' Compare IV estimators on a common cluster-robust SE (FLM Table 1)
#'
#' Returns OLS, 2SLS, JIVE and CJIVE for the same design, each reported with the
#' \emph{identical} just-identified cluster-robust IV sandwich standard error;
#' only the constructed instrument differs between rows.  This reproduces the
#' shape of Table 1 in Frandsen, Leslie and McIntyre (2025).
#'
#' @inheritParams cjive
#' @param y Outcome (numeric vector).
#' @param x Single endogenous regressor (numeric vector).
#' @param z Instruments (numeric matrix/vector or a grouping factor).
#'
#' @return A data frame with one row per estimator (in the order OLS, 2SLS,
#'   JIVE, CJIVE) and columns \code{estimator}, \code{coefficient}, \code{se},
#'   \code{statistic}, \code{p.value}, \code{conf.low}, \code{conf.high}.
#'
#' @details The constructed instruments are: OLS, the residualised \eqn{x}
#'   itself; 2SLS, the full-sample fit \eqn{Z\hat\pi}; JIVE, the leave-one-out
#'   fit \eqn{(\hat x - h x)/(1 - h)}; CJIVE, the leave-cluster-out block fit.
#'   The CJIVE row equals \code{cjive(..., method = "dense")} on the same design.
#'
#' @references
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
iv_compare <- function(y, x, z, cluster, controls = NULL, weights = NULL,
                       level = 0.95, intercept = TRUE) {
  d <- .prep_data(y, x, z, cluster, controls, weights, intercept)
  fs <- .first_stage(d$x, d$Z)

  # Each estimator differs only in the constructed instrument.
  phats <- list(
    OLS   = d$x,
    `2SLS` = .phat_2sls(fs),
    JIVE  = .phat_jive(fs, d$x),
    CJIVE = .phat_cjive(fs, d$x, d$Z, d$groups)
  )

  rows <- lapply(names(phats), function(nm) {
    inf <- .iv_inference(phats[[nm]], d$x, d$y, d$cluster, level)
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
