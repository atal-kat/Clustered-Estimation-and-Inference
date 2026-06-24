#' Cluster-jackknife instrumental variables estimation (CJIVE)
#'
#' Computes the cluster-jackknife IV estimator of Frandsen, Leslie and McIntyre
#' (2025) for a single endogenous regressor in a just-identified design, with
#' cluster-robust inference.  The first-stage value for each observation is
#' fitted from a regression that leaves out the observation's entire cluster,
#' which removes the many-instrument bias that survives clustering.
#'
#' @param y Outcome (numeric vector), or a two-sided formula \code{y ~ x | z}
#'   for the formula method (the bar separates the endogenous regressor from the
#'   instruments).
#' @param x Single endogenous regressor (numeric vector).
#' @param z Instruments: a numeric vector/matrix, or a factor/character grouping
#'   vector for a judge design (expanded internally to a dummy design with one
#'   reference level dropped, the intercept supplying the rest).
#' @param cluster Cluster identifiers (length n).  For the formula method a
#'   one-sided formula (\code{~ g}) or a column name is also accepted.
#' @param formula A formula \code{y ~ x | z1 + z2}.
#' @param data A data frame in which to evaluate the formula.
#' @param controls Optional covariates (FLM's \eqn{X}): a matrix or data frame,
#'   or a one-sided formula in the formula method.  May be rank deficient (fixed
#'   effects are allowed).  An intercept is added unless \code{intercept = FALSE}.
#' @param weights Optional strictly positive precision weights.
#' @param level Confidence level for the reported interval.
#' @param intercept Logical; partial out an intercept (default \code{TRUE}).
#' @param method One of \code{"auto"}, \code{"dense"}, \code{"leaveout_mean"}.
#'   \code{"auto"} and \code{"dense"} both use the dense Frisch-Waugh-Lovell
#'   block-jackknife and are the default.  \code{"leaveout_mean"} evaluates FLM's
#'   printed leave-cluster-out group-mean form and is available only for a
#'   grouping-factor \code{z} with intercept-only controls; it differs from the
#'   default by an O(1/n_g) intercept term and is never selected automatically.
#' @param ... Unused.
#'
#' @return An object of class \code{"cjive"}: a list with \code{coefficient},
#'   \code{se}, \code{statistic}, \code{p.value}, \code{conf.low},
#'   \code{conf.high}, \code{level}, the diagnostics \code{n}, \code{G}, \code{p},
#'   \code{path} (\code{"dense"} or \code{"leaveout_mean"}) and \code{maxlev}
#'   (the maximum within-cluster leverage \eqn{\max_g \lambda_{\max}(H_g)}, a
#'   conditioning diagnostic; \code{NA} on the mean path), and the \code{call}.
#'
#' @details
#' The estimator is the covariance ratio
#' \eqn{\hat\delta = \widehat{Cov}(Y,\hat p)/\widehat{Cov}(D,\hat p)} with the
#' cluster-jackknife constructed instrument \eqn{\hat p}.  Covariates are handled
#' by Frisch-Waugh-Lovell: \eqn{Y}, \eqn{D} and each instrument are residualised
#' on the covariates (with an intercept) once, up front, then the estimator runs
#' on the residuals.  This dense route is the single convention everywhere, so
#' \code{cjive()} and \code{\link{iv_compare}} return the identical CJIVE for any
#' design.  The leave-cluster-out fits are computed by a Woodbury block update
#' (one Cholesky of \eqn{Z'Z} plus a small solve per cluster), exact against the
#' brute-force definition, and collapsing to observation-level JIVE when every
#' cluster is a singleton.
#'
#' @references
#' Frandsen, B., Leslie, E. and McIntyre, S. (2025). Cluster Jackknife
#' Instrumental Variables Estimation. \emph{Review of Economics and Statistics}.
#'
#' @examples
#' set.seed(1)
#' G  <- 40; ng <- 6; n <- G * ng
#' cl <- rep(seq_len(G), each = ng)
#' j  <- factor(rep(rep(1:4, length.out = ng), G))   # judge identity
#' u  <- rnorm(G)[cl]
#' x  <- as.numeric(j) + u + rnorm(n)
#' y  <- 1.5 * x + u + rnorm(n)
#' fit <- cjive(y, x, j, cluster = cl)
#' print(fit)
#'
#' ## formula interface
#' dat <- data.frame(y = y, x = x, j = j, cl = cl)
#' cjive(y ~ x | j, data = dat, cluster = ~cl)
#'
#' @export
cjive <- function(y, ...) UseMethod("cjive")

#' @rdname cjive
#' @export
cjive.default <- function(y, x, z, cluster, controls = NULL, weights = NULL,
                          level = 0.95, intercept = TRUE,
                          method = c("auto", "dense", "leaveout_mean"), ...) {
  cl <- match.call()
  method <- match.arg(method)

  d <- .prep_data(y, x, z, cluster, controls, weights, intercept)

  if (method == "leaveout_mean") {
    if (!d$grouping)
      stop("method = \"leaveout_mean\" requires a grouping-factor `z`.", call. = FALSE)
    if (!is.null(controls))
      stop("method = \"leaveout_mean\" requires intercept-only controls.", call. = FALSE)

    phat <- .leaveout_mean(x, d$group, cluster, weights)
    # Inference in the same weighted, centred geometry as the dense route:
    # residualise outcome, regressor and constructed instrument on the intercept.
    po <- .partial_out(as.numeric(y), as.numeric(x), phat,
                       controls = NULL, weights = weights, intercept = intercept)
    inf <- .iv_inference(as.numeric(po$Z), po$x, po$y, d$cluster, level)
    path <- "leaveout_mean"
    maxlev <- NA_real_
  } else {
    fs <- .first_stage(d$x, d$Z)
    lo <- .leaveout_fit(d$x, d$Z, d$groups, fs)
    inf <- .iv_inference(lo$phat, d$x, d$y, d$cluster, level)
    path <- "dense"
    maxlev <- lo$maxlev
  }

  structure(c(inf[c("coefficient", "se", "statistic", "p.value",
                    "conf.low", "conf.high")],
              list(level = level, n = d$n, G = d$G, p = d$p,
                   path = path, maxlev = maxlev, call = cl)),
            class = "cjive")
}

#' @rdname cjive
#' @export
cjive.formula <- function(formula, data, cluster, controls = NULL,
                          weights = NULL, level = 0.95, intercept = TRUE,
                          method = c("auto", "dense", "leaveout_mean"), ...) {
  cl <- match.call()
  method <- match.arg(method)
  if (missing(data)) data <- environment(formula)

  parts <- .parse_iv_formula(formula, data)
  clval <- .eval_cluster(cluster, data, parent.frame())
  ctl <- .eval_controls(controls, data, parent.frame())
  wts <- if (is.null(weights)) NULL else .eval_side(weights, data, parent.frame())

  out <- cjive.default(parts$y, parts$x, parts$z, cluster = clval,
                       controls = ctl, weights = wts, level = level,
                       intercept = intercept, method = method)
  out$call <- cl
  out
}

# --- formula helpers -------------------------------------------------------

# Parse `y ~ x | z1 + z2` in base R: the LHS is the outcome, the part of the RHS
# left of `|` the single endogenous regressor, the part to its right the
# instruments.  Factors among the instruments expand to dummies (reference
# dropped); the intercept is partialled out downstream.
.parse_iv_formula <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be two-sided, of the form y ~ x | z.", call. = FALSE)
  rhs <- formula[[3L]]
  if (!(is.call(rhs) && identical(rhs[[1L]], as.name("|"))))
    stop("the right-hand side must be of the form x | z (endogenous | instruments).",
         call. = FALSE)

  y <- eval(formula[[2L]], data, environment(formula))
  x <- eval(rhs[[2L]], data, environment(formula))

  zform <- stats::reformulate(deparse(rhs[[3L]]))
  mm <- stats::model.matrix(zform, data)
  z <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  list(y = as.numeric(y), x = as.numeric(x), z = z)
}

.eval_side <- function(expr, data, env) {
  if (inherits(expr, "formula")) expr <- expr[[length(expr)]]
  eval(expr, data, env)
}

.eval_cluster <- function(cluster, data, env) {
  if (inherits(cluster, "formula")) return(eval(cluster[[length(cluster)]], data, env))
  if (is.character(cluster) && length(cluster) == 1L && cluster %in% names(data))
    return(data[[cluster]])
  cluster
}

.eval_controls <- function(controls, data, env) {
  if (is.null(controls)) return(NULL)
  if (inherits(controls, "formula")) {
    mm <- stats::model.matrix(controls, data)
    return(mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE])
  }
  controls
}
