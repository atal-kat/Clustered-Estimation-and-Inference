#' Tidiers for clusterIV objects
#'
#' Coherent \code{tidy()} and \code{glance()} methods for \code{cjive},
#' \code{cjar}, \code{cjscore} and \code{iv_infer} objects.  The methods return
#' plain data frames and are registered with \pkg{generics} when that optional
#' package is installed.  A conventional coefficient
#' table cannot represent every weak-ID-robust confidence-set topology
#' (disjoint or unbounded sets in particular); consumers needing the complete
#' set must read the \code{conf.set} list-column and the topology columns.
#' No optional package is required
#' at run time, and no adapter is claimed for packages that discard procedure
#' metadata.
#'
#' @details
#' A row represents exactly one inferential procedure.  The CJIVE row contains
#' its estimate, cluster-robust standard error, Wald statistic and p-value for
#' the zero null, and Wald interval.  Test rows contain no estimate or standard
#' error; their statistic and p-value refer to the explicit \code{null.value}
#' stored in the same row.  \code{tidy.iv_infer()} returns one row for each
#' component actually present, beginning with CJIVE.
#'
#' A confidence set is flattened into \code{conf.low}/\code{conf.high} only
#' when it is one finite, non-degenerate interval.  The authoritative endpoint
#' matrix is always retained in the ordinary list-column \code{conf.set}; the
#' columns \code{shape}, \code{n.components}, \code{unbounded.left}, and
#' \code{unbounded.right} describe its topology without relying on attributes.
#' In \code{glance.iv_infer()}, procedure-specific diagnostics use the
#' \code{cjar.*} and \code{cjscore.*} prefixes; a CJS strength statistic is
#' therefore never paired with a CJAR critical value or topology.  The
#' \code{variance.estimator} glance field is \code{NA} when no selected
#' jackknife test uses that choice.  Unknown or partially matched arguments are
#' errors.
#'
#' @param x A fitted \code{cjive}, \code{cjar}, \code{cjscore} or
#'   \code{iv_infer} object.
#' @param conf.int Include scalar \code{conf.low}/\code{conf.high} endpoints
#'   when the set is representable as one finite non-degenerate interval.
#' @param conf.level Confidence level for those endpoints.  Test inversion is
#'   recomputed from stored coefficients when necessary.
#' @param vcov,coef_rename Compatibility arguments used by downstream table
#'   tools.  \code{vcov} must be \code{NULL}; coefficient renaming is applied
#'   after extraction.
#' @param gof_map A downstream goodness-of-fit mapping, accepted and left for
#'   the calling tool to apply after extraction.
#' @param ... Must be empty.
#' @return A \code{data.frame}.  \code{tidy.iv_infer()} has one row per present
#'   procedure; the standalone tidiers and all glance methods have one row.
#' @seealso \code{\link{cjive}}, \code{\link{cjar}}, \code{\link{cjscore}},
#'   \code{\link{iv_infer}}.
#' @name clusterIV-tidiers
NULL

.conf_set_interval <- function(cs) {
  cs <- as.matrix(cs)
  if (nrow(cs) == 1L && ncol(cs) == 2L && all(is.finite(cs)) &&
      cs[1L, 1L] < cs[1L, 2L]) {
    unname(cs[1L, ])
  } else c(NA_real_, NA_real_)
}

.empty_conf_set <- function() {
  matrix(numeric(0), 0L, 2L,
         dimnames = list(NULL, c("lower", "upper")))
}

.tidy_row <- function(term, component, procedure, estimate, std.error,
                      statistic, p.value, null.value, conf_set, level,
                      shape, variance_estimator = NA_character_,
                      flatten = TRUE) {
  ci <- if (isTRUE(flatten) && !is.null(conf_set)) {
    .conf_set_interval(conf_set)
  } else c(NA_real_, NA_real_)
  tp <- if (is.null(conf_set)) {
    list(n_components = NA_integer_, unbounded_left = NA,
         unbounded_right = NA)
  } else .set_topology(conf_set)
  out <- data.frame(
    term = term,
    component = component,
    procedure = procedure,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    null.value = null.value,
    conf.low = ci[1L],
    conf.high = ci[2L],
    level = level,
    shape = shape,
    n.components = tp$n_components,
    unbounded.left = tp$unbounded_left,
    unbounded.right = tp$unbounded_right,
    variance.estimator = variance_estimator,
    stringsAsFactors = FALSE
  )
  out$conf.set <- I(list(conf_set))
  out
}

.tidy_cjive_row <- function(x, conf_set = NULL, level = x$level,
                            flatten = TRUE) {
  if (is.null(conf_set)) {
    conf_set <- matrix(c(x$conf.low, x$conf.high), 1L, 2L,
                       dimnames = list(NULL, c("lower", "upper")))
  }
  .tidy_row(
    term = .object_term(x), component = "cjive",
    procedure = "CJIVE/Wald", estimate = x$coefficient,
    std.error = x$se, statistic = x$statistic, p.value = x$p.value,
    null.value = 0, conf_set = conf_set, level = level, shape = "bounded",
    variance_estimator = "cluster-robust sandwich", flatten = flatten
  )
}

.tidy_test_row <- function(x, component, procedure, term = .object_term(x),
                           null_value = x$beta0,
                           variance_estimator = x$variance_estimator,
                           conf_set = x$conf_set, level = x$level,
                           shape = x$shape, flatten = TRUE) {
  .tidy_row(
    term, component, procedure, NA_real_, NA_real_, x$statistic, x$p.value,
    null_value, conf_set, level, shape, variance_estimator,
    flatten = flatten
  )
}

.tidy_conf_args <- function(conf.int, conf.level, fitted_level, fn) {
  .check_flag(conf.int, "conf.int")
  if (isTRUE(conf.int)) .check_level(conf.level)
  if (!isTRUE(conf.int)) conf.level <- fitted_level
  list(include = conf.int, level = conf.level)
}

# R partially matches named arguments before a method body runs. Inspect the
# original call so misspellings such as "conf.in" or "gof" cannot silently
# become supported adapter arguments.
.check_tidier_call <- function(cl, allowed, fn) {
  args <- as.list(cl)[-1L]
  nms <- names(args)
  if (is.null(nms)) return(invisible(NULL))
  bad <- unique(nms[nzchar(nms) & !(nms %in% allowed)])
  if (length(bad)) {
    stop(fn, "() received unsupported or partially matched argument name(s): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  invisible(NULL)
}

.tidy_adapter_args <- function(vcov, coef_rename) {
  if (!is.null(vcov)) {
    stop("custom `vcov` is not supported by clusterIV tidiers.", call. = FALSE)
  }
  if (!(is.logical(coef_rename) || is.function(coef_rename) ||
        is.character(coef_rename))) {
    stop("`coef_rename` must be logical, a function, or a character mapping.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' @rdname clusterIV-tidiers
tidy.cjive <- function(x, conf.int = TRUE, conf.level = x$level,
                       vcov = NULL, coef_rename = FALSE, ...) {
  .check_tidier_call(sys.call(),
                     c("x", "conf.int", "conf.level", "vcov", "coef_rename"),
                     "tidy.cjive")
  .check_dots(match.call(expand.dots = FALSE)$..., "tidy.cjive")
  .tidy_adapter_args(vcov, coef_rename)
  ca <- .tidy_conf_args(conf.int, conf.level, x$level, "tidy.cjive")
  cs <- if (identical(ca$level, x$level)) {
    matrix(c(x$conf.low, x$conf.high), 1L, 2L,
           dimnames = list(NULL, c("lower", "upper")))
  } else stats::confint(x, level = ca$level)
  .tidy_cjive_row(x, cs, ca$level, flatten = ca$include)
}

#' @rdname clusterIV-tidiers
glance.cjive <- function(x, gof_map = NULL, ...) {
  .check_tidier_call(sys.call(), c("x", "gof_map"), "glance.cjive")
  .check_dots(match.call(expand.dots = FALSE)$..., "glance.cjive")
  k_controls <- if (!is.null(x$k_controls) && length(x$k_controls) == 1L) {
    x$k_controls
  } else NA_integer_
  maxlev <- if (!is.null(x$maxlev) && length(x$maxlev) == 1L) {
    x$maxlev
  } else NA_real_
  path <- if (!is.null(x$path) && length(x$path) == 1L) {
    x$path
  } else NA_character_
  inference <- if (!is.null(x$inference) && length(x$inference) == 1L) {
    x$inference
  } else "asymptotic"
  data.frame(nobs = x$n, G = x$G,
             k = if (!is.null(x$k)) x$k else x$p,
             k_controls = k_controls, maxlev = maxlev, path = path,
             inference = inference, term = .object_term(x),
             stringsAsFactors = FALSE)
}

#' @rdname clusterIV-tidiers
tidy.cjar <- function(x, conf.int = TRUE, conf.level = x$level,
                      vcov = NULL, coef_rename = FALSE, ...) {
  .check_tidier_call(sys.call(),
                     c("x", "conf.int", "conf.level", "vcov", "coef_rename"),
                     "tidy.cjar")
  .check_dots(match.call(expand.dots = FALSE)$..., "tidy.cjar")
  .tidy_adapter_args(vcov, coef_rename)
  ca <- .tidy_conf_args(conf.int, conf.level, x$level, "tidy.cjar")
  cs <- if (identical(ca$level, x$level)) x$conf_set
        else stats::confint(x, level = ca$level)
  .tidy_test_row(x, "cjar", "cluster jackknife Anderson-Rubin",
                 variance_estimator = x$variance_estimator,
                 conf_set = cs, level = ca$level,
                 shape = .shape_from_set(cs), flatten = ca$include)
}

#' @rdname clusterIV-tidiers
glance.cjar <- function(x, gof_map = NULL, ...) {
  .check_tidier_call(sys.call(), c("x", "gof_map"), "glance.cjar")
  .check_dots(match.call(expand.dots = FALSE)$..., "glance.cjar")
  tp <- .set_topology(x$conf_set)
  data.frame(nobs = x$n, G = x$G, k = x$k, F_CJ = x$F_CJ, crit = x$crit,
             shape = x$shape, n.components = tp$n_components,
             unbounded.left = tp$unbounded_left,
             unbounded.right = tp$unbounded_right, bounded = x$bounded,
             maxlev = x$maxlev, level = x$level,
             calibration = x$calibration,
             variance.estimator = x$variance_estimator,
             stringsAsFactors = FALSE)
}

#' @rdname clusterIV-tidiers
tidy.cjscore <- function(x, conf.int = TRUE, conf.level = x$level,
                         vcov = NULL, coef_rename = FALSE, ...) {
  .check_tidier_call(sys.call(),
                     c("x", "conf.int", "conf.level", "vcov", "coef_rename"),
                     "tidy.cjscore")
  .check_dots(match.call(expand.dots = FALSE)$..., "tidy.cjscore")
  .tidy_adapter_args(vcov, coef_rename)
  ca <- .tidy_conf_args(conf.int, conf.level, x$level, "tidy.cjscore")
  cs <- if (identical(ca$level, x$level)) x$conf_set
        else stats::confint(x, level = ca$level)
  .tidy_test_row(x, "cjscore", "cluster jackknife score",
                 variance_estimator = x$variance_estimator,
                 conf_set = cs, level = ca$level,
                 shape = .shape_from_set(cs), flatten = ca$include)
}

#' @rdname clusterIV-tidiers
glance.cjscore <- function(x, gof_map = NULL, ...) {
  .check_tidier_call(sys.call(), c("x", "gof_map"), "glance.cjscore")
  .check_dots(match.call(expand.dots = FALSE)$..., "glance.cjscore")
  tp <- .set_topology(x$conf_set)
  data.frame(nobs = x$n, G = x$G, k = x$k, F_CJS = x$F_CJS,
             crit = x$crit, shape = x$shape,
             n.components = tp$n_components,
             unbounded.left = tp$unbounded_left,
             unbounded.right = tp$unbounded_right, bounded = x$bounded,
             maxlev = x$maxlev, level = x$level,
             variance.estimator = x$variance_estimator,
             stringsAsFactors = FALSE)
}

#' @rdname clusterIV-tidiers
tidy.iv_infer <- function(x, conf.int = TRUE, conf.level = x$level,
                          vcov = NULL, coef_rename = FALSE, ...) {
  .check_tidier_call(sys.call(),
                     c("x", "conf.int", "conf.level", "vcov", "coef_rename"),
                     "tidy.iv_infer")
  .check_dots(match.call(expand.dots = FALSE)$..., "tidy.iv_infer")
  .tidy_adapter_args(vcov, coef_rename)
  ca <- .tidy_conf_args(conf.int, conf.level, x$level, "tidy.iv_infer")
  term <- if (is.character(x$term) && length(x$term) == 1L) x$term
          else .object_term(x$cjive)
  cjive_set <- if (identical(ca$level, x$level)) {
    matrix(c(x$cjive$conf.low, x$cjive$conf.high), 1L, 2L,
           dimnames = list(NULL, c("lower", "upper")))
  } else stats::confint(x$cjive, level = ca$level)
  rows <- list(.tidy_cjive_row(x$cjive, cjive_set, ca$level,
                               flatten = ca$include))
  if (!is.null(x$cjar)) {
    cs <- if (identical(ca$level, x$level)) x$cjar$conf_set
          else stats::confint(x$cjar, level = ca$level)
    rows <- c(rows, list(.tidy_test_row(
      x$cjar, "cjar", "cluster jackknife Anderson-Rubin", term,
      variance_estimator = x$cjar$variance_estimator, conf_set = cs,
      level = ca$level, shape = .shape_from_set(cs), flatten = ca$include
    )))
  }
  if (!is.null(x$cjscore)) {
    cs <- if (identical(ca$level, x$level)) x$cjscore$conf_set
          else stats::confint(x$cjscore, level = ca$level)
    rows <- c(rows, list(.tidy_test_row(
      x$cjscore, "cjscore", "cluster jackknife score", term,
      variance_estimator = x$cjscore$variance_estimator, conf_set = cs,
      level = ca$level, shape = .shape_from_set(cs), flatten = ca$include
    )))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname clusterIV-tidiers
glance.iv_infer <- function(x, gof_map = NULL, ...) {
  .check_tidier_call(sys.call(), c("x", "gof_map"), "glance.iv_infer")
  .check_dots(match.call(expand.dots = FALSE)$..., "glance.iv_infer")
  ar <- x$cjar
  sc <- x$cjscore
  ar_tp <- if (!is.null(ar)) .set_topology(ar$conf_set) else
    list(n_components = NA_integer_, unbounded_left = NA,
         unbounded_right = NA)
  sc_tp <- if (!is.null(sc)) .set_topology(sc$conf_set) else
    list(n_components = NA_integer_, unbounded_left = NA,
         unbounded_right = NA)
  jackknife_used <- !is.null(ar) || !is.null(sc)
  jackknife_variance <- if (jackknife_used) {
    if (!is.null(x$variance) && length(x$variance) == 1L) x$variance
    else if (!is.null(ar)) ar$variance_estimator
    else sc$variance_estimator
  } else NA_character_
  data.frame(
    nobs = x$n, G = x$G, k = x$k,
    F_CJ = if (!is.null(ar)) ar$F_CJ else NA_real_,
    F_CJS = if (!is.null(sc)) sc$F_CJS else NA_real_,
    cjar.crit = if (!is.null(ar)) ar$crit else NA_real_,
    cjar.shape = if (!is.null(ar)) ar$shape else NA_character_,
    cjar.n.components = ar_tp$n_components,
    cjar.unbounded.left = ar_tp$unbounded_left,
    cjar.unbounded.right = ar_tp$unbounded_right,
    cjscore.crit = if (!is.null(sc)) sc$crit else NA_real_,
    cjscore.shape = if (!is.null(sc)) sc$shape else NA_character_,
    cjscore.n.components = sc_tp$n_components,
    cjscore.unbounded.left = sc_tp$unbounded_left,
    cjscore.unbounded.right = sc_tp$unbounded_right,
    maxlev = x$maxlev, k_controls = x$k_controls,
    variance.estimator = jackknife_variance,
    components = paste(c("cjive", x$tests), collapse = ","),
    stringsAsFactors = FALSE
  )
}
