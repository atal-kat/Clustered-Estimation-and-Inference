# Formula front end shared by every *.formula method (cjive, cjar, cjscore,
# iv_compare, iv_infer). One parser, one evaluation pass, one NA policy --
# the .default methods receive plain vectors/matrices and never see a formula.
#
# Two restricted layouts are accepted:
#
#   instruments-first:  y ~ x | z            or  y ~ x | z | fe
#   controls-inside:    y ~ exog | endo ~ inst  or  y ~ exog | fe | endo ~ inst
#
# Dispatch is unambiguous: `~` binds more loosely than `|`, so the controls-inside
# form parses as `(y ~ exog | fe | endo) ~ inst` -- the parsed formula's
# left-hand side is itself a `~` call, which the legacy grammar can never
# produce (its RHS contains no `~`).

# Flatten the left-associative `|` calls of an expression into a list of
# top-level parts: `a | b | c` -> list(a, b, c).
.split_bars <- function(e) {
  if (is.call(e) && identical(e[[1L]], as.name("|")))
    c(.split_bars(e[[2L]]), list(e[[3L]]))
  else list(e)
}

# Split the restricted additive grammar into individual terms. Subtracting the
# literal 1 is the conventional no-intercept marker; every other use of `-`
# remains intact and is rejected by .parse_additive() below.
.is_number <- function(e, value) {
  is.numeric(e) && length(e) == 1L && !is.na(e) && e == value
}

.is_minus_one <- function(e) {
  .is_number(e, -1) ||
    (is.call(e) && length(e) == 2L &&
       identical(e[[1L]], as.name("-")) && .is_number(e[[2L]], 1))
}

.split_additive <- function(e) {
  if (is.call(e) && length(e) == 3L &&
      identical(e[[1L]], as.name("+"))) {
    return(c(.split_additive(e[[2L]]), .split_additive(e[[3L]])))
  }
  if (is.call(e) && length(e) == 3L &&
      identical(e[[1L]], as.name("-")) && .is_number(e[[3L]], 1)) {
    return(c(.split_additive(e[[2L]]), list(quote(-1))))
  }
  list(e)
}

.join_plus <- function(x) {
  if (length(x) == 0L) return(NULL)
  if (length(x) == 1L) return(x[[1L]])
  Reduce(function(a, b) call("+", a, b), x)
}

# Parse one section of the deliberately small formula language. Model terms are
# bare variable names joined by `+`. Only the structural endogenous section in
# the instruments-first layout and the exogenous section in the controls-inside
# carry 0 / 1 / -1 intercept markers.
.parse_additive <- function(expr, what, allow_intercept = FALSE,
                            exactly_one = FALSE, allow_empty = FALSE) {
  pieces <- .split_additive(expr)
  vars <- list()
  drop_intercept <- FALSE
  keep_intercept <- FALSE

  for (piece in pieces) {
    if (is.name(piece) && !identical(piece, as.name("."))) {
      vars <- c(vars, list(piece))
    } else if (allow_intercept && .is_number(piece, 0)) {
      drop_intercept <- TRUE
    } else if (allow_intercept && .is_number(piece, 1)) {
      keep_intercept <- TRUE
    } else if (allow_intercept && .is_minus_one(piece)) {
      drop_intercept <- TRUE
    } else {
      stop("`", what, "` supports only bare variable names joined by `+`",
           if (allow_intercept) " and the intercept markers 0, 1 or -1" else "",
           "; arithmetic, subtraction and transformed terms are not supported.",
           call. = FALSE)
    }
  }

  if (drop_intercept && keep_intercept)
    stop("`", what, "` cannot combine an intercept with 0 or -1.",
         call. = FALSE)

  labels <- vapply(vars, as.character, character(1))
  if (anyDuplicated(labels))
    stop("`", what, "` contains a duplicated variable term.", call. = FALSE)

  if (exactly_one && length(vars) != 1L)
    stop("`", what, "` must specify exactly one bare variable; found ",
         length(vars), ".", call. = FALSE)
  if (!exactly_one && !allow_empty && length(vars) == 0L)
    stop("`", what, "` must specify at least one bare variable.",
         call. = FALSE)

  list(expr = .join_plus(vars), variables = vars,
       drop_intercept = drop_intercept)
}

.match_na_action <- function(na.action) {
  if (identical(na.action, stats::na.omit) || identical(na.action, "na.omit"))
    return("omit")
  if (identical(na.action, stats::na.fail) || identical(na.action, "na.fail"))
    return("fail")
  if (identical(na.action, stats::na.pass) || identical(na.action, "na.pass"))
    return("pass")
  stop("`na.action` must be stats::na.omit, na.fail or na.pass.",
       call. = FALSE)
}

# One-sided formula `~ expr` carrying `env` for variable lookup outside
# `data` (the environment of the user's formula, as in lm()).
.onesided <- function(expr, env) {
  f <- eval(call("~", expr))
  environment(f) <- env
  f
}

# Model frame of the terms in `expr`, evaluated in `data` then `env`, with
# NO row dropping (na.action = NULL): missing values are handled once,
# centrally, in .prep_formula.
.eval_frame <- function(expr, data, env) {
  stats::model.frame(.onesided(expr, env), data, na.action = NULL)
}

# The `cluster` / `weights` convention (all four forms work): a bare column
# name (evaluated in `data` first, then the caller), a one-sided formula
# `~ col`, a column name as a string, or a full vector. `expr` arrives
# unevaluated via substitute() in the *.formula methods.
.has_arithmetic <- function(expr) {
  if (!is.call(expr)) return(FALSE)
  op <- as.character(expr[[1L]])
  if (op %in% c("+", "-", "*", "/", "^", ":", "%%", "%/%")) return(TRUE)
  any(vapply(as.list(expr)[-1L], .has_arithmetic, logical(1)))
}

.eval_column <- function(expr, data, env, what) {
  if (is.null(expr) || identical(expr, quote(expr = ))) return(NULL)
  if (identical(what, "cluster") &&
      !(is.call(expr) && identical(expr[[1L]], as.name("~"))) &&
      .has_arithmetic(expr)) {
    stop("`cluster` must be one grouping variable; arithmetic expressions ",
         "that combine cluster codes are not supported.", call. = FALSE)
  }
  v <- eval(expr, data, env)
  if (is.null(v)) return(NULL)
  if (inherits(v, "formula")) {
    if (length(v) != 2L)
      stop(sprintf("`%s` must be a one-sided formula (~ col).", what),
           call. = FALSE)
    parsed <- .parse_additive(v[[2L]], paste0(what, " formula"),
                              exactly_one = TRUE)
    v <- eval(parsed$expr, data, environment(v))
  } else if (is.character(v) && length(v) == 1L) {
    if (!is.null(names(data)) && v %in% names(data)) {
      v <- data[[v]]
    } else {
      stop(sprintf("`%s` = \"%s\": no such column in `data`.", what, v),
           call. = FALSE)
    }
  }
  v
}

# Row-subset any component: vectors by [i], everything two-dimensional by
# [i, ]. Model frames keep their terms attribute (dropped by `[`), so the
# post-filter model.matrix call sees a genuine model frame.
.take_rows <- function(v, rows) {
  if (is.null(v)) return(NULL)
  if (is.data.frame(v) || is.matrix(v)) {
    tt <- attr(v, "terms")
    v <- v[rows, , drop = FALSE]
    if (!is.null(tt)) attr(v, "terms") <- tt
    v
  } else {
    v[rows]
  }
}

# Expand a model frame to its (intercept-dropped) model matrix, after
# droplevels(): row filtering can empty a factor level, and the design that
# reaches validation must be the post-filter one.
.frame_to_matrix <- function(mf, intercept = TRUE) {
  tt <- attr(mf, "terms")
  attr(tt, "intercept") <- as.integer(isTRUE(intercept))
  mf <- droplevels(mf)
  attr(mf, "terms") <- tt
  mm <- stats::model.matrix(tt, mf)
  mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
}

# The single formula front end. `cluster`, `weights` and `subset` arrive as
# unevaluated expressions (substitute() in the calling *.formula method);
# `controls`, `fixed_effects` and `na.action` arrive as values. Returns the
# pieces the .default methods take, plus `n_dropped` (rows removed by the
# missing-value policy) and `drop_intercept` (a parsed 0 or -1 marker).
#
# Evaluation order (deliberate, see ?cjive):
#   1. every component is evaluated to full length n, with no row dropping;
#   2. `subset` is applied;
#   3. complete.cases() across all components decides the NA drop, once;
#   4. only then are factors droplevel()ed and dummies expanded, so the
#      validated design is the post-filter one.
.prep_formula <- function(formula, data, cluster, controls = NULL,
                          fixed_effects = NULL, weights = NULL, subset = NULL,
                          na.action = stats::na.omit, intercept = TRUE,
                          env = parent.frame()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be two-sided: y ~ x | z or ",
         "y ~ exog | fe | endo ~ inst (controls-inside-formula IV layout).",
         call. = FALSE)
  na_action <- .match_na_action(na.action)
  fenv <- environment(formula)
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  ctl_expr <- NULL
  fe_expr <- NULL
  drop_intercept <- FALSE

  if (is.call(lhs) && identical(lhs[[1L]], as.name("~"))) {
    # ---- controls-inside layout: (y ~ exog | fe | endo) ~ inst ------------
    if (length(lhs) != 3L)
      stop("cannot parse the controls-inside-formula IV layout: the part before the IV `~` ",
           "must itself be two-sided (y ~ exog | fe | endo ~ inst).",
           call. = FALSE)
    if (length(.split_bars(rhs)) > 1L)
      stop("in the controls-inside-formula IV layout the IV part (endo ~ inst) must come last: ",
           "y ~ exog | fe | endo ~ inst, with no `|` after the second `~`.",
           call. = FALSE)
    parts <- .split_bars(lhs[[3L]])
    if (length(parts) == 1L)
      stop("the controls-inside-formula IV layout needs an exogenous part before the IV part: ",
           "write y ~ 1 | endo ~ inst for a model without controls.",
           call. = FALSE)
    if (length(parts) > 3L)
      stop("too many `|` parts: the controls-inside-formula IV layout is y ~ exog | endo ~ inst ",
           "or y ~ exog | fe | endo ~ inst.", call. = FALSE)
    yexpr <- lhs[[2L]]
    xparsed <- .parse_additive(parts[[length(parts)]], "endogenous regressor",
                               exactly_one = TRUE)
    xexpr <- xparsed$expr
    zexpr <- rhs
    if (length(parts) == 3L) fe_expr <- parts[[2L]]
    exog <- .parse_additive(parts[[1L]], "exogenous part",
                            allow_intercept = TRUE, allow_empty = TRUE)
    ctl_expr <- exog$expr
    drop_intercept <- exog$drop_intercept
  } else {
    # ---- instruments-first layout: y ~ x | z or y ~ x | z | fe -----------
    parts <- .split_bars(rhs)
    if (length(parts) < 2L)
      stop("the right-hand side must be of the form x | z ",
           "(endogenous | instruments).", call. = FALSE)
    if (length(parts) > 3L)
      stop("too many `|` parts: the instruments-first layout is y ~ x | z or ",
           "y ~ x | z | fe.", call. = FALSE)
    yexpr <- lhs
    xparsed <- .parse_additive(parts[[1L]], "endogenous regressor",
                               allow_intercept = TRUE, exactly_one = TRUE)
    xexpr <- xparsed$expr
    drop_intercept <- xparsed$drop_intercept
    zexpr <- parts[[2L]]
    if (length(parts) == 3L) fe_expr <- parts[[3L]]
  }

  yexpr <- .parse_additive(yexpr, "outcome", exactly_one = TRUE)$expr
  zexpr <- .parse_additive(zexpr, "instrument part")$expr
  if (!is.null(fe_expr))
    fe_expr <- .parse_additive(fe_expr, "fixed-effect part")$expr

  effective_intercept <- isTRUE(intercept) && !drop_intercept

  if (!is.null(ctl_expr) && !is.null(controls))
    stop("supply controls either in the formula (the exogenous part) or via ",
         "`controls`, not both.", call. = FALSE)
  if (!is.null(fe_expr) && !is.null(fixed_effects))
    stop("supply fixed effects either in the formula or via `fixed_effects`, ",
         "not both.", call. = FALSE)

  # ---- 1. evaluate every component to full length, no row dropping --------
  y <- eval(yexpr, data, fenv)
  x <- eval(xexpr, data, fenv)
  n <- length(y)
  if (length(x) != n)
    stop("`x` and `y` must have the same length.", call. = FALSE)

  # A grouping z -- a single instrument term evaluating to a factor or
  # character vector -- passes through untouched, so the group-spread
  # validation and method = "leaveout_mean" work from the formula interface.
  # Anything else (multiple terms, numeric columns, a matrix) is expanded via
  # model.matrix after the NA filter.
  zf <- .eval_frame(zexpr, data, fenv)
  if (nrow(zf) != n) stop("`z` must have n rows.", call. = FALSE)
  z_grouping <- ncol(zf) == 1L &&
    (is.factor(zf[[1L]]) || is.character(zf[[1L]]))

  ctl_frame <- NULL
  ctl_raw <- NULL
  if (!is.null(ctl_expr)) {
    ctl_frame <- .eval_frame(ctl_expr, data, fenv)
  } else if (!is.null(controls)) {
    if (inherits(controls, "formula")) {
      if (length(controls) != 2L)
        stop("`controls` must be a one-sided formula.", call. = FALSE)
      ctl_expr_arg <- .parse_additive(controls[[2L]], "controls")$expr
      ctl_frame <- .eval_frame(ctl_expr_arg, data, environment(controls))
    } else {
      ctl_raw <- controls
    }
  }
  if (!is.null(ctl_frame) && nrow(ctl_frame) != n)
    stop("`controls` must have n rows.", call. = FALSE)
  if (!is.null(ctl_raw) && NROW(ctl_raw) != n)
    stop("`controls` must have n rows.", call. = FALSE)

  fe_list <- NULL
  if (!is.null(fe_expr)) {
    fe_list <- lapply(.eval_frame(fe_expr, data, fenv), as.factor)
  } else if (!is.null(fixed_effects)) {
    if (inherits(fixed_effects, "formula")) {
      if (length(fixed_effects) != 2L)
        stop("`fixed_effects` must be a one-sided formula.", call. = FALSE)
      fe_expr_arg <- .parse_additive(fixed_effects[[2L]],
                                     "fixed_effects")$expr
      fe_list <- lapply(.eval_frame(fe_expr_arg, data,
                                    environment(fixed_effects)), as.factor)
    } else if (is.list(fixed_effects)) {
      fe_list <- lapply(fixed_effects, as.factor)
    } else {
      fe_list <- list(as.factor(fixed_effects))
    }
  }
  if (!is.null(fe_list) && any(lengths(fe_list) != n))
    stop("each fixed-effect factor must have length n.", call. = FALSE)

  clv <- .eval_column(cluster, data, env, "cluster")
  if (is.null(clv)) stop("`cluster` is required.", call. = FALSE)
  if (length(clv) != n) stop("`cluster` must have length n.", call. = FALSE)
  wtv <- .eval_column(weights, data, env, "weights")
  if (!is.null(wtv) && length(wtv) != n)
    stop("`weights` must have length n.", call. = FALSE)

  # ---- 2. subset -----------------------------------------------------------
  rows <- seq_len(n)
  if (!is.null(subset)) {
    ss <- eval(subset, data, env)
    if (!is.null(ss)) {
      if (is.logical(ss)) {
        if (length(ss) != n)
          stop("a logical `subset` must have length n.", call. = FALSE)
        rows <- which(ss & !is.na(ss))
      } else {
        rows <- rows[ss]
      }
    }
  }

  # ---- 3. one missing-value decision across all components ----------------
  comps <- list(.take_rows(y, rows), .take_rows(x, rows),
                .take_rows(zf, rows), .take_rows(clv, rows))
  for (extra in list(.take_rows(ctl_frame, rows), .take_rows(ctl_raw, rows),
                     .take_rows(wtv, rows))) {
    if (!is.null(extra)) comps <- c(comps, list(extra))
  }
  if (!is.null(fe_list)) {
    comps <- c(comps, lapply(fe_list, .take_rows, rows = rows))
  }
  cc <- do.call(stats::complete.cases, comps)

  n_dropped <- 0L
  if (!all(cc)) {
    if (na_action == "omit") {
      n_dropped <- sum(!cc)
      rows <- rows[cc]
    } else if (na_action == "fail") {
      stop(sum(!cc), " observation(s) have missing values ",
           "(na.action = na.fail).", call. = FALSE)
    } else {
      # keep the rows; downstream validation reports the NAs
    }
  }

  # ---- 4. filter once, then droplevels and expand --------------------------
  y <- .take_rows(y, rows)
  x <- .take_rows(x, rows)
  clv <- .take_rows(clv, rows)
  wtv <- .take_rows(wtv, rows)

  if (z_grouping) {
    z <- droplevels(as.factor(.take_rows(zf, rows)[[1L]]))
  } else {
    z <- .frame_to_matrix(.take_rows(zf, rows), effective_intercept)
  }

  ctl <- if (!is.null(ctl_frame))
           .frame_to_matrix(.take_rows(ctl_frame, rows), effective_intercept)
         else .take_rows(ctl_raw, rows)

  if (!is.null(fe_list)) {
    fe_list <- lapply(fe_list, function(f) droplevels(.take_rows(f, rows)))
  }

  list(y = y, x = x, term = as.character(xparsed$variables[[1L]]),
       z = z, cluster = clv, controls = ctl,
       fixed_effects = fe_list, weights = wtv,
       n_dropped = n_dropped, drop_intercept = drop_intercept)
}

# Shared body of every public *.formula method.  Each method evaluates its
# non-standard-evaluation pieces in its own frame -- substitute(cluster),
# substitute(weights), substitute(subset), parent.frame(), match.call() --
# and passes the results here together with `fit`, a closure that calls its
# .default method as `fit(prep, intercept)` with the resolved intercept
# flag.  Fitted-object interfaces get `call`, `term` and `n_dropped`
# stamped on the result and on any named `components`; iv_compare() returns
# a plain data frame and opts out via `stamp = FALSE`.  `data` may be a
# missing promise; it is only touched after `data_missing` says it is safe.
.formula_method <- function(formula, data, data_missing, cluster_expr,
                            weights_expr, subset_expr, controls,
                            fixed_effects, na.action, intercept, env, cl,
                            fit, stamp = TRUE, components = character(0)) {
  if (data_missing) data <- environment(formula)
  prep <- .prep_formula(formula, data, cluster = cluster_expr,
                        controls = controls, fixed_effects = fixed_effects,
                        weights = weights_expr, subset = subset_expr,
                        na.action = na.action, intercept = intercept,
                        env = env)
  if (prep$drop_intercept) intercept <- FALSE
  out <- fit(prep, intercept)
  if (stamp) {
    out$call <- cl
    out$term <- prep$term
    out$n_dropped <- prep$n_dropped
    for (nm in components) {
      if (!is.null(out[[nm]])) {
        out[[nm]]$call <- cl
        out[[nm]]$term <- prep$term
        out[[nm]]$n_dropped <- prep$n_dropped
      }
    }
  }
  out
}
