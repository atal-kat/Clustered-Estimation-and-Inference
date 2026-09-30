# Input-space and strict-validation oracles (base R, no testthat).
#
# These checks deliberately collect failures and stop once at the end.  That
# makes the oracle-first run useful: one run against the pre-repair package
# shows every validation gap instead of stopping at the first accepted bad
# design.

library(clusterIV)

failures <- character(0)

record <- function(cond, msg, detail = NULL) {
  if (isTRUE(cond)) {
    cat("PASS:", msg, "\n")
  } else {
    line <- if (is.null(detail)) msg else paste0(msg, " [", detail, "]")
    failures <<- c(failures, line)
    cat("FAIL:", line, "\n")
  }
  invisible(cond)
}

public <- list(
  cjive = cjive,
  iv_compare = iv_compare,
  cjar = cjar,
  cjscore = cjscore,
  iv_infer = iv_infer
)

run_call <- function(name, args) {
  tryCatch(
    list(value = suppressWarnings(do.call(public[[name]], args)), error = NULL),
    error = function(e) list(value = NULL, error = conditionMessage(e))
  )
}

error_summary <- function(results) {
  paste(vapply(names(results), function(nm) {
    msg <- results[[nm]]$error
    paste0(nm, "=", if (is.null(msg)) "ACCEPT" else msg)
  }, character(1)), collapse = " | ")
}

expect_errors <- function(label, args, which_api = base::names(public),
                          component, diagnosis, same = TRUE) {
  ans <- lapply(which_api, function(nm) run_call(nm, args))
  names(ans) <- which_api
  msgs <- vapply(ans, function(a) {
    if (is.null(a$error)) NA_character_ else a$error
  }, character(1))
  clear <- !is.na(msgs) & grepl(component, msgs, ignore.case = TRUE) &
    grepl(diagnosis, msgs, ignore.case = TRUE)
  record(all(clear), label, error_summary(ans))
  if (same && all(!is.na(msgs))) {
    record(length(unique(msgs)) == 1L,
           paste0(label, ": shared APIs use one error message"),
           paste(unique(msgs), collapse = " | "))
  }
  invisible(ans)
}

expect_accepts <- function(label, args, which_api = base::names(public)) {
  ans <- lapply(which_api, function(nm) run_call(nm, args))
  names(ans) <- which_api
  record(all(vapply(ans, function(a) is.null(a$error), logical(1))),
         label, error_summary(ans))
  invisible(ans)
}

# Recursive comparison for fitted objects.  Calls differ between factor,
# explicit-matrix and formula routes and are intentionally skipped.  Numeric
# fields, including confidence-set matrices, are compared elementwise.
objects_near <- function(a, b, tol = 1e-12) {
  if (is.numeric(a) && is.numeric(b)) {
    if (length(a) != length(b) || !identical(dim(a), dim(b))) return(FALSE)
    same <- (is.na(a) & is.na(b)) |
      (is.infinite(a) & is.infinite(b) & sign(a) == sign(b)) |
      (is.finite(a) & is.finite(b) & abs(a - b) <= tol)
    return(all(same))
  }
  if (is.list(a) && is.list(b)) {
    if (!identical(class(a), class(b)) || !identical(names(a), names(b)))
      return(FALSE)
    for (nm in setdiff(names(a), "call")) {
      if (!objects_near(a[[nm]], b[[nm]], tol = tol)) return(FALSE)
    }
    return(TRUE)
  }
  identical(a, b)
}

expect_equivalent <- function(label, lhs, rhs, tol = 1e-12,
                              which_api = base::names(public)) {
  for (nm in which_api) {
    a <- run_call(nm, lhs)
    b <- run_call(nm, rhs)
    good <- is.null(a$error) && is.null(b$error) &&
      objects_near(a$value, b$value, tol = tol)
    detail <- if (!is.null(a$error) || !is.null(b$error)) {
      paste0("lhs=", if (is.null(a$error)) "OK" else a$error,
             "; rhs=", if (is.null(b$error)) "OK" else b$error)
    } else {
      "returned fields differ"
    }
    record(good, paste0(label, " [", nm, "]"), if (good) NULL else detail)
  }
  invisible(NULL)
}

formula_args <- function(formula, data, cluster, ...) {
  c(list(formula), list(data = data, cluster = cluster), list(...))
}

# Shared, well-conditioned design.  Every factor level and cluster contains
# several observations, and every instrument group spans many clusters.
set.seed(20260717)
G <- 24L
ng <- 6L
n <- G * ng
cl <- rep(seq_len(G), each = ng)
w1 <- rnorm(n)
w2 <- rnorm(n)
z1 <- rnorm(n)
z2 <- rnorm(n)
u <- rnorm(G)[cl]
x <- 0.9 * z1 - 0.4 * z2 + 0.3 * w1 + u + rnorm(n)
y <- 0.7 * x + 0.5 * w1 - 0.2 * w2 + u + rnorm(n)
wt <- runif(n, 0.5, 2)
Z <- cbind(z1 = z1, z2 = z2)
base_args <- list(y = y, x = x, z = Z, cluster = cl,
                  controls = cbind(w1 = w1))

# === A: a dense control exactly absorbs the instrument =====================
expect_errors(
  "[A] instrument exactly absorbed by a dense control",
  list(y = y, x = x, z = w1, cluster = cl, controls = cbind(w1 = w1)),
  component = "instrument|`z`",
  diagnosis = "partial|residu|absorb|variation|rank|collinear|singular"
)

# === B: a fixed effect exactly absorbs the instrument ======================
fe <- factor(rep(seq_len(12L), each = n / 12L))
z_fe <- as.numeric(fe == levels(fe)[1L])
expect_errors(
  "[B] instrument exactly absorbed by fixed effects",
  list(y = y, x = x, z = z_fe, cluster = cl, fixed_effects = fe),
  component = "instrument|`z`",
  diagnosis = "partial|residu|absorb|variation|rank|collinear|singular"
)

for (s in c(2^-200, 2^200)) {
  expect_errors(
    paste0("[B/E] FE-absorbed instrument remains rejected at scale ",
           format(s, scientific = TRUE)),
    list(y = y, x = x, z = z_fe * s, cluster = cl, fixed_effects = fe),
    component = "instrument|`z`",
    diagnosis = "partial|residu|absorb|variation|rank|collinear|singular"
  )
}

# === C: the endogenous regressor has no residualized variation ============
expect_errors(
  "[C] endogenous regressor exactly absorbed by a dense control",
  list(y = y, x = w1, z = Z, cluster = cl, controls = cbind(w1 = w1)),
  component = "`x`|endogenous",
  diagnosis = "partial|residu|absorb|variation|constant|unidentified"
)
expect_errors(
  "[C] endogenous regressor exactly absorbed by fixed effects",
  list(y = y, x = z_fe, z = Z, cluster = cl, fixed_effects = fe),
  component = "`x`|endogenous",
  diagnosis = "partial|residu|absorb|variation|constant|unidentified"
)

# === D: duplicate and post-FWL-collinear instruments =======================
expect_errors(
  "[D] duplicate instrument columns",
  list(y = y, x = x, z = cbind(z1, z1), cluster = cl),
  component = "instrument|`z`",
  diagnosis = "rank|collinear|singular"
)
expect_errors(
  "[D] instrument columns become collinear after dense FWL",
  list(y = y, x = x, z = cbind(z1, z1 + w1), cluster = cl,
       controls = cbind(w1 = w1)),
  component = "instrument|`z`",
  diagnosis = "partial|residu|rank|collinear|singular"
)
expect_errors(
  "[D] instrument columns become collinear after FE absorption",
  list(y = y, x = x, z = cbind(z1, z1 + z_fe), cluster = cl,
       fixed_effects = fe),
  component = "instrument|`z`",
  diagnosis = "partial|residu|rank|collinear|singular"
)

# === E: scale-metamorphic acceptance decisions =============================
# The largest statistic is quartic in y/x, so these factors are extreme while
# keeping every exact polynomial term well inside double's finite range.
small <- 2^-200
large <- 2^200
z_tiny <- 2^-600
z_huge <- 2^400
valid_scaled <- list(
  base = base_args,
  y_small = within(base_args, y <- y * small),
  y_large = within(base_args, y <- y * large),
  x_small = within(base_args, x <- x * small),
  x_large = within(base_args, x <- x * large),
  z_small = within(base_args, z <- z * small),
  z_large = within(base_args, z <- z * large),
  z_underflow_scale = within(base_args, z <- z * z_tiny),
  z_mixed_column_units = within(
    base_args,
    z <- sweep(z, 2L, c(z_tiny, z_huge), "*")
  ),
  controls_small = within(base_args, controls <- controls * small),
  controls_large = within(base_args, controls <- controls * large)
)
for (nm in names(valid_scaled)) {
  expect_accepts(paste0("[E] valid design remains accepted: ", nm),
                 valid_scaled[[nm]])
}

absorbed_scales <- list(
  equal = c(z = 1, control = 1),
  small_large = c(z = small, control = large),
  large_small = c(z = large, control = small)
)
for (nm in names(absorbed_scales)) {
  s <- absorbed_scales[[nm]]
  expect_errors(
    paste0("[E] absorbed instrument remains rejected: ", nm),
    list(y = y, x = x, z = w1 * s[["z"]], cluster = cl,
         controls = cbind(w1 = w1 * s[["control"]])),
    component = "instrument|`z`",
    diagnosis = "partial|residu|absorb|variation|rank|collinear|singular"
  )
}

# Raw inputs can all be finite while the sqrt(weight) transformation
# overflows.  The front end must diagnose non-finite transformed/residualized
# data rather than passing it to a factorization or returning an estimate.
y_over <- rep(c(0.75, -0.75) * .Machine$double.xmax, length.out = n)
expect_errors(
  "[E] non-finite transformed design is rejected before factorization",
  list(y = y_over, x = x, z = Z, cluster = cl, weights = rep(4, n),
       intercept = FALSE),
  component = "finite|outcome|`y`|residu|transform",
  diagnosis = "finite|overflow|residu|transform",
  same = FALSE
)

# A column can contain only finite entries while its Euclidean norm lies
# outside the representable range.  The equilibrated Gram path must diagnose
# that scale explicitly before crossprod() turns it into an opaque Inf/NaN
# factorization failure.
z_norm_over <- numeric(n)
z_norm_over[1:2] <- c(.Machine$double.xmax, -.Machine$double.xmax)
expect_errors(
  "[E] finite instrument with overflowing column norm gets a scale remedy",
  list(y = y, x = x, z = z_norm_over, cluster = cl, intercept = FALSE),
  component = "instrument|`z`",
  diagnosis = "column norm|representable|rescale"
)

# === F/G: independent dummy oracles for factor instruments =================
inst <- factor(rep(c("a", "b", "c", "d"), length.out = n),
               levels = c("a", "b", "c", "d"))
# This construction is intentionally independent of model.matrix() and every
# package expansion helper.
D_full <- vapply(levels(inst), function(lev) as.numeric(inst == lev),
                 numeric(n))
colnames(D_full) <- paste0("inst_", levels(inst))
D_ref <- D_full[, -1L, drop = FALSE]

expect_equivalent(
  "[F] no-intercept factor equals independently built full dummies",
  list(y = y, x = x, z = inst, cluster = cl, intercept = FALSE),
  list(y = y, x = x, z = D_full, cluster = cl, intercept = FALSE)
)
expect_equivalent(
  "[G] intercept factor equals independently built reference dummies",
  list(y = y, x = x, z = inst, cluster = cl, intercept = TRUE),
  list(y = y, x = x, z = D_ref, cluster = cl, intercept = TRUE)
)

# Raw factor-level spread is not the dense-path rank condition.  With global
# intercept FWL, a dummy for a level confined to one cluster remains a nonzero
# constant outside that cluster.  Both an omitted reference level and an
# encoded non-reference level can therefore have full-rank leave-out Grams.
# The factor and its independently constructed dummy matrix must take the same
# path through every public estimator/test.
inst_ref_local <- character(n)
inst_ref_local[cl == 1L] <- "a"
inst_ref_local[cl != 1L] <- rep(c("b", "c", "d"),
                                length.out = sum(cl != 1L))
inst_ref_local <- factor(inst_ref_local, levels = c("a", "b", "c", "d"))
D_ref_local <- vapply(c("b", "c", "d"),
                      function(lev) as.numeric(inst_ref_local == lev),
                      numeric(n))
colnames(D_ref_local) <- paste0("inst_", c("b", "c", "d"))
expect_equivalent(
  "[G] cluster-local omitted reference equals explicit encoded dummies",
  list(y = y, x = x, z = inst_ref_local, cluster = cl, intercept = TRUE),
  list(y = y, x = x, z = D_ref_local, cluster = cl, intercept = TRUE)
)

dat_ref_local <- data.frame(y = y, x = x, inst = inst_ref_local, cl = cl)
expect_equivalent(
  "[G/J] cluster-local omitted reference survives the formula path",
  formula_args(y ~ x | inst, dat_ref_local, dat_ref_local$cl),
  list(y = y, x = x, z = D_ref_local, cluster = cl, intercept = TRUE)
)

inst_encoded_local <- character(n)
inst_encoded_local[cl == 1L] <- "d"
inst_encoded_local[cl != 1L] <- rep(c("a", "b", "c"),
                                    length.out = sum(cl != 1L))
inst_encoded_local <- factor(inst_encoded_local,
                             levels = c("a", "b", "c", "d"))
D_encoded_local <- vapply(c("b", "c", "d"),
                          function(lev) as.numeric(inst_encoded_local == lev),
                          numeric(n))
colnames(D_encoded_local) <- paste0("inst_", c("b", "c", "d"))
expect_equivalent(
  "[G] cluster-local encoded level equals explicit encoded dummies",
  list(y = y, x = x, z = inst_encoded_local, cluster = cl,
       intercept = TRUE),
  list(y = y, x = x, z = D_encoded_local, cluster = cl,
       intercept = TRUE)
)

# The group-mean shortcut has a genuinely stronger support condition than the
# dense encoded-dummy path, including for a level omitted by reference coding.
expect_errors(
  "[G] leaveout_mean rejects a group with no outside-cluster mass",
  list(y = y, x = x, z = inst_ref_local, cluster = cl,
       method = "leaveout_mean"),
  which_api = "cjive", component = "group",
  diagnosis = "undefined|entirely.*one cluster", same = FALSE
)

# Conversely, a one-level factor with an intercept encodes no instrument
# column at all and is unsupported before any estimator-specific work.
expect_errors(
  "[G] one-level intercept factor supplies no instrument columns",
  list(y = y, x = x, z = factor(rep("only", n)), cluster = cl,
       intercept = TRUE),
  component = "instrument|`z`|grouping factor",
  diagnosis = "no instrument columns|at least 2 levels"
)

# A genuinely rank-deficient leave-out design is still rejected by the paths
# that require a leave-cluster-out first-stage solve.  Here omitting cluster 1
# leaves only two distinct encoded rows for three instrument columns.  Factor
# coding and the identical explicit matrix must fail for the same mathematical
# reason; the standalone weak-ID tests do not require that solve and remain
# available (with their leverage advisory).
G_sing <- 8L
ng_sing <- 8L
n_sing <- G_sing * ng_sing
cl_sing <- rep(seq_len(G_sing), each = ng_sing)
inst_sing_chr <- character(n_sing)
inst_sing_chr[cl_sing == 1L] <- rep(c("c", "d"), each = ng_sing / 2L)
inst_sing_chr[cl_sing != 1L] <- rep(c("a", "b"),
                                    length.out = sum(cl_sing != 1L))
inst_sing <- factor(inst_sing_chr, levels = c("a", "b", "c", "d"))
D_sing <- vapply(c("b", "c", "d"),
                 function(lev) as.numeric(inst_sing == lev), numeric(n_sing))
set.seed(20260718)
x_sing <- as.numeric(inst_sing) + rnorm(n_sing)
y_sing <- 0.7 * x_sing + rnorm(n_sing)
sing_factor_args <- list(y = y_sing, x = x_sing, z = inst_sing,
                         cluster = cl_sing)
sing_matrix_args <- list(y = y_sing, x = x_sing, z = D_sing,
                         cluster = cl_sing)
point_paths <- c("cjive", "iv_compare", "iv_infer")
expect_errors(
  "[G] factor design with singular leave-out Gram is rejected",
  sing_factor_args, which_api = point_paths,
  component = "leave-cluster-out|instrument",
  diagnosis = "singular|no variation outside"
)
expect_errors(
  "[G] explicit copy of singular leave-out design is rejected identically",
  sing_matrix_args, which_api = point_paths,
  component = "leave-cluster-out|instrument",
  diagnosis = "singular|no variation outside"
)
expect_equivalent(
  "[G] standalone tests agree where no leave-out solve is required",
  sing_factor_args, sing_matrix_args, which_api = c("cjar", "cjscore")
)

# Coding is a package contract, not a function of the session's contrast
# options or whether the grouping factor happens to be ordered.
expand_z <- getFromNamespace(".expand_z", "clusterIV")
old_contrasts <- getOption("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
got_sum <- unname(expand_z(inst, intercept = TRUE)$Z)
got_ordered <- unname(expand_z(ordered(inst), intercept = TRUE)$Z)
options(contrasts = old_contrasts)
record(identical(got_sum, unname(D_ref)),
       "[G] intercept factor coding ignores ambient contrast options")
record(identical(got_ordered, unname(D_ref)),
       "[G] ordered factor still uses treatment/reference coding")

one_level <- factor(rep("only", n))
expect_equivalent(
  "[F] no-intercept one-level factor equals one full dummy",
  list(y = y, x = x, z = one_level, cluster = cl, intercept = FALSE),
  list(y = y, x = x, z = matrix(1, n, 1L), cluster = cl,
       intercept = FALSE)
)

dat_factor <- data.frame(y = y, x = x, inst = inst, cl = cl)
ref_no_intercept <- list(y = y, x = x, z = D_full, cluster = cl,
                         intercept = FALSE)
ref_intercept <- list(y = y, x = x, z = D_ref, cluster = cl,
                      intercept = TRUE)

no_intercept_formulas <- list(
  `legacy 0 + term` = y ~ 0 + x | inst,
  `legacy -1` = y ~ x - 1 | inst,
  `fixest-style 0` = y ~ 0 | x ~ inst,
  `fixest-style -1` = y ~ -1 | x ~ inst
)
for (nm in names(no_intercept_formulas)) {
  expect_equivalent(
    paste0("[F/J] formula intercept removal propagates to factor coding: ", nm),
    formula_args(no_intercept_formulas[[nm]], dat_factor, dat_factor$cl),
    ref_no_intercept
  )
}

intercept_formulas <- list(
  legacy = y ~ x | inst,
  `fixest-style` = y ~ 1 | x ~ inst
)
for (nm in names(intercept_formulas)) {
  expect_equivalent(
    paste0("[G] formula intercept propagates to factor coding: ", nm),
    formula_args(intercept_formulas[[nm]], dat_factor, dat_factor$cl),
    ref_intercept
  )
}

# A factor-valued data-frame control must obey the same no-intercept coding as
# the formula path: all levels are retained.  Reference coding without an
# intercept omits a real nuisance direction and changes the estimator.
ctrl_factor <- factor(rep(c("c1", "c2", "c3", "c4"), length.out = n))
ctrl_full <- vapply(levels(ctrl_factor),
                    function(lev) as.numeric(ctrl_factor == lev), numeric(n))
colnames(ctrl_full) <- paste0("ctrl_", levels(ctrl_factor))
dat_ctrl <- data.frame(y = y, x = x, z1 = z1, z2 = z2, cl = cl,
                       ctrl_factor = ctrl_factor)
expect_equivalent(
  "[F/G] no-intercept data-frame factor controls equal full dummies",
  list(y = y, x = x, z = Z, cluster = cl,
       controls = data.frame(ctrl_factor = ctrl_factor), intercept = FALSE),
  list(y = y, x = x, z = Z, cluster = cl,
       controls = ctrl_full, intercept = FALSE)
)
expect_equivalent(
  "[F/G/J] no-intercept data-frame factor controls equal formula controls",
  list(y = y, x = x, z = Z, cluster = cl,
       controls = data.frame(ctrl_factor = ctrl_factor), intercept = FALSE),
  formula_args(y ~ 0 + ctrl_factor | x ~ z1 + z2, dat_ctrl, dat_ctrl$cl)
)

# The object classes that report k must also expose the intended dummy count.
for (nm in c("cjive", "cjar", "cjscore", "iv_infer")) {
  got_full <- run_call(nm, list(y = y, x = x, z = inst, cluster = cl,
                                intercept = FALSE))
  got_ref <- run_call(nm, list(y = y, x = x, z = inst, cluster = cl,
                               intercept = TRUE))
  record(is.null(got_full$error) && got_full$value$k == nlevels(inst),
         paste0("[F] no-intercept factor reports all levels [", nm, "]"),
         if (is.null(got_full$error)) paste0("k=", got_full$value$k)
         else got_full$error)
  record(is.null(got_ref$error) && got_ref$value$k == nlevels(inst) - 1L,
         paste0("[G] intercept factor reports reference coding [", nm, "]"),
         if (is.null(got_ref$error)) paste0("k=", got_ref$value$k)
         else got_ref$error)
}

# === H: numeric type contracts =============================================
bad_vectors <- list(
  `factor y` = list(component = "y", value = factor(ifelse(y > 0, "hi", "lo"))),
  `character y` = list(component = "y", value = as.character(round(y, 6))),
  `factor x` = list(component = "x", value = factor(ifelse(x > 0, "hi", "lo"))),
  `character x` = list(component = "x", value = as.character(round(x, 6))),
  `factor weights` = list(component = "weights",
                          value = factor(rep(c("1", "2"), length.out = n))),
  `character weights` = list(component = "weights", value = as.character(wt))
)
for (nm in names(bad_vectors)) {
  spec <- bad_vectors[[nm]]
  args <- list(y = y, x = x, z = Z, cluster = cl)
  args[[spec$component]] <- spec$value
  expect_errors(
    paste0("[H] vector interface rejects ", nm), args,
    component = paste0("`", spec$component, "`|", spec$component),
    diagnosis = "numeric|number|type"
  )
}

dat_types <- data.frame(
  y = y, x = x, z1 = z1, z2 = z2, cl = cl,
  yf = factor(ifelse(y > 0, "hi", "lo")),
  xf = factor(ifelse(x > 0, "hi", "lo")),
  wf = factor(rep(c("1", "2"), length.out = n))
)
formula_type_cases <- list(
  `factor y` = list(formula = yf ~ x | z1 + z2, extra = list(), component = "y"),
  `factor x` = list(formula = y ~ xf | z1 + z2, extra = list(), component = "x"),
  `factor weights` = list(formula = y ~ x | z1 + z2,
                          extra = list(weights = dat_types$wf),
                          component = "weights")
)
for (case_nm in names(formula_type_cases)) {
  spec <- formula_type_cases[[case_nm]]
  args <- c(formula_args(spec$formula, dat_types, dat_types$cl), spec$extra)
  expect_errors(
    paste0("[H] formula interface rejects ", case_nm), args,
    component = paste0("`", spec$component, "`|", spec$component),
    diagnosis = "numeric|number|type"
  )
}

# Positivity and finiteness are part of the numeric weights contract.
bad_weights <- list(zero = replace(wt, 1L, 0),
                    negative = replace(wt, 1L, -1),
                    infinite = replace(wt, 1L, Inf))
for (nm in names(bad_weights)) {
  expect_errors(
    paste0("[H] all APIs reject ", nm, " weights"),
    list(y = y, x = x, z = Z, cluster = cl, weights = bad_weights[[nm]]),
    component = "weight",
    diagnosis = "positive|finite|numeric"
  )
}

# === I: scalar levels, nulls and logical flags =============================
bad_levels <- list(`0` = 0, `1` = 1, `-0.1` = -0.1, `1.1` = 1.1,
                   `NA` = NA_real_, `Inf` = Inf,
                   `length two` = c(0.9, 0.95))
for (nm in names(bad_levels)) {
  expect_errors(
    paste0("[I] invalid level rejected: ", nm),
    c(list(y = y, x = x, z = Z, cluster = cl),
      list(level = bad_levels[[nm]])),
    component = "level",
    diagnosis = "single|scalar|finite|between|numeric|number"
  )
}
expect_accepts(
  "[I] level = 0.4 is inside the required open interval (0, 1)",
  list(y = y, x = x, z = Z, cluster = cl, level = 0.4)
)

bad_beta0 <- list(`NA` = NA_real_, `Inf` = Inf, `length two` = c(0, 1),
                  character = "0", factor = factor("0"))
for (nm in names(bad_beta0)) {
  expect_errors(
    paste0("[I] invalid beta0 rejected: ", nm),
    c(list(y = y, x = x, z = Z, cluster = cl),
      list(beta0 = bad_beta0[[nm]])),
    which_api = c("cjar", "cjscore", "iv_infer"),
    component = "beta0",
    diagnosis = "single|scalar|finite|numeric|number"
  )
}

bad_intercepts <- list(`NA` = NA, `length two` = c(TRUE, FALSE), numeric = 1)
for (nm in names(bad_intercepts)) {
  expect_errors(
    paste0("[I] invalid intercept flag rejected: ", nm),
    c(list(y = y, x = x, z = Z, cluster = cl),
      list(intercept = bad_intercepts[[nm]])),
    component = "intercept",
    diagnosis = "logical|TRUE|FALSE|single|scalar"
  )
}

# === J: the formula language is strict and explicit ========================
dat_formula <- data.frame(y = y, x = x, x2 = x + rnorm(n), z1 = z1, z2 = z2,
                          w1 = w1, w2 = w2, fe = fe, cl1 = cl,
                          cl2 = rep(seq_len(12L), length.out = n))
bad_formulas <- list(
  `unsupported subtraction in endogenous part` = y ~ x - x2 | z1,
  `unsupported subtraction in instrument part` = y ~ x | z1 - z2,
  `multiple endogenous regressors (legacy)` = y ~ x + x2 | z1,
  `multiple endogenous regressors (fixest-style)` = y ~ 1 | x + x2 ~ z1,
  `missing instrument section` = y ~ x,
  `too many legacy sections` = y ~ x | z1 | fe | w1,
  `missing fixest-style exogenous section` = y ~ x ~ z1
)
for (nm in names(bad_formulas)) {
  expect_errors(
    paste0("[J] formula is rejected: ", nm),
    formula_args(bad_formulas[[nm]], dat_formula, dat_formula$cl1),
    component = "formula|fixest|exogenous|endogenous|instrument|part|section|operator|support|parse|`[xz]`",
    diagnosis = "support|exactly one|multiple|too many|must|need|form|operator|parse|section",
    same = FALSE
  )
}

multi_cluster <- ~ cl1 + cl2
expect_errors(
  "[J] multi-term cluster formula is rejected (multiway clustering unsupported)",
  formula_args(y ~ x | z1, dat_formula, multi_cluster),
  component = "cluster|multiway",
  diagnosis = "one|single|term|multiway|unsupported"
)

for (nm in names(public)) {
  ans <- tryCatch(
    list(value = suppressWarnings(public[[nm]](
      y ~ x | z1, data = dat_formula, cluster = cl1 + cl2
    )), error = NULL),
    error = function(e) list(value = NULL, error = conditionMessage(e))
  )
  record(!is.null(ans$error) &&
           grepl("cluster|arithmetic|expression|combine|multiway", ans$error,
                 ignore.case = TRUE),
         paste0("[J] arithmetic cluster expression is rejected [", nm, "]"),
         if (is.null(ans$error)) "ACCEPT" else ans$error)
}

# === K: no public fit silently ignores dots ================================
for (nm in names(public)) {
  args <- c(list(y = y, x = x, z = Z, cluster = cl),
            list(calbration = "normal"))
  ans <- run_call(nm, args)
  record(!is.null(ans$error) &&
           grepl("unused|unknown|argument|calbration|\\.\\.\\.", ans$error,
                 ignore.case = TRUE),
         paste0("[K] misspelled named argument rejected [", nm, "]"),
         if (is.null(ans$error)) "ACCEPT" else ans$error)
}

partial_misspellings <- list(
  cjive = list(leve = 0.9),
  iv_compare = list(inferenc = "t"),
  cjar = list(calibratio = "normal"),
  cjscore = list(varia = "plain"),
  iv_infer = list(inferenc = "t")
)
for (nm in names(partial_misspellings)) {
  args <- c(list(y = y, x = x, z = Z, cluster = cl),
            partial_misspellings[[nm]])
  ans <- run_call(nm, args)
  record(!is.null(ans$error) &&
           grepl("unused|unknown|argument|misspell|partial", ans$error,
                 ignore.case = TRUE),
         paste0("[K] partially matched misspelling is rejected [", nm, "]"),
         if (is.null(ans$error)) "ACCEPT" else ans$error)
}

# Exact-name validation belongs to the package's built-in default/formula
# methods.  It must not prevent third-party S3 methods from defining their own
# named arguments on these exported generics.
cjive.input_probe <- function(y, special, ...) special
iv_compare.input_probe <- function(y, special, ...) special
cjar.input_probe <- function(y, special, ...) special
cjscore.input_probe <- function(y, special, ...) special
iv_infer.input_probe <- function(y, special, ...) special
probe <- structure(1, class = "input_probe")
for (nm in names(public)) {
  ans <- tryCatch(public[[nm]](probe, special = 42), error = conditionMessage)
  record(identical(ans, 42),
         paste0("[K] third-party S3 arguments still dispatch [", nm, "]"),
         if (identical(ans, 42)) NULL else as.character(ans))
}

for (nm in names(public)) {
  args <- switch(nm,
    cjive = list(y = y, x = x, z = Z, cluster = cl, controls = NULL,
                 fixed_effects = NULL, weights = NULL, intercept = TRUE,
                 level = 0.95, method = "dense", inference = "asymptotic"),
    iv_compare = list(y = y, x = x, z = Z, cluster = cl, controls = NULL,
                      fixed_effects = NULL, weights = NULL, intercept = TRUE,
                      level = 0.95, inference = "asymptotic"),
    cjar = list(y = y, x = x, z = Z, cluster = cl, controls = NULL,
                fixed_effects = NULL, weights = NULL, intercept = TRUE,
                level = 0.95, beta0 = 0, calibration = "chisq",
                variance = "plain"),
    cjscore = list(y = y, x = x, z = Z, cluster = cl, controls = NULL,
                   fixed_effects = NULL, weights = NULL, intercept = TRUE,
                   level = 0.95, beta0 = 0, variance = "plain"),
    iv_infer = list(y = y, x = x, z = Z, cluster = cl, controls = NULL,
                    fixed_effects = NULL, weights = NULL, intercept = TRUE,
                    level = 0.95, beta0 = 0, calibration = "chisq",
                    inference = "asymptotic",
                    tests = c("cjar", "cjscore"),
                    variance = "plain")
  )
  args <- c(args, list(123))
  ans <- run_call(nm, args)
  record(!is.null(ans$error) &&
           grepl("unused|unnamed|argument|\\.\\.\\.", ans$error,
                 ignore.case = TRUE),
         paste0("[K] unnamed dots entry rejected [", nm, "]"),
         if (is.null(ans$error)) "ACCEPT" else ans$error)
}

for (nm in names(public)) {
  args <- formula_args(y ~ x | z1 + z2, dat_formula, dat_formula$cl1,
                       calbration = "normal")
  ans <- run_call(nm, args)
  record(!is.null(ans$error) &&
           grepl("unused|unknown|argument|calbration|\\.\\.\\.", ans$error,
                 ignore.case = TRUE),
         paste0("[K] formula method rejects misspelled argument [", nm, "]"),
         if (is.null(ans$error)) "ACCEPT" else ans$error)
}

cf_vector <- run_call(
  "iv_infer",
  list(y = y, x = x, z = Z, cluster = cl, variance = "crossfit",
       tests = c("cjar", "cjscore"))
)
record(is.null(cf_vector$error) &&
         identical(cf_vector$value$variance, "crossfit") &&
         identical(cf_vector$value$cjar$variance_estimator, "crossfit") &&
         identical(cf_vector$value$cjscore$variance_estimator, "crossfit"),
       "[K] iv_infer accepts explicit variance = 'crossfit'",
       cf_vector$error)
cf_formula <- run_call(
  "iv_infer",
  formula_args(y ~ x | z1 + z2, dat_formula, dat_formula$cl1,
               variance = "crossfit", tests = c("cjar", "cjscore"))
)
record(is.null(cf_formula$error) &&
         identical(cf_formula$value$variance, "crossfit") &&
         identical(cf_formula$value$cjar$variance_estimator, "crossfit") &&
         identical(cf_formula$value$cjscore$variance_estimator, "crossfit"),
       "[K] iv_infer.formula passes variance = 'crossfit'",
       cf_formula$error)

# === Estimator denominator: zero first stage, but nonconstant x =============
# In every cluster z'x is exactly zero.  Hence every leave-cluster-out first
# stage is zero although x has variation.  Point estimators must stop before
# dividing by zero; weak-ID-robust tests must remain available.
G0 <- 24L
cl0 <- rep(seq_len(G0), each = 4L)
z0 <- rep(c(-3, -1, 1, 3), G0)
x0 <- rep(c(1, -1, -1, 1), G0)
set.seed(20260718)
y0 <- rnorm(length(x0))

expect_errors(
  "[denominator] point-estimation APIs reject an exactly zero IV denominator",
  list(y = y0, x = x0, z = z0, cluster = cl0, intercept = FALSE),
  which_api = c("cjive", "iv_compare", "iv_infer"),
  component = "denominator|first.stage|identified|instrument",
  diagnosis = "zero|numerical|unidentified|undefined|variation|weak",
  same = FALSE
)
expect_accepts(
  "[denominator] weak-ID-robust CJAR/CJS remain available at zero first stage",
  list(y = y0, x = x0, z = z0, cluster = cl0, intercept = FALSE),
  which_api = c("cjar", "cjscore")
)

# The constructed instrument scales with x, so its raw product with x can
# underflow or overflow even though the coefficient and SE remain representable.
# Point inference must use an equivalent scaled calculation at that boundary.
point_fields <- function(name, value) {
  switch(name,
    cjive = c(value$coefficient, value$se),
    iv_compare = c(value$coefficient, value$se),
    iv_infer = c(value$cjive$coefficient, value$cjive$se)
  )
}
point_apis <- c("cjive", "iv_compare", "iv_infer")
point_base <- lapply(point_apis, function(nm) {
  # A CJIVE-only panel: this block audits point inference, and the jackknife
  # coefficient polynomials are deliberately unrepresentable at the extreme
  # scales probed below (their own guard is tested elsewhere).
  args <- if (nm == "iv_infer") c(base_args, list(tests = NULL)) else base_args
  run_call(nm, args)$value
})
names(point_base) <- point_apis
for (s_name in c("tiny", "huge")) {
  s <- if (s_name == "tiny") 2^-540 else 2^510
  for (nm in point_apis) {
    args <- within(base_args, x <- x * s)
    if (nm == "iv_infer") args["tests"] <- list(NULL)
    ans <- run_call(nm, args)
    good <- is.null(ans$error)
    if (good) {
      got <- point_fields(nm, ans$value)
      ref <- point_fields(nm, point_base[[nm]])
      good <- all(is.finite(got)) &&
        isTRUE(all.equal(got * abs(s), ref, tolerance = 1e-10))
    }
    record(good,
           paste0("[denominator/E] point inference is stable at ", s_name,
                  " x scale [", nm, "]"),
           if (!is.null(ans$error)) ans$error else "non-finite or rescaled mismatch")
  }
}

# The same unit-normalized sandwich must be equivariant to outcome units and
# invariant to a common positive scaling of precision weights.
for (s_name in c("tiny", "huge")) {
  s <- if (s_name == "tiny") 2^-540 else 2^510
  for (nm in point_apis) {
    args <- within(base_args, y <- y * s)
    if (nm == "iv_infer") args["tests"] <- list(NULL)
    ans <- run_call(nm, args)
    got <- if (is.null(ans$error)) point_fields(nm, ans$value) else NULL
    ref <- point_fields(nm, point_base[[nm]])
    good <- !is.null(got) && all(is.finite(got)) &&
      isTRUE(all.equal(got / abs(s), ref, tolerance = 1e-10))
    record(good,
           paste0("[sandwich/E] point inference is equivariant at ", s_name,
                  " y scale [", nm, "]"),
           if (!is.null(ans$error)) ans$error else "non-finite or rescaled mismatch")
  }
}

weighted_base_args <- c(base_args, list(weights = wt))
weighted_base <- lapply(point_apis, function(nm) {
  args <- weighted_base_args
  if (nm == "iv_infer") args["tests"] <- list(NULL)
  run_call(nm, args)$value
})
names(weighted_base) <- point_apis
for (s_name in c("tiny", "huge")) {
  s <- if (s_name == "tiny") 2^-900 else 2^900
  for (nm in point_apis) {
    args <- weighted_base_args
    args$weights <- args$weights * s
    if (nm == "iv_infer") args["tests"] <- list(NULL)
    ans <- run_call(nm, args)
    got <- if (is.null(ans$error)) point_fields(nm, ans$value) else NULL
    ref <- point_fields(nm, weighted_base[[nm]])
    good <- !is.null(got) && all(is.finite(got)) &&
      isTRUE(all.equal(got, ref, tolerance = 1e-10))
    record(good,
           paste0("[sandwich/E] point inference is invariant to ", s_name,
                  " common weight scale [", nm, "]"),
           if (!is.null(ans$error)) ans$error else "non-finite or scale mismatch")
  }
}

# Explicit zero-variance semantics: a nonzero estimate has an infinite test
# statistic, while coefficient == SE == 0 is an undefined 0/0 test (NA, not
# an accidental NaN).
iv_inf_internal <- getFromNamespace(".iv_inference", "clusterIV")
cl_zv <- rep(1:4, each = 2L)
p_zv <- rep(c(1, 2), 4L)
fit_zv_nonzero <- iv_inf_internal(p_zv, p_zv, 2 * p_zv, cl_zv)
record(fit_zv_nonzero$coefficient == 2 && fit_zv_nonzero$se == 0 &&
         is.infinite(fit_zv_nonzero$statistic) &&
         fit_zv_nonzero$statistic > 0 && fit_zv_nonzero$p.value == 0,
       "[sandwich] zero SE with nonzero coefficient gives an infinite statistic")
p_z0 <- rep(1, 8L)
y_z0 <- rep(c(1, -1), 4L)
fit_zv_zero <- iv_inf_internal(p_z0, p_z0, y_z0, cl_zv)
record(fit_zv_zero$coefficient == 0 && fit_zv_zero$se == 0 &&
         is.na(fit_zv_zero$statistic) && is.na(fit_zv_zero$p.value),
       "[sandwich] coefficient == SE == 0 is reported as an undefined test")

# === Algebraically exact first-stage fit: effective F convention ============
set.seed(20260719)
nF <- 120L
clF <- rep(seq_len(24L), each = 5L)
ZF <- matrix(rnorm(nF * 2L), nF, 2L)
xF <- drop(ZF %*% c(2, -1))
yF <- 0.5 * xF + rnorm(nF)
for (nm in c("cjive", "cjar", "cjscore", "iv_infer")) {
  ans <- run_call(nm, list(y = yF, x = xF, z = ZF, cluster = clF,
                           intercept = FALSE))
  good <- is.null(ans$error) && ans$value$F_eff > 1e20 &&
    ((is.finite(ans$value$F_eff) && is.finite(ans$value$K_eff) &&
      is.finite(ans$value$F_eff_crit)) ||
     (is.infinite(ans$value$F_eff) && is.na(ans$value$K_eff) &&
      is.na(ans$value$F_eff_crit)))
  record(good,
         paste0("[effective F] algebraically exact fit is overwhelmingly strong without a tolerance-forced Inf [",
                nm, "]"),
         if (!is.null(ans$error)) ans$error
         else paste0("F_eff=", format(ans$value$F_eff, digits = 8),
                     ", K_eff=", ans$value$K_eff,
                     ", crit=", ans$value$F_eff_crit))
}

eff_f_internal <- getFromNamespace(".eff_f", "clusterIV")
eff_00 <- eff_f_internal(matrix(0, 4L, 1L), 0, rep(0, 4L),
                         cluster = seq_len(4L))
record(all(is.na(unlist(eff_00))),
       "[effective F] genuine zero-over-zero case is NA throughout")
eff_inf <- eff_f_internal(matrix(1, 4L, 1L), 1, rep(0, 4L),
                          cluster = seq_len(4L))
record(is.infinite(eff_inf$F_eff) && eff_inf$F_eff > 0 &&
         is.na(eff_inf$K_eff) && is.na(eff_inf$F_eff_crit),
       "[effective F] nonzero-over-exact-zero case is Inf with undefined auxiliaries")

if (length(failures)) {
  stop(paste0("Input-space oracle failures (", length(failures), "):\n - ",
              paste(failures, collapse = "\n - ")),
       call. = FALSE)
}

cat("\nAll input-space and strict-validation tests passed.\n")
