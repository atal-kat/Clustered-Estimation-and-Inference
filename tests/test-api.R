# API-coherence tests (base R, no testthat), in the style of test-cjive.R:
# each block stops with an informative message on failure, so a failure fails
# `R CMD check`. Covers the WP2 gates: formula-dialect equivalence (G2), the
# grouping-factor fix (G3), subset/na.action (G4), and iv_infer() consistency
# with the standalone calls plus single-fire advisories (G5).

library(clusterIV)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}
errs_with <- function(expr, pattern) {
  msg <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  is.character(msg) && grepl(pattern, msg)
}

# Compare two fitted objects field by field at `tol`, skipping `call`.
# Numeric fields must match elementwise, treating NA == NA and matching
# infinities as equal; everything else must be identical().
fields_equal <- function(a, b, tol = 1e-12, skip = "call") {
  nms <- union(names(a), names(b))
  for (nm in setdiff(nms, skip)) {
    va <- a[[nm]]; vb <- b[[nm]]
    if (is.numeric(va) && is.numeric(vb)) {
      if (length(va) != length(vb)) return(nm)
      same <- (is.na(va) & is.na(vb)) |
        (is.infinite(va) & is.infinite(vb) & sign(va) == sign(vb)) |
        (!is.na(va) & !is.na(vb) & is.finite(va) & is.finite(vb) &
           abs(va - vb) < tol)
      if (!all(same)) return(nm)
    } else if (!identical(va, vb)) {
      return(nm)
    }
  }
  TRUE
}
ok_fields <- function(a, b, msg, tol = 1e-12, skip = "call") {
  r <- fields_equal(a, b, tol = tol, skip = skip)
  if (!isTRUE(r)) stop("FAILED: ", msg, " (field: ", r, ")", call. = FALSE)
  cat("PASS:", msg, "\n")
}

quiet <- function(expr) suppressWarnings(expr)

# --- a design with controls, fixed effects and a judge factor ---------------
set.seed(2026)
G <- 25L
sizes <- sample(5:12, G, replace = TRUE)
n <- sum(sizes)
cl <- rep(seq_len(G), sizes)
judge <- factor(sample(paste0("j", 1:7), n, replace = TRUE))
fe1 <- factor(sample(letters[1:9], n, replace = TRUE))
w1 <- rnorm(n); w2 <- rnorm(n)
u <- rnorm(G)[cl]
x <- 0.7 * as.numeric(judge) + 0.3 * w1 + u + rnorm(n)
y <- 0.5 * x + 0.4 * w1 - 0.2 * w2 + u + rnorm(n)
wt <- runif(n, 0.5, 2)
dat <- data.frame(y = y, x = x, judge = judge, fe1 = fe1,
                  w1 = w1, w2 = w2, cl = cl, wt = wt)

# === G2: dialect equivalence =================================================
# Vector interface, legacy formula (controls via the argument) and fixest
# formula (controls inside the formula) must agree on every returned field.
for (fn in list(list(f = cjive, tag = "cjive"),
                list(f = cjar, tag = "cjar"),
                list(f = cjscore, tag = "cjscore"),
                list(f = iv_infer, tag = "iv_infer"))) {
  f <- fn$f
  vec <- quiet(f(y, x, judge, cluster = cl, controls = cbind(w1, w2),
                 fixed_effects = fe1))
  leg <- quiet(f(y ~ x | judge | fe1, data = dat, cluster = cl,
                 controls = ~ w1 + w2))
  fxs <- quiet(f(y ~ w1 + w2 | fe1 | x ~ judge, data = dat, cluster = cl))
  if (fn$tag == "iv_infer") {
    for (part in c("cjive", "cjar", "cjscore")) {
      ok_fields(vec[[part]], leg[[part]],
                sprintf("[G2 %s$%s] legacy formula == vector (1e-12, all fields)",
                        fn$tag, part))
      ok_fields(vec[[part]], fxs[[part]],
                sprintf("[G2 %s$%s] fixest formula == vector (1e-12, all fields)",
                        fn$tag, part))
    }
    ok_fields(vec, leg, "[G2 iv_infer] legacy formula == vector (shared fields)",
              skip = c("call", "cjive", "cjar", "cjscore"))
    ok_fields(vec, fxs, "[G2 iv_infer] fixest formula == vector (shared fields)",
              skip = c("call", "cjive", "cjar", "cjscore"))
  } else {
    ok_fields(vec, leg,
              sprintf("[G2 %s] legacy formula == vector (1e-12, all fields)", fn$tag))
    ok_fields(vec, fxs,
              sprintf("[G2 %s] fixest formula == vector (1e-12, all fields)", fn$tag))
  }
}

# iv_compare through the same three routes (a data frame, compared column-wise).
tab_v <- quiet(iv_compare(y, x, judge, cluster = cl, controls = cbind(w1, w2),
                          fixed_effects = fe1))
tab_l <- quiet(iv_compare(y ~ x | judge | fe1, data = dat, cluster = cl,
                          controls = ~ w1 + w2))
tab_f <- quiet(iv_compare(y ~ w1 + w2 | fe1 | x ~ judge, data = dat, cluster = cl))
num <- vapply(tab_v, is.numeric, logical(1))
ok(identical(tab_v$estimator, tab_l$estimator) &&
   identical(tab_v$estimator, tab_f$estimator) &&
   max(abs(as.matrix(tab_v[, num]) - as.matrix(tab_l[, num]))) < 1e-12 &&
   max(abs(as.matrix(tab_v[, num]) - as.matrix(tab_f[, num]))) < 1e-12,
   "[G2 iv_compare] all three routes agree (1e-12)")

# The fixest exog codes: 1 = intercept only, 0 = no intercept.
fx1 <- quiet(cjive(y ~ 1 | x ~ judge, data = dat, cluster = cl))
v1 <- quiet(cjive(y, x, judge, cluster = cl))
ok_fields(v1, fx1, "[G2] fixest exog = 1 == intercept-only vector call")
fx0 <- quiet(cjive(y ~ 0 | x ~ judge, data = dat, cluster = cl))
v0 <- quiet(cjive(y, x, judge, cluster = cl, intercept = FALSE))
ok_fields(v0, fx0, "[G2] fixest exog = 0 == intercept = FALSE vector call")

# Every FE span contains the intercept.  A grouping instrument must therefore
# remain reference-coded when FE are present even if intercept = FALSE; full
# coding would be singular after absorption.  Pin the shared vector engine,
# the formula route, every fitted public API, and the raw nuisance count.
Zref_fe0 <- model.matrix(~ judge)[, -1L, drop = FALSE]
for (fn in list(list(f = cjive, tag = "cjive"),
                list(f = cjar, tag = "cjar"),
                list(f = cjscore, tag = "cjscore"))) {
  fac <- quiet(fn$f(y, x, judge, cluster = cl, fixed_effects = fe1,
                    intercept = FALSE))
  ref <- quiet(fn$f(y, x, Zref_fe0, cluster = cl, fixed_effects = fe1,
                    intercept = FALSE))
  ok_fields(fac, ref,
            sprintf("[G2 %s] FE + intercept=FALSE factor z uses reference coding",
                    fn$tag))
  ok(fac$k == ncol(Zref_fe0) && fac$k_controls == nlevels(fe1),
     sprintf("[G2 %s] FE factor-z k and nuisance count are non-redundant",
             fn$tag))
}

wf_fac <- quiet(iv_infer(y, x, judge, cluster = cl, fixed_effects = fe1,
                         intercept = FALSE))
wf_ref <- quiet(iv_infer(y, x, Zref_fe0, cluster = cl, fixed_effects = fe1,
                         intercept = FALSE))
for (part in c("cjive", "cjar", "cjscore")) {
  ok_fields(wf_fac[[part]], wf_ref[[part]],
            sprintf("[G2 iv_infer$%s] FE factor-z reference coding", part))
}
ok(wf_fac$k == ncol(Zref_fe0) && wf_fac$k_controls == nlevels(fe1),
   "[G2 iv_infer] FE factor-z k and nuisance count are non-redundant")

tab_fac <- quiet(iv_compare(y, x, judge, cluster = cl,
                            fixed_effects = fe1, intercept = FALSE))
tab_ref <- quiet(iv_compare(y, x, Zref_fe0, cluster = cl,
                            fixed_effects = fe1, intercept = FALSE))
num_fe0 <- vapply(tab_fac, is.numeric, logical(1))
ok(identical(tab_fac$estimator, tab_ref$estimator) &&
     max(abs(as.matrix(tab_fac[, num_fe0]) -
             as.matrix(tab_ref[, num_fe0]))) < 1e-12,
   "[G2 iv_compare] FE + intercept=FALSE factor z uses reference coding")

fx_fe0 <- quiet(cjive(y ~ 0 | fe1 | x ~ judge, data = dat, cluster = cl))
ok_fields(fx_fe0,
          quiet(cjive(y, x, judge, cluster = cl, fixed_effects = fe1,
                      intercept = FALSE)),
          "[G2 formula] FE + intercept=FALSE factor z equals vector route")

# === G3: the grouping factor survives the formula path ======================
lm_vec <- cjive(y, x, judge, cluster = cl, method = "leaveout_mean")
lm_leg <- cjive(y ~ x | judge, data = dat, cluster = cl,
                method = "leaveout_mean")
lm_fxs <- cjive(y ~ 1 | x ~ judge, data = dat, cluster = cl,
                method = "leaveout_mean")
ok_fields(lm_vec, lm_leg, "[G3] leaveout_mean via legacy formula == vector")
ok_fields(lm_vec, lm_fxs, "[G3] leaveout_mean via fixest formula == vector")

# A judge confined to one cluster remains legal for the dense encoded-dummy
# path when the actual leave-out Gram is full rank, but the group-mean shortcut
# is undefined.  That stricter requirement must survive both formula dialects.
jud_bad <- as.character(judge)
jud_bad[cl == 1L] <- "solo"
jud_bad[cl != 1L & jud_bad == "solo"] <- "j1"
dat_bad <- transform(dat, judge = factor(jud_bad))
ok(errs_with(cjive(y ~ x | judge, data = dat_bad, cluster = cl,
                   method = "leaveout_mean"),
             "entirely in one cluster"),
   "[G3] single-cluster judge: leaveout_mean error via legacy formula")
ok(errs_with(cjive(y ~ 1 | x ~ judge, data = dat_bad, cluster = cl,
                   method = "leaveout_mean"),
             "entirely in one cluster"),
   "[G3] single-cluster judge: leaveout_mean error via fixest formula")

# === G4: subset and na.action ================================================
# NAs scattered across y, a control, the cluster id and the weight vector.
dat_na <- dat
dat_na$y[c(4L, 40L)] <- NA
dat_na$w1[15L] <- NA
dat_na$cl[27L] <- NA
dat_na$wt[c(33L, 90L)] <- NA
keep <- stats::complete.cases(dat_na$y, dat_na$w1, dat_na$cl, dat_na$wt)

for (fn in list(list(f = cjive, tag = "cjive"),
                list(f = cjar, tag = "cjar"),
                list(f = cjscore, tag = "cjscore"))) {
  f <- fn$f
  fit_na <- quiet(f(y ~ x | judge, data = dat_na, cluster = cl,
                    controls = ~ w1, weights = wt))
  fit_man <- quiet(f(dat_na$y[keep], x[keep], droplevels(judge[keep]),
                     cluster = dat$cl[keep],
                     controls = cbind(w1 = dat_na$w1[keep]),
                     weights = dat$wt[keep]))
  ok_fields(fit_na, fit_man,
            sprintf("[G4 %s] na.omit formula == vector on filtered data (1e-12)",
                    fn$tag),
            skip = c("call", "n_dropped"))
  ok(fit_na$n_dropped == sum(!keep),
     sprintf("[G4 %s] n_dropped == %d", fn$tag, sum(!keep)))
  out <- capture.output(print(fit_na))
  ok(any(grepl("dropped due to missing values", out)),
     sprintf("[G4 %s] the drop is reported by print()", fn$tag))
}

# na.fail stops; na.pass leaves the NAs for the input validation.
ok(errs_with(cjive(y ~ x | judge, data = dat_na, cluster = cl,
                   na.action = stats::na.fail), "na.fail"),
   "[G4] na.action = na.fail stops on missing values")
ok(errs_with(cjive(y ~ x | judge, data = dat_na, cluster = cl,
                   na.action = stats::na.pass), "missing"),
   "[G4] na.action = na.pass leaves NAs to the input validation")

# subset: formula subset == manual row selection.
sub <- cl <= 20L
fit_sub <- cjive(y ~ x | judge, data = dat, cluster = cl, subset = cl <= 20L)
fit_man <- cjive(y[sub], x[sub], droplevels(judge[sub]), cluster = cl[sub])
ok_fields(fit_sub, fit_man, "[G4] subset == manual row selection (1e-12)")
ok(fit_sub$n_dropped == 0L, "[G4] subset alone drops nothing via na.action")

# Dropping rows can empty a judge level (must be droplevel()ed away, not
# turned into an all-zero dummy) and can confine a surviving judge to a
# single cluster.  The dense path is governed by its encoded leave-out Gram;
# the stricter leaveout_mean support check must use the post-filter design.
dat_lvl <- dat
lvl_gone <- dat_lvl$judge == "j7"                # j7 exists only in NA rows
dat_lvl$y[lvl_gone] <- NA
fit_lvl <- quiet(cjive(y ~ x | judge, data = dat_lvl, cluster = cl))
keep_lvl <- !lvl_gone
man_lvl <- quiet(cjive(dat$y[keep_lvl], x[keep_lvl],
                       droplevels(judge[keep_lvl]), cluster = cl[keep_lvl]))
ok_fields(fit_lvl, man_lvl,
          "[G4] emptied judge level is dropped, fit == manual (1e-12)",
          skip = c("call", "n_dropped"))
ok(fit_lvl$n_dropped == sum(lvl_gone) && fit_lvl$k == man_lvl$k,
   "[G4] emptied level: n_dropped and k reflect the post-filter design")

dat_conf <- dat
conf_gone <- dat_conf$judge == "j6" & dat_conf$cl != 1L   # confine j6 to cluster 1
dat_conf$y[conf_gone] <- NA
keep_conf <- !conf_gone
fit_conf <- quiet(cjive(y ~ x | judge, data = dat_conf, cluster = cl))
man_conf <- quiet(cjive(dat$y[keep_conf], x[keep_conf],
                        droplevels(judge[keep_conf]), cluster = cl[keep_conf]))
ok_fields(fit_conf, man_conf,
          "[G4] post-filter cluster-local level: dense formula == manual",
          skip = c("call", "n_dropped"))
ok(errs_with(quiet(cjive(y ~ x | judge, data = dat_conf, cluster = cl,
                         method = "leaveout_mean")),
             "entirely in one cluster"),
   "[G4] post-filter cluster-local level: leaveout_mean support error")

# === G5: iv_infer() consistency ==============================================
fit_all <- quiet(iv_infer(y, x, judge, cluster = cl, controls = cbind(w1, w2),
                          fixed_effects = fe1, beta0 = 0.4, level = 0.9,
                          calibration = "normal", inference = "t",
                          tests = c("cjar", "cjscore")))
sa_cjive <- quiet(cjive(y, x, judge, cluster = cl, controls = cbind(w1, w2),
                        fixed_effects = fe1, level = 0.9, inference = "t"))
sa_cjar <- quiet(cjar(y, x, judge, cluster = cl, controls = cbind(w1, w2),
                      fixed_effects = fe1, beta0 = 0.4, level = 0.9,
                      calibration = "normal"))
sa_cjs <- quiet(cjscore(y, x, judge, cluster = cl, controls = cbind(w1, w2),
                        fixed_effects = fe1, beta0 = 0.4, level = 0.9))
ok_fields(fit_all$cjive, sa_cjive, "[G5] iv_infer$cjive == standalone cjive (1e-12, all fields)")
ok_fields(fit_all$cjar, sa_cjar, "[G5] iv_infer$cjar == standalone cjar (1e-12, all fields)")
ok_fields(fit_all$cjscore, sa_cjs, "[G5] iv_infer$cjscore == standalone cjscore (1e-12, all fields)")
ok(fit_all$n == sa_cjive$n && fit_all$G == sa_cjive$G &&
   fit_all$k == sa_cjar$k && identical(fit_all$maxlev, sa_cjar$maxlev),
   "[G5] iv_infer shared diagnostics match the components")

# Existing methods keep working on the components.
ok(identical(coef(fit_all$cjive), coef(sa_cjive)) &&
   identical(confint(fit_all$cjar), confint(sa_cjar)) &&
   nobs(fit_all$cjscore) == n && nobs(fit_all) == n,
   "[G5] component methods (coef/confint/nobs) work unchanged")

# Each advisory fires exactly once per iv_infer() call. Design 1 triggers
# small-G and dominating-cluster; design 2 triggers the leverage advisory.
set.seed(7)
G5a <- 10L; n5a <- 80L
cl5a <- c(rep(1L, 30L), rep(2:G5a, length.out = n5a - 30L))
Z5a <- matrix(rnorm(n5a * 3L), n5a, 3L)
x5a <- drop(Z5a %*% c(1, -0.5, 0.4)) + rnorm(n5a)
y5a <- x5a + rnorm(n5a)
msgs <- character(0)
tmp <- withCallingHandlers(
  iv_infer(y5a, x5a, Z5a, cluster = cl5a),
  warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
  })
ok(sum(grepl("asymptotic in the number of clusters", msgs)) == 1L,
   "[G5] small-G advisory fires exactly once")
ok(sum(grepl("dominating cluster", msgs)) == 1L,
   "[G5] dominating-cluster advisory fires exactly once")

set.seed(88)
G5b <- 20L; ng5b <- 5L; n5b <- G5b * ng5b
cl5b <- rep(seq_len(G5b), each = ng5b)
# Nearly cluster-local instrument: dominant variation inside cluster 1, a
# whisper outside so the leave-cluster-out fit stays defined but the cluster
# leverage exceeds 0.99.
zloc <- rnorm(n5b) * 0.01
zloc[cl5b == 1L] <- c(1.5, -1.5, 1, -1, 0)
Z5b <- cbind(matrix(rnorm(n5b * 2L), n5b, 2L), zloc)
x5b <- drop(Z5b %*% c(1, -1, 0.5)) + rnorm(n5b)
y5b <- x5b + rnorm(n5b)
msgs <- character(0)
tmp <- withCallingHandlers(
  iv_infer(y5b, x5b, Z5b, cluster = cl5b),
  warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
  })
ok(sum(grepl("nearly spans the instrument space", msgs)) == 1L,
   "[G5] leverage advisory fires exactly once")

# The compact print renders the workflow summary.
out5 <- capture.output(print(fit_all))
ok(any(grepl("CJIVE/Wald (H0: beta = 0)", out5, fixed = TRUE)) &&
   any(grepl("CJAR (H0: beta = 0.4)", out5, fixed = TRUE)) &&
   any(grepl("CJS (H0: beta = 0.4)", out5, fixed = TRUE)) &&
   any(grepl("F_CJ", out5)) &&
   any(grepl("max within-cluster leverage", out5)),
   "[G5] print.iv_infer: estimate, set, CJS p-value, F_CJ, maxlev all render")

# === A8: confint() defaults to the fitted level ==============================
ar99 <- quiet(cjar(y, x, judge, cluster = cl, level = 0.99))
ok(identical(confint(ar99), ar99$conf_set),
   "[A8] confint.cjar defaults to the fitted level (0.99), not 0.95")
sc99 <- quiet(cjscore(y, x, judge, cluster = cl, level = 0.99))
ok(identical(confint(sc99), sc99$conf_set),
   "[A8] confint.cjscore defaults to the fitted level")
est99 <- quiet(cjive(y, x, judge, cluster = cl, level = 0.99))
ok(identical(unname(confint(est99)),
             unname(cbind(est99$conf.low, est99$conf.high))),
   "[A8] confint.cjive defaults to the fitted level")
ok(errs_with(confint(est99, level = 1.2), "strictly between"),
   "[A8] confint.cjive validates level")

# === W4: iv_infer() `tests` selector =========================================
# The panel accepts exactly "cjar" and "cjscore".  Selection must trim
# computation, every method must cope with each remaining subset, and the
# removed plain-panel selectors must fail with the unknown-selector error.

# `tests` trims computation; print and tidiers cope with absent components
fit_trim <- quiet(iv_infer(y, x, judge, cluster = cl, tests = "cjar"))
ok(is.null(fit_trim$cjscore) && !is.null(fit_trim$cjar) &&
     identical(fit_trim$tests, "cjar"),
   "[W4] tests = 'cjar': deselected fields are NULL, selector stored")
out_trim <- capture.output(print(fit_trim))
ok(!any(grepl("CJS (H0", out_trim, fixed = TRUE)) &&
     any(grepl("CJAR (H0: beta = 0)", out_trim, fixed = TRUE)),
   "[W4] print.iv_infer skips deselected rows")
ok(nrow(clusterIV:::tidy.iv_infer(fit_trim)) == 2L,
   "[W4] tidy.iv_infer returns CJIVE + CJAR when cjscore is deselected")

fit_min <- quiet(iv_infer(y, x, judge, cluster = cl, tests = "cjscore"))
ok(is.null(fit_min$cjar) && !is.null(fit_min$cjscore) &&
     identical(fit_min$tests, "cjscore"),
   "[W4] tests = 'cjscore': only the CJS row is computed")
out_min <- capture.output(print(fit_min))
ok(!any(grepl("CJAR (H0", out_min, fixed = TRUE)) &&
     !any(grepl("F_CJ =", out_min)) &&
     any(grepl("CJS (H0: beta = 0)", out_min, fixed = TRUE)),
   "[W4] print.iv_infer renders a cjar-free panel")
ok(nrow(clusterIV:::tidy.iv_infer(fit_min)) == 2L,
   "[W4] tidy.iv_infer returns CJIVE + CJS without CJAR")
ok(is.data.frame(clusterIV:::glance.iv_infer(fit_min)),
   "[W4] glance.iv_infer supports a panel without CJAR")

# Selector resolution: canonical order, deduplication, defaults, and the
# rejection of the removed plain-panel values.
fit_dup <- quiet(iv_infer(y, x, judge, cluster = cl,
                          tests = c("cjscore", "cjar", "cjscore")))
ok(identical(fit_dup$tests, c("cjar", "cjscore")),
   "[W4] selector values are deduplicated and stored in canonical order")
fit_default <- quiet(iv_infer(y, x, judge, cluster = cl))
ok(identical(fit_default$tests, c("cjar", "cjscore")),
   "[W4] omitted tests selects the recommended cjar + cjscore pair")
fit_none <- quiet(iv_infer(y, x, judge, cluster = cl, tests = NULL))
ok(is.null(fit_none$cjar) && is.null(fit_none$cjscore) &&
     identical(fit_none$tests, character(0)),
   "[W4] tests = NULL is a CJIVE-only panel")
ok(errs_with(plot(fit_none), "CJAR or CJS"),
   "[W4] plot.iv_infer requires an available plottable jackknife test")
ok(errs_with(iv_infer(y, x, judge, cluster = cl, tests = "wald"),
             "unknown"),
   "[W4] unknown tests value errors")
ok(errs_with(iv_infer(y, x, judge, cluster = cl, tests = "ar"),
             "unknown"),
   "[W4] removed plain selector 'ar' is an unknown tests value")
ok(errs_with(iv_infer(y, x, judge, cluster = cl, tests = "score"),
             "unknown"),
   "[W4] removed plain selector 'score' is an unknown tests value")

# k > G leaves the jackknife rows fully defined (that regime is what they
# are for).
set.seed(41)
G6 <- 5L; n6 <- 40L
cl6 <- rep(seq_len(G6), each = 8L)
Z6 <- matrix(rnorm(n6 * 6L), n6, 6L)
x6 <- drop(Z6 %*% rep(c(0.8, -0.5), 3L)) + rnorm(n6)
y6 <- x6 + rnorm(n6)
fit6 <- quiet(iv_infer(y6, x6, Z6, cluster = cl6))
ok(is.finite(fit6$cjar$statistic) && is.finite(fit6$cjscore$statistic),
   "[W4] k > G leaves the jackknife rows fully defined")

cat("\nAll API-coherence tests passed.\n")
