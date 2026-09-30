# Repair-item-4 oracle: publication-safe panel, tidier, print, and plot
# semantics.  Base R only.  Failures are collected so the pre-repair run
# records every unsafe behavior in one pass.

library(clusterIV)

.failures <- character(0)
check <- function(cond, msg, detail = NULL) {
  if (isTRUE(cond)) {
    cat("PASS:", msg, "\n")
  } else {
    line <- if (is.null(detail)) msg else paste0(msg, " [", detail, "]")
    .failures <<- c(.failures, line)
    cat("FAIL:", line, "\n")
  }
  invisible(cond)
}

capture_error <- function(expr) {
  tryCatch(list(value = expr, error = NULL),
           error = function(e) list(value = NULL, error = e))
}

same_num <- function(a, b, tol = 1e-12) {
  length(a) == length(b) &&
    all((is.na(a) & is.na(b)) |
          (is.infinite(a) & is.infinite(b) & sign(a) == sign(b)) |
          (is.finite(a) & is.finite(b) & abs(a - b) <= tol))
}

same_set <- function(a, b, tol = 1e-10) {
  identical(dim(a), dim(b)) &&
    all(is.infinite(a) == is.infinite(b)) &&
    same_num(a[is.finite(a)], b[is.finite(b)], tol)
}

quiet <- function(expr) suppressWarnings(expr)

# A small, well-conditioned design used throughout.
set.seed(2404)
G <- 24L
ng <- 5L
n <- G * ng
cl <- rep(seq_len(G), each = ng)
z1 <- rnorm(n)
z2 <- rnorm(n)
Z <- cbind(z1 = z1, z2 = z2)
u <- rnorm(G)[cl]
x <- 0.9 * z1 - 0.45 * z2 + u + rnorm(n)
y <- 0.65 * x + u + rnorm(n)
dat <- data.frame(y = y, x = x, z1 = z1, z2 = z2, cl = cl)

# -------------------------------------------------------------------------
# A. Every tidy row is one coherent inferential procedure.
# -------------------------------------------------------------------------
beta0 <- 1.25
panel <- quiet(iv_infer(y, x, Z, cluster = cl, beta0 = beta0,
                        tests = c("cjar", "cjscore")))
td <- capture_error(clusterIV:::tidy.iv_infer(panel))
needed <- c("term", "component", "procedure", "estimate", "std.error",
            "statistic", "p.value", "null.value", "conf.low", "conf.high",
            "shape", "n.components", "unbounded.left", "unbounded.right",
            "conf.set")
coherent <- is.null(td$error) && is.data.frame(td$value) &&
  nrow(td$value) == 3L && all(needed %in% names(td$value)) &&
  identical(td$value$component, c("cjive", "cjar", "cjscore"))
check(coherent,
      "tidy.iv_infer returns one explicitly labelled row per procedure",
      if (!is.null(td$error)) conditionMessage(td$error)
      else paste(names(td$value), collapse = ", "))

if (coherent) {
  wald <- td$value[td$value$component == "cjive", ]
  ar <- td$value[td$value$component == "cjar", ]
  sc <- td$value[td$value$component == "cjscore", ]
  check(same_num(wald$estimate, panel$cjive$coefficient) &&
          same_num(wald$std.error, panel$cjive$se) &&
          same_num(wald$statistic, panel$cjive$statistic) &&
          same_num(wald$p.value, panel$cjive$p.value) &&
          same_num(wald$conf.low, panel$cjive$conf.low) &&
          same_num(wald$conf.high, panel$cjive$conf.high) &&
          wald$null.value == 0,
        "CJIVE/Wald tidy row contains only CJIVE/Wald quantities")
  check(is.na(ar$estimate) && is.na(ar$std.error) &&
          same_num(ar$statistic, panel$cjar$statistic) &&
          same_num(ar$p.value, panel$cjar$p.value) &&
          ar$null.value == beta0 &&
          identical(ar$conf.set[[1L]], panel$cjar$conf_set),
        "CJAR tidy row is a test of its explicit stored null and retains its set")
  check(is.na(sc$estimate) && is.na(sc$std.error) &&
          same_num(sc$statistic, panel$cjscore$statistic) &&
          same_num(sc$p.value, panel$cjscore$p.value) &&
          sc$null.value == beta0 &&
          identical(sc$conf.set[[1L]], panel$cjscore$conf_set),
        "CJS tidy row is a test of its explicit stored null and retains its set")
  check(all(td$value$null.value[td$value$component != "cjive"] == beta0),
        "a nonzero beta0 cannot be mistaken for a zero-null table row")
}

td_ar <- capture_error(clusterIV:::tidy.cjar(panel$cjar))
td_sc <- capture_error(clusterIV:::tidy.cjscore(panel$cjscore))
check(is.null(td_ar$error) && td_ar$value$null.value == beta0 &&
        identical(td_ar$value$component, "cjar") &&
        identical(td_ar$value$conf.set[[1L]], panel$cjar$conf_set),
      "tidy.cjar exposes procedure, null, topology, and the complete set")
check(is.null(td_sc$error) && td_sc$value$null.value == beta0 &&
        identical(td_sc$value$component, "cjscore") &&
        identical(td_sc$value$conf.set[[1L]], panel$cjscore$conf_set),
      "tidy.cjscore exposes procedure, null, topology, and the complete set")
check(!is.null(capture_error(clusterIV:::tidy.cjar(panel$cjar,
                                                    typo = TRUE))$error),
      "tidy methods reject unknown arguments in ...")
dot_probes <- list(
  function() clusterIV:::tidy.cjive(panel$cjive, typo = TRUE),
  function() clusterIV:::tidy.cjar(panel$cjar, typo = TRUE),
  function() clusterIV:::tidy.cjscore(panel$cjscore, typo = TRUE),
  function() clusterIV:::tidy.iv_infer(panel, typo = TRUE),
  function() clusterIV:::glance.cjive(panel$cjive, typo = TRUE),
  function() clusterIV:::glance.cjar(panel$cjar, typo = TRUE),
  function() clusterIV:::glance.cjscore(panel$cjscore, typo = TRUE),
  function() clusterIV:::glance.iv_infer(panel, typo = TRUE)
)
dot_rejected <- vapply(dot_probes, function(f) {
  !is.null(capture_error(f())$error)
}, logical(1))
check(all(dot_rejected),
      "every tidy/glance method rejects unsupported arguments in ...")
partial_tidier_args <- list(
  capture_error(clusterIV:::tidy.cjar(panel$cjar, conf.in = FALSE))$error,
  capture_error(clusterIV:::tidy.cjar(panel$cjar, conf.lev = 0.9))$error,
  capture_error(clusterIV:::tidy.cjar(panel$cjar, coef = TRUE))$error,
  capture_error(clusterIV:::glance.cjar(panel$cjar, gof = NULL))$error
)
check(all(vapply(partial_tidier_args, inherits, logical(1), what = "error")),
      "tidy/glance reject partially matched compatibility argument names")

# -------------------------------------------------------------------------
# B. Term labels come from the endogenous regressor.
# -------------------------------------------------------------------------
`endogenous rate` <- x
term_vec <- quiet(cjive(y, `endogenous rate`, Z, cluster = cl))
term_subset1 <- quiet(cjive(y, dat[["x"]], Z, cluster = cl))
term_subset2 <- quiet(cjive(y, dat[, "x"], Z, cluster = cl))
term_ar <- quiet(cjar(y, `endogenous rate`, Z, cluster = cl))
term_sc <- quiet(cjscore(y, `endogenous rate`, Z, cluster = cl))
term_panel <- quiet(iv_infer(y, `endogenous rate`, Z, cluster = cl,
                             tests = c("cjar", "cjscore")))
term_legacy <- quiet(cjive(y ~ `endogenous rate` | z1 + z2,
                           data = transform(dat,
                             `endogenous rate` = x, check.names = FALSE),
                           cluster = cl))
term_fixest <- quiet(cjive(y ~ 1 | `endogenous rate` ~ z1 + z2,
                           data = transform(dat,
                             `endogenous rate` = x, check.names = FALSE),
                           cluster = cl))
check(identical(term_vec$term, "endogenous rate") &&
        identical(term_subset1$term, "x") &&
        identical(term_subset2$term, "x") &&
        identical(term_ar$term, "endogenous rate") &&
        identical(term_sc$term, "endogenous rate") &&
        identical(term_panel$term, "endogenous rate") &&
        all(clusterIV:::tidy.iv_infer(term_panel)$term ==
              "endogenous rate") &&
        identical(term_legacy$term, "endogenous rate") &&
        identical(term_fixest$term, "endogenous rate") &&
        identical(clusterIV:::tidy.cjive(term_vec)$term,
                  "endogenous rate"),
      "vector/subset calls and both formula dialects retain the endogenous name")

# -------------------------------------------------------------------------
# C. iv_infer exposes and shares the cross-fit variance computation.
# -------------------------------------------------------------------------
cf_panel <- capture_error(quiet(iv_infer(
  y, x, Z, cluster = cl, beta0 = beta0,
  tests = c("cjar", "cjscore"), variance = "crossfit")))
cf_ok <- is.null(cf_panel$error) &&
  identical(cf_panel$value$variance, "crossfit") &&
  identical(cf_panel$value$cjar$variance_estimator, "crossfit") &&
  identical(cf_panel$value$cjscore$variance_estimator, "crossfit")
check(cf_ok, "iv_infer accepts and records variance = 'crossfit'",
      if (!is.null(cf_panel$error)) conditionMessage(cf_panel$error))
if (cf_ok) {
  ar_cf <- quiet(cjar(y, x, Z, cluster = cl, beta0 = beta0,
                      variance = "crossfit"))
  sc_cf <- quiet(cjscore(y, x, Z, cluster = cl, beta0 = beta0,
                         variance = "crossfit"))
  check(same_num(cf_panel$value$cjar$statistic, ar_cf$statistic) &&
          same_num(cf_panel$value$cjar$coef_var, ar_cf$coef_var) &&
          same_set(cf_panel$value$cjar$conf_set, ar_cf$conf_set),
        "panel CJAR cross-fit component equals standalone cjar")
  check(same_num(cf_panel$value$cjscore$statistic, sc_cf$statistic) &&
          same_num(cf_panel$value$cjscore$coef_var, sc_cf$coef_var) &&
          same_set(cf_panel$value$cjscore$conf_set, sc_cf$conf_set),
        "panel CJS cross-fit component equals standalone cjscore")
  cf_print <- capture.output(print(cf_panel$value))
  check(any(grepl("variance estimator.*cross-fit", cf_print,
                  ignore.case = TRUE)),
        "print.iv_infer identifies the shared cross-fit variance estimator")
}

# A trace is a credible implementation gate: both selected jackknife tests
# must share one .crossfit_var() call.  This block runs only once the public
# cross-fit route above exists.
if (cf_ok) {
  ns <- asNamespace("clusterIV")
  crossfit_calls <- 0L
  trace(".crossfit_var", where = ns,
        tracer = quote(crossfit_calls <<- crossfit_calls + 1L), print = FALSE)
  traced <- capture_error(quiet(iv_infer(
    y, x, Z, cluster = cl, tests = c("cjar", "cjscore"),
    variance = "crossfit")))
  untrace(".crossfit_var", where = ns)
  check(is.null(traced$error) && crossfit_calls == 1L,
        "iv_infer computes cross-fit variance exactly once for CJAR + CJS",
        paste("calls =", crossfit_calls))

  crossfit_calls <- 0L
  trace(".crossfit_var", where = ns,
        tracer = quote(crossfit_calls <<- crossfit_calls + 1L), print = FALSE)
  skipped <- capture_error(quiet(iv_infer(
    y, x, Z, cluster = cl, tests = NULL, variance = "crossfit")))
  untrace(".crossfit_var", where = ns)
  check(is.null(skipped$error) && crossfit_calls == 0L,
        "iv_infer skips cross-fit work when no jackknife test is selected",
        paste("calls =", crossfit_calls))

  cf_formula <- quiet(iv_infer(y ~ x | z1 + z2, data = dat, cluster = cl,
                               beta0 = beta0,
                               tests = c("cjar", "cjscore"),
                               variance = "crossfit"))
  check(same_num(cf_formula$cjar$coef_var,
                 cf_panel$value$cjar$coef_var) &&
          same_num(cf_formula$cjscore$coef_var,
                   cf_panel$value$cjscore$coef_var),
        "formula and vector panels agree under variance = 'crossfit'")
}

# -------------------------------------------------------------------------
# D. tests has exact, canonical, CJIVE-only semantics.
# -------------------------------------------------------------------------
default_tests <- capture_error(quiet(iv_infer(y, x, Z, cluster = cl)))
default_formula_tests <- capture_error(quiet(iv_infer(
  y ~ x | z1 + z2, data = dat, cluster = cl)))
null_tests <- capture_error(quiet(iv_infer(y, x, Z, cluster = cl,
                                            tests = NULL)))
zero_tests <- capture_error(quiet(iv_infer(y, x, Z, cluster = cl,
                                            tests = character(0))))
is_cjive_only <- function(z) {
  is.null(z$error) && identical(z$value$tests, character(0)) &&
    !is.null(z$value$cjive) && is.null(z$value$cjar) &&
    is.null(z$value$cjscore)
}
check(is_cjive_only(null_tests), "tests = NULL returns a CJIVE-only panel")
check(is_cjive_only(zero_tests),
      "tests = character(0) returns a CJIVE-only panel")
check(is.null(default_tests$error) &&
        identical(default_tests$value$tests, c("cjar", "cjscore")) &&
        !is.null(default_tests$value$cjar) &&
        !is.null(default_tests$value$cjscore),
      "omitted tests selects the recommended CJAR/CJS pair")
check(is.null(default_formula_tests$error) &&
        identical(default_formula_tests$value$tests,
                  c("cjar", "cjscore")),
      "formula omitted tests selects the recommended CJAR/CJS pair")

canon <- capture_error(quiet(iv_infer(
  y, x, Z, cluster = cl,
  tests = c("cjscore", "cjar", "cjscore"))))
check(is.null(canon$error) &&
        identical(canon$value$tests, c("cjar", "cjscore")),
      "tests are deduplicated and stored in canonical order")
partial <- capture_error(quiet(iv_infer(y, x, Z, cluster = cl,
                                         tests = "cjs")))
check(is.null(partial$error) && identical(partial$value$tests, "cjscore"),
      "an unambiguous partial selector resolves to its canonical value")
bad_tests <- list(NA_character_, "", "cj", "wald", 1, TRUE,
                  "ar", "score")
for (v in bad_tests) {
  ans <- capture_error(quiet(iv_infer(y, x, Z, cluster = cl, tests = v)))
  check(!is.null(ans$error),
        paste0("tests rejects invalid selector: ", paste(v, collapse = ",")))
}
if (is_cjive_only(null_tests)) {
  po <- capture.output(print(null_tests$value))
  check(any(grepl("CJIVE", po)) &&
          !any(grepl("CJAR/CJS workflow", po, fixed = TRUE)),
        "CJIVE-only print header does not claim a CJAR/CJS workflow")
  check(nrow(clusterIV:::tidy.iv_infer(null_tests$value)) == 1L,
        "tidy.iv_infer supports a CJIVE-only object")
}

# -------------------------------------------------------------------------
# E. Compatibility: old positions and the historical p field keep working.
# -------------------------------------------------------------------------
old_named <- quiet(cjive(y, x, Z, cluster = cl, controls = NULL,
                         weights = NULL, level = 0.9, intercept = TRUE,
                         method = "dense"))
old_pos <- capture_error(quiet(cjive(y, x, Z, cl, NULL, NULL, 0.9, TRUE,
                                     "dense")))
check(is.null(old_pos$error) &&
        same_num(old_pos$value$coefficient, old_named$coefficient) &&
        identical(old_pos$value$level, 0.9),
      "historical positional cjive call retains its meaning",
      if (!is.null(old_pos$error)) conditionMessage(old_pos$error))
check(identical(old_named$p, old_named$k),
      "cjive restores p as a compatibility alias for k")
legacy_object <- old_named
legacy_object$k <- NULL
legacy_object$k_controls <- NULL
legacy_object$inference <- NULL
legacy_object$term <- NULL
legacy_glance <- capture_error(clusterIV:::glance.cjive(legacy_object))
check(is.null(legacy_glance$error) &&
        identical(legacy_glance$value$k, legacy_object$p) &&
        is.na(legacy_glance$value$k_controls) &&
        identical(legacy_glance$value$inference, "asymptotic") &&
        identical(legacy_glance$value$term, "x"),
      "glance.cjive remains usable on live-v0.1 serialized objects")

old_formula_named <- quiet(cjive(y ~ x | z1 + z2, data = dat,
                                  cluster = cl, controls = NULL,
                                  weights = NULL, level = 0.9,
                                  intercept = TRUE, method = "dense"))
old_formula_pos <- capture_error(quiet(cjive(
  y ~ x | z1 + z2, dat, cl, NULL, NULL, 0.9, TRUE, "dense")))
check(is.null(old_formula_pos$error) &&
        same_num(old_formula_pos$value$coefficient,
                 old_formula_named$coefficient) &&
        identical(old_formula_pos$value$level, 0.9),
      "historical positional cjive formula call retains its meaning",
      if (!is.null(old_formula_pos$error))
        conditionMessage(old_formula_pos$error))

cmp_named <- quiet(iv_compare(y, x, Z, cluster = cl, controls = NULL,
                              weights = NULL, level = 0.9, intercept = TRUE))
cmp_pos <- capture_error(quiet(iv_compare(y, x, Z, cl, NULL, NULL,
                                          0.9, TRUE)))
check(is.null(cmp_pos$error) &&
        same_num(cmp_pos$value$coefficient, cmp_named$coefficient),
      "historical positional iv_compare call retains its meaning",
      if (!is.null(cmp_pos$error)) conditionMessage(cmp_pos$error))
cmp_formula_named <- quiet(iv_compare(y ~ x | z1 + z2, data = dat,
                                      cluster = cl, controls = NULL,
                                      weights = NULL, level = 0.9,
                                      intercept = TRUE))
cmp_formula_pos <- capture_error(quiet(iv_compare(
  y ~ x | z1 + z2, dat, cl, NULL, NULL, 0.9, TRUE)))
check(is.null(cmp_formula_pos$error) &&
        same_num(cmp_formula_pos$value$coefficient,
                 cmp_formula_named$coefficient),
      "historical positional iv_compare formula call retains its meaning",
      if (!is.null(cmp_formula_pos$error))
        conditionMessage(cmp_formula_pos$error))

# -------------------------------------------------------------------------
# F. Print topology is read from conf_set, never from F diagnostics.
# -------------------------------------------------------------------------
bounded <- panel$cjar
bounded$F_CJ <- NA_real_
bounded_text <- capture.output(print(bounded))
bounded_summary_text <- capture.output(print(summary(bounded)))
check(!any(grepl("confidence set unbounded", bounded_text, fixed = TRUE)),
      "print.cjar does not declare a bounded set unbounded because F_CJ is NA")
check(!any(grepl("unbounded", bounded_summary_text, ignore.case = TRUE)),
      "summary.cjar does not infer unboundedness from an NA F diagnostic")

unbounded <- panel$cjar
unbounded$conf_set <- rbind(c(-Inf, -1), c(2, Inf))
colnames(unbounded$conf_set) <- c("lower", "upper")
unbounded$shape <- "two_rays"
unbounded$F_CJ <- Inf
unbounded_text <- capture.output(print(unbounded))
check(any(grepl("unbounded", unbounded_text, ignore.case = TRUE)),
      "print.cjar reports actual unbounded topology even when F_CJ is large")

three_part <- panel$cjscore
three_part$conf_set <- rbind(c(-Inf, -2), c(-0.5, 0.5), c(2, Inf))
colnames(three_part$conf_set) <- c("lower", "upper")
three_part$shape <- "two_rays"
three_part$F_CJS <- Inf
three_print <- capture.output(print(three_part))
three_summary <- capture.output(print(summary(three_part)))
three_tidy <- clusterIV:::tidy.cjscore(three_part)
check(sum(lengths(regmatches(three_print,
                             gregexpr(" U ", three_print,
                                      fixed = TRUE)))) >= 2L &&
        any(grepl("unbounded", three_summary, ignore.case = TRUE)),
      "CJS print/summary render both rays and every bounded middle component")
check(three_tidy$n.components == 3L &&
        three_tidy$unbounded.left && three_tidy$unbounded.right &&
        identical(three_tidy$conf.set[[1L]], three_part$conf_set),
      "tidy.cjscore retains complete three-component topology in columns")

# -------------------------------------------------------------------------
# G. Plot arguments, validation, selected-component behavior, and return.
# -------------------------------------------------------------------------
pdf(tempfile(fileext = ".pdf"))
plot_custom <- capture_error(withVisible(plot(
  panel$cjar, xlab = "effect", ylab = "tail probability",
  main = "custom", ylim = c(0, 1), col = "firebrick", lty = 2L)))
plot_bad_range <- capture_error(plot(panel$cjar, from = 1, to = 1))
plot_bad_n <- capture_error(plot(panel$cjar, n = 0))
plot_nonfinite <- capture_error(plot(panel$cjar, from = -1, to = Inf))
plot_vector_limit <- capture_error(plot(panel$cjar, from = c(-1, 0), to = 1))
plot_fractional_n <- capture_error(plot(panel$cjar, n = 2.5))
plot_bad_ylim <- capture_error(plot(panel$cjar, ylim = c(0, Inf)))
plot_xlim_collision <- capture_error(plot(panel$cjar, xlim = c(-1, 1)))
only_cjs <- quiet(iv_infer(y, x, Z, cluster = cl, tests = "cjscore"))
style_seen <- NULL
shade_input <- NULL
trace(".plot_pcurve", where = ns,
      tracer = quote(style_seen <<- list(main = main, col = col,
                                          lty = lty, lwd = lwd)),
      print = FALSE)
trace(".clip_regions", where = ns,
      tracer = quote(shade_input <<- conf_set), print = FALSE)
plot_cjs_panel <- capture_error(withVisible(plot(
  only_cjs, col = "purple", lty = 5L, lwd = 2.25)))
untrace(".clip_regions", where = ns)
untrace(".plot_pcurve", where = ns)
dev.off()
check(is.null(plot_custom$error) && !plot_custom$value$visible &&
        identical(plot_custom$value$value, panel$cjar),
      "plot accepts customization without duplicate matching and returns invisible(x)",
      if (!is.null(plot_custom$error)) conditionMessage(plot_custom$error))
check(!is.null(plot_bad_range$error), "plot rejects non-increasing from/to")
check(!is.null(plot_bad_n$error), "plot rejects n < 2 instead of drawing blank")
check(!is.null(plot_nonfinite$error) &&
        !is.null(plot_vector_limit$error) &&
        !is.null(plot_fractional_n$error) &&
        !is.null(plot_bad_ylim$error) &&
        !is.null(plot_xlim_collision$error),
      "plot rejects nonfinite/vector limits, fractional n, bad ylim, and xlim collisions")
check(is.null(plot_cjs_panel$error) && !plot_cjs_panel$value$visible &&
        identical(plot_cjs_panel$value$value, only_cjs),
      "plot.iv_infer coherently plots CJS when CJAR is deselected",
      if (!is.null(plot_cjs_panel$error)) conditionMessage(plot_cjs_panel$error))
check(is.null(plot_cjs_panel$error) &&
        identical(style_seen$main, "CJS p-value curve") &&
        identical(style_seen$col, "purple") &&
        identical(style_seen$lty, 5L) &&
        identical(style_seen$lwd, 2.25) &&
        identical(shade_input, only_cjs$cjscore$conf_set),
      "CJS-only title/styles are forwarded and shading uses only its conf_set")

if (length(.failures)) {
  stop(paste0("API-safety oracle failures (", length(.failures), "):\n - ",
              paste(.failures, collapse = "\n - ")), call. = FALSE)
}

cat("\nAll API-safety tests passed.\n")
