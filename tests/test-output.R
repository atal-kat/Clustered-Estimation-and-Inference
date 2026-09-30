# Tests for the WP3 output layer (base R, no testthat): tidy()/glance(),
# summary.cjar/.cjscore, and the p-value-curve plots. Graphics tests assert
# behaviour (no error, invisible return, shading == conf_set clipped to the
# range), never pixels.

library(clusterIV)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}
near <- function(a, b, tol) max(abs(a - b)) < tol
quiet <- function(expr) suppressWarnings(expr)

# --- fixtures ----------------------------------------------------------------
set.seed(42)
G <- 30L; ng <- 8L; n <- G * ng
cl <- rep(seq_len(G), each = ng)
judge <- factor(rep(1:6, length.out = n))
u <- rnorm(G)[cl]
x <- 0.6 * as.numeric(judge) + u + rnorm(n)
y <- 0.5 * x + u + rnorm(n)
x_w <- 0.02 * as.numeric(judge) + u + rnorm(n)     # weak first stage
y_w <- 0.5 * x_w + u + rnorm(n)

est <- cjive(y, x, judge, cluster = cl)
ar <- quiet(cjar(y, x, judge, cluster = cl))
sc <- quiet(cjscore(y, x, judge, cluster = cl))
wf <- quiet(iv_infer(y, x, judge, cluster = cl,
                     tests = c("cjar", "cjscore")))
ar_w <- quiet(cjar(y_w, x_w, judge, cluster = cl))  # unbounded set

# A valid public fit whose tail label alone is incomplete: two unbounded rays
# plus a bounded accepted middle component.
set.seed(64)
G_top <- 10L; ng_top <- 4L; n_top <- G_top * ng_top
cl_top <- rep(seq_len(G_top), each = ng_top)
Z_top <- matrix(rnorm(n_top * 3L), n_top, 3L)
u_top <- rnorm(G_top)[cl_top]
x_top <- drop(Z_top %*% rnorm(3L, sd = 0.12)) + u_top + rnorm(n_top)
y_top <- runif(1L, -1, 1) * x_top + u_top + rnorm(n_top)
sc_top <- quiet(cjscore(y_top, x_top, Z_top, cluster = cl_top,
                        intercept = FALSE))
ok(sc_top$shape == "two_rays" && sc_top$n_components == 3L &&
     sc_top$unbounded_left && sc_top$unbounded_right,
   "public fit retains two_rays as a tail label and stores complete topology")
td_top <- clusterIV:::tidy.cjscore(sc_top)
ok(td_top$n.components == 3L && td_top$unbounded.left &&
     td_top$unbounded.right && is.na(td_top$conf.low) &&
     is.na(td_top$conf.high) &&
     identical(td_top$conf.set[[1L]], sc_top$conf_set),
   "public three-component CJS topology survives the tidy schema")
top_text <- capture.output(print(sc_top))
top_set_line <- top_text[grepl("confidence set =", top_text, fixed = TRUE)]
ok(length(top_set_line) == 1L &&
     lengths(regmatches(top_set_line, gregexpr(" U ", top_set_line,
                                               fixed = TRUE))) == 2L,
   "print.cjscore renders both rays and the bounded middle component")

# === tidy() / glance() =======================================================
td <- clusterIV:::tidy.cjive(est)
ok(is.data.frame(td) && nrow(td) == 1L &&
   all(c("term", "component", "procedure", "estimate", "std.error",
         "statistic", "p.value", "null.value", "conf.low", "conf.high",
         "conf.set") %in% names(td)) && td$component == "cjive" &&
   td$null.value == 0 && td$estimate == est$coefficient &&
   td$conf.low == est$conf.low,
   "tidy.cjive: one coherent CJIVE/Wald row with an explicit zero null")
gl <- clusterIV:::glance.cjive(est)
ok(is.data.frame(gl) && nrow(gl) == 1L && gl$nobs == n && gl$k == est$k &&
   gl$path == "dense" && is.null(gl$F_eff),
   "glance.cjive: nobs/G/k/k_controls/maxlev/path/inference; no F_eff")

td_ar <- clusterIV:::tidy.cjar(ar)
ok(is.na(td_ar$estimate) && td_ar$statistic == ar$statistic &&
   td_ar$component == "cjar" && td_ar$null.value == ar$beta0 &&
   identical(td_ar$conf.set[[1L]], ar$conf_set) &&
   td_ar$conf.low == ar$conf_set[1, 1] && td_ar$conf.high == ar$conf_set[1, 2],
   "tidy.cjar: explicit test null and complete set; finite interval fills endpoints")
td_arw <- clusterIV:::tidy.cjar(ar_w)
ok(is.na(td_arw$conf.low) && is.na(td_arw$conf.high),
   "tidy.cjar: an unbounded set is never flattened into an interval (NA)")
deg <- structure(list(conf_set = cbind(lower = 1, upper = 1)), class = "cjar")
ok(all(is.na(clusterIV:::.conf_set_interval(deg$conf_set))),
   "tidy: a degenerate [b, b] set yields NA, not a zero-width interval")
two <- cbind(lower = c(-3, 1), upper = c(-1, 3))
ok(all(is.na(clusterIV:::.conf_set_interval(two))),
   "tidy: two disjoint intervals yield NA")
gl_ar <- clusterIV:::glance.cjar(ar)
ok(gl_ar$shape == "bounded" && gl_ar$bounded && gl_ar$F_CJ == ar$F_CJ &&
   gl_ar$calibration == "chisq" && is.null(gl_ar$F_eff),
   "glance.cjar: shape/bounded/F_CJ/crit/calibration; no F_eff")
gl_sc <- clusterIV:::glance.cjscore(sc)
ok(gl_sc$F_CJS == sc$F_CJS && is.null(gl_sc$calibration) && is.null(gl_sc$F_eff),
   "glance.cjscore: F_CJS, no calibration column, no F_eff")

td_wf <- clusterIV:::tidy.iv_infer(wf)
ok(nrow(td_wf) == 3L &&
   identical(td_wf$component, c("cjive", "cjar", "cjscore")),
   "tidy.iv_infer: one labelled row per present procedure")
wald_row <- td_wf[td_wf$component == "cjive", ]
ar_row <- td_wf[td_wf$component == "cjar", ]
ok(wald_row$estimate == wf$cjive$coefficient &&
   wald_row$std.error == wf$cjive$se &&
   wald_row$statistic == wf$cjive$statistic &&
   wald_row$p.value == wf$cjive$p.value &&
   wald_row$conf.low == wf$cjive$conf.low &&
   ar_row$statistic == wf$cjar$statistic && is.na(ar_row$estimate) &&
   identical(ar_row$conf.set[[1L]], wf$cjar$conf_set),
   "tidy.iv_infer: CJIVE/Wald and CJAR quantities never share a row")
ok(is.null(attr(td_wf, "conf.int.source")),
   "tidy.iv_infer semantics live in ordinary columns, not a fragile attribute")
gl_wf <- clusterIV:::glance.iv_infer(wf)
ok(gl_wf$F_CJ == wf$cjar$F_CJ &&
   gl_wf$cjar.shape == wf$cjar$shape &&
   gl_wf$cjscore.shape == wf$cjscore$shape &&
   gl_wf$cjar.crit == wf$cjar$crit &&
   gl_wf$cjscore.crit == wf$cjscore$crit &&
   is.null(gl_wf$F_eff),
   "glance.iv_infer: component-labelled CJAR/CJS diagnostics; no F_eff")
cjs_only <- quiet(iv_infer(y, x, judge, cluster = cl, tests = "cjscore"))
gl_cjs_only <- clusterIV:::glance.iv_infer(cjs_only)
ok(is.na(gl_cjs_only$F_CJ) && gl_cjs_only$F_CJS == cjs_only$cjscore$F_CJS &&
   is.na(gl_cjs_only$cjar.crit) &&
   gl_cjs_only$cjscore.crit == cjs_only$cjscore$crit &&
   gl_cjs_only$cjscore.shape == cjs_only$cjscore$shape,
   "glance.iv_infer retains CJS-only topology without a CJAR-labelled value")
cjive_only <- quiet(iv_infer(y, x, judge, cluster = cl, tests = NULL,
                             variance = "crossfit"))
ok(is.na(clusterIV:::glance.iv_infer(cjive_only)$variance.estimator),
   "glance.iv_infer does not report an unused crossfit variance on CJIVE-only panels")

# Conditional registration: with generics available, tidy()/glance() dispatch.
if (requireNamespace("generics", quietly = TRUE)) {
  ok(identical(generics::tidy(est), td) &&
     identical(generics::glance(wf), gl_wf),
     "generics::tidy()/glance() dispatch to the registered methods")
} else {
  cat("SKIP: generics not installed; direct calls tested above\n")
}
# Probes against UNDECLARED consumer packages (broom, modelsummary, texreg,
# fixest) live in dev/integration/consumer-probes.R, not here: package tests
# exercise only base R and declared Suggests.

# === summary() ===============================================================
out_ar <- capture.output(print(summary(ar)))
ok(any(grepl("-- summary", out_ar)) &&
   any(grepl("reading: one bounded accepted component", out_ar)) &&
   any(grepl("F_CJ  =", out_ar)) && any(grepl("F_eff =", out_ar)) &&
   any(grepl("K_eff", out_ar)) &&
   any(grepl("Advisories:", out_ar)) && any(grepl("none", out_ar)),
   "summary.cjar: set reading, strength block (F_CJ and F_eff side by side), advisories")
out_arw <- capture.output(print(summary(ar_w)))
ok(any(grepl("reading: unbounded|whole real line", out_arw)),
   "summary.cjar: unbounded topology gets an endpoint-matrix reading")
out_sc <- capture.output(print(summary(sc)))
ok(any(grepl("Cluster jackknife score test", out_sc)) &&
   any(grepl("F_CJS\\^2 =", out_sc)) && any(grepl("F_eff =", out_sc)),
   "summary.cjscore: strength block present")
out_est <- capture.output(print(summary(est)))
ok(sum(grepl("F_eff =", out_est)) == 1L &&
   any(grepl("cjar\\(\\) fit", out_est)),
   "summary.cjive: effective F prints once with the cjar()/iv_infer() pointer")
# Advisory lines fire in the summary of a small-G, dominating-cluster design.
cl_dom <- ifelse(cl <= 12L, 1L, cl)   # G = 19 (< 20) and cluster 1 holds 40%
ar_dom <- quiet(cjar(y, x, judge, cluster = cl_dom))
out_dom <- capture.output(print(summary(ar_dom)))
ok(any(grepl("asymptotic in the number of clusters", out_dom)) &&
   any(grepl("dominating cluster", out_dom)),
   "summary.cjar: advisory lines listed for a small-G dominating-cluster design")

# === plot(): every shape runs without error and returns invisibly ============
# Synthetic cjar objects from the frozen coefficient tuples of test-cjar.R
# (Tests A/B/C and I), inverted by the package's own inverter.
synth_cjar <- function(nc, wc, k, crit, level = 0.95, calibration = "chisq",
                       beta0 = 0) {
  inv <- clusterIV:::.cjar_invert(nc, wc, k, crit)
  structure(list(coef_num = nc, coef_var = wc, k = k, crit = crit,
                 level = level, calibration = calibration, beta0 = beta0,
                 conf_set = inv$conf_set, shape = inv$shape),
            class = "cjar")
}
shapes <- list(
  bounded        = ar,
  two_intervals  = synth_cjar(c(1, 0, 1), c(0.1, 0, 3, 0, 0.9), 1, 1),
  ray            = synth_cjar(c(0, 1, 0), c(1, 0, 0.1, 0, 0), 1, 1),
  two_rays       = synth_cjar(c(1, 0, 0), c(0.01, 0, 1, 0, 0), 1, 1),
  whole_line     = synth_cjar(c(-1, 0, 0), c(1, 0, 0.1, 0, 0), 1, 1),
  empty          = synth_cjar(c(1, 0, 1), c(1e-4, 0, 1e-4, 0, 1e-4), 5, 1.96),
  degenerate_pt  = synth_cjar(c(4, 2, 1), c(1, 0, 0, 0, 0), 4, 1.5)
)
pdf(tempfile(fileext = ".pdf"))
for (nm in names(shapes)) {
  r <- withVisible(plot(shapes[[nm]]))
  ok(!r$visible && identical(r$value, shapes[[nm]]),
     sprintf("plot.cjar returns invisible(x) [shape: %s]", nm))
}
r <- withVisible(plot(sc))
ok(!r$visible && identical(r$value, sc),
   "plot.cjscore returns invisible(x)")
r <- withVisible(plot(quiet(cjscore(y_w, x_w, judge, cluster = cl))))
ok(!r$visible, "plot.cjscore runs on an unbounded-set fit")
r <- withVisible(plot(wf))
ok(!r$visible && identical(r$value, wf),
   "plot.iv_infer returns invisible(x)")

# Separate changes of outcome and regressor units can leave every stored
# coefficient finite while overflowing the raw w0/w4 or v0/v2 ratio used for
# an unbounded set's default range.  The range calculation must use the same
# scale-safe principle as point evaluation and inversion.
ar_scaled <- quiet(cjar(y_w * 1e75, x_w * 1e-75, judge, cluster = cl))
sc_scaled <- quiet(cjscore(y_w * 1e80, x_w * 1e-80, judge, cluster = cl))
wf_ar_scaled <- quiet(iv_infer(y_w * 1e75, x_w * 1e-75, judge,
                               cluster = cl))
wf_sc_scaled <- quiet(iv_infer(y_w * 1e80, x_w * 1e-80, judge,
                               cluster = cl, tests = "cjscore"))
ok(ar_scaled$shape == "whole_line" &&
     is.infinite(ar_scaled$coef_var[1L] / ar_scaled$coef_var[5L]) &&
     is.finite(clusterIV:::.plot_root_ratio(
       ar_scaled$coef_var[1L], ar_scaled$coef_var[5L], 0.25)),
   "CJAR whole-line default scale survives an overflowing coefficient ratio")
ok(sc_scaled$shape == "whole_line" &&
     is.infinite(sc_scaled$coef_var[1L] / sc_scaled$coef_var[3L]) &&
     is.finite(clusterIV:::.plot_root_ratio(
       sc_scaled$coef_var[1L], sc_scaled$coef_var[3L], 0.5)),
   "CJS whole-line default scale survives an overflowing coefficient ratio")
scaled_plots <- list(cjar = ar_scaled, cjscore = sc_scaled,
                     iv_infer_cjar = wf_ar_scaled,
                     iv_infer_cjscore = wf_sc_scaled)
for (nm in names(scaled_plots)) {
  obj <- scaled_plots[[nm]]
  r <- withVisible(plot(obj))
  ok(!r$visible && identical(r$value, obj),
     paste0("scale-safe unbounded default plot returns invisible(x) [", nm,
            "]"))
}

# Ordinary designs retain the exact direct-arithmetic result.
ok(identical(clusterIV:::.plot_root_ratio(ar$coef_var[1L],
                                          ar$coef_var[5L], 0.25),
             (ar$coef_var[1L] / ar$coef_var[5L])^0.25),
   "plot scale helper preserves the ordinary CJAR arithmetic path")
plot(ar, from = -1, to = 2)
plot(ar, n = 32L)
cat("PASS: plot.cjar honours from/to and n overrides\n")
dev.off()

# === the shading helper: regions == conf_set clipped to the range ============
cases <- list(
  list(cs = rbind(c(0.4, 0.7)), rng = c(0, 1),      exp = rbind(c(0.4, 0.7))),
  list(cs = rbind(c(0.4, 0.7)), rng = c(0.5, 0.6),  exp = rbind(c(0.5, 0.6))),
  list(cs = rbind(c(0.4, 0.7)), rng = c(0.8, 2),    exp = matrix(numeric(0), 0, 2)),
  list(cs = rbind(c(-Inf, -0.3), c(1.2, Inf)), rng = c(-2, 3),
       exp = rbind(c(-2, -0.3), c(1.2, 3))),
  list(cs = rbind(c(-Inf, Inf)), rng = c(-5, 5),    exp = rbind(c(-5, 5))),
  list(cs = matrix(numeric(0), 0, 2), rng = c(-1, 1), exp = matrix(numeric(0), 0, 2)),
  list(cs = rbind(c(-3, -1), c(1, 1), c(2, 4)), rng = c(-2, 3),
       exp = rbind(c(-2, -1), c(1, 1), c(2, 3)))
)
for (i in seq_along(cases)) {
  got <- clusterIV:::.clip_regions(cases[[i]]$cs, cases[[i]]$rng[1], cases[[i]]$rng[2])
  ok(nrow(got) == nrow(cases[[i]]$exp) &&
     (nrow(got) == 0L || near(unname(got), unname(cases[[i]]$exp), 1e-12)),
     sprintf(".clip_regions case %d: regions == conf_set clipped to range", i))
}
# On a real fit: the shaded regions returned by the drawing core equal the
# stored set whenever the range covers it.
pdf(tempfile(fileext = ".pdf"))
regs <- clusterIV:::.plot_pcurve(ar,
  list(CJAR = function(b) clusterIV:::.pval_curve_cjar(ar, b)),
  from = ar$conf_set[1, 1] - 1, to = ar$conf_set[1, 2] + 1, n = 64L,
  main = "")
dev.off()
ok(near(unname(regs), unname(ar$conf_set), 1e-12),
   "drawing core: shaded regions are exactly x$conf_set (range covers the set)")

cat("\nAll output-layer tests passed.\n")
