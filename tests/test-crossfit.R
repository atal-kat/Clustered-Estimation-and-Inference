# Correctness tests for variance = "crossfit" on cjar() and cjscore()
# (Ligtenberg 2025, Section 5.4), in the style of test-cjar.R: base R, no
# testthat, each block stops with an informative message on failure. The
# brute-force reference is a literal dense transcription of the paper's
# eq. (7) written inline below -- oracle E4-O1 of
# dev/extraction/ligtenberg-crossfit-extraction.md (dense half executed
# pre-implementation in dev/oracle-e4-crossfit-dense.R). Also here: the G0
# regression (plain path unchanged, oracle E4-O4), the conservative
# non-positive-variance conventions on the cross-fit path (Section 3 of the
# extraction doc), and the frozen Mikusheva-Sun non-identity (oracle E4-O6)
# that keeps the D4 labelling honest.

library(clusterIV)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}
near <- function(a, b, tol) max(abs(a - b)) < tol
relnear <- function(a, b, tol) max(abs(a - b) / pmax(1, abs(b))) < tol
errs_msg <- function(expr) {
  tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
}
cjar_q <- function(...) suppressWarnings(cjar(...))
cjs_q <- function(...) suppressWarnings(cjscore(...))

Qpoly <- function(nc, b) nc[1] - nc[2] * b + nc[3] * b^2
Vpoly <- function(wc, b) wc[1] + wc[2]*b + wc[3]*b^2 + wc[4]*b^3 + wc[5]*b^4
in_set <- function(b, cs) nrow(cs) > 0 && any(b >= cs[, 1] & b <= cs[, 2])
set_near <- function(A, B, tol) {
  nrow(A) == nrow(B) && all(is.infinite(A) == is.infinite(B)) &&
    (sum(is.finite(A)) == 0L ||
       max(abs(A[is.finite(A)] - B[is.finite(B)])) < tol)
}

# Grid agreement for the cross-fit acceptance rule: b is accepted iff
# Vhat(b) <= 0 (conservative convention) or Q(b) <= crit * sqrt(k Vhat(b)).
grid_agrees_cf <- function(nc, wc, k, crit, cs, lo = -8, hi = 8,
                           npts = 4001L) {
  fe <- cs[is.finite(cs)]
  if (length(fe)) { lo <- min(lo, min(fe) - 2); hi <- max(hi, max(fe) + 2) }
  bs <- seq(lo, hi, length.out = npts)
  step <- bs[2L] - bs[1L]
  if (length(fe)) {
    bs <- bs[vapply(bs, function(b) min(abs(b - fe)) > step, logical(1))]
  }
  V <- Vpoly(wc, bs)
  acc <- V <= 0 | Qpoly(nc, bs) <= crit * sqrt(k * pmax(V, 0))
  all(acc == vapply(bs, in_set, logical(1), cs = cs))
}

# Literal dense transcription of eq. (7), chain reading of the triple sum
# (extraction doc notes N-1..N-4): explicit solve(crossprod(Z[-idx, ]), ...)
# leave-out projections, no whitening, no downdates. Returns Vhat_AR_CF(b)
# and n * Vhat_S_CF(b) -- exactly the scales stored in coef_var by cjar()
# and cjscore().
brute_cf <- function(y, x, Z, cl, b) {
  idx <- split(seq_along(y), droplevels(as.factor(cl)))
  G <- length(idx)
  n <- length(y)
  k <- ncol(Z)
  e <- y - b * x
  P <- Z %*% solve(crossprod(Z), t(Z))
  lto <- function(drop_idx) {
    Zm <- Z[-drop_idx, , drop = FALSE]
    drop(Z %*% solve(crossprod(Zm), crossprod(Zm, e[-drop_idx])))
  }
  etp <- lapply(seq_len(G), function(g) vector("list", G))
  for (g in seq_len(G - 1L)) for (h in (g + 1L):G) {
    v <- lto(c(idx[[g]], idx[[h]]))
    etp[[g]][[h]] <- v
    etp[[h]][[g]] <- v
  }
  et1 <- lapply(seq_len(G), function(h) lto(idx[[h]]))

  V_AR <- 0; T2 <- 0
  for (g in seq_len(G)) for (h in seq_len(G)) if (g != h) {
    ig <- idx[[g]]; ih <- idx[[h]]
    Pgh <- P[ig, ih, drop = FALSE]
    dg <- (e - etp[[g]][[h]])[ig]
    dh <- (e - etp[[g]][[h]])[ih]
    V_AR <- V_AR + drop(t(dg) %*% Pgh %*% e[ih]) *
                   drop(t(dh) %*% P[ih, ig, drop = FALSE] %*% e[ig])
    T2 <- T2 + drop(t(x[ig]) %*% Pgh %*% e[ih]) *
               drop(t(dg) %*% Pgh %*% x[ih])
  }
  T1 <- 0
  for (h in seq_len(G)) {
    ih <- idx[[h]]
    dh <- (e - et1[[h]])[ih]
    for (g in seq_len(G)) {
      if (g == h) next
      ig <- idx[[g]]
      left <- drop(t(x[ig]) %*% P[ig, ih, drop = FALSE] %*% dh)
      for (j in seq_len(G)) {
        if (j == h) next
        ij <- idx[[j]]
        T1 <- T1 + left * drop(t(e[ih]) %*% P[ih, ij, drop = FALSE] %*% x[ij])
      }
    }
  }
  list(V_AR = 2 / k * V_AR, VSn = T1 + T2)
}

# --- designs -----------------------------------------------------------------
# A: judge design, G = 12 clusters of 5, 6 judges (k = 5).
set.seed(2306)
G_A <- 12L; ng_A <- 5L; n_A <- G_A * ng_A
cl_A <- rep(seq_len(G_A), each = ng_A)
jud_A <- factor(rep(1:6, length.out = n_A))
u_A <- rnorm(G_A)[cl_A]
x_A <- 0.8 * as.numeric(jud_A) + u_A + rnorm(n_A)
y_A <- 0.7 * x_A + u_A + rnorm(n_A)
Zj_A <- stats::model.matrix(~jud_A)[, -1L, drop = FALSE]

# B: unbalanced clusters, continuous matrix z.
set.seed(78)
G_B <- 14L
sz_B <- sample(2:7, G_B, replace = TRUE)
n_B <- sum(sz_B)
cl_B <- rep(seq_len(G_B), sz_B)
u_B <- rnorm(G_B)[cl_B]
Z_B <- matrix(rnorm(n_B * 3L), n_B, 3L)
x_B <- drop(Z_B %*% c(1, -0.6, 0.4)) + u_B + rnorm(n_B)
y_B <- 0.9 * x_B + u_B + rnorm(n_B)

# S: singleton clusters (independent data).
set.seed(32)
n_S <- 30L; k_S <- 3L
Z_S <- matrix(rnorm(n_S * k_S), n_S, k_S)
v_S <- rnorm(n_S)
x_S <- drop(Z_S %*% c(0.8, 0.5, -0.4)) + v_S + rnorm(n_S)
y_S <- 0.6 * x_S + v_S + rnorm(n_S)

# weights and fixed effects for design B
set.seed(58)
fe_B <- factor(sample(letters[1:4], n_B, replace = TRUE))
w_B <- runif(n_B, 0.5, 2)

# === Test 1 (ORACLE E4-O1): kernel == literal eq. (7) transcription =========
betas1 <- c(-2, -0.5, 0, 0.7, 3)
specs <- list(
  list(tag = "judge A", fit_args = list(y_A, x_A, jud_A, cluster = cl_A),
       po_args = list(y_A, x_A, Zj_A, NULL, NULL, TRUE), cl = cl_A),
  list(tag = "matrix B", fit_args = list(y_B, x_B, Z_B, cluster = cl_B),
       po_args = list(y_B, x_B, Z_B, NULL, NULL, TRUE), cl = cl_B),
  list(tag = "singleton S",
       fit_args = list(y_S, x_S, Z_S, cluster = seq_len(n_S)),
       po_args = list(y_S, x_S, Z_S, NULL, NULL, TRUE), cl = seq_len(n_S)),
  list(tag = "weights + FE",
       fit_args = list(y_B, x_B, Z_B, cluster = cl_B, weights = w_B,
                       fixed_effects = fe_B),
       po_args = list(y_B, x_B, Z_B, NULL, w_B, TRUE, list(fe_B)), cl = cl_B))
for (sp in specs) {
  ar <- do.call(cjar_q, c(sp$fit_args, list(variance = "crossfit")))
  sc <- do.call(cjs_q, c(sp$fit_args, list(variance = "crossfit")))
  po <- do.call(clusterIV:::.partial_out, sp$po_args)
  for (b in betas1) {
    bf <- brute_cf(po$y, po$x, po$Z, sp$cl, b)
    ok(relnear(Vpoly(ar$coef_var, b), bf$V_AR, 1e-10),
       sprintf("[%s] Vhat_AR_CF(%.1f) == dense eq. (7) transcription (1e-10 rel)",
               sp$tag, b))
    ok(relnear(sc$coef_var[1] + sc$coef_var[2]*b + sc$coef_var[3]*b^2, bf$VSn,
               1e-10),
       sprintf("[%s] n Vhat_S_CF(%.1f) == dense eq. (7) transcription (1e-10 rel)",
               sp$tag, b))
  }
  # statistic and p-value at beta0 = 0 from the same brute-force objects
  bf0 <- brute_cf(po$y, po$x, po$Z, sp$cl, 0)
  T0 <- Qpoly(ar$coef_num, 0) / sqrt(ar$k * bf0$V_AR)
  ok(relnear(ar$statistic, T0, 1e-10),
     sprintf("[%s] crossfit T(0) == brute force (1e-10 rel)", sp$tag))
  ok(relnear(ar$p.value,
             pchisq(ar$k + sqrt(2 * ar$k) * T0, df = ar$k, lower.tail = FALSE),
             1e-10),
     sprintf("[%s] crossfit p-value == brute force (1e-10 rel)", sp$tag))
  S0 <- Qpoly(c(sc$coef_score[1], sc$coef_score[2], 0), 0)  # s0
  ok(relnear(sc$statistic, S0^2 / bf0$VSn, 1e-10),
     sprintf("[%s] crossfit LM(0) == brute force (1e-10 rel)", sp$tag))
  # statistic coefficients are untouched by the variance choice
  ar_p <- do.call(cjar_q, sp$fit_args)
  sc_p <- do.call(cjs_q, sp$fit_args)
  ok(identical(ar$coef_num, ar_p$coef_num) &&
     identical(sc$coef_score, sc_p$coef_score),
     sprintf("[%s] numerator/score coefficients identical across variance =", sp$tag))
  # inversion consistency at beta0 (NA statistic means accepted)
  ok(isTRUE((is.na(ar$statistic) || ar$statistic <= ar$crit) ==
              in_set(ar$beta0, ar$conf_set)),
     sprintf("[%s] crossfit AR: beta0 classified consistently with the set", sp$tag))
  ok(isTRUE((is.na(sc$statistic) || sc$statistic <= sc$crit) ==
              in_set(sc$beta0, sc$conf_set)),
     sprintf("[%s] crossfit LM: beta0 classified consistently with the set", sp$tag))
  ok(grid_agrees_cf(ar$coef_num, ar$coef_var, ar$k, ar$crit, ar$conf_set),
     sprintf("[%s] crossfit AR conf_set == fine grid inversion", sp$tag))
}

# Supplying the already-computed whitened design is an efficiency-only path:
# it must be coefficient-identical to .crossfit_var() doing its own whitening.
po_share <- clusterIV:::.partial_out(y_B, x_B, Z_B, NULL, NULL, TRUE)
R_share <- clusterIV:::.gram_chol(po_share$Z)
Zt_share <- clusterIV:::.whiten(po_share$Z, R_share)
cf_own <- clusterIV:::.crossfit_var(po_share$y, po_share$x, po_share$Z, cl_B,
                                    R = R_share)
cf_shared <- clusterIV:::.crossfit_var(po_share$y, po_share$x, po_share$Z,
                                       cl_B, R = R_share, Ztil = Zt_share)
ok(identical(cf_shared$w, cf_own$w) && identical(cf_shared$v, cf_own$v),
   "shared Ztil crossfit path is coefficient-identical to internal whitening")

# A leave-two-out Gram can remain formally nonsingular while being too close
# to rank loss for its inverse, and hence the cross-fit variance, to be a
# reproducible double-precision quantity.  The guard is spectral: it must make
# the same decision after either an orthogonal rotation or a nonsingular,
# non-orthogonal reparameterisation of the instruments.  A well-conditioned
# companion design also checks that ordinary outputs are not changed by the
# guard.
set.seed(817)
G_cond <- 6L; ng_cond <- 5L; n_cond <- G_cond * ng_cond
cl_cond <- rep(seq_len(G_cond), each = ng_cond)
z1_cond <- rnorm(n_cond)
tail_cond <- rnorm(n_cond)
z2_cond <- ifelse(cl_cond <= 2L, rnorm(n_cond), 0)
x_cond <- 0.4 * z1_cond - 0.2 * z2_cond + rnorm(n_cond)
y_cond <- 0.8 * x_cond + rnorm(n_cond)
theta_cond <- 0.37
Z_maps <- list(
  identity = diag(2L),
  orthogonal = matrix(c(cos(theta_cond), -sin(theta_cond),
                        sin(theta_cond),  cos(theta_cond)), 2L),
  nonorthogonal = matrix(c(1.2, -0.4, 0.3, 0.9), 2L)
)

Z_safe <- cbind(z1_cond, z2_cond + 1e-2 * tail_cond)
cf_safe <- lapply(Z_maps, function(A)
  clusterIV:::.crossfit_var(y_cond, x_cond, Z_safe %*% A, cl_cond))
ok(all(vapply(cf_safe[-1L], function(z)
       relnear(c(z$w, z$v), c(cf_safe[[1L]]$w, cf_safe[[1L]]$v), 1e-10),
       logical(1))),
   "well-conditioned crossfit coefficients are invariant to orthogonal and non-orthogonal Z transforms")

Z_bad <- cbind(z1_cond, z2_cond + 1e-7 * tail_cond)
cond_msg <- vapply(Z_maps, function(A)
  errs_msg(cjar_q(y_cond, x_cond, Z_bad %*% A, cluster = cl_cond,
                  intercept = FALSE, variance = "crossfit")), character(1))
ok(all(nzchar(cond_msg)) && length(unique(unname(cond_msg))) == 1L &&
     all(grepl("spectral reciprocal condition number", cond_msg,
               fixed = TRUE)),
   "near deletion singularity fails closed identically in every instrument basis")

# === Test 2 (ORACLE E4-O4 / G0): plain path unchanged ========================
# Frozen values produced by the shipped plain path (commit 9624d16 state,
# verified bit-identical through the crossfit change on the reference
# machine; 1e-12 here for cross-platform arithmetic).
set.seed(2306)
G_F <- 15L; ng_F <- 6L; n_F <- G_F * ng_F
cl_F <- rep(seq_len(G_F), each = ng_F)
jud_F <- factor(rep(1:6, length.out = n_F))
u_F <- rnorm(G_F)[cl_F]
x_F <- 0.8 * as.numeric(jud_F) + u_F + rnorm(n_F)
y_F <- 0.7 * x_F + u_F + rnorm(n_F)
arF <- cjar_q(y_F, x_F, jud_F, cluster = cl_F, beta0 = 0.7)
scF <- cjs_q(y_F, x_F, jud_F, cluster = cl_F, beta0 = 0.7)
ok(relnear(c(arF$statistic, arF$p.value),
           c(0.27809263952029956, 0.31813067735592687), 1e-12) &&
   relnear(arF$coef_var,
           c(22.700760773249506, -107.30220471508599, 205.76817936016133,
             -178.82818305727335, 60.823694965710018), 1e-12) &&
   relnear(as.numeric(arF$conf_set),
           c(0.54589321637897903, 0.90634837758600784), 1e-12),
   "G0: plain cjar statistic, coef_var and conf_set match frozen values (1e-12)")
ok(relnear(c(scF$statistic, scF$p.value),
           c(0.089106228157538325, 0.76531649107762445), 1e-12) &&
   relnear(scF$coef_var,
           c(1258.7555858622814, -3122.4503577136093, 2103.1038937719618),
           1e-12) &&
   relnear(as.numeric(scF$conf_set),
           c(0.5668606754024581, 0.84979306465244786), 1e-12),
   "G0: plain cjscore statistic, coef_var and conf_set match frozen values (1e-12)")
# default == explicit variance = "plain", field by field (call differs)
arF2 <- cjar_q(y_F, x_F, jud_F, cluster = cl_F, beta0 = 0.7,
               variance = "plain")
scF2 <- cjs_q(y_F, x_F, jud_F, cluster = cl_F, beta0 = 0.7,
              variance = "plain")
same_fields <- function(a, b) {
  nm <- setdiff(names(a), "call")
  identical(names(a), names(b)) && all(vapply(nm, function(f)
    identical(a[[f]], b[[f]]), logical(1)))
}
ok(same_fields(arF, arF2) && same_fields(scF, scF2),
   "G0: default fit identical to variance = \"plain\", field by field")
ok(identical(arF$variance_estimator, "plain") &&
   identical(scF$variance_estimator, "plain"),
   "variance_estimator recorded as \"plain\" on default fits")

# === Test 3: conservative conventions on the cross-fit path ==================
# (i) Inverter, w4 < 0: Vhat = 1 - 0.01 b^4 is negative for |b| > 100^(1/4),
# so both tails are accepted by the convention; with Q = 1 + b^2, k = 4,
# c_a = 1.5 the middle acceptance interval solves 1.09 b^4 + 2 b^2 - 8 = 0.
# Three rows, two infinite endpoints (the Dufour two-ray signature), F_CJ NA.
bv <- 100^0.25
bm <- sqrt((-2 + sqrt(4 + 4 * 1.09 * 8)) / (2 * 1.09))
inv3 <- clusterIV:::.cjar_invert(c(1, 0, 1), c(1, 0, 0, 0, -0.01), 4, 1.5,
                                 nonpos_var = "accept")
ok(nrow(inv3$conf_set) == 3L && inv3$shape == "two_rays" && is.na(inv3$F_CJ),
   "accept convention, w4 < 0: 3 rows, two_rays, F_CJ NA")
ok(set_near(inv3$conf_set, rbind(c(-Inf, -bv), c(-bm, bm), c(bv, Inf)), 1e-6),
   "accept convention, w4 < 0: endpoints == analytic values")
ok(grid_agrees_cf(c(1, 0, 1), c(1, 0, 0, 0, -0.01), 4, 1.5, inv3$conf_set),
   "accept convention, w4 < 0: grid agreement")
# (ii) Inverter, interior Vhat < 0 window and no boundary roots of h at all:
# Q = 4 + b^2 rejects wherever Vhat > 0, so the set is exactly the window
# where Vhat = -0.5 + b^2 + b^4 <= 0 -- pooled variance roots only.
b0 <- sqrt((-1 + sqrt(3)) / 2)
inv3b <- clusterIV:::.cjar_invert(c(4, 0, 1), c(-0.5, 0, 1, 0, 1), 1, 1,
                                  nonpos_var = "accept")
ok(inv3b$shape == "bounded" && nrow(inv3b$conf_set) == 1L &&
   near(inv3b$conf_set, cbind(-b0, b0), 1e-6),
   "accept convention: set == the Vhat <= 0 window (analytic endpoints)")
ok(grid_agrees_cf(c(4, 0, 1), c(-0.5, 0, 1, 0, 1), 1, 1, inv3b$conf_set),
   "accept convention, Vhat window: grid agreement")
# Same coefficients under the plain convention.  The plain path mirrors the
# point-statistic degenerate rule exactly (V0 <= 0 => T := 0, so the point
# is accepted iff 0 <= c_a): with c_a = 1 the non-positive-variance window
# is accepted -- the set equals the window, matching what the statistic
# itself decides at every beta inside it.  (Previously the inverter
# returned an empty set here while the statistic accepted, an internal
# contradiction.)
inv3c <- clusterIV:::.cjar_invert(c(4, 0, 1), c(-0.5, 0, 1, 0, 1), 1, 1)
ok(inv3c$shape == "bounded" && nrow(inv3c$conf_set) == 1L &&
   near(inv3c$conf_set, cbind(-b0, b0), 1e-6),
   "zero convention, c_a >= 0: set == the T := 0 window (statistic-consistent)")
zero_grid_agrees <- function(nc, wc, k, crit, cs, lo = -8, hi = 8,
                             npts = 4001L) {
  fe <- cs[is.finite(cs)]
  if (length(fe)) { lo <- min(lo, min(fe) - 2); hi <- max(hi, max(fe) + 2) }
  bs <- seq(lo, hi, length.out = npts)
  step <- bs[2L] - bs[1L]
  if (length(fe)) {
    bs <- bs[vapply(bs, function(b) min(abs(b - fe)) > step, logical(1))]
  }
  V <- Vpoly(wc, bs)
  acc <- ifelse(V <= 0, crit >= 0,
                Qpoly(nc, bs) <= crit * sqrt(k * pmax(V, 0)))
  all(acc == vapply(bs, in_set, logical(1), cs = cs))
}
ok(zero_grid_agrees(c(4, 0, 1), c(-0.5, 0, 1, 0, 1), 1, 1, inv3c$conf_set),
   "zero convention, c_a >= 0: grid agreement with the statistic rule")
# The two conventions remain genuinely different objects: at a negative
# critical value the plain rule rejects the T := 0 window (0 > c_a) and the
# set is empty, while the cross-fit accept convention keeps the window
# regardless of the sign of c_a.
inv3d <- clusterIV:::.cjar_invert(c(4, 0, 1), c(-0.5, 0, 1, 0, 1), 1, -0.5)
ok(inv3d$shape == "empty" && nrow(inv3d$conf_set) == 0L,
   "zero convention, c_a < 0: empty set (contrast)")
inv3e <- clusterIV:::.cjar_invert(c(4, 0, 1), c(-0.5, 0, 1, 0, 1), 1, -0.5,
                                  nonpos_var = "accept")
ok(nrow(inv3e$conf_set) == 1L && near(inv3e$conf_set, cbind(-b0, b0), 1e-6),
   "accept convention, c_a < 0: window kept (contrast)")
# (iii) Statistic convention at the object level, mocked cross-fit
# coefficients with Vhat(beta0) < 0 (E4-O7 "construct or mock"): statistic
# NA, p-value 1, beta0 accepted, print note rendered.
d3 <- clusterIV:::.prep_data(y_B, x_B, Z_B, cl_B, NULL, NULL, TRUE)
R3 <- clusterIV:::.gram_chol(d3$Z)
cs3 <- clusterIV:::.cluster_sums(d3$y, d3$x, d3$Z, d3$cluster, R = R3)
fs3 <- clusterIV:::.first_stage(d3$x, d3$Z, R = R3)
ml3 <- clusterIV:::.cluster_leverage(fs3$Ztil, d3$groups, d3$k)
ef3 <- clusterIV:::.eff_f(fs3$Ztil, fs3$t, fs3$e, d3$cluster)
mock <- clusterIV:::.cjar_build(d3, cs3, ml3, ef3, beta0 = 0, level = 0.95,
                                calibration = "chisq", call = quote(mock()),
                                variance = "crossfit",
                                cf = list(w = c(-1, 0, 0, 0, 1e-6),
                                          v = c(0, 0, 0), k = d3$k, G = d3$G))
ok(is.na(mock$statistic) && mock$p.value == 1,
   "crossfit Vhat(beta0) <= 0: statistic NA, p-value 1")
ok(in_set(0, mock$conf_set),
   "crossfit Vhat(beta0) <= 0: beta0 is accepted (in the set)")
out_mock <- capture.output(print(mock))
ok(any(grepl("undefined and beta0 is accepted", out_mock)),
   "crossfit Vhat(beta0) <= 0: print renders the conservative-convention note")
summary_mock <- capture.output(print(summary(mock)))
ok(any(grepl("undefined and beta0 is accepted", summary_mock)),
   "crossfit Vhat(beta0) <= 0: summary renders the conservative-convention note")
panel_mock <- suppressWarnings(iv_infer(
  y_B, x_B, Z_B, cluster = cl_B, tests = c("cjar", "cjscore")))
panel_mock$cjar <- mock
panel_mock$cjscore$statistic <- NA_real_
panel_mock$cjscore$p.value <- 1
panel_mock$cjscore$variance <- -1
panel_mock$cjscore$variance_estimator <- "crossfit"
panel_mock$variance <- "crossfit"
panel_mock_text <- capture.output(print(panel_mock))
ok(sum(grepl("undefined and beta0 is accepted", panel_mock_text,
             fixed = TRUE)) == 2L,
   "iv_infer print explains both conservative non-positive-variance rows")

# === Test 4 (ORACLE E4-O6): NOT the Mikusheva-Sun cross-fit variance ========
# Frozen non-identity, inverted expectation: at singleton clusters the
# cross-fit variance is Ligtenberg's leave-two-out construction, which is
# NOT Mikusheva & Sun's Phihat_2 (their eq. 4, transcribed literally below).
# Both are unbiased; they are not algebraically equal, and no setting of the
# package produces Phihat_2 -- this test keeps the D4 labelling honest.
set.seed(14)
n4 <- 14L; k4 <- 3L
Z4 <- matrix(rnorm(n4 * k4), n4, k4)
x4 <- drop(Z4 %*% c(1, -0.5, 0.7)) + rnorm(n4)
y4 <- 0.4 * x4 + rnorm(n4)
ar4 <- cjar_q(y4, x4, Z4, cluster = seq_len(n4), variance = "crossfit")
po4 <- clusterIV:::.partial_out(y4, x4, Z4, NULL, NULL, TRUE)
P4 <- po4$Z %*% solve(crossprod(po4$Z), t(po4$Z))
M4 <- diag(n4) - P4
b4 <- 0.4
e4 <- po4$y - b4 * po4$x
Me4 <- drop(M4 %*% e4)
Phi2 <- 0
for (i in seq_len(n4)) for (j in seq_len(n4)) if (i != j) {
  Phi2 <- Phi2 + P4[i, j]^2 / (M4[i, i] * M4[j, j] + M4[i, j]^2) *
    (e4[i] * Me4[i]) * (e4[j] * Me4[j])
}
Phi2 <- 2 / k4 * Phi2
gap4 <- abs(Vpoly(ar4$coef_var, b4) - Phi2) / abs(Phi2)
ok(gap4 > 1e-6,
   sprintf("E4-O6: crossfit at singletons differs from MS Phihat_2 (gap %.2g)", gap4))

# === Test 5: recording, labels, interfaces, guards ==========================
arX <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, variance = "crossfit")
scX <- cjs_q(y_A, x_A, jud_A, cluster = cl_A, variance = "crossfit")
ok(identical(arX$variance_estimator, "crossfit") &&
   identical(scX$variance_estimator, "crossfit"),
   "variance_estimator recorded as \"crossfit\"")
outX <- capture.output({ print(arX); print(scX) })
ok(sum(grepl("variance estimator: cross-fit (Ligtenberg 2025, Section 5.4)",
             outX, fixed = TRUE)) == 2L,
   "print: cross-fit variance line rendered by both tests")
# honest labels at G == n: name the construction, claim neither identity
arS <- cjar_q(y_S, x_S, Z_S, cluster = seq_len(n_S), variance = "crossfit")
scS <- cjs_q(y_S, x_S, Z_S, cluster = seq_len(n_S), variance = "crossfit")
outS <- capture.output({ print(arS); print(scS) })
ok(sum(grepl("cross-fit variance, Ligtenberg 2025 construction", outS)) == 2L &&
   !any(grepl("Mikusheva-Sun statistic", outS)) &&
   !any(grepl("Matsushita-Otsu", outS)),
   "G == n crossfit print labels: construction named, neither identity claimed")
# formula interface passes variance through
datA <- data.frame(y = y_A, x = x_A, judge = jud_A, cl = cl_A)
arf <- cjar_q(y ~ x | judge, data = datA, cluster = ~cl,
              variance = "crossfit")
ok(identical(arf$statistic, arX$statistic) &&
   identical(arf$coef_var, arX$coef_var) &&
   set_near(arf$conf_set, arX$conf_set, 1e-12),
   "formula interface == default interface under variance = \"crossfit\"")
# confint at a new level re-inverts under the accept convention
ci90 <- confint(arX, level = 0.90)
dir90 <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, variance = "crossfit",
                level = 0.90)$conf_set
ok(set_near(ci90, dir90, 1e-6),
   "crossfit confint at a new level == direct refit at that level")
# guard: fewer than 3 clusters cannot support the leave-two-out construction
msg5 <- errs_msg(cjar_q(y_B, x_B, Z_B, cluster = rep(1:2, length.out = n_B),
                        variance = "crossfit"))
ok(nzchar(msg5) && grepl("at least 3 clusters", msg5),
   "G = 2 stops with the at-least-3-clusters message")

cat("\nAll crossfit tests passed.\n")
