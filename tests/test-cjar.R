# Correctness tests for cjar() (base R, no testthat), in the style of
# test-cjive.R: each block stops with an informative message on failure, so a
# failure fails `R CMD check`. Brute-force references are written inline.

library(clusterIV)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}
near <- function(a, b, tol) max(abs(a - b)) < tol
relnear <- function(a, b, tol) max(abs(a - b) / pmax(1, abs(b))) < tol
errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
cjar_q <- function(...) suppressWarnings(cjar(...))   # quiet small-G advisories

# Polynomial evaluation from stored coefficients (the object contract).
Qpoly <- function(nc, b) nc[1] - nc[2] * b + nc[3] * b^2
Vpoly <- function(wc, b) wc[1] + wc[2]*b + wc[3]*b^2 + wc[4]*b^3 + wc[5]*b^4
hpoly <- function(nc, wc, k, crit, b)
  Qpoly(nc, b) - crit * sqrt(k * pmax(Vpoly(wc, b), 0))

in_set <- function(b, cs) nrow(cs) > 0 && any(b >= cs[, 1] & b <= cs[, 2])
set_near <- function(A, B, tol) {
  nrow(A) == nrow(B) && all(is.infinite(A) == is.infinite(B)) &&
    (sum(is.finite(A)) == 0L ||
       max(abs(A[is.finite(A)] - B[is.finite(B)])) < tol)
}

# Grid inversion reference: membership by the sign of h must match interval
# membership everywhere except within one grid step of each finite endpoint.
# A fixed point count keeps the step proportional to the range, so designs
# whose endpoints sit far out (the near-degenerate quartic) stay cheap.
grid_agrees <- function(nc, wc, k, crit, cs, lo = -6, hi = 6, npts = 4001L) {
  fe <- cs[is.finite(cs)]
  if (length(fe)) { lo <- min(lo, min(fe) - 2); hi <- max(hi, max(fe) + 2) }
  bs <- seq(lo, hi, length.out = npts)
  step <- bs[2L] - bs[1L]
  if (length(fe)) {
    bs <- bs[vapply(bs, function(b) min(abs(b - fe)) > step, logical(1))]
  }
  hb <- hpoly(nc, wc, k, crit, bs)
  ms <- vapply(bs, in_set, logical(1), cs = cs)
  all((hb <= 0) == ms)
}
grid_agrees_fit <- function(fit, ...)
  grid_agrees(fit$coef_num, fit$coef_var, fit$k, fit$crit, fit$conf_set, ...)

# Dense brute force on partialled data: Pddot explicitly, and
# Vhat = 2/k sum_{g != h} c_gh(beta)^2 by a double loop over cluster pairs.
brute_QV <- function(y, x, Z, cl, beta) {
  k <- ncol(Z)
  cl <- droplevels(as.factor(cl))
  e <- y - beta * x
  P <- Z %*% solve(crossprod(Z), t(Z))
  Pdd <- P
  for (g in levels(cl)) {
    ig <- which(cl == g)
    Pdd[ig, ig] <- 0
  }
  V <- 0
  for (g in levels(cl)) for (h in levels(cl)) if (g != h) {
    ig <- which(cl == g); ih <- which(cl == h)
    cgh <- drop(t(e[ig]) %*% P[ig, ih, drop = FALSE] %*% e[ih])
    V <- V + cgh^2
  }
  list(Q = drop(t(e) %*% Pdd %*% e), V = 2 / k * V)
}

# --- designs -----------------------------------------------------------------
# A: judge design, n = 90, G = 15, k = 5 (6 judges, one reference dropped).
set.seed(2306)
G_A <- 15L; ng_A <- 6L; n_A <- G_A * ng_A
cl_A <- rep(seq_len(G_A), each = ng_A)
jud_A <- factor(rep(1:6, length.out = n_A))
u_A <- rnorm(G_A)[cl_A]
x_A <- 0.8 * as.numeric(jud_A) + u_A + rnorm(n_A)
y_A <- 0.7 * x_A + u_A + rnorm(n_A)
Zj_A <- stats::model.matrix(~jud_A)[, -1L, drop = FALSE]

# B: unbalanced clusters, continuous matrix z.
set.seed(77)
G_B <- 24L
sz_B <- sample(2:9, G_B, replace = TRUE)
n_B <- sum(sz_B)
cl_B <- rep(seq_len(G_B), sz_B)
u_B <- rnorm(G_B)[cl_B]
Z_B <- matrix(rnorm(n_B * 3L), n_B, 3L)
x_B <- drop(Z_B %*% c(1, -0.6, 0.4)) + u_B + rnorm(n_B)
y_B <- 0.9 * x_B + u_B + rnorm(n_B)

# Weak variant of A: first stage scaled toward zero.
set.seed(41)
x_Aw <- 0.02 * as.numeric(jud_A) + u_A + rnorm(n_A)
y_Aw <- 0.5 * x_Aw + u_A + rnorm(n_A)

# === Test 1: brute force vs kernel ==========================================
betas1 <- c(-2, -0.5, 0, 0.7, 3)
for (des in list(list(y = y_A, x = x_A, Z = Zj_A, cl = cl_A, tag = "judge A"),
                 list(y = y_B, x = x_B, Z = Z_B, cl = cl_B, tag = "matrix B"))) {
  po <- clusterIV:::.partial_out(des$y, des$x, des$Z, NULL, NULL, TRUE)
  cf <- clusterIV:::.cjar_coefs(clusterIV:::.cluster_sums(po$y, po$x, po$Z, des$cl))
  for (b in betas1) {
    bf <- brute_QV(po$y, po$x, po$Z, des$cl, b)
    ok(relnear(Qpoly(cf$n, b), bf$Q, 1e-10),
       sprintf("[%s] Q(%.1f) polynomial == dense brute force (1e-10 rel)", des$tag, b))
    ok(relnear(Vpoly(cf$w, b), bf$V, 1e-10),
       sprintf("[%s] Vhat(%.1f) polynomial == 2/k sum c_gh^2 brute force (1e-10 rel)", des$tag, b))
  }
}

fit1 <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, beta0 = 0.7)
po1 <- clusterIV:::.partial_out(y_A, x_A, Zj_A, NULL, NULL, TRUE)
bf1 <- brute_QV(po1$y, po1$x, po1$Z, cl_A, 0.7)
T1 <- bf1$Q / sqrt(fit1$k * bf1$V)
ok(relnear(fit1$statistic, T1, 1e-10), "T(beta0) == brute force (1e-10 rel)")
ok(relnear(fit1$p.value,
           pchisq(fit1$k + sqrt(2 * fit1$k) * T1, df = fit1$k, lower.tail = FALSE),
           1e-10), "chisq p-value == brute force (1e-10 rel)")
fit1n <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, beta0 = 0.7,
                calibration = "normal")
ok(identical(fit1n$statistic, fit1$statistic) &&
   near(fit1n$crit, qnorm(0.95), 1e-12) &&
   relnear(fit1n$p.value, pnorm(T1, lower.tail = FALSE), 1e-10),
   "normal calibration: same T, z critical value, normal p-value")
ok(near(fit1$crit, (qchisq(0.95, fit1$k) - fit1$k) / sqrt(2 * fit1$k), 1e-12),
   "chisq calibration critical value = (qchisq - k)/sqrt(2k)")

# === Test 2: inversion vs grid; F_CJ boundedness criterion ==================
fitA <- cjar_q(y_A, x_A, jud_A, cluster = cl_A)
fitB <- cjar_q(y_B, x_B, Z_B, cluster = cl_B)
fitW <- cjar_q(y_Aw, x_Aw, jud_A, cluster = cl_A)

ok(fitA$shape == "bounded" && fitA$F_CJ > fitA$crit,
   "strong judge design: bounded set, F_CJ > crit")
ok(fitB$shape == "bounded" && fitB$F_CJ > fitB$crit,
   "strong continuous design: bounded set, F_CJ > crit")
ok(fitW$shape %in% c("two_rays", "whole_line") && fitW$F_CJ < fitW$crit,
   "weak design: unbounded set, F_CJ < crit")

for (f in list(fitA, fitB, fitW)) {
  ok(grid_agrees_fit(f),
     sprintf("conf_set == fine grid inversion (shape %s)", f$shape))
  # isTRUE() so an NA F_CJ fails this assertion loudly instead of opaquely.
  ok(isTRUE((f$F_CJ > f$crit) == all(is.finite(f$conf_set))),
     sprintf("F_CJ > crit iff set bounded (shape %s)", f$shape))
  ok(f$bounded == all(is.finite(f$conf_set)), "bounded flag matches endpoints")
  # 5.2 sanity: beta0 classified consistently with the set.
  ok((f$statistic <= f$crit) == in_set(f$beta0, f$conf_set),
     "T(beta0) <= crit iff beta0 in conf_set")
  # 5.2 sanity: if min Q < 0 the set is non-empty.
  if (f$coef_num[3] > 0 &&
      f$coef_num[1] - f$coef_num[2]^2 / (4 * f$coef_num[3]) < 0) {
    ok(nrow(f$conf_set) > 0L, "min Q < 0 implies a non-empty set")
  }
}

# === Test 2b: literal oracle and unsquared boundary on a covariate design ===
# C: weighted + dense control + fixed effect + factor instrument.  The dense
# equation oracle runs on the package-prepared data, so the whole covariate
# layer (sqrt-weight transform, FWL, FE absorption, reference coding) sits
# upstream of the equation check.
set.seed(93)
G_C <- 18L; ng_C <- 6L; n_C <- G_C * ng_C
cl_C <- rep(seq_len(G_C), each = ng_C)
jud_C <- factor(rep(1:5, length.out = n_C))
fe_C <- factor(rep_len(1:9, n_C))
w_C <- 0.5 + (seq_len(n_C) %% 7L) / 4
c_C <- sin(seq_len(n_C) * 0.21)
u_C <- rnorm(G_C)[cl_C]
x_C <- 0.9 * as.numeric(jud_C) + 0.5 * c_C + as.numeric(fe_C) / 3 +
  u_C + rnorm(n_C)
y_C <- 0.6 * x_C + 0.4 * c_C + as.numeric(fe_C) / 5 + u_C + rnorm(n_C)
d_C <- clusterIV:::.prep_data(y_C, x_C, jud_C, cl_C, cbind(c_C), w_C, TRUE,
                              fixed_effects = fe_C)
cf_C <- clusterIV:::.cjar_coefs(
  clusterIV:::.cluster_sums(d_C$y, d_C$x, d_C$Z, d_C$cluster))
for (b in betas1) {
  bf <- brute_QV(d_C$y, d_C$x, d_C$Z, d_C$cluster, b)
  ok(relnear(Qpoly(cf_C$n, b), bf$Q, 1e-10),
     sprintf("[covariate C] Q(%.1f) polynomial == dense brute force (1e-10 rel)", b))
  ok(relnear(Vpoly(cf_C$w, b), bf$V, 1e-10),
     sprintf("[covariate C] Vhat(%.1f) polynomial == 2/k sum c_gh^2 brute force (1e-10 rel)", b))
}
fitC <- cjar_q(y_C, x_C, jud_C, cluster = cl_C, controls = cbind(c_C),
               weights = w_C, fixed_effects = fe_C)
ok(grid_agrees_fit(fitC),
   sprintf("[covariate C] conf_set == fine grid inversion (shape %s)",
           fitC$shape))

# Every finite CJAR endpoint must solve the UNSQUARED boundary equation
# Q(b) = crit * sqrt(k * V(b)) or be a variance root: checking only the
# squared quartic would also accept sign-flipped spurious roots.
cjar_unsq <- function(f, b) {
  Q <- Qpoly(f$coef_num, b)
  V <- Vpoly(f$coef_var, b)
  vscale <- sum(abs(f$coef_var) * abs(b)^(0:4))
  vr <- abs(V) / max(vscale, .Machine$double.xmin)
  tr <- if (V > 0) {
    rhs <- f$crit * sqrt(f$k * V)
    abs(Q - rhs) / max(1, abs(Q) + abs(rhs))
  } else Inf
  min(vr, tr)
}
for (f in list(fitA, fitB, fitW, fitC)) {
  for (b in unique(as.numeric(f$conf_set[is.finite(f$conf_set)]))) {
    ok(cjar_unsq(f, b) <= 1e-8,
       sprintf("unsquared boundary residual at endpoint %.6g (shape %s)",
               b, f$shape))
  }
}

# === Test 3 (ORACLE R1): singletons == Mikusheva-Sun statistic, plain var ===
# Named oracle R1 (dev/extraction/mikusheva-sun-extraction.md, Section 5).
# With all clusters of size 1, Pddot is the diagonal-zeroed projection
# (C = P - D) and cjar() must equal a literal transcription of Mikusheva &
# Sun (2022, eq. 2) studentised by the plain ("naive") variance
# Phihat_1 = (2/k) sum_{i != j} (P_ij e_i e_j)^2 -- the "Mikusheva-Sun
# statistic with plain variance". (Never labelled Crudu-Mellace-Sandor:
# their test uses a different centring and is not this reduction.)
# Statistic and p-value at 1e-10; the numerator quadratic and the variance
# quartic are pinned pointwise at 5 betas, which determine both.
set.seed(31)
n3 <- 80L; k3 <- 4L
Z3 <- matrix(rnorm(n3 * k3), n3, k3)
v3 <- rnorm(n3)
x3 <- drop(Z3 %*% c(0.6, 0.4, -0.5, 0.3)) + v3 + rnorm(n3)
y3 <- 0.9 * x3 + v3 + rnorm(n3)
fit3 <- cjar_q(y3, x3, Z3, cluster = seq_len(n3), beta0 = 0.3)

po3 <- clusterIV:::.partial_out(y3, x3, Z3, NULL, NULL, TRUE)
P3 <- po3$Z %*% solve(crossprod(po3$Z), t(po3$Z))
diag(P3) <- 0                                   # C = P - D, MS eq. 2
jack_Q <- function(b) {
  e <- po3$y - b * po3$x
  drop(t(e) %*% P3 %*% e)
}
jack_V1 <- function(b) {                        # Phihat_1, MS Section 4.1
  e <- po3$y - b * po3$x
  (2 / k3) * sum((P3 * tcrossprod(e))^2)
}
jack_T <- function(b) jack_Q(b) / sqrt(k3 * jack_V1(b))
for (b in c(-2, -0.5, 0.3, 1, 3)) {
  ok(relnear(Qpoly(fit3$coef_num, b), jack_Q(b), 1e-10),
     sprintf("R1: Q(%.1f) == MS eq. 2 numerator e'(P - D)e (1e-10 rel)", b))
  ok(relnear(Vpoly(fit3$coef_var, b), jack_V1(b), 1e-10),
     sprintf("R1: Vhat(%.1f) == Phihat_1 plain variance (1e-10 rel)", b))
}
ok(relnear(fit3$statistic, jack_T(0.3), 1e-10),
   "R1: T(beta0) == Mikusheva-Sun statistic, plain variance (1e-10 rel)")
ok(relnear(fit3$p.value,
           pchisq(k3 + sqrt(2 * k3) * jack_T(0.3), df = k3, lower.tail = FALSE),
           1e-10), "R1: p-value == shifted-scaled chisq of the MS statistic")
bs3 <- seq(-4, 5, by = 0.01)
fe3 <- fit3$conf_set[is.finite(fit3$conf_set)]
bs3 <- bs3[vapply(bs3, function(b) min(abs(b - fe3)) > 0.01, logical(1))]
memb3 <- vapply(bs3, function(b) jack_T(b) <= fit3$crit, logical(1))
ok(all(memb3 == vapply(bs3, in_set, logical(1), cs = fit3$conf_set)),
   "R1: conf_set == grid-inverted jackknife AR")
# D4 honest label: at G == n the print header says independent data and
# names the reduction; the "cluster jackknife" wording must not appear.
out3 <- capture.output(print(fit3))
ok(any(grepl("Jackknife Anderson-Rubin test (independent data; Mikusheva-Sun statistic, plain variance)",
             out3, fixed = TRUE)) &&
   !any(grepl("Cluster jackknife", out3)),
   "G == n print label: independent-data jackknife AR header")

# === Test 3b: KS null calibration at singleton clusters =====================
# Mirrors the CJS null-calibration check in tests/test-cjscore.R: independent
# homoskedastic data, fixed k, statistic evaluated at the true beta. Under
# the null k + sqrt(2k) T(beta_true) -> chisq_k (Mikusheva & Sun 2022, the
# small-K regime after their Lemma 1) -- exactly the "chisq" calibration.
set.seed(20260716)
R_ks3 <- 2000L
n_K3 <- 200L; k_K3 <- 3L
t_ks3 <- numeric(R_ks3)
for (r in seq_len(R_ks3)) {
  Z_K3 <- matrix(rnorm(n_K3 * k_K3), n_K3, k_K3)
  v_K3 <- rnorm(n_K3)
  x_K3 <- drop(Z_K3 %*% c(1, 0.6, -0.5)) + v_K3
  y_K3 <- 0.5 * x_K3 + 0.7 * v_K3 + rnorm(n_K3)
  po_K3 <- clusterIV:::.partial_out(y_K3, x_K3, Z_K3, NULL, NULL, TRUE)
  cf_K3 <- clusterIV:::.cjar_coefs(
    clusterIV:::.cluster_sums(po_K3$y, po_K3$x, po_K3$Z, seq_len(n_K3)))
  t_ks3[r] <- Qpoly(cf_K3$n, 0.5) / sqrt(k_K3 * Vpoly(cf_K3$w, 0.5))
}
ks3 <- suppressWarnings(stats::ks.test(k_K3 + sqrt(2 * k_K3) * t_ks3,
                                       stats::pchisq, df = k_K3))
rej3 <- mean(t_ks3 > (qchisq(0.95, k_K3) - k_K3) / sqrt(2 * k_K3))
cat(sprintf("INFO: singleton KS D = %.4f (p = %.3f), null rejection at 5%% = %.3f\n",
            ks3$statistic, ks3$p.value, rej3))
ok(ks3$statistic < 0.05,
   "singleton null calibration: KS distance to chisq_k < 0.05")
ok(rej3 > 0.03 && rej3 < 0.07,
   "singleton null calibration: 5% one-sided rejection within [0.03, 0.07]")

# === Test 4: invariances ====================================================
# (a) affine: y -> y + c*x shifts every endpoint by +c, same shape.
fit4a <- cjar_q(y_A + 0.8 * x_A, x_A, jud_A, cluster = cl_A)
ok(fit4a$shape == fitA$shape &&
   set_near(fit4a$conf_set, fitA$conf_set + 0.8, 1e-6),
   "affine invariance: endpoints shift by +c, shape unchanged")
# (b) scale: x -> s*x scales endpoints by 1/s; T, F_CJ, p unchanged.
fit4b <- cjar_q(y_A, 3 * x_A, jud_A, cluster = cl_A)
ok(set_near(fit4b$conf_set, fitA$conf_set / 3, 1e-6) &&
   relnear(fit4b$statistic, fitA$statistic, 1e-10) &&
   relnear(fit4b$F_CJ, fitA$F_CJ, 1e-10) &&
   relnear(fit4b$p.value, fitA$p.value, 1e-10),
   "scale invariance: endpoints / s; T, F_CJ, p invariant")
# (c) Vhat(beta) >= 0 on a wide grid.
for (f in list(fitA, fitB, fitW, fit3)) {
  Vb <- Vpoly(f$coef_var, seq(-50, 50, by = 0.25))
  ok(min(Vb) >= -1e-8 * max(abs(Vb)), "Vhat >= 0 on a wide beta grid")
}
# (d) consistency at probe betas away from the boundary.
for (f in list(fitA, fitB, fitW)) {
  fe <- f$conf_set[is.finite(f$conf_set)]
  probes <- seq(-4, 4, by = 0.37)
  if (length(fe)) {
    probes <- probes[vapply(probes, function(b) min(abs(b - fe)) > 1e-3,
                            logical(1))]
  }
  hb <- hpoly(f$coef_num, f$coef_var, f$k, f$crit, probes)
  ok(all((hb <= 0) == vapply(probes, in_set, logical(1), cs = f$conf_set)),
     "probe betas: T(b) <= crit iff b in conf_set")
}

# === Test 5: covariate paths ================================================
# Dense controls == brute-force manual partialling.
set.seed(55)
W5 <- cbind(rnorm(n_B), rnorm(n_B))
fit5c <- cjar_q(y_B, x_B, Z_B, cluster = cl_B, controls = W5)
C5 <- cbind(1, W5)
res5 <- function(v) v - C5 %*% solve(crossprod(C5), crossprod(C5, v))
fit5m <- cjar_q(drop(res5(y_B)), drop(res5(x_B)), res5(Z_B), cluster = cl_B,
                intercept = FALSE)
ok(relnear(fit5c$statistic, fit5m$statistic, 1e-8) &&
   relnear(fit5c$p.value, fit5m$p.value, 1e-8) &&
   set_near(fit5c$conf_set, fit5m$conf_set, 1e-6),
   "dense controls == brute-force manual partialling")

# FE path == equivalent dummy-controls path.
set.seed(56)
f5 <- factor(sample(letters[1:8], n_B, replace = TRUE))
fit5f <- cjar_q(y_B, x_B, Z_B, cluster = cl_B, fixed_effects = f5)
fit5d <- cjar_q(y_B, x_B, Z_B, cluster = cl_B,
                controls = stats::model.matrix(~f5)[, -1L, drop = FALSE])
ok(relnear(fit5f$statistic, fit5d$statistic, 1e-8) &&
   relnear(fit5f$p.value, fit5d$p.value, 1e-8) &&
   set_near(fit5f$conf_set, fit5d$conf_set, 1e-6),
   "fixed-effects path == dense dummy-controls path")

# Weights == brute force on the sqrt(w)-transformed data.
set.seed(57)
w5 <- runif(n_B, 0.5, 2)
fit5w <- cjar_q(y_B, x_B, Z_B, cluster = cl_B, weights = w5)
rw5 <- sqrt(w5)
Cw5 <- matrix(rw5, ncol = 1L)
resw <- function(v) v - Cw5 %*% solve(crossprod(Cw5), crossprod(Cw5, v))
fit5t <- cjar_q(drop(resw(y_B * rw5)), drop(resw(x_B * rw5)), resw(Z_B * rw5),
                cluster = cl_B, intercept = FALSE)
ok(relnear(fit5w$statistic, fit5t$statistic, 1e-8) &&
   relnear(fit5w$p.value, fit5t$p.value, 1e-8) &&
   set_near(fit5w$conf_set, fit5t$conf_set, 1e-6),
   "weights == brute force on the rw-transformed model")

# === Test 6: formula interface ==============================================
datA <- data.frame(y = y_A, x = x_A, judge = jud_A, cl = cl_A)
f6 <- cjar_q(y ~ x | judge, data = datA, cluster = ~cl)
ok(identical(f6$statistic, fitA$statistic) &&
   identical(f6$p.value, fitA$p.value) &&
   set_near(f6$conf_set, fitA$conf_set, 1e-12),
   "formula interface == default interface (judge design)")
datB <- data.frame(y = y_B, x = x_B, z1 = Z_B[, 1], z2 = Z_B[, 2],
                   z3 = Z_B[, 3], f5 = f5, cl = cl_B)
f6b <- cjar_q(y ~ x | z1 + z2 + z3 | f5, data = datB, cluster = ~cl)
ok(identical(f6b$statistic, fit5f$statistic) &&
   set_near(f6b$conf_set, fit5f$conf_set, 1e-12),
   "three-part formula y ~ x | z | fe == fixed_effects argument")

# === Test 7: error paths and advisories =====================================
msg7 <- tryCatch(cjar_q(y_B, x_B, cbind(Z_B, Z_B[, 1]), cluster = cl_B),
                 error = function(e) conditionMessage(e))
ok(is.character(msg7) && grepl("singular", msg7),
   "collinear instruments hit the singular-Gram message")
nested_z <- factor(cl_A)
nested_D <- vapply(levels(nested_z)[-1L],
                   function(lev) as.numeric(nested_z == lev), numeric(n_A))
nested_factor <- cjar_q(y_A, x_A, nested_z, cluster = cl_A)
nested_matrix <- cjar_q(y_A, x_A, nested_D, cluster = cl_A)
ok(identical(nested_factor$statistic, nested_matrix$statistic) &&
   identical(nested_factor$p.value, nested_matrix$p.value) &&
   set_near(nested_factor$conf_set, nested_matrix$conf_set, 1e-12),
   "cluster-local factor levels equal their explicit encoded dummies")
fit04 <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, level = 0.4)
ok(identical(fit04$level, 0.4) && is.finite(fit04$crit),
   "accepts level = 0.4 (the contract is the open interval (0, 1))")
ok(errs(cjar(y_A, x_A, jud_A, cluster = cl_A, beta0 = c(0, 1))),
   "stops on non-scalar beta0")
ok(errs(cjar(y_A, x_A, jud_A, cluster = rep(1L, n_A))),
   "stops on a single cluster")
wmsg <- tryCatch(cjar(y_A, x_A, jud_A, cluster = cl_A),
                 warning = function(w) conditionMessage(w))
ok(is.character(wmsg) && grepl("asymptotic in the number of clusters", wmsg),
   "small-G advisory fires once for G < 20")
cl_dom <- ifelse(cl_B <= 8L, 1L, cl_B)
msgs <- character(0)
fit_dom <- withCallingHandlers(cjar(y_B, x_B, Z_B, cluster = cl_dom),
  warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
  })
ok(any(grepl("dominating cluster", msgs)), "dominating-cluster advisory fires")

# === Test 8: degenerate quartic (F_CJ on the boundary) ======================
# Nudge n2 so that F_CJ sits within ~1e-6 of the critical value; at eps = 0
# the leading quartic coefficient vanishes exactly and the trimmed cubic is
# the correct object. The inversion must return and match the grid.
k8 <- fitW$k; wc8 <- fitW$coef_var; crit8 <- fitW$crit
for (eps in c(3e-7, -3e-7)) {
  nc8 <- fitW$coef_num
  nc8[3] <- crit8 * sqrt(k8 * wc8[5]) * (1 + eps)
  inv8 <- clusterIV:::.cjar_invert(nc8, wc8, k8, crit8)
  ok(grid_agrees(nc8, wc8, k8, crit8, inv8$conf_set, lo = -20, hi = 20),
     sprintf("degenerate quartic (eps = %g): inversion matches grid", eps))
}
# At eps = 0 whether the leading coefficient cancels to exactly zero is
# platform arithmetic (exact on R 4.3/4.4, a few ulps off on R 4.5), so
# the exact case and +-1..4 ulp nudges are all checked.  A leading
# coefficient at roundoff puts the far root near 1/roundoff, where neither
# this grid nor the inversion is a reference; the check is therefore on
# the moderate region: grid membership on [-20, 20] and the moderate
# endpoint against a direct root of h at 1e-10.
base8 <- crit8 * sqrt(k8 * wc8[5])
nc8 <- fitW$coef_num
nc8[3] <- base8
h8 <- function(b) hpoly(nc8, wc8, k8, crit8, b)
b8 <- uniroot(h8, c(0.4, 0.6), tol = 1e-14)$root
bm8 <- seq(-20, 20, length.out = 20001L)
bm8 <- bm8[abs(bm8 - b8) > 1e-6]
for (u in -4:4) {
  nc8[3] <- base8 * (1 + u * .Machine$double.eps / 2)
  cs8 <- clusterIV:::.cjar_invert(nc8, wc8, k8, crit8)$conf_set
  fe8 <- cs8[is.finite(cs8) & abs(cs8) <= 20]
  ok(length(fe8) == 1L && abs(fe8 - b8) < 1e-10 &&
       all((h8(bm8) <= 0) == vapply(bm8, in_set, logical(1), cs = cs8)),
     sprintf("degenerate quartic (eps = 0, %+d ulp): moderate region exact", u))
}

# === Test 9: methods ========================================================
ok(nobs(fitA) == n_A, "nobs.cjar returns n")
ok(identical(confint(fitA), fitA$conf_set),
   "confint at the fitted level returns the stored set")
ci90 <- confint(fitA, level = 0.90)
dir90 <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, level = 0.90)$conf_set
ok(set_near(ci90, dir90, 1e-6),
   "confint at a new level == direct refit at that level")
ci30 <- suppressWarnings(confint(fitA, level = 0.30))
dir30 <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, level = 0.30)$conf_set
ok(set_near(ci30, dir30, 1e-6),
   "confint accepts level = 0.3 and agrees with a direct refit")
out9 <- capture.output({ print(fitA); print(fitW) })
ok(any(grepl("Cluster jackknife Anderson-Rubin", out9)) &&
   any(grepl("unbounded", out9)) && any(grepl("F_CJ =", out9, fixed = TRUE)),
   "print method: header, endpoint-derived topology and F_CJ diagnostic render")

# === Test A: empty set ======================================================
# Q(b) = 1 + b^2 > 0 everywhere and a tiny variance: since
# sqrt(1 + b^2 + b^4) <= 1 + b^2, the boundary c*sqrt(kV) stays below Q for
# every beta and the set is empty. Empty occurs only in the bounded regime.
ncA <- c(1, 0, 1); wcA <- c(1e-4, 0, 1e-4, 0, 1e-4)
invA <- clusterIV:::.cjar_invert(ncA, wcA, 5, 1.96)
ok(invA$shape == "empty" && nrow(invA$conf_set) == 0L,
   "synthetic empty set: shape and 0-row matrix")
ok(grid_agrees(ncA, wcA, 5, 1.96, invA$conf_set),
   "synthetic empty set: grid agreement")
ok(isTRUE(invA$F_CJ > 1.96), "empty set occurs in the bounded regime (F_CJ > crit)")

# === Test B: two disjoint bounded intervals =================================
# Constructed inversely: with Q(b) = 1 + b^2 and W(b) = c^2 k V(b) - Q(b)^2
# = -0.1 (b^2 - 9)(b^2 - 1), acceptance (W >= 0, Q > 0) is exactly
# [-3, -1] U [1, 3]; V = Q^2 + W = 0.9 b^4 + 3 b^2 + 0.1 >= 0 is a valid
# variance polynomial and the tails are rejected (n2 = 1 > sqrt(w4)).
ncB <- c(1, 0, 1); wcB <- c(0.1, 0, 3, 0, 0.9)
invB <- clusterIV:::.cjar_invert(ncB, wcB, 1, 1)
ok(invB$shape == "bounded" && nrow(invB$conf_set) == 2L &&
   all(is.finite(invB$conf_set)),
   "two-interval quartic: 2 rows, all finite, shape bounded")
ok(near(invB$conf_set, rbind(c(-3, -1), c(1, 3)), 1e-6),
   "two-interval quartic: endpoints == analytic values (-3, -1, 1, 3)")
ok(grid_agrees(ncB, wcB, 1, 1, invB$conf_set),
   "two-interval quartic: grid agreement")

# === Test C: w4 = 0 degenerate family =======================================
# w4 = 0 forces n2 = 0 (all d_gh vanish), so Q is linear and the tail signs
# can differ; F_CJ is NA by convention in this family.
# (i) one tail accepted -> single ray [-1/sqrt(0.9), Inf).
ncC1 <- c(0, 1, 0); wcC1 <- c(1, 0, 0.1, 0, 0)
invC1 <- clusterIV:::.cjar_invert(ncC1, wcC1, 1, 1)
ok(invC1$shape == "ray" && nrow(invC1$conf_set) == 1L &&
   is.finite(invC1$conf_set[1, 1]) && is.infinite(invC1$conf_set[1, 2]),
   "w4 = 0, differing tails: single ray with one finite endpoint")
ok(near(invC1$conf_set[1, 1], -1 / sqrt(0.9), 1e-6),
   "w4 = 0 ray: finite endpoint == analytic value")
ok(is.na(invC1$F_CJ), "w4 = 0 ray: F_CJ is NA")
ok(grid_agrees(ncC1, wcC1, 1, 1, invC1$conf_set), "w4 = 0 ray: grid agreement")
# (ii) both tails accepted, middle rejected -> two rays.
ncC2 <- c(1, 0, 0); wcC2 <- c(0.01, 0, 1, 0, 0)
invC2 <- clusterIV:::.cjar_invert(ncC2, wcC2, 1, 1)
ok(invC2$shape == "two_rays" && is.na(invC2$F_CJ),
   "w4 = 0, both tails accepted: two rays, F_CJ NA")
ok(grid_agrees(ncC2, wcC2, 1, 1, invC2$conf_set),
   "w4 = 0 two rays: grid agreement")
# (iii) h < 0 everywhere -> whole line (also exercises the spurious
# double root at Q = -c*sqrt(kV)).
ncC3 <- c(-1, 0, 0); wcC3 <- c(1, 0, 0.1, 0, 0)
invC3 <- clusterIV:::.cjar_invert(ncC3, wcC3, 1, 1)
ok(invC3$shape == "whole_line" && is.na(invC3$F_CJ),
   "w4 = 0, always accepted: whole line, F_CJ NA")
ok(grid_agrees(ncC3, wcC3, 1, 1, invC3$conf_set),
   "w4 = 0 whole line: grid agreement")

# === Test D: combined weights + fixed effects ===============================
fitDW <- cjar_q(y_B, x_B, Z_B, cluster = cl_B, weights = w5,
                fixed_effects = f5)
DD <- stats::model.matrix(~f5)[, -1L, drop = FALSE]
CwD <- cbind(rw5, DD * rw5)
resD <- function(v) v - CwD %*% solve(crossprod(CwD), crossprod(CwD, v))
fitDT <- cjar_q(drop(resD(y_B * rw5)), drop(resD(x_B * rw5)),
                resD(Z_B * rw5), cluster = cl_B, intercept = FALSE)
ok(relnear(fitDW$statistic, fitDT$statistic, 1e-8) &&
   relnear(fitDW$p.value, fitDT$p.value, 1e-8) &&
   set_near(fitDW$conf_set, fitDT$conf_set, 1e-6),
   "weights + fixed effects == brute force on rw-transformed dummy route")

# === Test E: normal calibration through the inversion =======================
fitE <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, calibration = "normal")
ok(near(fitE$crit, qnorm(0.95), 1e-12), "normal calibration: crit == qnorm(level)")
ok(grid_agrees_fit(fitE), "normal calibration: conf_set == grid inversion")
ciE <- confint(fitE, level = 0.90)
dirE <- cjar_q(y_A, x_A, jud_A, cluster = cl_A, calibration = "normal",
               level = 0.90)$conf_set
ok(set_near(ciE, dirE, 1e-6),
   "normal calibration: confint at a new level == direct refit")

# === Test F: near-0.5 level, negative critical value ========================
# With k = 5, F_chisq5(5) > 0.51, so level = 0.51 passes the guard yet gives
# crit < 0; the inversion must survive and the F_CJ criterion must hold.
msgsF <- character(0)
fitF <- withCallingHandlers(
  cjar(y_A, x_A, jud_A, cluster = cl_A, level = 0.51),
  warning = function(w) {
    msgsF <<- c(msgsF, conditionMessage(w)); invokeRestart("muffleWarning")
  })
ok(any(grepl("non-positive critical value", msgsF)),
   "level = 0.51 advisory (non-positive critical value) fires")
ok(fitF$crit < 0, "level = 0.51 with k = 5 yields crit < 0")
ok(grid_agrees_fit(fitF), "crit < 0: conf_set == grid inversion")
ok(isTRUE((fitF$F_CJ > fitF$crit) == all(is.finite(fitF$conf_set))),
   "crit < 0: F_CJ > crit iff set bounded")
ok(isTRUE((fitF$statistic <= fitF$crit) == in_set(fitF$beta0, fitF$conf_set)),
   "crit < 0: beta0 classified consistently with the set")

# === Test G: maxlev warning =================================================
# One instrument lives (mean-zero) entirely inside cluster 1, so that cluster
# nearly spans the instrument space: maxlev > 0.99 and the A2(ii) advisory
# must fire.
set.seed(88)
G_G <- 20L; ng_G <- 5L; n_G <- G_G * ng_G
cl_G <- rep(seq_len(G_G), each = ng_G)
zloc <- numeric(n_G)
zloc[cl_G == 1L] <- c(1.5, -1.5, 1, -1, 0)
Z_G <- cbind(matrix(rnorm(n_G * 2L), n_G, 2L), zloc)
x_G <- drop(Z_G %*% c(1, -1, 0.5)) + rnorm(n_G)
y_G <- x_G + rnorm(n_G)
msgsG <- character(0)
fitG <- withCallingHandlers(cjar(y_G, x_G, Z_G, cluster = cl_G),
  warning = function(w) {
    msgsG <<- c(msgsG, conditionMessage(w)); invokeRestart("muffleWarning")
  })
ok(fitG$maxlev > 0.99, "constructed design reaches maxlev > 0.99")
ok(any(grepl("nearly spans the instrument space", msgsG)),
   "maxlev advisory (Assumption A2(ii)) fires")

# === Test H: zero instrument columns ========================================
# A single-level factor z expands to a zero-column Z; the front end must stop
# before chol() sees a 0 x 0 Gram matrix.
msgH <- tryCatch(cjar_q(y_A, x_A, factor(rep("only", n_A)), cluster = cl_A),
                 error = function(e) conditionMessage(e))
ok(is.character(msgH) && grepl("no instrument columns", msgH),
   "single-level factor z stops: no instrument columns")

# === Test I: isolated accepted points (numerical-kernel unit tests) =========
# These exercise .cjar_invert on synthetic coefficient tuples; whether such a
# tuple is reachable from data is a separate question.
# (i) Tangency from above: constant variance w = w0, upward parabola Q with
# min Q = c_a * sqrt(k * w0) exactly (all quantities exact in doubles):
# k = 4, w0 = 1, c_a = 1.5 -> boundary at 3; Q = 4 - 2b + b^2 has min 3 at
# b* = n1 / (2 n2) = 1. h(b) = (b - 1)^2 touches zero from above at b*, an
# isolated accepted point.
invI1 <- clusterIV:::.cjar_invert(c(4, 2, 1), c(1, 0, 0, 0, 0), 4, 1.5)
ok(invI1$shape == "bounded" && nrow(invI1$conf_set) == 1L &&
   near(invI1$conf_set, cbind(1, 1), 1e-8),
   "tangency from above: single degenerate row [b*, b*], shape bounded")
# (ii) Near-miss control: min Q strictly above the boundary -> empty set,
# no spurious point introduced.
invI2 <- clusterIV:::.cjar_invert(c(4 + 1e-3, 2, 1), c(1, 0, 0, 0, 0), 4, 1.5)
ok(invI2$shape == "empty" && nrow(invI2$conf_set) == 0L,
   "near-miss tangency: empty set, no spurious degenerate point")
# (iii) Regression guard: a non-degenerate bounded geometry is unchanged by
# the degenerate-point fix (frozen endpoints from the pre-fix implementation).
ok(near(fitA$conf_set,
        cbind(0.54589321637897903, 0.90634837758600784), 1e-8),
   "fitA conf_set identical to its pre-fix value (frozen endpoints)")

cat("\nAll cjar tests passed.\n")
