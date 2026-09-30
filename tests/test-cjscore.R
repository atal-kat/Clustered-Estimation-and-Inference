# Correctness tests for cjscore() (base R, no testthat), in the style of
# test-cjar.R. The oracle below was written and committed BEFORE the fast
# path (protocol: dense-oracle-first, with the derived identities
# spot-checked in dev/spot-check-cjs.R); it is never edited to match the
# fast path. If the two disagree, the fast path is wrong until proven
# otherwise against the paper.

library(clusterIV)

# ============================================================================
# ORACLE -- frozen. Literal transcription of Ligtenberg (2025,
# arXiv:2306.08559v3), Section 4, for p = 1:
#   statistic  S_CLJ(beta) = X' Pddot eps(beta) / sqrt(n)          [Section 4]
#   variance   Vhat^S(beta) = (1/n) [ X' Pddot B_{eps eps'} Pddot X
#                + sum_{g != h} X_[g]' P_[g,h] eps_[h]
#                              * eps_[g]' P_[g,h] X_[h] ]          [eq. (6)]
# The dense n x n projector is formed explicitly, its diagonal cluster
# blocks are zeroed explicitly, B_{eps eps'} is assembled block by block,
# and the cross term is a literal double loop over (g, h) with BOTH factors
# on the SAME [g, h] block, exactly as printed in the paper.
# ============================================================================
.cjs_brute <- function(y, x, Z, cl, beta) {
  n <- length(y)
  cl <- droplevels(as.factor(cl))
  idx <- split(seq_len(n), cl)
  G <- length(idx)
  eps <- y - beta * x
  P <- Z %*% solve(crossprod(Z), t(Z))
  Pdd <- P
  for (g in seq_len(G)) Pdd[idx[[g]], idx[[g]]] <- 0
  S <- drop(t(x) %*% Pdd %*% eps) / sqrt(n)
  # term 1: X' Pddot B_{eps(beta) eps(beta)'} Pddot X
  Bee <- matrix(0, n, n)
  for (g in seq_len(G)) {
    ig <- idx[[g]]
    Bee[ig, ig] <- eps[ig] %*% t(eps[ig])
  }
  t1 <- drop(t(x) %*% Pdd %*% Bee %*% Pdd %*% x)
  # term 2: sum_{g != h} (x_g' P_[g,h] eps_h) * (eps_g' P_[g,h] x_h)
  t2 <- 0
  for (g in seq_len(G)) for (h in seq_len(G)) if (g != h) {
    ig <- idx[[g]]; ih <- idx[[h]]
    Pgh <- P[ig, ih, drop = FALSE]
    t2 <- t2 + drop(t(x[ig]) %*% Pgh %*% eps[ih]) *
               drop(t(eps[ig]) %*% Pgh %*% x[ih])
  }
  V <- (t1 + t2) / n
  LM <- if (V > 0) S^2 / V else if (S == 0 && V == 0) 0 else Inf
  list(S = S, V = V, LM = LM)
}
# ======================= end of frozen oracle ===============================
# NOTE: the oracle's LM line predates the conservative non-positive-variance
# convention (Vhat <= 0 now accepts, statistic NA); its S and V remain the
# gates. LM comparisons against the oracle are made only where V > 0.

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}
near <- function(a, b, tol) max(abs(a - b)) < tol
relnear <- function(a, b, tol) max(abs(a - b) / pmax(1, abs(b))) < tol
errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
cjs_q <- function(...) suppressWarnings(cjscore(...))

Spoly <- function(sc, b) sc[1] - sc[2] * b
Vpoly <- function(vc, b) vc[1] + vc[2] * b + vc[3] * b^2
in_set <- function(b, cs) nrow(cs) > 0 && any(b >= cs[, 1] & b <= cs[, 2])
set_near <- function(A, B, tol) {
  nrow(A) == nrow(B) && all(is.infinite(A) == is.infinite(B)) &&
    (sum(is.finite(A)) == 0L ||
       max(abs(A[is.finite(A)] - B[is.finite(B)])) < tol)
}
# Dense grid inversion under the conservative non-positive-variance
# convention: reject b iff v(b) > 0 AND A(b) > 0. This oracle is the gate for
# the analytic inversion; it was rewritten with the convention change (the
# pre-change kernel disagrees with it wherever v(b) <= 0).
grid_agrees_s <- function(sc, vc, crit, cs, lo = -8, hi = 8, step = 0.01) {
  fe <- cs[is.finite(cs)]
  if (length(fe)) { lo <- min(lo, min(fe) - 2); hi <- max(hi, max(fe) + 2) }
  bs <- seq(lo, hi, by = step)
  if (length(fe)) {
    bs <- bs[vapply(bs, function(b) min(abs(b - fe)) > step, logical(1))]
  }
  vb <- Vpoly(vc, bs)
  Ab <- Spoly(sc, bs)^2 - crit * vb
  rej <- vb > 0 & Ab > 0
  ms <- vapply(bs, in_set, logical(1), cs = cs)
  all((!rej) == ms)
}
grid_agrees_fit <- function(f, ...)
  grid_agrees_s(f$coef_score, f$coef_var, f$crit, f$conf_set, ...)

# --- designs -----------------------------------------------------------------
# A: judge design, n = 90, G = 15 balanced, k = 5.
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

# Weak variant of A.
set.seed(41)
x_Aw <- 0.02 * as.numeric(jud_A) + u_A + rnorm(n_A)
y_Aw <- 0.5 * x_Aw + u_A + rnorm(n_A)

# Singletons: every cluster of size 1.
set.seed(31)
n_S <- 60L; k_S <- 4L
Z_S <- matrix(rnorm(n_S * k_S), n_S, k_S)
v_S <- rnorm(n_S)
x_S <- drop(Z_S %*% c(0.6, 0.4, -0.5, 0.3)) + v_S + rnorm(n_S)
y_S <- 0.9 * x_S + v_S + rnorm(n_S)

# === T-S1: oracle agreement battery (1e-10 relative) ========================
# The fast-path coefficient vectors are pinned by matching the oracle's S and
# V pointwise at four betas (a linear/quadratic is determined by 2/3 points).
battery <- list(
  list(tag = "judge A (balanced, k = 5)",
       y = y_A, x = x_A, Z = Zj_A, cl = cl_A),
  list(tag = "matrix B (unbalanced, k = 3)",
       y = y_B, x = x_B, Z = Z_B, cl = cl_B),
  list(tag = "weak A (first stage ~ 0)",
       y = y_Aw, x = x_Aw, Z = Zj_A, cl = cl_A),
  # At singleton clusters the frozen oracle IS a literal transcription of
  # Matsushita & Otsu (2024, eqs. 3-4): every diagonal cluster block is a
  # diagonal entry, so Pdd = P*, Bee = Sigma_0, and the pair term is their
  # second sum -- identity R3, confirmed term by term in
  # dev/extraction/matsushita-otsu-extraction.md, Section 3. This battery
  # row is therefore the named R3 oracle: cjscore(plain) at singletons ==
  # the Matsushita-Otsu jackknife LM test.
  list(tag = "singleton clusters (n = G = 60) [R3: Matsushita-Otsu JLM]",
       y = y_S, x = x_S, Z = Z_S, cl = seq_len(n_S))
)
betasT <- c(-2, 0, 0.7, 3)
for (des in battery) {
  po <- clusterIV:::.partial_out(des$y, des$x, des$Z, NULL, NULL, TRUE)
  cf <- clusterIV:::.cjs_coefs(clusterIV:::.cluster_sums(po$y, po$x, po$Z, des$cl))
  nn <- length(des$y)
  for (b in betasT) {
    br <- .cjs_brute(po$y, po$x, po$Z, des$cl, b)
    ok(relnear(Spoly(cf$s, b) / sqrt(nn), br$S, 1e-10),
       sprintf("[%s] S(%.1f) == oracle", des$tag, b))
    ok(relnear(Vpoly(cf$v, b) / nn, br$V, 1e-10),
       sprintf("[%s] Vhat(%.1f) == oracle", des$tag, b))
  }
  fit <- cjs_q(des$y, des$x, if (identical(des$Z, Zj_A)) jud_A else des$Z,
               cluster = des$cl, beta0 = 0.7)
  br <- .cjs_brute(po$y, po$x, po$Z, des$cl, 0.7)
  ok(relnear(fit$statistic, br$LM, 1e-10) &&
     relnear(fit$score, br$S, 1e-10) &&
     relnear(fit$variance, br$V, 1e-10),
     sprintf("[%s] cjscore() statistic/score/variance == oracle", des$tag))
  ok(relnear(fit$p.value, pchisq(br$LM, 1, lower.tail = FALSE), 1e-10),
     sprintf("[%s] p-value == chisq_1 of oracle LM", des$tag))
}

# Weights: fast path == oracle on the sqrt(w)-transformed model.
set.seed(57)
w5 <- runif(n_B, 0.5, 2)
rw5 <- sqrt(w5)
Cw5 <- matrix(rw5, ncol = 1L)
resw <- function(v) v - Cw5 %*% solve(crossprod(Cw5), crossprod(Cw5, v))
fitW <- cjs_q(y_B, x_B, Z_B, cluster = cl_B, weights = w5, beta0 = 0.4)
brW <- .cjs_brute(drop(resw(y_B * rw5)), drop(resw(x_B * rw5)),
                  resw(Z_B * rw5), cl_B, 0.4)
ok(relnear(fitW$statistic, brW$LM, 1e-10) &&
   relnear(fitW$variance, brW$V, 1e-10),
   "[weights] cjscore() == oracle on rw-transformed data")

# Fixed effects: fast path == oracle on the dummy-partialled data.
set.seed(56)
f5 <- factor(sample(letters[1:8], n_B, replace = TRUE))
C5 <- cbind(1, stats::model.matrix(~f5)[, -1L, drop = FALSE])
res5 <- function(v) v - C5 %*% solve(crossprod(C5), crossprod(C5, v))
fitF <- cjs_q(y_B, x_B, Z_B, cluster = cl_B, fixed_effects = f5, beta0 = 0.4)
brF <- .cjs_brute(drop(res5(y_B)), drop(res5(x_B)), res5(Z_B), cl_B, 0.4)
ok(relnear(fitF$statistic, brF$LM, 1e-10) &&
   relnear(fitF$variance, brF$V, 1e-10),
   "[fixed effects] cjscore() == oracle on dummy-partialled data")

# Weights + fixed effects combined (T-S1 requirement).
Cwf <- cbind(rw5, stats::model.matrix(~f5)[, -1L, drop = FALSE] * rw5)
reswf <- function(v) v - Cwf %*% solve(crossprod(Cwf), crossprod(Cwf, v))
fitWF <- cjs_q(y_B, x_B, Z_B, cluster = cl_B, weights = w5,
               fixed_effects = f5, beta0 = 0.4)
brWF <- .cjs_brute(drop(reswf(y_B * rw5)), drop(reswf(x_B * rw5)),
                   reswf(Z_B * rw5), cl_B, 0.4)
ok(relnear(fitWF$statistic, brWF$LM, 1e-10) &&
   relnear(fitWF$variance, brWF$V, 1e-10),
   "[weights + fixed effects] cjscore() == oracle, combined path")

# Full covariate design: weighted + dense control + fixed effect + factor
# instrument (mirrors test-cjar.R Test 2b).  The oracle runs on the
# package-prepared data, so every covariate-handling layer sits upstream of
# the equation check.
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
cf_C <- clusterIV:::.cjs_coefs(
  clusterIV:::.cluster_sums(d_C$y, d_C$x, d_C$Z, d_C$cluster))
for (b in betasT) {
  br <- .cjs_brute(d_C$y, d_C$x, d_C$Z, d_C$cluster, b)
  ok(relnear(Spoly(cf_C$s, b) / sqrt(n_C), br$S, 1e-10),
     sprintf("[covariate C] S(%.1f) == oracle on prepared data", b))
  ok(relnear(Vpoly(cf_C$v, b) / n_C, br$V, 1e-10),
     sprintf("[covariate C] Vhat(%.1f) == oracle on prepared data", b))
}
fitC <- cjs_q(y_C, x_C, jud_C, cluster = cl_C, controls = cbind(c_C),
              weights = w_C, fixed_effects = fe_C, beta0 = 0.7)
brC <- .cjs_brute(d_C$y, d_C$x, d_C$Z, d_C$cluster, 0.7)
ok(relnear(fitC$statistic, brC$LM, 1e-10) &&
   relnear(fitC$score, brC$S, 1e-10) &&
   relnear(fitC$variance, brC$V, 1e-10),
   "[covariate C] cjscore() statistic/score/variance == oracle")

# === T-S2: analytic inversion with prescribed roots =========================
# A(beta) = (s0 - s1 b)^2 - crit v(b) constructed so that A = (b-1)(b-3):
# sc = (0, 2), crit = 4, v = (4b^2 - (b^2-4b+3))/4 -> vc = (-3/4, 1, 3/4).
# Under the conservative rule the acceptance set is {A <= 0} U {v <= 0}:
# [1, 3] plus the closed window between the roots of v.
vrP <- sort((-1 + c(-1, 1) * sqrt(3.25)) / 1.5)
invP <- clusterIV:::.cjs_invert(c(0, 2), c(-0.75, 1, 0.75), 4)
ok(invP$shape == "bounded" && nrow(invP$conf_set) == 2L &&
   near(invP$conf_set, rbind(vrP, c(1, 3)), 1e-8),
   "prescribed roots: acceptance [v-window] U [1, 3] recovered to 1e-8")
ok(grid_agrees_s(c(0, 2), c(-0.75, 1, 0.75), 4, invP$conf_set),
   "prescribed roots: grid agreement")
# Mirror construction with A = -(b-1)(b-3): vc = (3/4, -1, 5/4) -> two rays.
invP2 <- clusterIV:::.cjs_invert(c(0, 2), c(0.75, -1, 1.25), 4)
ok(invP2$shape == "two_rays" &&
   set_near(invP2$conf_set, rbind(c(-Inf, 1), c(3, Inf)), 1e-8),
   "prescribed roots: complement of (1, 3) recovered to 1e-8")
# Dataset-based endpoint recovery: endpoints solve the stored quadratic.
fitA <- cjs_q(y_A, x_A, jud_A, cluster = cl_A)
Aq <- c(fitA$coef_score[1]^2 - fitA$crit * fitA$coef_var[1],
        -2 * prod(fitA$coef_score) - fitA$crit * fitA$coef_var[2],
        fitA$coef_score[2]^2 - fitA$crit * fitA$coef_var[3])
rts <- sort(Re(polyroot(Aq)))
ok(near(fitA$conf_set, matrix(rts, 1L, 2L), 1e-8),
   "dataset endpoints == independent quadratic-formula roots (1e-8)")

# === T-S3: shape classification, isTRUE()-hardened, NA-safe =================
fitB <- cjs_q(y_B, x_B, Z_B, cluster = cl_B)
fitWk <- cjs_q(y_Aw, x_Aw, jud_A, cluster = cl_A)
ok(isTRUE(fitA$shape == "bounded" && fitA$F_CJS^2 > fitA$crit),
   "strong judge design: bounded, F_CJS^2 > crit")
ok(isTRUE(fitB$shape == "bounded" && fitB$F_CJS^2 > fitB$crit),
   "strong continuous design: bounded, F_CJS^2 > crit")
ok(isTRUE(fitWk$shape %in% c("two_rays", "whole_line") &&
          fitWk$F_CJS^2 < fitWk$crit),
   "weak design: unbounded, F_CJS^2 < crit")
for (f in list(fitA, fitB, fitWk, fitW, fitF, fitWF)) {
  ok(grid_agrees_fit(f), sprintf("conf_set == grid (shape %s)", f$shape))
  ok(isTRUE((f$F_CJS^2 > f$crit) == all(is.finite(f$conf_set))),
     "F_CJS^2 > crit iff set bounded (NA-safe)")
  # The invariant the pre-change code violated: beta0 is in the set iff the
  # statistic is undefined (NA, accepted) or does not exceed the critical
  # value.
  ok(isTRUE((is.na(f$statistic) || f$statistic <= f$crit) ==
              in_set(f$beta0, f$conf_set)),
     "beta0 classified consistently with the set (NA accepts)")
}
# Synthetic degenerate family (numerical-kernel unit tests; whether such a
# tuple is reachable from data is a separate question):
# v = -1 + 0.01 b^2 <= 0 exactly on [-10, 10] while A = 1 + 0.99 b^2 > 0
# everywhere: the conservative rule accepts exactly the closed interval where
# the variance estimate is non-positive.
invE <- clusterIV:::.cjs_invert(c(0, 1), c(-1, 0, 0.01), 1)
ok(isTRUE(invE$shape == "bounded" && nrow(invE$conf_set) == 1L &&
          near(invE$conf_set, cbind(-10, 10), 1e-8) && invE$F_CJS^2 > 1),
   "v <= 0 window: conf_set == the closed interval where v <= 0")
ok(grid_agrees_s(c(0, 1), c(-1, 0, 0.01), 1, invE$conf_set),
   "v <= 0 window: grid agreement")
invL <- clusterIV:::.cjs_invert(c(0, 0), c(1, 0, 1), 4)
ok(isTRUE(invL$shape == "whole_line" && invL$F_CJS == 0),
   "S identically 0 with positive variance: whole line")
invR <- clusterIV:::.cjs_invert(c(1, 1), c(0, 0, 1), 1)
ok(isTRUE(invR$shape == "ray" && nrow(invR$conf_set) == 2L &&
          near(invR$conf_set[1L, ], c(0, 0), 1e-8) &&
          is.infinite(invR$conf_set[2L, 2L]) &&
          abs(invR$conf_set[2L, 1L] - 0.5) < 1e-8),
   "degenerate leading coefficient: isolated point {0} (v = 0 there) plus ray [0.5, Inf)")
ok(grid_agrees_s(c(1, 1), c(0, 0, 1), 1, invR$conf_set),
   "ray: grid agreement")
invN <- clusterIV:::.cjs_invert(c(1, 1), c(1, 0, 0), 1)
ok(isTRUE(is.na(invN$F_CJS) && invN$shape == "bounded" &&
          near(invN$conf_set, cbind(0, 2), 1e-8)),
   "v2 = 0: F_CJS is NA, bounded set [0, 2] still exact")
# Tangency from above (twin of the cjar kernel test): A = (b - 2)^2 with
# A2 > 0 and v(2) = 1 > 0, from sc = (1, 0), vc = (-3, 4, -1), crit = 1
# (v2 < 0 is legal at the kernel level; .cjs_coefs guarantees v2 >= 0).
# v > 0 only on (1, 3), so the rule accepts (-Inf, 1] U {2} U [3, Inf):
# the tangency point 2 survives as a single degenerate row.
invT <- clusterIV:::.cjs_invert(c(1, 0), c(-3, 4, -1), 1)
ok(isTRUE(invT$shape == "two_rays" && nrow(invT$conf_set) == 3L &&
          is.infinite(invT$conf_set[1L, 1L]) &&
          abs(invT$conf_set[1L, 2L] - 1) < 1e-8 &&
          near(invT$conf_set[2L, ], c(2, 2), 1e-6) &&
          abs(invT$conf_set[3L, 1L] - 3) < 1e-8 &&
          is.infinite(invT$conf_set[3L, 2L])),
   "tangency from above: degenerate row [2, 2] kept between the rays")
ok(grid_agrees_s(c(1, 0), c(-3, 4, -1), 1, invT$conf_set),
   "tangency from above: grid agreement")

# === T-S4: cross-statistic coherence ========================================
est4 <- cjive(y_A, x_A, jud_A, cluster = cl_A)
ar4 <- suppressWarnings(cjar(y_A, x_A, jud_A, cluster = cl_A))
ok(est4$conf.low <= est4$coefficient && est4$coefficient <= est4$conf.high &&
   in_set(est4$coefficient, ar4$conf_set) &&
   in_set(est4$coefficient, fitA$conf_set),
   "CJIVE point estimate covered by Wald, CJAR and CJS 95% sets")
bogus <- est4$coefficient + 10
ar_rej <- suppressWarnings(cjar(y_A, x_A, jud_A, cluster = cl_A,
                                beta0 = bogus))
sc_rej <- cjs_q(y_A, x_A, jud_A, cluster = cl_A, beta0 = bogus)
ok(ar_rej$p.value < 0.01 && sc_rej$p.value < 0.01,
   "CJAR and CJS agree on rejecting a grossly false null")

# === T-S5: the non-positive-variance path is exercised ======================
# Seed 6 (found by search) yields min_beta Vhat(beta) < 0 on a small weak
# design; at the variance minimiser the LM statistic is undefined and the
# conservative convention must ACCEPT, in the statistic and in the inversion
# alike.
set.seed(6)
G_N <- 5L; sz_N <- sample(3:6, G_N, replace = TRUE); n_N <- sum(sz_N)
cl_N <- rep(seq_len(G_N), sz_N)
Z_N <- matrix(rnorm(n_N * 3L), n_N, 3L)
u_N <- rnorm(G_N)[cl_N]
x_N <- 0.05 * Z_N[, 1L] + u_N + rnorm(n_N)
y_N <- u_N + rnorm(n_N)
po_N <- clusterIV:::.partial_out(y_N, x_N, Z_N, NULL, NULL, TRUE)
cf_N <- clusterIV:::.cjs_coefs(clusterIV:::.cluster_sums(po_N$y, po_N$x, po_N$Z, cl_N))
bstar <- -cf_N$v[2L] / (2 * cf_N$v[3L])
ok(Vpoly(cf_N$v, bstar) < 0, "seed-6 design: Vhat(beta*) < 0 is reached")
fit_N <- cjs_q(y_N, x_N, Z_N, cluster = cl_N, beta0 = bstar)
br_N <- .cjs_brute(po_N$y, po_N$x, po_N$Z, cl_N, bstar)
ok(relnear(fit_N$variance, br_N$V, 1e-10) && fit_N$variance < 0,
   "Vhat(beta0) < 0 confirmed against the oracle's V")
ok(is.na(fit_N$statistic) && fit_N$p.value == 1,
   "Vhat(beta0) < 0: statistic NA, p = 1 (conservative convention)")
ok(isTRUE(in_set(bstar, fit_N$conf_set)),
   "Vhat(beta0) < 0: beta0 contained in the confidence set")
ok(isTRUE((is.na(fit_N$statistic) || fit_N$statistic <= fit_N$crit) ==
            in_set(fit_N$beta0, fit_N$conf_set)),
   "non-positive-variance fit: statistic/set consistency invariant")
ok(grid_agrees_fit(fit_N), "non-positive-variance fit: conf_set == grid")

# === Invariances =============================================================
fit0 <- cjs_q(y_A, x_A, jud_A, cluster = cl_A)   # beta0 = 0
# y -> c y: endpoints scale by c; LM at 0 and F_CJS invariant.
fitYc <- cjs_q(2.5 * y_A, x_A, jud_A, cluster = cl_A)
ok(set_near(fitYc$conf_set, fit0$conf_set * 2.5, 1e-6) &&
   relnear(fitYc$statistic, fit0$statistic, 1e-10) &&
   relnear(fitYc$F_CJS, fit0$F_CJS, 1e-10),
   "invariance y -> c y: endpoints x c; LM(0), F_CJS invariant")
# x -> c x: endpoints scale by 1/c; LM at 0 and F_CJS invariant.
fitXc <- cjs_q(y_A, 3 * x_A, jud_A, cluster = cl_A)
ok(set_near(fitXc$conf_set, fit0$conf_set / 3, 1e-6) &&
   relnear(fitXc$statistic, fit0$statistic, 1e-10) &&
   relnear(fitXc$F_CJS, fit0$F_CJS, 1e-10),
   "invariance x -> c x: endpoints / c; LM(0), F_CJS invariant")
# y -> y + c x: endpoints shift by +c; LM at beta0 + c matches LM at beta0.
fitSh <- cjs_q(y_A + 0.8 * x_A, x_A, jud_A, cluster = cl_A, beta0 = 0.8)
ok(set_near(fitSh$conf_set, fit0$conf_set + 0.8, 1e-6) &&
   relnear(fitSh$statistic, fit0$statistic, 1e-10),
   "invariance y -> y + c x: endpoints + c; LM(beta0 + c) invariant")
# Z -> Z A, A nonsingular: everything invariant to roundoff.
set.seed(9)
Amix <- matrix(rnorm(25), 5L, 5L) + 5 * diag(5L)
fitZA <- cjs_q(y_A, x_A, Zj_A %*% Amix, cluster = cl_A)
ok(relnear(fitZA$statistic, fit0$statistic, 1e-8) &&
   relnear(fitZA$F_CJS, fit0$F_CJS, 1e-8) &&
   set_near(fitZA$conf_set, fit0$conf_set, 1e-6),
   "invariance Z -> Z A (nonsingular): statistic, F_CJS, set unchanged")
# Cluster relabelling and full row permutation.
fitRel <- cjs_q(y_A, x_A, jud_A, cluster = paste0("c", 99 - cl_A))
prm <- sample(n_A)
fitPrm <- cjs_q(y_A[prm], x_A[prm], jud_A[prm], cluster = cl_A[prm])
ok(relnear(fitRel$statistic, fit0$statistic, 1e-12) &&
   set_near(fitRel$conf_set, fit0$conf_set, 1e-10) &&
   relnear(fitPrm$statistic, fit0$statistic, 1e-10) &&
   set_near(fitPrm$conf_set, fit0$conf_set, 1e-8),
   "invariance: cluster relabelling and row reordering")

# === Null calibration: KS sanity on R = 2000 draws ==========================
# Strong instruments, correctly specified null: LM(beta_true) ~ chisq_1.
set.seed(20250614)
R_ks <- 2000L
G_K <- 40L; ng_K <- 5L; n_K <- G_K * ng_K
cl_K <- rep(seq_len(G_K), each = ng_K)
lm_ks <- numeric(R_ks)
for (r in seq_len(R_ks)) {
  Z_K <- matrix(rnorm(n_K * 3L), n_K, 3L)
  u_K <- rnorm(G_K)[cl_K]
  v_K <- rnorm(n_K)
  x_K <- drop(Z_K %*% c(1, 0.6, -0.5)) + u_K + v_K
  y_K <- 0.5 * x_K + u_K + 0.7 * v_K + rnorm(n_K)
  po_K <- clusterIV:::.partial_out(y_K, x_K, Z_K, NULL, NULL, TRUE)
  cf_K <- clusterIV:::.cjs_coefs(clusterIV:::.cluster_sums(po_K$y, po_K$x, po_K$Z, cl_K))
  lm_ks[r] <- Spoly(cf_K$s, 0.5)^2 / Vpoly(cf_K$v, 0.5)
}
ks <- suppressWarnings(stats::ks.test(lm_ks, stats::pchisq, df = 1))
rej_ks <- mean(lm_ks > qchisq(0.95, 1))
cat(sprintf("INFO: KS D = %.4f (p = %.3f), null rejection at 5%% = %.3f\n",
            ks$statistic, ks$p.value, rej_ks))
ok(ks$statistic < 0.05, "null calibration: KS distance to chisq_1 < 0.05")
ok(rej_ks > 0.03 && rej_ks < 0.07,
   "null calibration: 5% rejection rate within [0.03, 0.07]")

# === Formula interface, methods, errors =====================================
datA <- data.frame(y = y_A, x = x_A, judge = jud_A, cl = cl_A)
f_A <- cjs_q(y ~ x | judge, data = datA, cluster = ~cl)
ok(identical(f_A$statistic, fit0$statistic) &&
   set_near(f_A$conf_set, fit0$conf_set, 1e-12),
   "formula interface == default interface")
datB <- data.frame(y = y_B, x = x_B, z1 = Z_B[, 1], z2 = Z_B[, 2],
                   z3 = Z_B[, 3], f5 = f5, cl = cl_B)
f_B <- cjs_q(y ~ x | z1 + z2 + z3 | f5, data = datB, cluster = ~cl,
             beta0 = 0.4)
ok(identical(f_B$statistic, fitF$statistic) &&
   set_near(f_B$conf_set, fitF$conf_set, 1e-12),
   "three-part formula y ~ x | z | fe == fixed_effects argument")
ok(nobs(fitA) == n_A, "nobs.cjscore returns n")
ok(identical(confint(fitA), fitA$conf_set),
   "confint at the fitted level returns the stored set")
ci90 <- confint(fitA, level = 0.90)
dir90 <- cjs_q(y_A, x_A, jud_A, cluster = cl_A, level = 0.90)$conf_set
ok(set_near(ci90, dir90, 1e-6),
   "confint at a new level == direct refit at that level")
ci30 <- confint(fitA, level = 0.30)
dir30 <- cjs_q(y_A, x_A, jud_A, cluster = cl_A, level = 0.30)$conf_set
ok(set_near(ci30, dir30, 1e-6),
   "confint accepts level = 0.3 and agrees with a direct refit")
out_pr <- capture.output({ print(fitA); print(fitWk); print(fit_N) })
ok(any(grepl("Cluster jackknife score", out_pr)) &&
   any(grepl("unbounded", out_pr)) &&
   any(grepl("spurious regions", out_pr)) &&
   any(grepl("undefined and beta0 is accepted (conservative convention",
             out_pr, fixed = TRUE)),
   "print method: header, endpoint-derived topology, spurious-region and conservative notes render")
ok(any(grepl("F_CJS^2 = ", out_pr, fixed = TRUE)) &&
   any(grepl("first-stage LM statistic", out_pr, fixed = TRUE)),
   "print method: F_CJS^2 rendered as the first-stage LM statistic")
# D4 honest label: at G == n the print header names the Matsushita-Otsu
# reduction; the "cluster jackknife" wording must not appear.
fit_S <- cjs_q(y_S, x_S, Z_S, cluster = seq_len(n_S))
out_S <- capture.output(print(fit_S))
ok(any(grepl("Jackknife score test (independent data; Matsushita-Otsu 2024)",
             out_S, fixed = TRUE)) &&
   !any(grepl("Cluster jackknife", out_S)),
   "G == n print label: Matsushita-Otsu independent-data header")
# The boundedness criterion is sign-blind (F_CJS^2 vs crit), so a bounded set
# with a negative F_CJS must not print the unbounded flag.
fit_neg <- fitA
fit_neg$F_CJS <- -abs(fit_neg$F_CJS)
out_neg <- capture.output(print(fit_neg))
ok(any(grepl("first-stage LM statistic", out_neg, fixed = TRUE)) &&
   !any(grepl("confidence set unbounded", out_neg, fixed = TRUE)),
   "bounded set with F_CJS < 0: no contradictory unbounded flag")
fit04 <- cjs_q(y_A, x_A, jud_A, cluster = cl_A, level = 0.4)
ok(identical(fit04$level, 0.4) && is.finite(fit04$crit),
   "accepts level = 0.4 (the contract is the open interval (0, 1))")
ok(errs(cjscore(y_A, x_A, jud_A, cluster = cl_A, beta0 = c(0, 1))),
   "stops on non-scalar beta0")
ok(errs(cjscore(y_A, x_A, jud_A, cluster = rep(1L, n_A))),
   "stops on a single cluster")
msgS <- tryCatch(cjs_q(y_B, x_B, cbind(Z_B, Z_B[, 1]), cluster = cl_B),
                 error = function(e) conditionMessage(e))
ok(is.character(msgS) && grepl("singular", msgS),
   "collinear instruments hit the singular-Gram message")
wm <- tryCatch(cjscore(y_A, x_A, jud_A, cluster = cl_A),
               warning = function(w) conditionMessage(w))
ok(is.character(wm) && grepl("asymptotic in the number of clusters", wm),
   "small-G advisory fires for G < 20")

cat("\nAll cjscore tests passed.\n")
