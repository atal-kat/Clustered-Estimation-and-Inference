# Correctness tests for the Montiel Olea-Pflueger effective F (base R, no
# testthat), in the style of the other test files. Oracles O1-O4 of the WP3
# spec; the math reference is dev/extraction/mop-effective-f-extraction.md.

library(clusterIV)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}
near <- function(a, b, tol) max(abs(a - b)) < tol
quiet <- function(expr) suppressWarnings(expr)

# === O1: the paper's Table 1 (published frozen oracle) =======================
# 5% Patnaik critical values for the simplified TSLS procedure, K_eff = 1..30.
tab_x10 <- c(23.11, 19.29, 17.67, 16.72, 16.08, 15.62, 15.26, 14.97, 14.73, 14.53,
             14.36, 14.21, 14.08, 13.96, 13.86, 13.77, 13.68, 13.60, 13.53, 13.46,
             13.40, 13.35, 13.29, 13.24, 13.20, 13.15, 13.11, 13.07, 13.04, 13.00)
tab_x5  <- c(15.06, 12.17, 10.95, 10.23,  9.75,  9.40,  9.14,  8.92,  8.74,  8.59,
              8.47,  8.36,  8.26,  8.17,  8.10,  8.03,  7.96,  7.91,  7.85,  7.80,
              7.76,  7.72,  7.68,  7.64,  7.61,  7.57,  7.54,  7.51,  7.49,  7.46)
got_x10 <- vapply(1:30, clusterIV:::.eff_f_crit, numeric(1), x = 10)
got_x5  <- vapply(1:30, clusterIV:::.eff_f_crit, numeric(1), x = 5)
ok(max(abs(got_x10 - tab_x10)) < 0.01,
   "O1: Table 1, x = 10 (tau = 10%), K_eff = 1..30 reproduced to +-0.01")
ok(max(abs(got_x5 - tab_x5)) < 0.01,
   "O1: Table 1, x = 5 (tau = 20%), K_eff = 1..30 reproduced to +-0.01")
ok(abs(got_x10[1] - 23.11) < 0.01,
   "O1: K_eff = 1, x = 10 gives 23.11 (the 'effective F > 23.1' rule of thumb)")

# === a design generator ======================================================
make_design <- function(sizes, k, seed) {
  set.seed(seed)
  n <- sum(sizes)
  cl <- rep(seq_along(sizes), sizes)
  Z <- matrix(rnorm(n * k), n, k)
  u <- rnorm(length(sizes))[cl]
  x <- drop(Z %*% rep(c(1, -0.6), length.out = k)) + u + rnorm(n)
  y <- 0.7 * x + u + rnorm(n)
  list(y = y, x = x, Z = Z, cl = cl, n = n)
}

eff_of <- function(fit) c(fit$F_eff, fit$K_eff, fit$F_eff_crit)

# === O2: invariance to reparametrisation of Z ================================
# For any non-singular M, Z -> Z M leaves the column span (hence P_Z and the
# first-stage residuals) unchanged and rotates Ztil orthogonally, which
# preserves tr(W2), tr(W2'W2) and the max eigenvalue. F_eff, K_eff and
# F_eff_crit must be unchanged to 1e-10. An implementation that forgets to
# orthonormalise the instruments fails this and passes almost everything else.
d2 <- make_design(sample(3:9, 20, replace = TRUE), 4L, 101)
base2 <- eff_of(quiet(cjar(d2$y, d2$x, d2$Z, cluster = d2$cl)))
set.seed(7)
Ms <- list(
  `random M`            = matrix(rnorm(16), 4, 4),
  `column permutation`  = diag(4)[, c(3, 1, 4, 2)],
  `diagonal rescaling`  = diag(c(0.01, 5, 100, 0.3))
)
for (nm in names(Ms)) {
  alt <- eff_of(quiet(cjar(d2$y, d2$x, d2$Z %*% Ms[[nm]], cluster = d2$cl)))
  ok(near(base2, alt, 1e-10),
     sprintf("O2: F_eff/K_eff/crit invariant under Z %%*%% M (%s), 1e-10", nm))
}

# === O3: just-identified case ================================================
# With K = 1 the effective F coincides with the robust first-stage F: the
# square of the cluster-robust t statistic on the single instrument, with no
# small-sample correction (the paper's What_2 carries none). Hand-rolled
# cluster-robust vcov on the intercept-partialled first stage.
d3 <- make_design(sample(4:10, 25, replace = TRUE), 1L, 202)
fit3 <- quiet(cjar(d3$y, d3$x, d3$Z, cluster = d3$cl))
zc <- d3$Z[, 1L] - mean(d3$Z[, 1L])          # intercept partialled out
xc <- d3$x - mean(d3$x)
b1 <- sum(zc * xc) / sum(zc^2)               # first-stage OLS slope
v <- xc - b1 * zc                            # first-stage residuals
meat <- sum(tapply(zc * v, d3$cl, sum)^2)    # clustered meat, no correction
t2 <- b1^2 / (meat / sum(zc^2)^2)            # squared cluster-robust t
ok(near(fit3$F_eff, t2, 1e-10),
   "O3: K = 1: F_eff == squared cluster-robust first-stage t (1e-10)")
ok(near(fit3$K_eff, 1, 1e-10), "O3: K = 1 gives K_eff = 1 exactly")
ok(near(fit3$F_eff_crit, 23.11, 0.01), "O3: K = 1 critical value is 23.11")

# === O4: dense brute force (literal eq. 16) ==================================
# Explicit orthonormalisation Zstar = sqrt(n) Z chol(Z'Z)^-1 (so that
# Zstar'Zstar / n == I_k, asserted), What_2 built by a per-cluster for loop
# of outer products, F_eff = x'P_Z x / tr(What_2). Compared to the fast
# whitened path at 1e-10 on balanced, unbalanced, singleton, weighted and
# fixed-effect-absorbed designs.
brute_eff_f <- function(x, Z, cl) {
  n <- length(x)
  Zstar <- sqrt(n) * Z %*% solve(chol(crossprod(Z)))
  stopifnot(max(abs(crossprod(Zstar) / n - diag(ncol(Z)))) < 1e-10)
  pihat <- solve(crossprod(Zstar), crossprod(Zstar, x))
  xPx <- drop(t(x) %*% Zstar %*% pihat)
  v <- x - drop(Zstar %*% pihat)
  W2 <- matrix(0, ncol(Z), ncol(Z))
  for (g in unique(cl)) {
    ig <- which(cl == g)
    sg <- crossprod(Zstar[ig, , drop = FALSE], v[ig])   # k x 1
    W2 <- W2 + sg %*% t(sg)
  }
  W2 <- W2 / n
  x10 <- 10
  lmax <- max(eigen(W2, symmetric = TRUE, only.values = TRUE)$values)
  trW2 <- sum(diag(W2))
  K_eff <- trW2^2 * (1 + 2 * x10) /
    (sum(diag(crossprod(W2))) + 2 * x10 * trW2 * lmax)
  c(F_eff = xPx / trW2, K_eff = K_eff,
    crit = qchisq(0.95, df = K_eff, ncp = x10 * K_eff) / K_eff)
}

designs <- list(
  `balanced (n_g = 6, k = 4)`   = make_design(rep(6L, 18), 4L, 301),
  `unbalanced (n_g in 1..12)`   = make_design(c(1L, 2L, 4L, 7L, 12L, 1L, 9L, 3L, 6L, 10L), 3L, 302),
  `singletons (n = G = 50)`     = make_design(rep(1L, 50), 3L, 303)
)
for (nm in names(designs)) {
  d <- designs[[nm]]
  fit <- quiet(cjar(d$y, d$x, d$Z, cluster = d$cl))
  po <- clusterIV:::.partial_out(d$y, d$x, d$Z, NULL, NULL, TRUE)
  bf <- brute_eff_f(po$x, po$Z, d$cl)
  ok(near(eff_of(fit), unname(bf), 1e-10),
     sprintf("O4: fast whitened path == literal eq. 16 brute force [%s] (1e-10)", nm))
}

# The implementation forms the smaller of Sg'Sg and Sg Sg'.  Their nonzero
# eigenvalues and squared Frobenius traces agree, so the G x G branch used when
# k > G must still equal the literal k x k oracle.
d_wide <- make_design(rep(10L, 8L), 20L, 3031)
po_wide <- clusterIV:::.partial_out(d_wide$y, d_wide$x, d_wide$Z,
                                    NULL, NULL, TRUE)
fs_wide <- clusterIV:::.first_stage(po_wide$x, po_wide$Z)
fast_wide <- clusterIV:::.eff_f(fs_wide$Ztil, fs_wide$t, fs_wide$e,
                                d_wide$cl)
bf_wide <- brute_eff_f(po_wide$x, po_wide$Z, d_wide$cl)
ok(near(unlist(fast_wide), unname(bf_wide), 1e-10),
   "O4: k > G smaller-Gram branch == literal k x k brute force (1e-10)")

# F_eff and K_eff are homogeneous of degree zero in the first-stage scale.
# These scales force the old quadratic/quartic raw sums to underflow/overflow.
d_scale <- make_design(rep(6L, 18L), 4L, 3032)
po_scale <- clusterIV:::.partial_out(d_scale$y, d_scale$x, d_scale$Z,
                                     NULL, NULL, TRUE)
fs_scale <- clusterIV:::.first_stage(po_scale$x, po_scale$Z)
eff_scale <- clusterIV:::.eff_f(fs_scale$Ztil, fs_scale$t, fs_scale$e,
                                d_scale$cl)
for (s in c(2^-600, 2^500)) {
  alt <- clusterIV:::.eff_f(fs_scale$Ztil, fs_scale$t * s,
                            fs_scale$e * s, d_scale$cl)
  ok(near(unlist(alt), unlist(eff_scale), 1e-10),
     sprintf("O4: effective F is invariant at first-stage scale %.3g", s))
}

# Weighted design: the statistic is the effective F of the sqrt(w)-transformed
# model, so the brute force runs on the transformed, intercept-partialled data.
dw <- make_design(sample(3:8, 22, replace = TRUE), 3L, 304)
set.seed(305)
w <- runif(dw$n, 0.5, 2)
fit_w <- quiet(cjar(dw$y, dw$x, dw$Z, cluster = dw$cl, weights = w))
po_w <- clusterIV:::.partial_out(dw$y, dw$x, dw$Z, NULL, w, TRUE)
bf_w <- brute_eff_f(po_w$x, po_w$Z, dw$cl)
ok(near(eff_of(fit_w), unname(bf_w), 1e-10),
   "O4: weighted design == brute force on the rw-transformed data (1e-10)")

# Fixed-effect-absorbed design: brute force on the FE-partialled data.
df <- make_design(sample(4:9, 20, replace = TRUE), 3L, 306)
set.seed(307)
fe <- factor(sample(letters[1:7], df$n, replace = TRUE))
fit_f <- quiet(cjar(df$y, df$x, df$Z, cluster = df$cl, fixed_effects = fe))
po_f <- clusterIV:::.partial_out(df$y, df$x, df$Z, NULL, NULL, TRUE,
                                 fixed_effects = list(fe))
bf_f <- brute_eff_f(po_f$x, po_f$Z, df$cl)
ok(near(eff_of(fit_f), unname(bf_f), 1e-10),
   "O4: fixed-effect-absorbed design == brute force on the FE-partialled data (1e-10)")

# === consistency across the four object classes ==============================
d5 <- make_design(sample(4:9, 24, replace = TRUE), 4L, 401)
f_est <- quiet(cjive(d5$y, d5$x, d5$Z, cluster = d5$cl))
f_ar <- quiet(cjar(d5$y, d5$x, d5$Z, cluster = d5$cl))
f_sc <- quiet(cjscore(d5$y, d5$x, d5$Z, cluster = d5$cl))
f_all <- quiet(iv_infer(d5$y, d5$x, d5$Z, cluster = d5$cl))
ok(identical(eff_of(f_est), eff_of(f_ar)) &&
   identical(eff_of(f_ar), eff_of(f_sc)) &&
   identical(eff_of(f_sc), eff_of(f_all)) &&
   identical(eff_of(f_all), eff_of(f_all$cjive)),
   "F_eff/K_eff/crit identical across cjive, cjar, cjscore, iv_infer and components")

# The leaveout_mean path has no whitened first stage: F_eff is NA there and
# print() stays silent about it.
set.seed(402)
jm <- factor(sample(1:5, 200, replace = TRUE))
clm <- rep(1:25, each = 8)
xm <- as.numeric(jm) + rnorm(200)
ym <- xm + rnorm(200)
fm <- cjive(ym, xm, jm, cluster = clm, method = "leaveout_mean")
ok(is.na(fm$F_eff) && is.na(fm$K_eff) && is.na(fm$F_eff_crit),
   "leaveout_mean path: F_eff/K_eff/crit are NA")
ok(!any(grepl("effective F", capture.output(print(fm)))),
   "leaveout_mean path: print() omits the effective-F line")

# print() shows F_eff beside F_CJ (both lines present, F_CJ first).
out <- capture.output(print(f_ar))
i_fcj <- grep("^  F_CJ = ", out)
i_eff <- grep("effective F \\(Montiel Olea-Pflueger\\)", out)
ok(length(i_fcj) == 1L && length(i_eff) == 1L && i_eff == i_fcj + 1L,
   "print.cjar: effective F appears directly beside F_CJ, never instead of it")

cat("\nAll effective-F tests passed.\n")
