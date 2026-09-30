# Correctness tests for the clusterIV package (base R, no testthat).
# Each block stops with an informative message on failure, so a failure fails
# `R CMD check`.

library(clusterIV)

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}
near <- function(a, b, tol) max(abs(a - b)) < tol

# --- data with a clustered error structure ---------------------------------
set.seed(123)
G <- 60L
sizes <- sample(3:8, G, replace = TRUE)
n <- sum(sizes)
cl <- rep(seq_len(G), sizes)
ucl <- rnorm(G)[cl]                              # common within-cluster shock
J <- 10L
judge <- factor(sample(seq_len(J), n, replace = TRUE))
x <- as.numeric(judge) + ucl + rnorm(n)
y <- 1.0 * x + ucl + rnorm(n)
Zc <- matrix(rnorm(n * 4L), n, 4L)               # continuous instruments
xc <- Zc %*% c(1, -1, 0.5, 0.25) + ucl + rnorm(n)
yc <- 1.0 * xc + ucl + rnorm(n)

# Residualised (intercept-only) quantities for the judge design, via the same
# front end the package uses.
Zj <- stats::model.matrix(~judge)[, -1L, drop = FALSE]
po <- clusterIV:::.partial_out(y, x, Zj, controls = NULL, weights = NULL, intercept = TRUE)
yr <- po$y; xr <- po$x; Zr <- po$Z
groups <- split(seq_len(n), factor(cl))

# Brute-force FLM leave-cluster-out fitted values on the residualised data.
brute_cjive <- function(D, Z, groups) {
  p <- numeric(length(D))
  for (idx in groups) {
    keep <- setdiff(seq_len(length(D)), idx)
    pi_g <- solve(crossprod(Z[keep, , drop = FALSE]),
                  crossprod(Z[keep, , drop = FALSE], D[keep]))
    p[idx] <- Z[idx, , drop = FALSE] %*% pi_g
  }
  p
}

# === Test 1: dense CJIVE == brute-force FLM definition ======================
fs <- clusterIV:::.first_stage(xr, Zr)
phat_pkg <- clusterIV:::.leaveout_fit(xr, fs$Ztil, fs$t, groups, ncol(Zr))$phat
phat_bf <- brute_cjive(xr, Zr, groups)
ok(near(phat_pkg, phat_bf, 1e-10), "dense CJIVE p-hat vector == brute force (1e-10)")

# iv_compare() does not report maxlev.  Its comparison-only kernel skips that
# eigendecomposition but must return the same fitted values as the diagnostic
# path and the literal dense leave-out definition.
lo_fast <- clusterIV:::.leaveout_fit(xr, fs$Ztil, fs$t, groups, ncol(Zr),
                                     leverage = FALSE)
ok(near(lo_fast$phat, phat_bf, 1e-10) && is.na(lo_fast$maxlev),
   "comparison-only leave-out path skips leverage and matches brute force (1e-10)")

beta_pkg <- cjive(y, x, judge, cluster = cl)$coefficient
beta_bf <- sum(phat_bf * yr) / sum(phat_bf * xr)
ok(near(beta_pkg, beta_bf, 1e-10), "dense CJIVE coefficient == brute force (1e-10)")

# === Test 2: singleton clusters: CJIVE == improved JIVE ====================
# Both sides run on the residualised (FWL-partialled) instruments, so the
# singleton limit is the improved JIVE (IJIVE; Ackerberg-Devereux 2009), the
# leave-one-out fit -- not the original AIK (1999) JIVE.
sing <- split(seq_len(n), factor(seq_len(n)))
phat_sing <- clusterIV:::.leaveout_fit(xr, fs$Ztil, fs$t, sing, ncol(Zr))$phat
phat_jive <- clusterIV:::.phat_jive(fs, xr)
ok(near(phat_sing, phat_jive, 1e-10),
   "singleton-cluster CJIVE == improved JIVE (1e-10)")

# === Test 3: cjive() default == iv_compare() CJIVE row =====================
fit3 <- cjive(y, x, judge, cluster = cl)
tab3 <- iv_compare(y, x, judge, cluster = cl)
row3 <- tab3[tab3$estimator == "CJIVE", ]
ok(near(fit3$coefficient, row3$coefficient, 1e-10) &&
   near(fit3$se, row3$se, 1e-10),
   "cjive() == iv_compare() CJIVE row (1e-10)")

# === Test 4: leaveout_mean == brute-force leave-cluster-out group mean ======
brute_mean <- function(x, group, cluster) {
  vapply(seq_along(x), function(i) {
    sel <- group == group[i] & cluster != cluster[i]
    mean(x[sel])
  }, numeric(1))
}
phat_mean_pkg <- clusterIV:::.leaveout_mean(x, judge, cl, weights = NULL)
phat_mean_bf <- brute_mean(x, judge, cl)
ok(near(phat_mean_pkg, phat_mean_bf, 1e-8),
   "leaveout_mean p-hat == brute-force group mean (1e-8)")

fit_dense <- cjive(y, x, judge, cluster = cl, method = "dense")
fit_mean <- cjive(y, x, judge, cluster = cl, method = "leaveout_mean")
gap <- abs(fit_dense$coefficient - fit_mean$coefficient)
cat(sprintf("INFO: dense-vs-mean coefficient gap = %.3g (intercept term)\n", gap))
ok(gap < 5e-3, "dense-vs-mean gap is the small intercept term (< 5e-3)")

# === Test 5: formula interface == matrix interface =========================
datj <- data.frame(y = y, x = x, judge = judge, cl = cl)
f_judge <- cjive(y ~ x | judge, data = datj, cluster = ~cl)
m_judge <- cjive(y, x, judge, cluster = cl)
ok(near(f_judge$coefficient, m_judge$coefficient, 1e-10) &&
   near(f_judge$se, m_judge$se, 1e-10),
   "formula == matrix interface, judge design (1e-10)")

datc <- data.frame(y = yc, x = xc, z1 = Zc[, 1], z2 = Zc[, 2],
                   z3 = Zc[, 3], z4 = Zc[, 4], cl = cl)
f_cont <- cjive(y ~ x | z1 + z2 + z3 + z4, data = datc, cluster = ~cl)
m_cont <- cjive(yc, xc, Zc, cluster = cl)
ok(near(f_cont$coefficient, m_cont$coefficient, 1e-10) &&
   near(f_cont$se, m_cont$se, 1e-10),
   "formula == matrix interface, continuous z (1e-10)")

# === Test 6: iv_compare internal checks ====================================
tab6 <- iv_compare(yc, xc, Zc, cluster = cl)
ok(identical(tab6$estimator, c("OLS", "2SLS", "JIVE", "CJIVE")),
   "iv_compare returns the four estimators in order")
# Hand OLS on residualised (FWL) data: slope of yr on xr.
poc <- clusterIV:::.partial_out(yc, xc, Zc, NULL, NULL, TRUE)
ols_hand <- sum(poc$x * poc$y) / sum(poc$x * poc$x)
ok(near(tab6$coefficient[1], ols_hand, 1e-10), "iv_compare OLS row == hand OLS (1e-10)")

# === Test 7: SE finite/positive; weighted path matches independent recompute =
ok(is.finite(fit3$se) && fit3$se > 0, "SE is finite and positive")

set.seed(7)
Gw <- 25L; sw <- sample(2:5, Gw, replace = TRUE); nw <- sum(sw)
clw <- rep(seq_len(Gw), sw)
uw <- rnorm(Gw)[clw]
Zw <- matrix(rnorm(nw * 3L), nw, 3L)
xw <- Zw %*% c(1, 0.5, -0.5) + uw + rnorm(nw)
yw <- xw + uw + rnorm(nw)
w <- runif(nw, 0.5, 2)

fit_w <- cjive(yw, xw, Zw, cluster = clw, weights = w)

# Independent weighted recomputation in the sqrt(w) geometry.
rw <- sqrt(w)
M1 <- cbind(rw)                                  # weighted intercept
resid <- function(v) v - M1 %*% solve(crossprod(M1), crossprod(M1, v))
Yt <- resid(yw * rw); Xt <- resid(xw * rw); Zt <- resid(Zw * rw)
gw <- split(seq_len(nw), factor(clw))
phat_w <- brute_cjive(as.numeric(Xt), Zt, gw)
beta_w <- sum(phat_w * Yt) / sum(phat_w * Xt)
den_w <- sum(phat_w * Xt)
s_w <- phat_w * (Yt - beta_w * Xt)
Sg_w <- tapply(s_w, factor(clw), sum)
se_w <- sqrt(sum(Sg_w^2) * Gw / (Gw - 1)) / abs(den_w)
ok(near(fit_w$coefficient, beta_w, 1e-10) && near(fit_w$se, se_w, 1e-10),
   "weighted CJIVE == independent weighted recomputation (1e-10)")

# A nominally invertible projection complement below the numerical rank
# frontier is not silently solved.  Test both the diagnostic and the fast
# comparison-only branches on a directly controlled 2 x 2 block.
Z_frontier <- diag(c(sqrt(1 - .Machine$double.eps), 0.5), 2L)
frontier_error <- function(leverage) inherits(tryCatch(
  clusterIV:::.cluster_block(Z_frontier, 2L, c(1, 1), leverage = leverage),
  error = function(e) e), "error")
ok(frontier_error(TRUE) && frontier_error(FALSE),
   "leave-out kernel rejects a numerically rank-deficient downdate")

# === Test 8: input validation (clear, specific errors) =====================
errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
ok(errs(cjive(y, x, judge, cluster = rep(1L, n))), "stops on < 2 clusters")
ok(errs(cjive(y, x, judge, cluster = cl, weights = replace(rep(1, n), 1, -1))),
   "stops on non-positive weights")
ok(errs(cjive(y[-1], x, judge, cluster = cl)), "stops on length mismatch")
ok(errs(cjive(replace(y, 1, NA), x, judge, cluster = cl)), "stops on NA in y")
ok(errs(cjive(y, x, judge, cluster = replace(cl, 1, NA))), "stops on NA in cluster")
nested_dense <- suppressWarnings(cjive(y, x, factor(cl), cluster = cl))
ok(is.finite(nested_dense$coefficient) && nested_dense$maxlev < 1,
   "dense path uses the full-rank encoded design when groups equal clusters")
ok(errs(cjive(y, x, factor(cl), cluster = cl, method = "leaveout_mean")),
   "leaveout_mean stops when each instrument group lies in one cluster")

# === Test 9: agreement with FLM/McIntyre Stata `cjive` =====================
# Literal translation of the Mata core (leverage trick + ivregress 2sls,
# noconstant, vce(cluster)). The point estimate is exact; the SE matches because
# covariates are partialled out, so the final model has K = 1 parameter and
# Stata's q_c = (N-1)/(N-K) * G/(G-1) reduces to G/(G-1).
stata_cjive_ref <- function(y, d, Z, Wexog, cluster) {
  N <- length(y)
  Cx <- if (is.null(Wexog)) matrix(1, N, 1L) else cbind(1, Wexog)
  rr <- function(v) v - Cx %*% solve(crossprod(Cx), crossprod(Cx, v))
  yres <- as.numeric(rr(y)); dres <- as.numeric(rr(d))
  Zres <- apply(Z, 2, function(z) as.numeric(rr(z)))
  ZTZI <- solve(crossprod(Zres))
  Pi <- ZTZI %*% crossprod(Zres, dres)
  Lev <- numeric(N)
  for (g in unique(cluster)) {
    ix <- which(cluster == g); ng <- length(ix)
    AZC <- Zres[ix, , drop = FALSE]
    AHC <- AZC %*% ZTZI %*% t(AZC)
    A <- AZC %*% Pi - AHC %*% dres[ix]
    Lev[ix] <- solve(diag(ng) - AHC) %*% A
  }
  beta <- sum(Lev * yres) / sum(Lev * dres)
  S <- tapply(Lev * (yres - beta * dres), cluster, sum)
  M <- length(unique(cluster))
  se <- sqrt((M / (M - 1)) * sum(S^2)) / abs(sum(Lev * dres))
  list(beta = beta, se = se)
}

set.seed(99)
Wq <- cbind(rnorm(n), rnorm(n))                  # exogenous covariates
dq <- as.numeric(judge) + Wq %*% c(.4, -.3) + ucl + rnorm(n)
yq <- 1.3 * dq + Wq %*% c(1, 1) + ucl + rnorm(n)
ref <- stata_cjive_ref(yq, dq, model.matrix(~judge)[, -1, drop = FALSE], Wq, cl)
fitq <- cjive(yq, dq, judge, cluster = cl, controls = Wq)
ok(near(fitq$coefficient, ref$beta, 1e-10) && near(fitq$se, ref$se, 1e-10),
   "matches FLM Stata cjive (coefficient and SE, 1e-10)")

# === Test 10: single FE == dense dummy controls ============================
set.seed(11)
f1 <- factor(sample(letters[1:12], n, replace = TRUE))
fit_fe1 <- cjive(yc, xc, Zc, cluster = cl, fixed_effects = f1)
fit_dd1 <- cjive(yc, xc, Zc, cluster = cl,
                 controls = model.matrix(~f1)[, -1, drop = FALSE])
ok(near(fit_fe1$coefficient, fit_dd1$coefficient, 1e-8) &&
   near(fit_fe1$se, fit_dd1$se, 1e-8),
   "single FE == dense dummy route (coef and se, 1e-8)")

# === Test 11: two crossed FE dimensions == dense two-FE dummy route ========
set.seed(12)
f2 <- factor(sample(LETTERS[1:9], n, replace = TRUE))   # crossed with f1
fit_fe2 <- cjive(yc, xc, Zc, cluster = cl, fixed_effects = list(f1, f2))
fit_dd2 <- cjive(yc, xc, Zc, cluster = cl,
                 controls = cbind(model.matrix(~f1)[, -1, drop = FALSE],
                                  model.matrix(~f2)[, -1, drop = FALSE]))
ok(near(fit_fe2$coefficient, fit_dd2$coefficient, 1e-8) &&
   near(fit_fe2$se, fit_dd2$se, 1e-8),
   "two crossed FE == dense two-FE dummy route (coef and se, 1e-8)")

# === Test 12: weighted FE case == weighted dense dummy route ===============
set.seed(13)
w12 <- runif(n, 0.5, 2)
fit_few <- cjive(yc, xc, Zc, cluster = cl, weights = w12,
                 fixed_effects = list(f1, f2))
fit_ddw <- cjive(yc, xc, Zc, cluster = cl, weights = w12,
                 controls = cbind(model.matrix(~f1)[, -1, drop = FALSE],
                                  model.matrix(~f2)[, -1, drop = FALSE]))
ok(near(fit_few$coefficient, fit_ddw$coefficient, 1e-8) &&
   near(fit_few$se, fit_ddw$se, 1e-8),
   "weighted two-FE == weighted dense dummy route (coef and se, 1e-8)")

# === Test 13: three-part formula == matrix interface ========================
datf <- data.frame(y = yc, x = xc, z1 = Zc[, 1], z2 = Zc[, 2],
                   z3 = Zc[, 3], z4 = Zc[, 4], f1 = f1, f2 = f2, cl = cl)
f_fe <- cjive(y ~ x | z1 + z2 + z3 + z4 | f1 + f2, data = datf, cluster = ~cl)
m_fe <- cjive(yc, xc, Zc, cluster = cl, fixed_effects = list(f1, f2))
ok(near(f_fe$coefficient, m_fe$coefficient, 1e-10) &&
   near(f_fe$se, m_fe$se, 1e-10),
   "three-part formula == matrix interface (1e-10)")
a_fe <- cjive(y ~ x | z1 + z2 + z3 + z4, data = datf, cluster = ~cl,
              fixed_effects = ~ f1 + f2)
ok(near(a_fe$coefficient, m_fe$coefficient, 1e-10) &&
   near(a_fe$se, m_fe$se, 1e-10),
   "fixed_effects = ~ f1 + f2 argument == matrix interface (1e-10)")

# === Test 14: FE + additional dense controls == all-dense route =============
set.seed(14)
W14 <- cbind(rnorm(n), rnorm(n))
fit_mix <- cjive(yc, xc, Zc, cluster = cl, controls = W14,
                 fixed_effects = list(f1, f2))
fit_all <- cjive(yc, xc, Zc, cluster = cl,
                 controls = cbind(W14,
                                  model.matrix(~f1)[, -1, drop = FALSE],
                                  model.matrix(~f2)[, -1, drop = FALSE]))
ok(near(fit_mix$coefficient, fit_all$coefficient, 1e-8) &&
   near(fit_mix$se, fit_all$se, 1e-8),
   "FE + dense controls == all-dense route (coef and se, 1e-8)")

# === Test 15: FE input validation ===========================================
ok(errs(cjive(yc, xc, Zc, cluster = cl,
              fixed_effects = replace(as.character(f1), 1, NA))),
   "stops on NA in fixed effects")
ok(errs(cjive(yc, xc, Zc, cluster = cl, fixed_effects = f1[-1])),
   "stops on wrong-length fixed effects")
# Non-convergence message path: a deterministic crossed design and one solver
# iteration cannot meet the deliberately tight certificate.
i15 <- seq_len(n)
f15a <- factor((i15 - 1L) %% 13L)
f15b <- factor(((i15 - 1L) %/% 7L + 3L * i15) %% 11L)
m15 <- cbind(sin(i15 * 0.137) + cos(i15 * 0.043))
msg <- tryCatch(clusterIV:::.demean(m15, list(f15a, f15b),
                                    tol = 1e-14, maxit = 1L),
                error = function(e) conditionMessage(e))
ok(is.character(msg) && grepl("converge|certif", msg, ignore.case = TRUE) &&
     grepl("sweep|iteration", msg, ignore.case = TRUE) &&
     grepl("change|orthogon|residual", msg, ignore.case = TRUE) &&
     grepl("tol|tolerance", msg, ignore.case = TRUE),
   "non-convergence stop() message is reachable and informative")

# === Test 16: inference dispatcher ==========================================
fit_asym <- cjive(y, x, judge, cluster = cl, inference = "asymptotic")
ok(identical(fit_asym$coefficient, fit3$coefficient) &&
   identical(fit_asym$se, fit3$se) &&
   identical(fit_asym$conf.low, fit3$conf.low) &&
   identical(fit_asym$p.value, fit3$p.value),
   "inference = \"asymptotic\" is bit-identical to the default")
zc16 <- qnorm(1 - (1 - fit_asym$level) / 2)   # mirrors the v0.1.0 construction
ok(identical(fit_asym$conf.low, fit_asym$coefficient - zc16 * fit_asym$se) &&
   identical(fit_asym$p.value, 2 * pnorm(-abs(fit_asym$statistic))),
   "asymptotic branch reproduces the v0.1.0 normal-quantile construction")
fit_t <- cjive(y, x, judge, cluster = cl, inference = "t")
ok(identical(fit_t$coefficient, fit_asym$coefficient) &&
   identical(fit_t$se, fit_asym$se),
   "inference = \"t\" leaves coefficient and se unchanged")
width_t <- fit_t$conf.high - fit_t$conf.low
width_a <- fit_asym$conf.high - fit_asym$conf.low
ok(width_t > width_a, "t(G-1) interval is wider")
ok(near(width_t / width_a, qt(0.975, G - 1) / qnorm(0.975), 1e-12),
   "CI width ratio == qt(.975, G-1)/qnorm(.975) (1e-12)")

# === Test 17: print labels follow the inference dispatcher ==================
out_a <- capture.output({ print(fit_asym); print(summary(fit_asym)) })
out_t <- capture.output({ print(fit_t); print(summary(fit_t)) })
ok(any(grepl("  z = ", out_a, fixed = TRUE)) &&
   !any(grepl("  t = ", out_a, fixed = TRUE)) &&
   any(grepl("z value", out_a, fixed = TRUE)) &&
   any(grepl("Pr(>|z|)", out_a, fixed = TRUE)) &&
   !any(grepl("t value", out_a, fixed = TRUE)) &&
   !any(grepl("Pr(>|t|)", out_a, fixed = TRUE)),
   "asymptotic fit prints z labels only (header and coefficient table)")
ok(any(grepl("  t = ", out_t, fixed = TRUE)) &&
   !any(grepl("  z = ", out_t, fixed = TRUE)) &&
   any(grepl("t value", out_t, fixed = TRUE)) &&
   any(grepl("Pr(>|t|)", out_t, fixed = TRUE)) &&
   !any(grepl("z value", out_t, fixed = TRUE)) &&
   !any(grepl("Pr(>|z|)", out_t, fixed = TRUE)),
   "t fit prints t labels only (header and coefficient table)")

# === Test 18: JIVE guard against unit leverage ==============================
# The third instrument column is an indicator for observation 1, so h_1 = 1
# exactly and JIVE's division by 1 - h is undefined; the whole comparison
# table must fail with the leverage message (the CJIVE row is degenerate for
# the same reason).
set.seed(21)
n18 <- 40L; cl18 <- rep(1:10, each = 4L)
Z18 <- cbind(rnorm(n18), rnorm(n18), c(1, rep(0, n18 - 1L)))
x18 <- drop(Z18 %*% c(1, -0.5, 2)) + rnorm(n18)
y18 <- x18 + rnorm(n18)
msg18 <- tryCatch(iv_compare(y18, x18, Z18, cluster = cl18, intercept = FALSE),
                  error = function(e) conditionMessage(e))
ok(is.character(msg18) && grepl("leverage is numerically 1", msg18),
   "iv_compare stops when an observation's leverage is 1 (JIVE undefined)")

cat("\nAll clusterIV tests passed.\n")
