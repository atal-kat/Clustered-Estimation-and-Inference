# Correctness tests for the cjive package (base R, no testthat).
# Each block stops with an informative message on failure, so a failure fails
# `R CMD check`.

library(cjive)

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
po <- cjive:::.partial_out(y, x, Zj, controls = NULL, weights = NULL, intercept = TRUE)
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
fs <- cjive:::.first_stage(xr, Zr)
phat_pkg <- cjive:::.leaveout_fit(xr, Zr, groups, fs)$phat
phat_bf <- brute_cjive(xr, Zr, groups)
ok(near(phat_pkg, phat_bf, 1e-8), "dense CJIVE p-hat vector == brute force (1e-8)")

beta_pkg <- cjive(y, x, judge, cluster = cl)$coefficient
beta_bf <- sum(phat_bf * yr) / sum(phat_bf * xr)
ok(near(beta_pkg, beta_bf, 1e-8), "dense CJIVE coefficient == brute force (1e-8)")

# === Test 2: singleton clusters: CJIVE == JIVE =============================
sing <- split(seq_len(n), factor(seq_len(n)))
phat_sing <- cjive:::.leaveout_fit(xr, Zr, sing, fs)$phat
phat_jive <- cjive:::.phat_jive(fs, xr)
ok(near(phat_sing, phat_jive, 1e-8), "singleton-cluster CJIVE == JIVE (1e-8)")

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
phat_mean_pkg <- cjive:::.leaveout_mean(x, judge, cl, weights = NULL)
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
poc <- cjive:::.partial_out(yc, xc, Zc, NULL, NULL, TRUE)
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
ok(near(fit_w$coefficient, beta_w, 1e-8) && near(fit_w$se, se_w, 1e-8),
   "weighted CJIVE == independent weighted recomputation (1e-8)")

# === Test 8: input validation (clear, specific errors) =====================
errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
ok(errs(cjive(y, x, judge, cluster = rep(1L, n))), "stops on < 2 clusters")
ok(errs(cjive(y, x, judge, cluster = cl, weights = replace(rep(1, n), 1, -1))),
   "stops on non-positive weights")
ok(errs(cjive(y[-1], x, judge, cluster = cl)), "stops on length mismatch")
ok(errs(cjive(replace(y, 1, NA), x, judge, cluster = cl)), "stops on NA in y")
ok(errs(cjive(y, x, judge, cluster = replace(cl, 1, NA))), "stops on NA in cluster")
ok(errs(cjive(y, x, factor(cl), cluster = cl)),
   "stops when each instrument group lies in a single cluster")

cat("\nAll cjive tests passed.\n")
