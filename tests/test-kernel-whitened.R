# Equivalence tests for the whitened-coordinates leave-cluster-out kernel
# (base R, no testthat), in the style of the other test files. Folded in
# from the retired dev/spot-check-whitened.R.
#
# The pre-refactor .first_stage, .leaveout_fit and .cjar_maxlev bodies are
# frozen below as brute-force reference implementations, copied verbatim from
# R/internals.R as of commit 9624d16 (the state the refactor replaced); they
# are never edited to match the fast path. For designs covering singleton
# clusters, n_g < k, n_g = k, n_g >> k, unbalanced sizes, weights and
# absorbed fixed effects, the shipped kernel
# (.gram_chol / .whiten / .leaveout_fit / .cluster_leverage) must agree with
# the frozen one at 1e-10 relative, and .cluster_leverage() must return
# exactly .leaveout_fit()$maxlev.

library(clusterIV)

rel <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))
ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE)
  cat("PASS:", msg, "\n")
}

# --- frozen reference implementations (pre-refactor, verbatim) ---------------

ref_first_stage <- function(D, Z) {
  Q <- crossprod(Z)
  ch <- tryCatch(chol(Q), error = function(e) NULL)
  if (is.null(ch)) {
    stop("instrument Gram matrix is singular after partialling out covariates: ",
         "the instruments are collinear (possibly with the controls).",
         call. = FALSE)
  }
  Qinv <- chol2inv(ch)
  pihat <- Qinv %*% crossprod(Z, D)
  Xhat <- as.numeric(Z %*% pihat)
  ZQinv <- Z %*% Qinv
  list(Q = Q, Qinv = Qinv, pihat = pihat, Xhat = Xhat,
       e = D - Xhat, ZQinv = ZQinv, h = rowSums(ZQinv * Z), ch = ch)
}

ref_leaveout_fit <- function(D, Z, groups, fs = NULL) {
  if (is.null(fs)) fs <- ref_first_stage(D, Z)
  phat <- fs$Xhat
  maxlev <- 0

  for (idx in groups) {
    ng <- length(idx)
    Zg <- Z[idx, , drop = FALSE]
    Hg <- fs$ZQinv[idx, , drop = FALSE] %*% t(Zg)
    M <- diag(ng) - Hg
    v <- tryCatch(solve(M, fs$e[idx]), error = function(err) NULL)
    if (is.null(v)) {
      stop("leave-cluster-out fit is undefined for a cluster (I - H_g is ",
           "singular): an instrument has no variation outside that cluster. ",
           sprintf("The offending cluster has %d observation(s).", ng),
           call. = FALSE)
    }
    phat[idx] <- fs$Xhat[idx] - as.numeric(Hg %*% v)
    lev <- if (ng == 1L) Hg[1L, 1L]
           else max(eigen(Hg, symmetric = TRUE, only.values = TRUE)$values)
    if (lev > maxlev) maxlev <- lev
  }
  list(phat = phat, maxlev = maxlev)
}

ref_cjar_maxlev <- function(Z, groups, R) {
  maxlev <- 0
  for (idx in groups) {
    Zt <- backsolve(R, t(Z[idx, , drop = FALSE]), transpose = TRUE)  # k x n_g
    lev <- if (min(dim(Zt)) == 1L) sum(Zt^2)
           else svd(Zt, nu = 0, nv = 0)$d[1L]^2
    if (lev > maxlev) maxlev <- lev
  }
  maxlev
}

# --- designs -----------------------------------------------------------------
# Each design supplies raw (y-free) inputs; partialling out runs through the
# package's own .partial_out so the kernel sees exactly what cjive()/cjar()
# feed it. `sizes` gives the cluster sizes; k the instrument columns.

make_design <- function(sizes, k, seed, weights = FALSE, fe = FALSE) {
  set.seed(seed)
  n <- sum(sizes)
  cl <- rep(seq_along(sizes), sizes)
  Z <- matrix(rnorm(n * k), n, k)
  u <- rnorm(length(sizes))[cl]
  x <- drop(Z %*% rep(c(1, -0.6), length.out = k)) + u + rnorm(n)
  w <- if (weights) runif(n, 0.5, 3) else NULL
  f <- if (fe) list(factor(sample(1:5, n, replace = TRUE))) else NULL
  po <- clusterIV:::.partial_out(rnorm(n), x, Z, controls = NULL,
                                 weights = w, intercept = TRUE,
                                 fixed_effects = f)
  list(D = po$x, Z = po$Z, groups = split(seq_len(n), cl))
}

designs <- list(
  `singletons (n = G = 40, k = 4)`   = make_design(rep(1L, 40), 4L, 11),
  `n_g < k (n_g = 3, k = 6)`         = make_design(rep(3L, 12), 6L, 12),
  `n_g = k (n_g = 5, k = 5)`         = make_design(rep(5L, 10), 5L, 13),
  `n_g >> k (n_g = 40, k = 3)`       = make_design(rep(40L, 6), 3L, 14),
  `unbalanced (n_g in 1..12, k = 4)` = make_design(c(1L, 2L, 4L, 7L, 12L, 1L, 9L, 3L), 4L, 15),
  `weights (unbalanced, k = 4)`      = make_design(c(2L, 5L, 8L, 3L, 11L, 6L), 4L, 16, weights = TRUE),
  `fixed effects (n_g = 6, k = 4)`   = make_design(rep(6L, 9), 4L, 17, fe = TRUE)
)

# --- comparison --------------------------------------------------------------

rows <- lapply(names(designs), function(tag) {
  d <- designs[[tag]]
  k <- ncol(d$Z)

  fs_ref  <- ref_first_stage(d$D, d$Z)
  lo_ref  <- ref_leaveout_fit(d$D, d$Z, d$groups, fs_ref)
  lev_ref <- ref_cjar_maxlev(d$Z, d$groups, fs_ref$ch)

  R    <- clusterIV:::.gram_chol(d$Z)
  Ztil <- clusterIV:::.whiten(d$Z, R)
  tt   <- crossprod(Ztil, d$D)
  lo_new  <- clusterIV:::.leaveout_fit(d$D, Ztil, tt, d$groups, k)
  lev_new <- clusterIV:::.cluster_leverage(Ztil, d$groups, k)

  e_phat <- rel(lo_new$phat, lo_ref$phat)
  e_lev  <- rel(lo_new$maxlev, lo_ref$maxlev)
  e_lev2 <- rel(lev_new, lev_ref)

  ok(e_phat <= 1e-10, sprintf("[%s] phat == frozen kernel (%.2e rel)", tag, e_phat))
  ok(e_lev  <= 1e-10, sprintf("[%s] maxlev == frozen .leaveout_fit (%.2e rel)", tag, e_lev))
  ok(e_lev2 <= 1e-10, sprintf("[%s] maxlev == frozen .cjar_maxlev (%.2e rel)", tag, e_lev2))
  ok(identical(lev_new, lo_new$maxlev),
     sprintf("[%s] .cluster_leverage() identical to .leaveout_fit()$maxlev", tag))

  data.frame(design = tag, phat_rel = e_phat, maxlev_rel = e_lev,
             maxlev_vs_cjar_rel = e_lev2, identical_lev = TRUE)
})

cat("\n")
print(do.call(rbind, rows), row.names = FALSE, digits = 3)
cat("\nAll spot-checks passed.\n")
