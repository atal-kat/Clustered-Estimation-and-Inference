# Numerical-safety oracle for high-dimensional fixed-effect absorption.
# Base R only.  The dense reference below explicitly constructs the weighted
# dummy design and never calls .demean() or uses an alternating-projection
# result.

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

stable_norm <- function(x) {
  ax <- abs(as.numeric(x))
  if (!length(ax)) return(0)
  m <- max(ax)
  if (!is.finite(m) || m == 0) return(m)
  m * sqrt(sum((ax / m)^2))
}

col_relative_error <- function(got, ref, scale) {
  got <- as.matrix(got)
  ref <- as.matrix(ref)
  scale <- as.matrix(scale)
  vapply(seq_len(ncol(got)), function(j) {
    den <- max(stable_norm(scale[, j]), stable_norm(ref[, j]))
    if (den == 0) {
      if (stable_norm(got[, j] - ref[, j]) == 0) 0 else Inf
    } else {
      stable_norm(got[, j] - ref[, j]) / den
    }
  }, numeric(1))
}

# Explicit treatment/reference coding, constructed without model.matrix().
# The intercept plus L-1 columns for every factor spans the joint FE space;
# harmless nesting, disconnected components and duplicated dimensions may
# make the coefficient vector unidentified, but lm.fit's residual is unique.
dense_weighted_fwl <- function(M, fe_list, weights = NULL) {
  M <- as.matrix(M)
  n <- nrow(M)
  w <- if (is.null(weights)) rep.int(1, n) else as.numeric(weights)
  rw <- sqrt(w)
  D <- matrix(rw, n, 1L)
  for (f in fe_list) {
    f <- droplevels(as.factor(f))
    lev <- levels(f)
    if (length(lev) > 1L) {
      block <- vapply(lev[-1L], function(z) as.numeric(f == z), numeric(n))
      D <- cbind(D, block * rw)
    }
  }
  out <- stats::lm.fit(D, M, tol = 1e-12)$residuals
  matrix(out, nrow = n, ncol = ncol(M), dimnames = dimnames(M))
}

dummy_controls <- function(fe_list) {
  n <- length(fe_list[[1L]])
  out <- NULL
  for (f in fe_list) {
    f <- droplevels(as.factor(f))
    lev <- levels(f)
    if (length(lev) > 1L) {
      block <- vapply(lev[-1L], function(z) as.numeric(f == z), numeric(n))
      out <- if (is.null(out)) block else cbind(out, block)
    }
  }
  out
}

capture_error <- function(expr) {
  tryCatch(list(value = expr, error = NULL),
           error = function(e) list(value = NULL,
                                     error = conditionMessage(e)))
}

fe_correction <- function(R, f, weights = NULL) {
  R <- as.matrix(R)
  n <- nrow(R)
  w <- if (is.null(weights)) rep.int(1, n) else as.numeric(weights)
  rw <- sqrt(w)
  f <- droplevels(as.factor(f))
  g <- as.integer(f)
  sw <- as.numeric(rowsum(w, g, reorder = TRUE))
  mu <- rowsum(R * rw, g, reorder = TRUE) / sw
  rw * mu[g, , drop = FALSE]
}

# Independent weighted normal-equation check.  For a transformed residual R,
# every FE level must have rowsum(sqrt(w) * R) equal to zero.  Dividing by the
# square root of the level weight makes the Euclidean score norm comparable
# across highly unequal levels; this is computed directly rather than through
# the production projector or fe_correction().
fe_group_score <- function(R, f, weights = NULL) {
  R <- as.matrix(R)
  n <- nrow(R)
  w <- if (is.null(weights)) rep.int(1, n) else as.numeric(weights)
  rw <- sqrt(w)
  f <- droplevels(as.factor(f))
  g <- as.integer(f)
  sw <- as.numeric(rowsum(w, g, reorder = TRUE))
  sweep(rowsum(R * rw, g, reorder = TRUE), 1L, sqrt(sw), "/")
}

assert_dense_case <- function(label, raw, fe, weights = NULL,
                              tol = 2e-9) {
  rw <- if (is.null(weights)) rep.int(1, nrow(raw)) else sqrt(weights)
  M <- as.matrix(raw) * rw
  oracle <- dense_weighted_fwl(M, fe, weights)
  ans <- capture_error(clusterIV:::.demean(M, fe, weights = weights))
  check(is.null(ans$error), paste(label, "returns a certified result"),
        ans$error)
  if (!is.null(ans$error)) return(invisible(NULL))
  err <- col_relative_error(ans$value, oracle, M)
  check(max(err) <= tol, paste(label, "agrees with dense weighted FWL"),
        paste("max relative column error =", format(max(err), digits = 8)))

  correction <- vapply(fe, function(f) {
    C <- fe_correction(ans$value, f, weights)
    max(col_relative_error(C, matrix(0, nrow(C), ncol(C)), M))
  }, numeric(1))
  score <- vapply(fe, function(f) {
    S <- fe_group_score(ans$value, f, weights)
    max(col_relative_error(S, matrix(0, nrow(S), ncol(S)), M))
  }, numeric(1))
  check(max(correction) <= max(tol, 5e-10),
        paste(label, "has small correction for every FE dimension"),
        paste("max normalized correction =",
              format(max(correction), digits = 8)))
  check(max(score) <= max(tol, 5e-10),
        paste(label, "has small weighted group score for every FE dimension"),
        paste("max normalized group score =", format(max(score), digits = 8)))
  invisible(ans$value)
}

make_columns <- function(n, p = 6L, offset = 0) {
  i <- seq_len(n)
  vapply(seq_len(p), function(j) {
    sin(i * (0.017 * j + 0.003 * offset)) +
      cos(i * (0.031 + 0.007 * j)) +
      sin((i + j * offset) * 0.0013 * (j + 1L))
  }, numeric(n))
}

# -------------------------------------------------------------------------
# A. Frozen one-swap design: the instrument is exactly in the joint FE span.
# -------------------------------------------------------------------------
n_np <- 4000L
f1_np <- factor(rep(1:2, each = n_np / 2L))
f2_np <- f1_np
f2_np[c(1L, n_np / 2L + 1L)] <- f2_np[c(n_np / 2L + 1L, 1L)]
z_np <- as.numeric(f1_np) + 0.37 * as.numeric(f2_np)
M_np <- cbind(z = z_np)
oracle_np <- dense_weighted_fwl(M_np, list(f1_np, f2_np))

np12 <- capture_error(clusterIV:::.demean(M_np, list(f1_np, f2_np)))
np21 <- capture_error(clusterIV:::.demean(M_np, list(f2_np, f1_np)))
clear_stop <- function(x) {
  !is.null(x$error) &&
    grepl("fixed.effect|demean|absorp", x$error, ignore.case = TRUE) &&
    grepl("converg|certif|condition", x$error, ignore.case = TRUE)
}
check(is.null(np12$error) || clear_stop(np12),
      "one-swap order (f1, f2) either certifies or stops clearly", np12$error)
check(is.null(np21$error) || clear_stop(np21),
      "one-swap order (f2, f1) either certifies or stops clearly", np21$error)
check((is.null(np12$error) && is.null(np21$error)) ||
        (!is.null(np12$error) && !is.null(np21$error)),
      "one-swap FE orders make the same return-or-stop decision",
      paste("(f1,f2) =", if (is.null(np12$error)) "RETURN" else "STOP",
            "; (f2,f1) =", if (is.null(np21$error)) "RETURN" else "STOP"))
if (is.null(np12$error)) {
  e12 <- col_relative_error(np12$value, oracle_np, M_np)
  check(max(e12) <= 1e-10,
        "one-swap order (f1, f2) matches the dense projection",
        paste("relative error =", format(max(e12), digits = 9),
              "; residual norm =", format(stable_norm(np12$value), digits = 9),
              "; dense norm =", format(stable_norm(oracle_np), digits = 9)))
}
if (is.null(np21$error)) {
  e21 <- col_relative_error(np21$value, oracle_np, M_np)
  check(max(e21) <= 1e-10,
        "one-swap order (f2, f1) matches the dense projection",
        paste("relative error =", format(max(e21), digits = 9),
              "; residual norm =", format(stable_norm(np21$value), digits = 9)))
}
if (is.null(np12$error) && is.null(np21$error)) {
  eo <- col_relative_error(np12$value, np21$value, M_np)
  check(max(eo) <= 1e-10, "one-swap result is FE-order invariant",
        paste("relative order discrepancy =", format(max(eo), digits = 9)))
}

# -------------------------------------------------------------------------
# B. Frozen public-API design.  Every datum is stored by a closed-form rule;
# no result depends on replaying an RNG stream.  z_np is exactly absorbed.
# -------------------------------------------------------------------------
i_np <- seq_len(n_np)
cl_np <- rep(seq_len(100L), length.out = n_np)
x_np <- 0.8 * z_np + sin(i_np * 0.17) + cos(i_np * 0.031)
y_np <- 1.4 * x_np + cos(i_np * 0.11) + sin(i_np * 0.071)
public <- list(cjive = cjive, iv_compare = iv_compare, cjar = cjar,
               cjscore = cjscore, iv_infer = iv_infer)
absorbed <- lapply(public, function(fn) capture_error(suppressWarnings(
  fn(y_np, x_np, z_np, cluster = cl_np,
     fixed_effects = list(f1_np, f2_np))
)))
absorbed_msg <- vapply(absorbed, function(z) {
  if (is.null(z$error)) "ACCEPT" else z$error
}, character(1))
absorbed_ok <- vapply(absorbed, function(z) {
  if (is.null(z$error)) return(FALSE)
  instrument <- grepl("instrument|`z`", z$error, ignore.case = TRUE) &&
    grepl("absorb|variation|rank|singular|collinear|partial|residu",
          z$error, ignore.case = TRUE)
  solver <- grepl("fixed.effect|demean|absorp", z$error, ignore.case = TRUE) &&
    grepl("converg|certif|condition", z$error, ignore.case = TRUE)
  instrument || solver
}, logical(1))
check(all(absorbed_ok),
      "all five public APIs reject the exactly absorbed one-swap instrument",
      paste(names(absorbed_msg), absorbed_msg, sep = "=", collapse = " | "))
if (all(absorbed_ok)) {
  check(length(unique(absorbed_msg)) == 1L,
        "all five public APIs use one shared absorbed-instrument diagnosis",
        paste(unique(absorbed_msg), collapse = " | "))
}

# A weighted version makes the relevant FE-span direction nearly null under
# S = P1 + P2.  Backward group scores alone are deliberately insufficient:
# the transformed instrument has unit norm but its individual FE corrections
# are O(sqrt(e)).  Its exact dense-FWL residual is algebraically zero because
# it is a known linear combination of the transformed FE dummy columns.
n_wn <- 40L
e_wn <- 1e-200
f1_wn <- factor(rep(c("A", "B"), each = n_wn / 2L))
f2_wn <- f1_wn
f2_wn[c(1L, 21L)] <- f2_wn[c(21L, 1L)]
w_wn <- rep.int(1, n_wn)
w_wn[c(1L, 21L)] <- e_wn
rw_wn <- sqrt(w_wn)
Mz_wn <- numeric(n_wn)
Mz_wn[c(1L, 21L)] <- c(1, -1) / sqrt(2)

wn12 <- capture_error(clusterIV:::.demean(
  cbind(z = Mz_wn), list(f1_wn, f2_wn), weights = w_wn))
wn21 <- capture_error(clusterIV:::.demean(
  cbind(z = Mz_wn), list(f2_wn, f1_wn), weights = w_wn))
wn_clear <- function(z) is.null(z$error) || clear_stop(z)
check(wn_clear(wn12) && wn_clear(wn21),
      "weighted near-null FE-span direction either projects or stops clearly",
      paste(wn12$error, wn21$error, sep = " | "))
check((is.null(wn12$error) && is.null(wn21$error)) ||
        (!is.null(wn12$error) && !is.null(wn21$error)),
      "weighted near-null FE orders make the same return-or-stop decision")
if (is.null(wn12$error) && is.null(wn21$error)) {
  wn_err <- c(stable_norm(wn12$value), stable_norm(wn21$value)) /
    stable_norm(Mz_wn)
  check(max(wn_err) <= 1e-10,
        "weighted near-null FE-span direction is removed in both orders",
        paste("relative residuals =", paste(format(wn_err, digits = 8),
                                             collapse = ", ")))
}

# Freeze the same geometry at every public entry point.  k1 and k2 are exact
# joint-FE-orthogonal transformed-coordinate components; z is exactly in the
# joint transformed-dummy span.  Correct behavior is shared absorbed-z
# rejection, or a shared solver conditioning stop before first-stage work.
k1_wn <- k2_wn <- numeric(n_wn)
k1_wn[c(2L, 3L, 22L, 23L)] <- c(1, -1, 0.7, -0.7)
k2_wn[c(4L, 5L, 24L, 25L)] <- c(0.8, -0.8, 1.1, -1.1)
Mx_wn <- k1_wn + 0.8 * Mz_wn
My_wn <- 1.5 * Mx_wn + k2_wn
weighted_args <- list(
  y = My_wn / rw_wn,
  x = Mx_wn / rw_wn,
  z = Mz_wn / rw_wn,
  cluster = factor(seq_len(n_wn)),
  weights = w_wn,
  fixed_effects = list(f1_wn, f2_wn)
)
weighted_absorbed <- lapply(public, function(fn) capture_error(
  suppressWarnings(do.call(fn, weighted_args))))
weighted_msg <- vapply(weighted_absorbed, function(z) {
  if (is.null(z$error)) "ACCEPT" else z$error
}, character(1))
weighted_ok <- vapply(weighted_absorbed, function(z) {
  if (is.null(z$error)) return(FALSE)
  grepl("instrument|`z`|fixed.effect|demean|absorp", z$error,
        ignore.case = TRUE) &&
    grepl("variation|rank|singular|collinear|partial|residu|converg|certif|condition",
          z$error, ignore.case = TRUE)
}, logical(1))
check(all(weighted_ok),
      "all five public APIs reject the weighted near-null absorbed instrument",
      paste(names(weighted_msg), weighted_msg, sep = "=", collapse = " | "))
if (all(weighted_ok)) {
  check(length(unique(weighted_msg)) == 1L,
        "all five public APIs share the weighted near-null diagnosis",
        paste(unique(weighted_msg), collapse = " | "))
}

# A well-conditioned FE-span component must not mask a near-null FE-span
# component in the same column.  Both are known algebraically to be in the
# transformed dummy span, so the exact residual is zero.  A backward-only
# certificate accepts the slow component after the fast one is removed; the
# mandatory normalized refinement must instead remove it or stop on condition.
e_mix <- 1e-30
f1_mix <- factor(c("a", "a", "b", "b", "c", "c", "d", "d"))
f2_mix <- factor(c("f", "e", "e", "f", "g", "h", "g", "h"))
w_mix <- c(e_mix, 1, e_mix, 1, 1, 1, 1, 1)
slow_mix <- c(1 / sqrt(2), 0, -1 / sqrt(2), 0, 0, 0, 0, 0)
fast_mix <- c(0, 0, 0, 0, 1, -1, 1, -1) / 2
M_mix <- cbind(z = (slow_mix + fast_mix) / sqrt(2))
mix12 <- capture_error(clusterIV:::.demean(
  M_mix, list(f1_mix, f2_mix), weights = w_mix))
mix21 <- capture_error(clusterIV:::.demean(
  M_mix, list(f2_mix, f1_mix), weights = w_mix))
check(wn_clear(mix12) && wn_clear(mix21),
      "mixed fast/near-null FE-span direction either projects or stops clearly",
      paste(mix12$error, mix21$error, sep = " | "))
check((is.null(mix12$error) && is.null(mix21$error)) ||
        (!is.null(mix12$error) && !is.null(mix21$error)),
      "mixed-spectrum FE orders make the same return-or-stop decision")
if (is.null(mix12$error) && is.null(mix21$error)) {
  mix_err <- c(stable_norm(mix12$value), stable_norm(mix21$value)) /
    stable_norm(M_mix)
  check(max(mix_err) <= 1e-10,
        "mixed-spectrum FE-span direction is removed in both orders",
        paste("relative residuals =", paste(format(mix_err, digits = 8),
                                             collapse = ", ")))
}

# Extreme but finite asymmetric weights can erase a weak incidence link by
# cancellation even when a later normalized correction looks well
# conditioned.  This five-edge bipartite graph is a tree, so its joint dummy
# design spans all five observation rows.  M_tree is also constructed
# explicitly from FE coefficients and therefore has exact residual zero.
f1_tree <- factor(c(1L, 5L, 3L, 5L, 3L))
f2_tree <- factor(c(4L, 6L, 4L, 4L, 5L))
w_tree <- 2^c(1000, 700, 700, 700, 1000)
coef_tree <- c(0, 1, -1, 0, 1, -1)
M_tree <- cbind(z = sqrt(w_tree) *
                  (coef_tree[c(1L, 5L, 3L, 5L, 3L)] +
                   coef_tree[c(4L, 6L, 4L, 4L, 5L)]))
tree12 <- capture_error(clusterIV:::.demean(
  M_tree, list(f1_tree, f2_tree), weights = w_tree))
tree21 <- capture_error(clusterIV:::.demean(
  M_tree, list(f2_tree, f1_tree), weights = w_tree))
check(wn_clear(tree12) && wn_clear(tree21),
      "asymmetric extreme-weight FE tree either projects or stops clearly",
      paste(tree12$error, tree21$error, sep = " | "))
check((is.null(tree12$error) && is.null(tree21$error)) ||
        (!is.null(tree12$error) && !is.null(tree21$error)),
      "asymmetric extreme-weight FE tree is order-consistent")
if (is.null(tree12$error) && is.null(tree21$error)) {
  tree_err <- c(stable_norm(tree12$value), stable_norm(tree21$value)) /
    stable_norm(M_tree)
  check(max(tree_err) <= 1e-10,
        "asymmetric extreme-weight FE tree is removed in both orders",
        paste("relative residuals =", paste(format(tree_err, digits = 8),
                                             collapse = ", ")))
}

# The sharper boundary case has a minimum normalized loading just above
# 1e-12, so a loading-only cutoff misses it.  The relevant principal-angle
# gap is quadratic in that loading.  This tree again has full row rank, hence
# every transformed-coordinate column is in the joint FE span.
w_edge <- 2^c(0, -79, -79, -79, 0)
M_edge <- cbind(z = c(1, -2, 3, -4, 5))
edge12 <- capture_error(clusterIV:::.demean(
  M_edge, list(f1_tree, f2_tree), weights = w_edge))
edge21 <- capture_error(clusterIV:::.demean(
  M_edge, list(f2_tree, f1_tree), weights = w_edge))
check(wn_clear(edge12) && wn_clear(edge21),
      "quadratic-gap FE tree either projects or stops clearly",
      paste(edge12$error, edge21$error, sep = " | "))
check((is.null(edge12$error) && is.null(edge21$error)) ||
        (!is.null(edge12$error) && !is.null(edge21$error)),
      "quadratic-gap FE tree is order-consistent")
if (is.null(edge12$error) && is.null(edge21$error)) {
  edge_err <- c(stable_norm(edge12$value), stable_norm(edge21$value)) /
    stable_norm(M_edge)
  check(max(edge_err) <= 1e-10,
        "quadratic-gap full-rank FE tree removes every input column",
        paste("relative residuals =", paste(format(edge_err, digits = 8),
                                             collapse = ", ")))
}

# A global bottleneck need not contain any nearly identical pair of individual
# dummy columns.  Two balanced K2,2 incidence components joined by one tiny
# edge have a near-common direction formed from component-wide contrasts.
# For every positive bridge weight the graph is connected and its joint dummy
# design has full row rank, so the bridge-row unit vector is exactly absorbed.
f1_bridge <- factor(c(1, 1, 2, 2, 3, 3, 4, 4, 1))
f2_bridge <- factor(c(1, 2, 1, 2, 3, 4, 3, 4, 3))
w_bridge <- c(rep(1, 8), 2^-80)
M_bridge <- cbind(z = c(rep(0, 8), 1))
bridge12 <- capture_error(clusterIV:::.demean(
  M_bridge, list(f1_bridge, f2_bridge), weights = w_bridge))
bridge21 <- capture_error(clusterIV:::.demean(
  M_bridge, list(f2_bridge, f1_bridge), weights = w_bridge))
check(wn_clear(bridge12) && wn_clear(bridge21),
      "global weak-bridge FE graph either projects or stops clearly",
      paste(bridge12$error, bridge21$error, sep = " | "))
check((is.null(bridge12$error) && is.null(bridge21$error)) ||
        (!is.null(bridge12$error) && !is.null(bridge21$error)),
      "global weak-bridge FE graph is order-consistent")
if (is.null(bridge12$error) && is.null(bridge21$error)) {
  bridge_err <- c(stable_norm(bridge12$value),
                  stable_norm(bridge21$value)) / stable_norm(M_bridge)
  check(max(bridge_err) <= 1e-10,
        "global weak-bridge full-rank FE graph removes the bridge column",
        paste("relative residuals =", paste(format(bridge_err, digits = 8),
                                             collapse = ", ")))
}

# A three-factor frozen batch guards against factor-order-dependent loss of a
# weak weighted direction.  Both columns are exact sums of transformed dummy
# columns.  A geometry stop must occur before the batch arithmetic can make
# one FE order appear certifiable and another singular.
w_ext3 <- 2^c(300, 100, -100, -100, 20, -300,
              -700, 0, -100, -20, 300, -1000)
f1_ext3 <- factor(c(4, 4, 1, 3, 1, 3, 1, 4, 2, 3, 2, 3))
f2_ext3 <- factor(c(2, 3, 7, 5, 1, 4, 5, 7, 6, 6, 1, 6))
f3_ext3 <- factor(c(4, 4, 2, 1, 6, 3, 5, 1, 4, 2, 4, 4))
v1_ext3 <- c(1.906875180807351, -1.570624981517808,
             -0.7052092221727121, 0.7139262980820532)
v2_ext3 <- c(-0.9427022884070128, -0.6086698336335642,
             1.01897725614122, -0.5321938788002789,
             -1.836863572690253, 0.8039131800193664,
             0.08880404880057666)
v3_ext3 <- c(-0.02949646371690473, -0.8214431807712392,
             1.111992468617097, 2.247735501753645,
             0.8590537810800808, -0.7096421616154794)
pattern_ext3 <- function(v) seq_along(v) %% 3 - 1
M_ext3 <- cbind(
  z1 = sqrt(w_ext3) * (v1_ext3[f1_ext3] + v2_ext3[f2_ext3] +
                         v3_ext3[f3_ext3]),
  z2 = sqrt(w_ext3) * (pattern_ext3(v1_ext3)[f1_ext3] +
                         pattern_ext3(v2_ext3)[f2_ext3] +
                         pattern_ext3(v3_ext3)[f3_ext3]))
ext3_orders <- list(
  forward = list(f1_ext3, f2_ext3, f3_ext3),
  reverse = list(f3_ext3, f2_ext3, f1_ext3),
  rotated = list(f2_ext3, f3_ext3, f1_ext3))
ext3 <- lapply(ext3_orders, function(ff) capture_error(
  clusterIV:::.demean(M_ext3, ff, weights = w_ext3)))
check(all(vapply(ext3, wn_clear, logical(1))),
      "extreme-weight three-factor batch either projects or stops clearly",
      paste(vapply(ext3, function(z) if (is.null(z$error)) "RETURN" else
        z$error, character(1)), collapse = " | "))
ext3_return <- vapply(ext3, function(z) is.null(z$error), logical(1))
check(length(unique(ext3_return)) == 1L,
      "extreme-weight three-factor batch is FE-order consistent",
      paste(names(ext3_return), ifelse(ext3_return, "RETURN", "STOP"),
            collapse = ", "))
if (all(ext3_return)) {
  ext3_err <- vapply(ext3, function(z) {
    max(col_relative_error(z$value,
                           matrix(0, nrow(M_ext3), ncol(M_ext3)), M_ext3))
  }, numeric(1))
  check(max(ext3_err) <= 1e-10,
        "extreme-weight three-factor FE-span batch is removed in every order",
        paste("relative residuals =", paste(format(ext3_err, digits = 8),
                                             collapse = ", ")))
}

# Canonical internal FE ordering prevents finite-precision Krylov behavior
# from turning a correct return in one two-factor order into a maxit stop in
# the reverse order.  This column is an explicit transformed FE combination.
w_canonical <- 2^c(-40, -70, -30, 0, -50, 0)
f1_canonical <- factor(c(3, 1, 1, 2, 3, 2))
f2_canonical <- factor(c(2, 1, 2, 3, 3, 4))
M_canonical <- cbind(z = sqrt(w_canonical) * c(0, -2, 0, 0, -1, 0))
canonical12 <- capture_error(clusterIV:::.demean(
  M_canonical, list(f1_canonical, f2_canonical), weights = w_canonical))
canonical21 <- capture_error(clusterIV:::.demean(
  M_canonical, list(f2_canonical, f1_canonical), weights = w_canonical))
check((is.null(canonical12$error) && is.null(canonical21$error)) ||
        (!is.null(canonical12$error) && !is.null(canonical21$error)),
      "canonical FE ordering makes the maxit decision order-invariant",
      paste(canonical12$error, canonical21$error, sep = " | "))
if (is.null(canonical12$error) && is.null(canonical21$error)) {
  canonical_err <- c(stable_norm(canonical12$value),
                     stable_norm(canonical21$value)) /
    stable_norm(M_canonical)
  check(max(canonical_err) <= 1e-10,
        "canonical-order FE-span column is removed in both user orders",
        paste("relative residuals =",
              paste(format(canonical_err, digits = 8), collapse = ", ")))
}

# Three FE spaces can have a genuinely multi-space near-null direction even
# when no pair of individual dummy columns is nearly identical.  The observed
# Krylov Rayleigh quotient must therefore be large enough to support the
# requested forward tolerance; a small positive value is not sufficient.
w_global3 <- 2^c(-22, -18, 0, -38, -2, -10,
                 -32, 0, -36, -18, -8, -22)
f1_global3 <- factor(c(2, 1, 2, 4, 1, 4, 2, 3, 3, 4, 1, 1))
f2_global3 <- factor(c(1, 1, 3, 1, 1, 2, 1, 1, 4, 2, 4, 4))
f3_global3 <- factor(c(3, 2, 4, 4, 3, 4, 2, 4, 1, 4, 1, 4))
M_global3 <- cbind(z = sqrt(w_global3) *
                     c(0, -3, 0, 2, 0, 0, -3, 0, 0, 0, 0, 0))
global3_orders <- list(
  forward = list(f1_global3, f2_global3, f3_global3),
  reverse = list(f3_global3, f2_global3, f1_global3),
  rotated = list(f2_global3, f3_global3, f1_global3))
global3 <- lapply(global3_orders, function(ff) capture_error(
  clusterIV:::.demean(M_global3, ff, weights = w_global3)))
check(length(unique(vapply(global3, function(z) is.null(z$error),
                           logical(1)))) == 1L,
      "global three-space near-null design is FE-order consistent")
if (all(vapply(global3, function(z) is.null(z$error), logical(1)))) {
  global3_err <- vapply(global3, function(z) {
    stable_norm(z$value) / stable_norm(M_global3)
  }, numeric(1))
  check(max(global3_err) <= 1e-10,
        "global three-space exact FE combination is removed",
        paste("relative residuals =",
              paste(format(global3_err, digits = 8), collapse = ", ")))
} else {
  check(all(vapply(global3, clear_stop, logical(1))),
        "global three-space near-null design stops with conditioning diagnosis",
        paste(vapply(global3, function(z) z$error, character(1)),
              collapse = " | "))
}

# A second full-rank three-FE oracle freezes the case where a backward
# certificate near 1e-25 formerly coexisted with a 3.16e-10 forward error.
ex_full3 <- c(23, 6, 15, -13, 4, 23, -32, 1, 30, -34, -15, 10, 11,
              -35, 31, 2, 4, 13, -7, -9, -21, 3, 22, 9, -2, -37, -29)
w_full3 <- 2^ex_full3
f1_full3 <- factor(c(5, 1, 4, 5, 1, 4, 1, 5, 1, 2, 5, 3, 5, 3,
                     5, 1, 3, 2, 5, 3, 5, 4, 5, 3, 1, 1, 2))
f2_full3 <- factor(c(4, 4, 2, 1, 5, 4, 4, 4, 5, 5, 5, 3, 4, 4,
                     2, 1, 3, 5, 4, 4, 2, 5, 4, 4, 5, 5, 4))
f3_full3 <- factor(c(3, 5, 2, 6, 2, 7, 1, 1, 6, 1, 1, 5, 2, 2,
                     1, 2, 7, 6, 3, 5, 2, 5, 3, 4, 6, 2, 1))
b1_full3 <- c(0.601533408314009, 0.721753444985399,
              -1.9293177801191, -0.492478739985114, 0.477853616230073)
b2_full3 <- c(-2.09712485415776, 0.64074511437132,
              0.146537329553128, 0.532475135907491,
              -0.428063528154498)
b3_full3 <- c(-0.747928533805026, 0.177969336296416,
              0.348134903003661, -1.17802597991584,
              0.414538689049682, 0.657730142759932,
              -0.28242991949567)
M_full3 <- cbind(z = sqrt(w_full3) *
                   (b1_full3[f1_full3] + b2_full3[f2_full3] +
                      b3_full3[f3_full3]))
full3 <- lapply(list(list(f1_full3, f2_full3, f3_full3),
                     list(f3_full3, f1_full3, f2_full3)), function(ff) {
  capture_error(clusterIV:::.demean(M_full3, ff, weights = w_full3))
})
check((all(vapply(full3, function(z) is.null(z$error), logical(1))) ||
        all(vapply(full3, function(z) !is.null(z$error), logical(1)))),
      "full-rank three-FE near-null design is FE-order consistent")
if (all(vapply(full3, function(z) is.null(z$error), logical(1)))) {
  full3_err <- vapply(full3, function(z) stable_norm(z$value) /
                        stable_norm(M_full3), numeric(1))
  check(max(full3_err) <= 1e-10,
        "full-rank three-FE exact combination meets forward tolerance",
        paste("relative residuals =",
              paste(format(full3_err, digits = 8), collapse = ", ")))
} else {
  check(all(vapply(full3, clear_stop, logical(1))),
        "full-rank three-FE near-null design stops with conditioning diagnosis")
}

# -------------------------------------------------------------------------
# C. Ordinary crossed, nested, redundant, disconnected and singleton designs.
# -------------------------------------------------------------------------
n_c <- 360L
i_c <- seq_len(n_c)
raw_c <- make_columns(n_c, 6L, 3L)
colnames(raw_c) <- c("y", "x", "z1", "z2", "c1", "c2")
f1_c <- factor((i_c - 1L) %% 12L)
f2_c <- factor(((i_c - 1L) %/% 7L + 3L * i_c) %% 11L)
f3_c <- factor(((i_c - 1L) %/% 13L + 5L * i_c) %% 9L)
f4_c <- factor(((i_c - 1L) %/% 17L + 7L * i_c) %% 8L)
assert_dense_case("ordinary crossed two-way FE", raw_c, list(f1_c, f2_c))
assert_dense_case("ordinary crossed three-way FE", raw_c,
                  list(f1_c, f2_c, f3_c))

parent_c <- factor(rep(seq_len(8L), length.out = n_c))
child_c <- interaction(parent_c, factor((i_c - 1L) %% 5L), drop = TRUE)
assert_dense_case("two-way nested FE", raw_c, list(parent_c, child_c))
assert_dense_case("two-way nested FE in reverse order", raw_c,
                  list(child_c, parent_c))
assert_dense_case("two identical FE dimensions", raw_c,
                  list(parent_c, parent_c))
parent_relabel_c <- factor(paste0("label-", 9L - as.integer(parent_c)))
assert_dense_case("two relabelled-equivalent FE dimensions", raw_c,
                  list(parent_c, parent_relabel_c))
assert_dense_case("nested and duplicated FE", raw_c,
                  list(parent_c, child_c, parent_c))

component_c <- factor(rep(seq_len(3L), each = n_c / 3L))
left_c <- interaction(component_c, factor((i_c - 1L) %% 5L), drop = TRUE)
right_c <- interaction(component_c,
                       factor(((i_c - 1L) %/% 5L) %% 4L), drop = TRUE)
assert_dense_case("disconnected FE incidence components", raw_c,
                  list(left_c, right_c))

unequal_c <- factor(c("singleton-a", "singleton-b", "singleton-c",
                      rep("large-a", 137L), rep("large-b", 89L),
                      rep("large-c", n_c - 229L)))
partner_c <- factor((i_c^2 + 3L * i_c) %% 17L)
assert_dense_case("singletons and unequal FE group sizes", raw_c,
                  list(unequal_c, partner_c))

# One FE is one exact weighted projection.
w_one <- 2^seq(-12, 12, length.out = n_c)
M_one <- raw_c * sqrt(w_one)
one_prod <- clusterIV:::.demean(M_one, list(f1_c), weights = w_one)
one_ref <- dense_weighted_fwl(M_one, list(f1_c), weights = w_one)
check(max(col_relative_error(one_prod, one_ref, M_one)) <= 5e-13,
      "one FE dimension is an exact one-sweep weighted projection")

# -------------------------------------------------------------------------
# D. Weighted geometry and common rescaling of all weights.
# -------------------------------------------------------------------------
w_c <- 2^seq(-20, 20, length.out = n_c)
weighted_base <- assert_dense_case(
  "heterogeneous weighted three-way FE", raw_c,
  list(f1_c, f2_c, f3_c), weights = w_c, tol = 3e-9)
w_big <- w_c * 2^100
weighted_big <- assert_dense_case(
  "common-rescaled weighted three-way FE", raw_c,
  list(f1_c, f2_c, f3_c), weights = w_big, tol = 3e-9)
if (!is.null(weighted_base) && !is.null(weighted_big)) {
  ew <- col_relative_error(weighted_big / sqrt(2^100), weighted_base,
                           raw_c * sqrt(w_c))
  check(max(ew) <= 3e-9,
        "common positive rescaling of weights leaves the projection unchanged",
        paste("max relative error =", format(max(ew), digits = 8)))
}

# Highly unbalanced total group weight: one level receives almost all mass.
w_unbalanced <- rep(2^-30, n_c)
w_unbalanced[f1_c == levels(f1_c)[1L]] <- 2^30
assert_dense_case("highly unbalanced FE group weights", raw_c,
                  list(f1_c, f2_c), weights = w_unbalanced, tol = 4e-9)

# -------------------------------------------------------------------------
# E. Independent column scaling, including the requested binary exponents.
# -------------------------------------------------------------------------
scale_fe <- list(f1_c, f2_c, f3_c)
M_scale <- raw_c * sqrt(w_c)
scale_base <- capture_error(clusterIV:::.demean(M_scale, scale_fe,
                                                weights = w_c))
check(is.null(scale_base$error), "column-scaling base design is accepted",
      scale_base$error)
if (is.null(scale_base$error)) {
  scales <- c(2^-200, 2^-50, 2^50, 2^200)
  for (j in seq_len(ncol(M_scale))) for (s in scales) {
    Ms <- M_scale
    Ms[, j] <- Ms[, j] * s
    ans <- capture_error(clusterIV:::.demean(Ms, scale_fe, weights = w_c))
    label <- sprintf("column %s scaling by 2^%d", colnames(M_scale)[j],
                     round(log2(s)))
    check(is.null(ans$error), paste(label, "preserves acceptance"), ans$error)
    if (is.null(ans$error)) {
      ans$value[, j] <- ans$value[, j] / s
      er <- col_relative_error(ans$value, scale_base$value, M_scale)
      check(max(er) <= 4e-9, paste(label, "is equivariant"),
            paste("max relative error =", format(max(er), digits = 8)))
    }
  }
}

zero_prod <- assert_dense_case(
  "weighted three-way FE with an exact zero column",
  cbind(raw_c, zero = 0), scale_fe, weights = w_c, tol = 3e-9)
if (!is.null(zero_prod)) {
  check(all(zero_prod[, "zero"] == 0),
        "an exact zero column remains exactly zero")
}

# -------------------------------------------------------------------------
# F. FE-order invariance for two, three and four dimensions.
# -------------------------------------------------------------------------
order_sets <- list(
  two = list(list(f1_c, f2_c), list(f2_c, f1_c)),
  three = list(list(f1_c, f2_c, f3_c), list(f3_c, f1_c, f2_c),
               list(f2_c, f3_c, f1_c)),
  four = list(list(f1_c, f2_c, f3_c, f4_c),
              list(f4_c, f3_c, f2_c, f1_c),
              list(f2_c, f4_c, f1_c, f3_c))
)
for (nm in names(order_sets)) {
  canonical <- dense_weighted_fwl(M_scale, order_sets[[nm]][[1L]], w_c)
  results <- lapply(order_sets[[nm]], function(ff) {
    capture_error(clusterIV:::.demean(M_scale, ff, weights = w_c))
  })
  check(all(vapply(results, function(z) is.null(z$error), logical(1))),
        paste(nm, "FE permutations all certify"),
        paste(vapply(results, function(z) if (is.null(z$error)) "OK" else z$error,
                     character(1)), collapse = " | "))
  if (all(vapply(results, function(z) is.null(z$error), logical(1)))) {
    ee <- vapply(results, function(z)
      max(col_relative_error(z$value, canonical, M_scale)), numeric(1))
    pair <- 0
    for (a in seq_along(results)) for (b in seq_along(results)) {
      pair <- max(pair, max(col_relative_error(results[[a]]$value,
                                                results[[b]]$value, M_scale)))
    }
    check(max(ee) <= 4e-9 && pair <= 4e-9,
          paste(nm, "FE permutations agree with dense FWL and one another"),
          paste("dense error =", format(max(ee), digits = 8),
                "; pairwise error =", format(pair, digits = 8)))
  }
}

# -------------------------------------------------------------------------
# G. Idempotence and direct per-dimension correction checks.
# -------------------------------------------------------------------------
R_id <- capture_error(clusterIV:::.demean(M_scale, scale_fe, weights = w_c))
if (is.null(R_id$error)) {
  RR_id <- capture_error(clusterIV:::.demean(R_id$value, scale_fe,
                                             weights = w_c))
  check(is.null(RR_id$error), "demeaning a certified residual certifies again",
        RR_id$error)
  if (is.null(RR_id$error)) {
    eid <- col_relative_error(RR_id$value, R_id$value, M_scale)
    check(max(eid) <= 5e-10, "certified HDFE projection is idempotent",
          paste("max relative error =", format(max(eid), digits = 8)))
  }
  corr <- vapply(scale_fe, function(f) {
    C <- fe_correction(R_id$value, f, w_c)
    max(col_relative_error(C, matrix(0, nrow(C), ncol(C)), M_scale))
  }, numeric(1))
  check(max(corr) <= 5e-10,
        "every individual FE correction is below the certificate tolerance",
        paste("max correction =", format(max(corr), digits = 8)))
  score <- vapply(scale_fe, function(f) {
    S <- fe_group_score(R_id$value, f, w_c)
    max(col_relative_error(S, matrix(0, nrow(S), ncol(S)), M_scale))
  }, numeric(1))
  check(max(score) <= 5e-10,
        "every weighted within-level group score is below tolerance",
        paste("max group score =", format(max(score), digits = 8)))
}

# -------------------------------------------------------------------------
# H. Iteration cap and internal argument validation.
# -------------------------------------------------------------------------
cap <- capture_error(clusterIV:::.demean(raw_c,
                                         list(f1_c, f2_c, f3_c, f4_c),
                                         tol = 1e-14, maxit = 1L))
cap_ok <- !is.null(cap$error) &&
  grepl("converg|certif", cap$error, ignore.case = TRUE) &&
  grepl("1", cap$error) &&
  grepl("diagnostic|change|orthogon|residual", cap$error, ignore.case = TRUE) &&
  grepl("tol|tolerance|1e-14", cap$error, ignore.case = TRUE)
check(cap_ok,
      "iteration-cap error reports sweeps, final diagnostic and tolerance",
      if (is.null(cap$error)) "ACCEPT" else cap$error)

bad_internal <- list(
  `tol zero` = quote(clusterIV:::.demean(raw_c, list(f1_c), tol = 0)),
  `tol negative` = quote(clusterIV:::.demean(raw_c, list(f1_c), tol = -1e-8)),
  `tol missing` = quote(clusterIV:::.demean(raw_c, list(f1_c), tol = NA_real_)),
  `tol character` = quote(clusterIV:::.demean(raw_c, list(f1_c), tol = "1e-8")),
  `tol non-scalar` = quote(clusterIV:::.demean(raw_c, list(f1_c),
                                               tol = c(1e-8, 1e-9))),
  `tol non-finite` = quote(clusterIV:::.demean(raw_c, list(f1_c), tol = Inf)),
  `maxit zero` = quote(clusterIV:::.demean(raw_c, list(f1_c), maxit = 0L)),
  `maxit negative` = quote(clusterIV:::.demean(raw_c, list(f1_c), maxit = -1L)),
  `maxit missing` = quote(clusterIV:::.demean(raw_c, list(f1_c), maxit = NA_real_)),
  `maxit non-finite` = quote(clusterIV:::.demean(raw_c, list(f1_c), maxit = Inf)),
  `maxit character` = quote(clusterIV:::.demean(raw_c, list(f1_c), maxit = "2")),
  `maxit non-scalar` = quote(clusterIV:::.demean(raw_c, list(f1_c),
                                                 maxit = c(1L, 2L))),
  `maxit fractional` = quote(clusterIV:::.demean(raw_c, list(f1_c),
                                                 maxit = 1.5)),
  `non-finite M` = quote(clusterIV:::.demean(replace(raw_c, 1L, Inf),
                                             list(f1_c))),
  `zero-column M` = quote(clusterIV:::.demean(matrix(numeric(), n_c, 0L),
                                              list(f1_c))),
  `wrong FE length` = quote(clusterIV:::.demean(raw_c, list(f1_c[-1L]))),
  `missing FE` = quote(clusterIV:::.demean(raw_c,
                                           list(replace(f1_c, 1L, NA)))),
  `non-finite weights` = quote(clusterIV:::.demean(raw_c, list(f1_c),
                                                   weights = replace(w_c, 1L, Inf)))
)
for (nm in names(bad_internal)) {
  ans <- capture_error(eval(bad_internal[[nm]]))
  check(!is.null(ans$error), paste("internal validation rejects", nm),
        if (is.null(ans$error)) "ACCEPT" else ans$error)
}

# -------------------------------------------------------------------------
# I. Downstream HDFE versus explicit dense-dummy equivalence.
# -------------------------------------------------------------------------
n_d <- 240L
i_d <- seq_len(n_d)
cl_d <- rep(seq_len(30L), each = 8L)
fe1_d <- factor((i_d - 1L) %% 8L)
fe2_d <- factor(((i_d - 1L) %/% 8L + 3L * i_d) %% 7L)
fe_d <- list(fe1_d, fe2_d)
Z_d <- cbind(sin(i_d * 0.13) + cos(i_d * 0.031),
             cos(i_d * 0.17) + sin(i_d * 0.047))
c_d <- sin(i_d * 0.071) + cos(i_d * 0.019)
u_d <- sin(cl_d * 0.37)
x_d <- drop(Z_d %*% c(0.9, -0.45)) + 0.3 * c_d + u_d +
  sin(i_d * 0.29)
y_d <- 0.65 * x_d + 0.4 * c_d + u_d + cos(i_d * 0.23)
w_d <- 0.4 + (1 + (i_d %% 17L))^2 / 30
D_d <- dummy_controls(fe_d)

h_args <- list(y = y_d, x = x_d, z = Z_d, cluster = cl_d,
               controls = cbind(c_d), fixed_effects = fe_d, weights = w_d)
d_args <- list(y = y_d, x = x_d, z = Z_d, cluster = cl_d,
               controls = cbind(c_d, D_d), weights = w_d)

num_near <- function(a, b, tol) {
  length(a) == length(b) && identical(dim(a), dim(b)) &&
    all((is.na(a) & is.na(b)) |
          (is.infinite(a) & is.infinite(b) & sign(a) == sign(b)) |
          (is.finite(a) & is.finite(b) &
             abs(a - b) <= tol * pmax(1, abs(a), abs(b))))
}
fields_near <- function(a, b, fields, tol = 2e-7) {
  all(vapply(fields, function(nm) {
    x <- a[[nm]]
    y <- b[[nm]]
    if (is.numeric(x) && is.numeric(y)) num_near(x, y, tol)
    else identical(x, y)
  }, logical(1)))
}

hd_cj <- suppressWarnings(do.call(cjive, h_args))
dd_cj <- suppressWarnings(do.call(cjive, d_args))
check(fields_near(hd_cj, dd_cj,
                  c("coefficient", "se", "statistic", "p.value",
                    "conf.low", "conf.high", "maxlev", "F_eff", "K_eff",
                    "F_eff_crit"), 2e-8),
      "cjive HDFE route matches the explicit dense-dummy route")

hd_tab <- suppressWarnings(do.call(iv_compare, h_args))
dd_tab <- suppressWarnings(do.call(iv_compare, d_args))
check(identical(hd_tab$estimator, dd_tab$estimator) &&
        num_near(as.matrix(hd_tab[, -1L]), as.matrix(dd_tab[, -1L]), 2e-8),
      "iv_compare HDFE route matches the explicit dense-dummy route")

hd_ar <- suppressWarnings(do.call(cjar, h_args))
dd_ar <- suppressWarnings(do.call(cjar, d_args))
check(fields_near(hd_ar, dd_ar,
                  c("statistic", "p.value", "crit", "conf_set", "shape",
                    "bounded", "F_CJ", "maxlev", "F_eff", "K_eff",
                    "F_eff_crit", "coef_num", "coef_var"), 2e-7),
      "cjar HDFE statistics, coefficients and set match dense dummies")

hd_sc <- suppressWarnings(do.call(cjscore, h_args))
dd_sc <- suppressWarnings(do.call(cjscore, d_args))
check(fields_near(hd_sc, dd_sc,
                  c("statistic", "p.value", "score", "variance", "crit",
                    "conf_set", "shape", "bounded", "F_CJS", "maxlev",
                    "F_eff", "K_eff", "F_eff_crit", "coef_score",
                    "coef_var"), 2e-7),
      "cjscore HDFE statistics, coefficients and set match dense dummies")

all_panel_tests <- c("cjar", "cjscore")
hd_panel <- suppressWarnings(do.call(iv_infer,
                                     c(h_args, list(tests = all_panel_tests))))
dd_panel <- suppressWarnings(do.call(iv_infer,
                                     c(d_args, list(tests = all_panel_tests))))
panel_ok <- fields_near(hd_panel$cjive, dd_panel$cjive,
                        c("coefficient", "se", "statistic", "p.value",
                          "conf.low", "conf.high", "F_eff"), 2e-8) &&
  fields_near(hd_panel$cjar, dd_panel$cjar,
              c("statistic", "p.value", "conf_set", "shape", "F_CJ",
                "coef_num", "coef_var"), 2e-7) &&
  fields_near(hd_panel$cjscore, dd_panel$cjscore,
              c("statistic", "p.value", "conf_set", "shape", "F_CJS",
                "coef_score", "coef_var"), 2e-7)
check(panel_ok,
      "iv_infer HDFE estimate and jackknife tests match dense dummies")

# A dense control exactly in the joint FE span must be treated as absorbed,
# not promoted from projection noise into a rank-one QR regressor.  The check
# is weighted and uses extreme control units so an absolute cutoff cannot pass.
fe_span_d <- as.numeric(fe1_d) + 0.37 * as.numeric(fe2_d)
for (s in c(2^-200, 2^200)) {
  r_args <- h_args
  r_args$controls <- cbind(c_d, fe_span_d * s)
  red_cj <- suppressWarnings(do.call(cjive, r_args))
  check(fields_near(red_cj, hd_cj,
                    c("coefficient", "se", "statistic", "p.value",
                      "conf.low", "conf.high", "maxlev", "F_eff", "K_eff",
                      "F_eff_crit"), 2e-10),
        sprintf("joint-FE-span control at scale 2^%d is numerically redundant",
                round(log2(s))))
}

r_args <- h_args
r_args$controls <- cbind(c_d, fe_span_d * 2^200)
red_tab <- suppressWarnings(do.call(iv_compare, r_args))
check(identical(red_tab$estimator, hd_tab$estimator) &&
        num_near(as.matrix(red_tab[, -1L]), as.matrix(hd_tab[, -1L]), 2e-10),
      "iv_compare drops a control absorbed by the joint FE span")

red_ar <- suppressWarnings(do.call(cjar, r_args))
check(fields_near(red_ar, hd_ar,
                  c("statistic", "p.value", "crit", "conf_set", "shape",
                    "bounded", "F_CJ", "maxlev", "F_eff", "K_eff",
                    "F_eff_crit", "coef_num", "coef_var"), 2e-9),
      "cjar drops a control absorbed by the joint FE span")

red_sc <- suppressWarnings(do.call(cjscore, r_args))
check(fields_near(red_sc, hd_sc,
                  c("statistic", "p.value", "score", "variance", "crit",
                    "conf_set", "shape", "bounded", "F_CJS", "maxlev",
                    "F_eff", "K_eff", "F_eff_crit", "coef_score",
                    "coef_var"), 2e-9),
      "cjscore drops a control absorbed by the joint FE span")

red_panel <- suppressWarnings(do.call(
  iv_infer, c(r_args, list(tests = all_panel_tests))))
red_panel_ok <- fields_near(red_panel$cjive, hd_panel$cjive,
                            c("coefficient", "se", "statistic", "p.value",
                              "conf.low", "conf.high", "F_eff"), 2e-10) &&
  fields_near(red_panel$cjar, hd_panel$cjar,
              c("statistic", "p.value", "conf_set", "shape", "F_CJ",
                "coef_num", "coef_var"), 2e-9) &&
  fields_near(red_panel$cjscore, hd_panel$cjscore,
              c("statistic", "p.value", "conf_set", "shape", "F_CJS",
                "coef_score", "coef_var"), 2e-9)
check(red_panel_ok,
      "iv_infer drops a control absorbed by the joint FE span in every row")

# A linear combination may be absorbed while neither individual control is.
# After FE removal c1 and c2 span only residual(v); the unit-scaled pivoted QR
# must therefore agree with supplying v alone, independent of their raw scale.
v_joint_d <- sin(i_d * 0.113) + cos(i_d * 0.057)
c1_joint_d <- v_joint_d + fe_span_d
c2_joint_d <- v_joint_d - fe_span_d
joint_args <- h_args
joint_args$controls <- cbind(c_d, c1_joint_d * 2^80,
                             c2_joint_d * 2^-80)
one_args <- h_args
one_args$controls <- cbind(c_d, v_joint_d)
joint_fit <- suppressWarnings(do.call(cjive, joint_args))
one_fit <- suppressWarnings(do.call(cjive, one_args))
check(fields_near(joint_fit, one_fit,
                  c("coefficient", "se", "statistic", "p.value",
                    "conf.low", "conf.high", "maxlev", "F_eff", "K_eff",
                    "F_eff_crit"), 2e-9),
      "pivoted QR handles a jointly absorbed control combination across units")

# -------------------------------------------------------------------------
# J. Cluster fixed effects and cluster-nested fixed effects: the regimes
#    the CJAR/CJS theory covers unchanged (?cjar, Ligtenberg 2025 Section
#    5.3).  On exactly these designs the HDFE route must equal the
#    dense-dummy FWL route for the estimate and both jackknife tests.
# -------------------------------------------------------------------------
clfe_d <- factor(cl_d)
# Two sub-units per cluster; every nest level lies wholly inside one cluster,
# so the nested dummies alone span the cluster dummies (nesting dedup).
nest_d <- factor((cl_d - 1L) * 2L + ((i_d - 1L) %/% 4L) %% 2L)
cluster_fe_cases <- list(
  `cluster FE` = list(fe = list(clfe_d),
                      dummies = dummy_controls(list(clfe_d))),
  `cluster-nested FE` = list(fe = list(clfe_d, nest_d),
                             dummies = dummy_controls(list(nest_d)))
)
for (nm in names(cluster_fe_cases)) {
  cse <- cluster_fe_cases[[nm]]
  hj <- list(y = y_d, x = x_d, z = Z_d, cluster = cl_d,
             controls = cbind(c_d), fixed_effects = cse$fe, weights = w_d)
  dj <- list(y = y_d, x = x_d, z = Z_d, cluster = cl_d,
             controls = cbind(c_d, cse$dummies), weights = w_d)
  cj_h <- suppressWarnings(do.call(cjive, hj))
  cj_dd <- suppressWarnings(do.call(cjive, dj))
  check(fields_near(cj_h, cj_dd,
                    c("coefficient", "se", "statistic", "p.value",
                      "conf.low", "conf.high", "maxlev", "F_eff", "K_eff",
                      "F_eff_crit"), 2e-10),
        paste("cjive:", nm, "HDFE route equals dense dummies"))
  ar_h <- suppressWarnings(do.call(cjar, hj))
  ar_dd <- suppressWarnings(do.call(cjar, dj))
  check(fields_near(ar_h, ar_dd,
                    c("statistic", "p.value", "crit", "conf_set", "shape",
                      "bounded", "F_CJ", "maxlev", "coef_num", "coef_var"),
                    2e-10),
        paste("cjar:", nm, "HDFE route equals dense dummies"))
  sc_h <- suppressWarnings(do.call(cjscore, hj))
  sc_dd <- suppressWarnings(do.call(cjscore, dj))
  check(fields_near(sc_h, sc_dd,
                    c("statistic", "p.value", "score", "variance", "crit",
                      "conf_set", "shape", "bounded", "F_CJS", "maxlev",
                      "coef_score", "coef_var"), 2e-10),
        paste("cjscore:", nm, "HDFE route equals dense dummies"))
}

if (length(.failures)) {
  stop(paste0("HDFE safety oracle failures (", length(.failures), "):\n - ",
              paste(.failures, collapse = "\n - ")), call. = FALSE)
}

cat("\nAll HDFE numerical-safety tests passed.\n")
