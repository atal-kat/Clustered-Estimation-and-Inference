# Focused numerical-safety oracles for confidence-set inversion (base R).
#
# This file was created before the repair under R/.  Every failing design is
# stored literally or by a closed-form construction: no expected value depends
# on replaying an unstored random-data generator.  The direct-acceptance
# helpers below are deliberately independent of the production inverters and
# of .assemble_set().

library(clusterIV)

.failures <- character(0)
check <- function(cond, msg, detail = NULL) {
  if (isTRUE(cond)) {
    cat("PASS:", msg, "\n")
  } else {
    .failures <<- c(.failures, msg)
    cat("FAIL:", msg, "\n")
    if (!is.null(detail)) cat("      ", detail, "\n", sep = "")
  }
  invisible(cond)
}

relnear <- function(a, b, tol = 1e-12) {
  length(a) == length(b) &&
    all((is.na(a) & is.na(b)) |
          (is.infinite(a) & is.infinite(b) & sign(a) == sign(b)) |
          (is.finite(a) & is.finite(b) &
             abs(a - b) <= tol * pmax(1, abs(a), abs(b))))
}

set_near <- function(A, B, tol = 1e-10) {
  identical(dim(A), dim(B)) &&
    all(is.infinite(A) == is.infinite(B)) &&
    relnear(A[is.finite(A)], B[is.finite(B)], tol)
}

in_set <- function(beta, set) {
  nrow(set) > 0L && any(beta >= set[, 1L] & beta <= set[, 2L])
}

# Affine image of a union of closed intervals.  Each row is re-ordered before
# the rows themselves are sorted, so negative changes of units are covered.
map_set <- function(set, mult = 1, shift = 0) {
  if (!nrow(set)) return(set)
  ans <- t(apply(set, 1L, function(z) sort(mult * z + shift)))
  ans <- ans[order(ans[, 1L], ans[, 2L]), , drop = FALSE]
  colnames(ans) <- c("lower", "upper")
  ans
}

# Universal interval/gap/tail oracle.  `accept` is always a direct predicate,
# never a production inverter.  Boundary points may be zeros of either the
# test equation or the variance equation, hence boundary_residual returns the
# smaller appropriately scaled residual.
audit_set <- function(set, accept, boundary_residual, beta0 = NULL,
                      beta0_accept = NULL, label = "set", tol = 1e-8) {
  ordered <- nrow(set) == 0L ||
    (all(set[, 1L] <= set[, 2L]) &&
       (nrow(set) == 1L || all(set[-nrow(set), 2L] < set[-1L, 1L])))
  check(ordered, paste(label, "has ordered, non-overlapping components"))

  ep <- sort(unique(as.numeric(set[is.finite(set)])))
  if (length(ep)) {
    er <- vapply(ep, boundary_residual, numeric(1))
    check(all(is.finite(er) & er <= tol),
          paste(label, "finite endpoints have small direct residual"),
          paste("max residual", format(max(er), digits = 8)))
    check(all(vapply(ep, function(b)
      accept(b) || boundary_residual(b) <= tol, logical(1))),
          paste(label, "finite endpoints obey the closed-set convention"))
  }

  if (nrow(set)) for (i in seq_len(nrow(set))) {
    lo <- set[i, 1L]; hi <- set[i, 2L]
    probe <- if (lo == hi) lo
             else if (!is.finite(lo) && !is.finite(hi)) 0
             else if (!is.finite(lo)) hi - max(1, abs(hi))
             else if (!is.finite(hi)) lo + max(1, abs(lo))
             else .5 * lo + .5 * hi
    check(accept(probe), sprintf("%s component %d has an accepted probe",
                                 label, i))
  }
  if (nrow(set) > 1L) for (i in seq_len(nrow(set) - 1L)) {
    probe <- .5 * set[i, 2L] + .5 * set[i + 1L, 1L]
    check(!accept(probe), sprintf("%s omitted gap %d has a rejected probe",
                                  label, i))
  }

  if (length(ep)) {
    tail <- c(ep[1L] - max(1, abs(ep[1L])),
              ep[length(ep)] + max(1, abs(ep[length(ep)])) )
  } else {
    tail <- c(-1, 1)
  }
  for (b in tail) {
    check(in_set(b, set) == accept(b),
          sprintf("%s direct and returned-set tail membership agree", label))
  }

  # Carefully scaled one-sided neighbours catch close-root swaps without a
  # production grid.  Skip a side when it cannot be separated from its next
  # endpoint in double precision.
  if (length(ep)) for (j in seq_along(ep)) {
    gap <- min(c(if (j > 1L) ep[j] - ep[j - 1L] else Inf,
                 if (j < length(ep)) ep[j + 1L] - ep[j] else Inf))
    h <- min(sqrt(.Machine$double.eps) * max(1, abs(ep[j])), gap / 4)
    if (is.finite(h) && h > 0 && ep[j] - h < ep[j] && ep[j] + h > ep[j]) {
      for (b in c(ep[j] - h, ep[j] + h)) {
        check(in_set(b, set) == accept(b),
              sprintf("%s neighbour membership agrees near endpoint %d",
                      label, j))
      }
    }
  }

  if (!is.null(beta0)) {
    check(in_set(beta0, set) == isTRUE(beta0_accept),
          paste(label, "beta0 membership agrees with its stored point test"))
  }
  invisible(TRUE)
}

# Direct cross-fit-CJAR acceptance from the stored coefficient contract.
# No root finder or set assembler enters this predicate.
cjar_accept <- function(beta, nc, wc, k, crit, nonpos_var = FALSE) {
  Q <- nc[1L] - nc[2L] * beta + nc[3L] * beta^2
  V <- wc[1L] + wc[2L] * beta + wc[3L] * beta^2 +
    wc[4L] * beta^3 + wc[5L] * beta^4
  if (nonpos_var && V <= 0) return(TRUE)
  V > 0 && Q <= crit * sqrt(k * V)
}

cjar_boundary_residual <- function(beta, nc, wc, k, crit) {
  jp <- 0:4
  Q <- nc[1L] - nc[2L] * beta + nc[3L] * beta^2
  V <- sum(wc * beta^jp)
  qscale <- sum(abs(c(nc[1L], -nc[2L], nc[3L])) * abs(beta)^(0:2))
  vscale <- sum(abs(wc) * abs(beta)^jp)
  vr <- abs(V) / max(vscale, .Machine$double.xmin)
  tr <- if (V > 0) {
    rhs <- crit * sqrt(k * V)
    abs(Q - rhs) / max(qscale + abs(rhs), .Machine$double.xmin)
  } else Inf
  min(vr, tr)
}

# Direct CJS acceptance from its two empirical polynomials.
cjs_accept <- function(beta, sc, vc, crit) {
  S <- sc[1L] - sc[2L] * beta
  V <- vc[1L] + vc[2L] * beta + vc[3L] * beta^2
  V <= 0 || S^2 <= crit * V
}


cjs_boundary_residual <- function(beta, sc, vc, crit) {
  S <- sc[1L] - sc[2L] * beta
  V <- vc[1L] + vc[2L] * beta + vc[3L] * beta^2
  A <- S^2 - crit * V
  ss <- abs(sc[1L]) + abs(sc[2L] * beta)
  vs <- sum(abs(vc) * abs(beta)^(0:2))
  ar <- abs(A) / max(ss^2 + abs(crit) * vs, .Machine$double.xmin)
  vr <- abs(V) / max(vs, .Machine$double.xmin)
  min(ar, vr)
}

# ---------------------------------------------------------------------------
# Frozen public design: n = 8, G = 4, k = 2, intercept = FALSE.  It exposes
# the signed-leading-crossfit failure while also freezing all four public test
# paths at an ordinary, finite null.
# ---------------------------------------------------------------------------
y_cf <- c(0.0718093708707297, 0.899504774473639, -0.450540547480756,
          -0.227023913379762, 0.640353180167933, 1.29289226957527,
          1.33131391377535, -1.26642849427646)
x_cf <- c(-0.685020471279976, 1.1945418925033, -0.368431300969284,
          1.55466058814939, 1.39812904409569, 1.09368838671946,
          -0.034207047192113, -1.45455095866683)
Z_cf <- structure(c(
  -0.626453810742332, 0.183643324222082, -0.835628612410047,
  1.59528080213779, 0.32950777181536, -0.820468384118015,
  0.487429052428485, 0.738324705129217, 0.575781351653492,
  -0.305388387156356, 1.51178116845085, 0.389843236411431,
  -0.621240580541804, -2.2146998871775, 1.12493091814311,
  -0.0449336090152309), dim = c(8L, 2L))
cl_cf <- rep(seq_len(4L), each = 2L)
b0_cf <- 1.135
quiet <- function(expr) suppressWarnings(expr)

ar_plain <- quiet(cjar(y_cf, x_cf, Z_cf, cluster = cl_cf,
                       intercept = FALSE, beta0 = b0_cf))
ar_cf <- quiet(cjar(y_cf, x_cf, Z_cf, cluster = cl_cf,
                    intercept = FALSE, beta0 = b0_cf,
                    variance = "crossfit"))
sc_plain <- quiet(cjscore(y_cf, x_cf, Z_cf, cluster = cl_cf,
                          intercept = FALSE, beta0 = b0_cf))
sc_cf <- quiet(cjscore(y_cf, x_cf, Z_cf, cluster = cl_cf,
                       intercept = FALSE, beta0 = b0_cf,
                       variance = "crossfit"))

# Point statistics, p-values and coefficient kernels are outside the repair.
# The frozen values below are a regression guard, compared at 1e-10 relative
# (the package-wide oracle gate) rather than bit for bit: chol/backsolve/eigen
# are not guaranteed identical across BLAS/platforms, so an identical() pin
# would ERROR on CRAN's test farm on a legitimate ulp-level difference. The
# exact bit-level freeze lives in the private dev regression scripts.
check(relnear(ar_plain$statistic, 0x1.3fa099e6a09cep-4, 1e-10) &&
        relnear(ar_plain$p.value, 0x1.5c6e1e6e84cc5p-2, 1e-10) &&
        relnear(ar_plain$coef_num,
                c(-0x1.58ed8ea3edc9ep-1, -0x1.08b2b90233887p+0,
                  -0x1.32af38801188p-2), 1e-10) &&
        relnear(ar_plain$coef_var,
                c(0x1.5910aa5615faep-1, -0x1.25add35715fap+0,
                  0x1.3c4b26367601p+0, -0x1.4cf80c5f2239cp+0,
                  0x1.3582a349dcf6p+0), 1e-10) &&
        set_near(ar_plain$conf_set,
                 rbind(c(-Inf, Inf)), 1e-10),
      "frozen plain cjar() statistic, p-value, coefficients and ordinary set")
check(relnear(ar_cf$statistic, 0x1.af7e4d13fb089p-4, 1e-10) &&
        relnear(ar_cf$p.value, 0x1.530afdac713cp-2, 1e-10) &&
        relnear(ar_cf$coef_num,
                c(-0x1.58ed8ea3edc9ep-1, -0x1.08b2b90233887p+0,
                  -0x1.32af38801188p-2), 1e-10) &&
        relnear(ar_cf$coef_var,
                c(0x1.284c6355fdc3ep+0, -0x1.f6b98653647fcp-2,
                  0x1.1c9e9995505bp-6, 0x1.e09ff6ccebba8p-4,
                  -0x1.00b6128e81304p-3), 1e-10) &&
        set_near(ar_cf$conf_set,
                 rbind(c(-Inf, 1.64230745651949),
                       c(1.64575185389322, Inf)), 1e-10),
      "frozen crossfit cjar() statistic, p-value, coefficients and ordinary set")
check(relnear(sc_plain$statistic, 0x1.989eef3b068efp-6, 1e-10) &&
        relnear(sc_plain$p.value, 0x1.bfc09201d39e3p-1, 1e-10) &&
        relnear(sc_plain$coef_score,
                c(-0x1.08b2b90233887p-1, -0x1.32af38801188p-2), 1e-10) &&
        relnear(sc_plain$coef_var,
                c(0x1.c3f257700199ep-1, -0x1.09f2d14d10657p+1,
                  0x1.0f81ea47d9b1p+1), 1e-10) &&
        set_near(sc_plain$conf_set,
                 rbind(c(-Inf, Inf)), 1e-10),
      "frozen plain cjscore() statistic, p-value, coefficients and ordinary set")
check(relnear(sc_cf$statistic, 0x1.5caf249ec024p-5, 1e-10) &&
        relnear(sc_cf$p.value, 0x1.ac50141c0ed42p-1, 1e-10) &&
        relnear(sc_cf$coef_score,
                c(-0x1.08b2b90233887p-1, -0x1.32af38801188p-2), 1e-10) &&
        relnear(sc_cf$coef_var,
                c(0x1.5a607a6895531p-1, -0x1.2100571cb0196p+1,
                  0x1.04974e1fbfcf2p+1), 1e-10) &&
        set_near(sc_cf$conf_set,
                 rbind(c(-Inf, Inf)), 1e-10),
      "frozen crossfit cjscore() statistic, p-value, coefficients and ordinary set")

# A. Crossfit CJAR: the base set has two very close finite boundaries.  After
# x -> a*x the point statistic is invariant and the set must divide by a.
a_cf <- 1e8
ar_cf_scaled <- quiet(cjar(y_cf, a_cf * x_cf, Z_cf, cluster = cl_cf,
                           intercept = FALSE, beta0 = b0_cf / a_cf,
                           variance = "crossfit"))
check(relnear(ar_cf_scaled$statistic, ar_cf$statistic, 1e-12) &&
        relnear(ar_cf_scaled$p.value, ar_cf$p.value, 1e-12),
      "crossfit CJAR point statistic and p-value are invariant to x units")
check(set_near(ar_cf_scaled$conf_set * a_cf, ar_cf$conf_set, 1e-8),
      "crossfit CJAR confidence set is scale-equivariant",
      paste("base shape", ar_cf$shape, "scaled shape", ar_cf_scaled$shape))
for (b in c(1.64230745651949, mean(c(1.64230745651949,
                                     1.64575185389322)),
            1.64575185389322)) {
  check(in_set(b, ar_cf$conf_set) ==
          cjar_accept(b, ar_cf$coef_num, ar_cf$coef_var,
                      ar_cf$k, ar_cf$crit, nonpos_var = TRUE),
        sprintf("crossfit CJAR direct membership at beta = %.15g", b))
}

# The same public design across ordinary and extreme finite units.  Normalise
# transformed sets back to the base beta coordinate before comparing them.
x_scales <- c(1e-16, -1e-16, 1e-8, -1e-8,
              1e8, -1e8, 1e16, -1e16)
for (a in x_scales) {
  fit <- quiet(cjar(y_cf, a * x_cf, Z_cf, cluster = cl_cf,
                    intercept = FALSE, beta0 = b0_cf / a,
                    variance = "crossfit"))
  tag <- sprintf("crossfit CJAR x units a = %s", format(a, scientific = TRUE))
  check(relnear(fit$statistic, ar_cf$statistic, 2e-11) &&
          relnear(fit$p.value, ar_cf$p.value, 2e-11),
        paste(tag, "preserve the point test"))
  check(set_near(map_set(fit$conf_set, a), ar_cf$conf_set, 2e-8),
        paste(tag, "preserve the confidence set in normalised units"))
}

y_scales <- x_scales
for (a in y_scales) {
  fit <- quiet(cjar(a * y_cf, x_cf, Z_cf, cluster = cl_cf,
                    intercept = FALSE, beta0 = a * b0_cf,
                    variance = "crossfit"))
  tag <- sprintf("crossfit CJAR y units a = %s", format(a, scientific = TRUE))
  check(relnear(fit$statistic, ar_cf$statistic, 2e-11) &&
          relnear(fit$p.value, ar_cf$p.value, 2e-11),
        paste(tag, "preserve the point test"))
  check(set_near(map_set(fit$conf_set, 1 / a), ar_cf$conf_set, 2e-8),
        paste(tag, "preserve the confidence set in normalised units"))
}

for (cc in c(0.75, 100, -100)) {
  fit <- quiet(cjar(y_cf + cc * x_cf, x_cf, Z_cf, cluster = cl_cf,
                    intercept = FALSE, beta0 = b0_cf + cc,
                    variance = "crossfit"))
  check(relnear(fit$statistic, ar_cf$statistic, 2e-7) &&
          set_near(map_set(fit$conf_set, 1, -cc), ar_cf$conf_set, 2e-7),
        sprintf("crossfit CJAR affine shift c = %g is equivariant", cc))
}

M_cf <- matrix(c(2, -0.4, 0.3, 1.5), 2L, 2L)
fit_zm <- quiet(cjar(y_cf, x_cf, Z_cf %*% M_cf, cluster = cl_cf,
                     intercept = FALSE, beta0 = b0_cf,
                     variance = "crossfit"))
check(relnear(fit_zm$statistic, ar_cf$statistic, 2e-11) &&
        relnear(fit_zm$p.value, ar_cf$p.value, 2e-11) &&
        set_near(fit_zm$conf_set, ar_cf$conf_set, 2e-9),
      "crossfit CJAR is invariant to a nonsingular instrument transformation")

perm_cf <- c(8L, 1L, 6L, 3L, 5L, 2L, 7L, 4L)
relab_cf <- c("north", "east", "south", "west")
fit_perm <- quiet(cjar(y_cf[perm_cf], x_cf[perm_cf], Z_cf[perm_cf, ],
                       cluster = relab_cf[cl_cf][perm_cf],
                       intercept = FALSE, beta0 = b0_cf,
                       variance = "crossfit"))
check(relnear(fit_perm$statistic, ar_cf$statistic, 2e-11) &&
        relnear(fit_perm$p.value, ar_cf$p.value, 2e-11) &&
        set_near(fit_perm$conf_set, ar_cf$conf_set, 2e-9),
      "crossfit CJAR is invariant to row permutation and cluster relabeling")

audit_set(
  ar_cf$conf_set,
  function(b) cjar_accept(b, ar_cf$coef_num, ar_cf$coef_var,
                          ar_cf$k, ar_cf$crit, TRUE),
  function(b) cjar_boundary_residual(b, ar_cf$coef_num, ar_cf$coef_var,
                                     ar_cf$k, ar_cf$crit),
  beta0 = b0_cf,
  beta0_accept = ar_cf$statistic <= ar_cf$crit,
  label = "crossfit CJAR base set")

# Powers near 2^+/-500 are exercised on a degree-degenerate but non-trivial
# CJAR tuple, where all transformed coefficients and roots remain finite.
nc_pow <- c(2, 1, 0)
wc_pow <- c(1, 0, 0, 0, 0)
pow_base <- clusterIV:::.cjar_invert(nc_pow, wc_pow, 1, 1, "accept")
check(set_near(pow_base$conf_set, cbind(1, Inf), 1e-14),
      "degree-degenerate CJAR power-of-two base tuple is [1, Inf)")
for (a in c(2^-500, -2^-500, 2^500, -2^500)) {
  inv <- clusterIV:::.cjar_invert(c(2, a, 0), wc_pow, 1, 1, "accept")
  check(set_near(map_set(inv$conf_set, a), pow_base$conf_set, 2e-12),
        sprintf("degree-degenerate CJAR is equivariant at a = %s",
                format(a, scientific = TRUE)))
}

# Public crossfit-CJS metamorphic design.  Every datum is a closed-form
# function of its row index; no random-number generator is involved.  Its
# confidence set has two bounded components, so each transformation exercises
# endpoint construction rather than only a whole-line special case.
n_cfs <- 30L
i_cfs <- seq_len(n_cfs)
cl_cfs <- rep(seq_len(10L), each = 3L)
Z_cfs <- cbind(
  sin(i_cfs * .171) + cos(i_cfs * .031),
  cos(i_cfs * .232) + sin(i_cfs * .071),
  sin(i_cfs * .11 + .07) + cos(i_cfs * .291))
x_cfs <- drop(Z_cfs %*% c(.8, -.4, .3)) +
  sin(i_cfs * .3707) + cos(i_cfs * .13)
y_cfs <- .7 * x_cfs + drop(Z_cfs %*% c(.2, .1, -.15)) +
  cos(i_cfs * .4109) + sin(i_cfs * .19)
b0_cfs <- .2
sc_cfs <- quiet(cjscore(y_cfs, x_cfs, Z_cfs, cluster = cl_cfs,
                         intercept = FALSE, beta0 = b0_cfs,
                         variance = "crossfit"))
check(nrow(sc_cfs$conf_set) == 2L && all(is.finite(sc_cfs$conf_set)),
      "public crossfit CJS metamorphic design has two finite components")
audit_set(
  sc_cfs$conf_set,
  function(b) cjs_accept(b, sc_cfs$coef_score, sc_cfs$coef_var, sc_cfs$crit),
  function(b) cjs_boundary_residual(
    b, sc_cfs$coef_score, sc_cfs$coef_var, sc_cfs$crit),
  beta0 = b0_cfs,
  beta0_accept = sc_cfs$variance <= 0 || sc_cfs$statistic <= sc_cfs$crit,
  label = "public crossfit CJS base set", tol = 2e-8)

for (a in c(1e-16, -1e-16, 1e-8, -1e-8,
            1e8, -1e8, 1e16, -1e16)) {
  fit_x <- quiet(cjscore(y_cfs, a * x_cfs, Z_cfs, cluster = cl_cfs,
                          intercept = FALSE, beta0 = b0_cfs / a,
                          variance = "crossfit"))
  fit_y <- quiet(cjscore(a * y_cfs, x_cfs, Z_cfs, cluster = cl_cfs,
                          intercept = FALSE, beta0 = a * b0_cfs,
                          variance = "crossfit"))
  tag <- format(a, scientific = TRUE)
  check(relnear(fit_x$statistic, sc_cfs$statistic, 2e-10) &&
          relnear(fit_x$p.value, sc_cfs$p.value, 2e-10) &&
          set_near(map_set(fit_x$conf_set, a), sc_cfs$conf_set, 2e-6),
        paste("public crossfit CJS x-unit equivariance at a =", tag))
  check(relnear(fit_y$statistic, sc_cfs$statistic, 2e-10) &&
          relnear(fit_y$p.value, sc_cfs$p.value, 2e-10) &&
          set_near(map_set(fit_y$conf_set, 1 / a), sc_cfs$conf_set, 2e-6),
        paste("public crossfit CJS y-unit equivariance at a =", tag))
}

for (cc in c(.75, 100, -100)) {
  fit <- quiet(cjscore(y_cfs + cc * x_cfs, x_cfs, Z_cfs,
                        cluster = cl_cfs, intercept = FALSE,
                        beta0 = b0_cfs + cc, variance = "crossfit"))
  check(relnear(fit$statistic, sc_cfs$statistic, 1e-8) &&
          relnear(fit$p.value, sc_cfs$p.value, 1e-8) &&
          set_near(map_set(fit$conf_set, 1, -cc),
                   sc_cfs$conf_set, 2e-7),
        sprintf("public crossfit CJS affine shift c = %g is equivariant", cc))
}

M_cfs <- matrix(c(2, -.4, .3, 1.5, .2, -.1, .4, .2, 1.3), 3L, 3L)
fit_cfs_zm <- quiet(cjscore(y_cfs, x_cfs, Z_cfs %*% M_cfs,
                             cluster = cl_cfs, intercept = FALSE,
                             beta0 = b0_cfs, variance = "crossfit"))
check(relnear(fit_cfs_zm$statistic, sc_cfs$statistic, 2e-10) &&
        relnear(fit_cfs_zm$p.value, sc_cfs$p.value, 2e-10) &&
        set_near(fit_cfs_zm$conf_set, sc_cfs$conf_set, 2e-6),
      "public crossfit CJS is invariant to nonsingular ZM")

perm_cfs <- c(30L, 1L, 17L, 4L, 25L, 8L, 12L, 6L, 21L, 3L,
              28L, 2L, 15L, 10L, 23L, 19L, 5L, 27L, 11L, 16L,
              7L, 24L, 14L, 29L, 9L, 18L, 26L, 13L, 22L, 20L)
labs_cfs <- paste0("cluster-", c(8L, 2L, 10L, 1L, 9L,
                                  4L, 6L, 3L, 7L, 5L))
fit_cfs_perm <- quiet(cjscore(
  y_cfs[perm_cfs], x_cfs[perm_cfs], Z_cfs[perm_cfs, ],
  cluster = labs_cfs[cl_cfs][perm_cfs], intercept = FALSE,
  beta0 = b0_cfs, variance = "crossfit"))
check(relnear(fit_cfs_perm$statistic, sc_cfs$statistic, 2e-10) &&
        relnear(fit_cfs_perm$p.value, sc_cfs$p.value, 2e-10) &&
        set_near(fit_cfs_perm$conf_set, sc_cfs$conf_set, 2e-6),
      "public crossfit CJS is invariant to rows and cluster labels")

# B. Degree-degenerate CJS: the tiny linear score coefficient creates a real
# far endpoint, not a ray.  The companion ordinary-scale tuple is frozen too.
inv_cjs0 <- clusterIV:::.cjs_invert(c(1, 1), c(1, 0, 0), 1)
inv_cjs1 <- clusterIV:::.cjs_invert(c(1, 1e-16), c(1, 0, 0), 1)
check(set_near(inv_cjs0$conf_set, cbind(0, 2), 1e-14),
      "CJS ordinary linear/constant tuple remains [0, 2]")
check(set_near(inv_cjs1$conf_set, cbind(0, 2e16), 1e-12),
      "CJS preserves a finite far root when the linear coefficient is tiny",
      paste("returned", paste(inv_cjs1$conf_set, collapse = ", ")))
check(cjs_accept(1e16, c(1, 1e-16), c(1, 0, 0), 1) &&
        !cjs_accept(3e16, c(1, 1e-16), c(1, 0, 0), 1),
      "CJS direct oracle separates the far interval from its rejected tail")

# Near cancellation in the leading coefficient is mathematical information,
# not interpolation noise.  These two stored-double tuples share
# A(b) = -2*b + delta*b^2 with delta = 1 - (1 - 1e-15), hence the finite far
# boundary 2/delta must survive in both low-degree inverters.
d_cancel <- 1e-15
delta_cancel <- 1 - (1 - d_cancel)
far_cancel <- 2 / delta_cancel
inv_cancel_cjs <- clusterIV:::.cjs_invert(
  c(1, 1), c(1, 0, 1 - d_cancel), 1)
check(set_near(inv_cancel_cjs$conf_set, cbind(0, far_cancel), 2e-12) &&
        !cjs_accept(3e15, c(1, 1), c(1, 0, 1 - d_cancel), 1),
      "CJS retains a genuine far root after leading-coefficient cancellation")
inv_cancel_cjar <- clusterIV:::.cjar_invert(
  c(-1, -1, 0), c(1, 0, 1 - d_cancel, 0, 0), 1, 1, "accept")
check(set_near(inv_cancel_cjar$conf_set, cbind(-Inf, far_cancel), 2e-12) &&
        !cjar_accept(3e15, c(-1, -1, 0),
                     c(1, 0, 1 - d_cancel, 0, 0), 1, 1, TRUE),
      "CJAR retains a genuine far root after leading-coefficient cancellation")

# Every feasible empirical-degree combination of A = S^2 - crit*v and v,
# plus zero, repeated/coincident-root, non-positive-variance and isolated-point
# cases.  The tuples are literal and their sets are audited only through the
# direct inequalities above.
cjs_cases <- list(
  quadratic_quadratic = list(s = c(1, 1), v = c(.2, .1, .3), crit = 1),
  quadratic_linear = list(s = c(1, 1), v = c(1, .2, 0), crit = 1),
  linear_quadratic = list(s = c(1, 1), v = c(.5, 0, 1), crit = 1),
  linear_linear = list(s = c(1, 0), v = c(.5, .2, 0), crit = 1),
  constant_quadratic = list(s = c(1, 1), v = c(.5, -2, 1), crit = 1),
  constant_linear = list(s = c(1, 0), v = c(.5, .2, 0), crit = 0),
  zero_A = list(s = c(1, 1), v = c(1, -2, 1), crit = 1),
  zero_v = list(s = c(1, 1), v = c(0, 0, 0), crit = 1),
  repeated_A_root = list(s = c(0, -1), v = c(-1, 2, 0), crit = 1),
  coincident_A_v_root = list(s = c(0, -1), v = c(0, 1, 0), crit = 1),
  nonpositive_variance_tails = list(s = c(2, 0),
                                    v = c(1, 0, -1), crit = 1),
  isolated_accepted_point = list(s = c(0, -1),
                                 v = c(-1, 2, 0), crit = 1)
)
for (nm in names(cjs_cases)) {
  z <- cjs_cases[[nm]]
  inv <- clusterIV:::.cjs_invert(z$s, z$v, z$crit)
  audit_set(
    inv$conf_set,
    function(b) cjs_accept(b, z$s, z$v, z$crit),
    function(b) cjs_boundary_residual(b, z$s, z$v, z$crit),
    label = paste("CJS", nm), tol = 2e-8)
}

# The isolated repeated root at beta = 1 must survive as [1,1], separately
# from the non-positive-variance ray (-Inf, 0.5].
iso <- clusterIV:::.cjs_invert(c(0, -1), c(-1, 2, 0), 1)
check(set_near(iso$conf_set, rbind(c(-Inf, .5), c(1, 1)), 2e-12),
      "CJS preserves an isolated accepted repeated root")

# Exact changes of x units include roots near both ends of double precision.
cjs_pow_base <- clusterIV:::.cjs_invert(c(2, 1), c(1, 0, 0), 1)
check(set_near(cjs_pow_base$conf_set, cbind(1, 3), 1e-14),
      "CJS power-of-two base tuple is [1, 3]")
for (a in c(1e-16, -1e-16, 1e-8, -1e-8,
            1e8, -1e8, 1e16, -1e16,
            2^-500, -2^-500, 2^500, -2^500)) {
  inv <- clusterIV:::.cjs_invert(c(2, a), c(1, 0, 0), 1)
  check(set_near(map_set(inv$conf_set, a), cjs_pow_base$conf_set, 2e-12),
        sprintf("CJS x-unit equivariance at a = %s",
                format(a, scientific = TRUE)))
}

# Shared assembler: topology and closed endpoints are tested independently of
# any polynomial inverter.  The rays-plus-middle example intentionally does
# not assert `shape`: the current public vocabulary cannot describe it
# truthfully, and adding a new value requires the maintainer's API decision.
assembly_cases <- list(
  empty = list(r = numeric(0), f = function(x) FALSE, l = FALSE, u = FALSE,
               set = matrix(numeric(0), 0L, 2L)),
  whole = list(r = numeric(0), f = function(x) TRUE, l = TRUE, u = TRUE,
               set = cbind(-Inf, Inf)),
  bounded = list(r = c(-1, 1), f = function(x) abs(x) <= 1,
                 l = FALSE, u = FALSE, set = cbind(-1, 1)),
  multiple = list(r = c(-3, -2, 1, 2),
                  f = function(x) (x >= -3 && x <= -2) ||
                    (x >= 1 && x <= 2),
                  l = FALSE, u = FALSE,
                  set = rbind(c(-3, -2), c(1, 2))),
  ray = list(r = 0, f = function(x) x >= 0, l = FALSE, u = TRUE,
             set = cbind(0, Inf)),
  two_rays = list(r = c(-1, 1), f = function(x) abs(x) >= 1,
                  l = TRUE, u = TRUE,
                  set = rbind(c(-Inf, -1), c(1, Inf))),
  rays_plus_middle = list(
    r = c(-3, -2, -1, 1),
    f = function(x) x <= -3 || (x >= -2 && x <= -1) || x >= 1,
    l = TRUE, u = TRUE,
    set = rbind(c(-Inf, -3), c(-2, -1), c(1, Inf))),
  isolated = list(r = 0, f = function(x) x == 0, l = FALSE, u = FALSE,
                  set = cbind(0, 0)),
  repeated = list(r = c(-1, -1, 1), f = function(x) abs(x) <= 1,
                  l = FALSE, u = FALSE, set = cbind(-1, 1)),
  narrow_gap = list(r = c(0, 1e-12),
                    f = function(x) x <= 0 || x >= 1e-12,
                    l = TRUE, u = TRUE,
                    set = rbind(c(-Inf, 0), c(1e-12, Inf)))
)
for (nm in names(assembly_cases)) {
  z <- assembly_cases[[nm]]
  got <- clusterIV:::.assemble_set(z$r, z$f, z$l, z$u,
                                   root_accept = z$f)$conf_set
  check(set_near(got, z$set, 1e-14),
        paste("shared set assembler handles", nm))
}

# === Regression: G = 2 plain CJAR sets (boundary on {Vhat = 0}) ============
# At G = 2 the plain variance is a single square, T(b) = sign(c_12(b)), and
# the acceptance boundary sits exactly on the zeros of Vhat.  The external
# audit of 2026-07-18 found the plain inverter dropped those boundary points
# (variance roots were pooled only on the cross-fit path), returning an
# empty set against a non-empty pointwise acceptance region for every level
# with c_a < 1.  Variance roots are now pooled on both paths and numerically
# non-positive Vhat mirrors the point-statistic rule (T := 0, accepted iff
# 0 <= c_a).
set.seed(107)
cl_g2 <- sort(rep(1:2, length.out = 40))
Z_g2 <- matrix(rnorm(80), 40, 2)
u_g2 <- rnorm(2)[cl_g2]
x_g2 <- drop(Z_g2 %*% rep(0.7, 2)) + u_g2 + rnorm(40)
y_g2 <- 0.5 * x_g2 + u_g2 + rnorm(40)
for (lv_g2 in c(0.5, 0.8, 0.95)) {
  f_g2 <- quiet(cjar(y_g2, x_g2, Z_g2, cluster = cl_g2, level = lv_g2))
  bs_g2 <- seq(-8, 8, by = 0.01)
  agree <- vapply(bs_g2, function(b) {
    ep <- as.numeric(f_g2$conf_set); ep <- ep[is.finite(ep)]
    if (length(ep) && min(abs(b - ep)) < 1e-6) return(TRUE)
    Q <- f_g2$coef_num[1] - f_g2$coef_num[2] * b + f_g2$coef_num[3] * b^2
    V <- sum(f_g2$coef_var * b^(0:4))
    acc <- if (V <= 0) f_g2$crit >= 0 else Q / sqrt(f_g2$k * V) <= f_g2$crit
    inn <- nrow(f_g2$conf_set) > 0 &&
      any(b >= f_g2$conf_set[, 1] & b <= f_g2$conf_set[, 2])
    acc == inn
  }, logical(1))
  check(all(agree),
        sprintf("G = 2 plain CJAR set matches its own statistic pointwise (level %.2f)", lv_g2))
}
f_g2a <- quiet(cjar(y_g2, x_g2, Z_g2, cluster = cl_g2, level = 0.8))
check(f_g2a$shape == "bounded" && nrow(f_g2a$conf_set) == 1L,
      "G = 2, level 0.8: one bounded interval (was: empty)")
f_g2b <- quiet(cjar(y_g2, x_g2, Z_g2, cluster = cl_g2, level = 0.95))
check(f_g2b$shape == "whole_line" && nrow(f_g2b$conf_set) == 1L,
      "G = 2, level 0.95: whole line in one row (no split at variance roots)")

# === Regression: coefficient-kernel scale guards ===========================
# The variance coefficients are quartic in the data scale.  Overflow used to
# crash with a bare NaN comparison; underflow silently returned an all-zero
# variance polynomial, i.e. statistic 0 / p ~ 0.4 next to an empty set.
# Both directions must now be one informative error naming the remedy, and
# the wide representable band must be untouched.
errs_with <- function(expr) {
  tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
}
set.seed(101)
cl_sc <- sort(rep(1:20, length.out = 160))
Z_sc <- matrix(rnorm(640), 160, 4)
u_sc <- rnorm(20)[cl_sc]
x_sc <- drop(Z_sc %*% rep(0.7, 4)) + u_sc + rnorm(160)
y_sc <- 0.5 * x_sc + u_sc + rnorm(160)
f_sc0 <- quiet(cjar(y_sc, x_sc, Z_sc, cluster = cl_sc))
for (s_sc in c(1e80, 1e-90)) {
  msg_a <- errs_with(quiet(cjar(y_sc * s_sc, x_sc * s_sc, Z_sc, cluster = cl_sc)))
  msg_s <- errs_with(quiet(cjscore(y_sc * s_sc, x_sc * s_sc, Z_sc, cluster = cl_sc)))
  check(grepl("not representable", msg_a) && grepl("scale-invariant", msg_a),
        sprintf("cjar at joint scale %g: informative scale error", s_sc))
  check(grepl("not representable", msg_s),
        sprintf("cjscore at joint scale %g: informative scale error", s_sc))
}
msg_cf <- errs_with(quiet(cjar(y_sc * 1e80, x_sc * 1e80, Z_sc, cluster = cl_sc,
                               variance = "crossfit")))
check(grepl("not representable", msg_cf),
      "cross-fit at joint scale 1e80: informative scale error")
f_sc1 <- quiet(cjar(y_sc * 1e60, x_sc * 1e60, Z_sc, cluster = cl_sc))
check(set_near(f_sc1$conf_set, f_sc0$conf_set, 1e-8) &&
        abs(f_sc1$statistic - f_sc0$statistic) < 1e-9,
      "joint scale 1e60 stays on the unchanged raw path (same set/statistic)")
f_keff <- quiet(cjive(y_sc, x_sc * 1e80, Z_sc, cluster = cl_sc))
f_keff0 <- quiet(cjive(y_sc, x_sc, Z_sc, cluster = cl_sc))
check(is.finite(f_keff$K_eff) &&
        abs(f_keff$K_eff - f_keff0$K_eff) < 1e-8 * f_keff0$K_eff,
      "K_eff survives x * 1e80 via the normalized fallback (scale-invariant)")

# === Regression: extreme finite point nulls and p-value curves =============
# The coefficient polynomials are representable on this ordinary design, but
# direct beta^4 evaluation is not.  Before the signed-log point helper, CJAR
# silently returned T = 0 at beta = 1e100 and crashed at 1e150; CJS and both
# plot curves became NaN once beta^2 overflowed.  Inversion itself was already
# scale-conditioned, so point membership must agree with the stored set too.
for (var_ext in c("plain", "crossfit")) {
  for (b_ext in c(-1e100, 1e100, -1e160, 1e160)) {
    ar_ext <- quiet(cjar(y_sc, x_sc, Z_sc, cluster = cl_sc, beta0 = b_ext,
                         variance = var_ext))
    sc_ext <- quiet(cjscore(y_sc, x_sc, Z_sc, cluster = cl_sc, beta0 = b_ext,
                            variance = var_ext))

    check(!is.nan(ar_ext$statistic) && is.finite(ar_ext$p.value) &&
            relnear(clusterIV:::.pval_curve_cjar(ar_ext, b_ext),
                    ar_ext$p.value, 2e-12),
          sprintf("%s CJAR extreme beta %g: stable point and curve agree",
                  var_ext, b_ext))
    ar_accept <- is.na(ar_ext$statistic) || ar_ext$statistic <= ar_ext$crit
    check(identical(ar_accept, in_set(b_ext, ar_ext$conf_set)),
          sprintf("%s CJAR extreme beta %g: point membership agrees with set",
                  var_ext, b_ext))

    check(!is.nan(sc_ext$statistic) && is.finite(sc_ext$p.value) &&
            relnear(clusterIV:::.pval_curve_cjs(sc_ext, b_ext),
                    sc_ext$p.value, 2e-12),
          sprintf("%s CJS extreme beta %g: stable point and curve agree",
                  var_ext, b_ext))
    sc_accept <- is.na(sc_ext$statistic) || sc_ext$statistic <= sc_ext$crit
    check(identical(sc_accept, in_set(b_ext, sc_ext$conf_set)),
          sprintf("%s CJS extreme beta %g: point membership agrees with set",
                  var_ext, b_ext))

    # With a positive leading variance coefficient, the tails converge to the
    # stored leading-coefficient diagnostics at this scale.
    if (ar_ext$coef_var[5L] > 0) {
      check(relnear(ar_ext$statistic, ar_ext$F_CJ, 2e-12),
            sprintf("%s CJAR extreme beta %g: statistic reaches F_CJ tail",
                    var_ext, b_ext))
    }
    if (sc_ext$coef_var[3L] > 0) {
      check(relnear(sc_ext$statistic, sc_ext$F_CJS^2, 2e-12),
            sprintf("%s CJS extreme beta %g: LM reaches F_CJS^2 tail",
                    var_ext, b_ext))
    }
  }
}

# Degree-degenerate tuples are especially important: multiplying a zero
# quartic coefficient by beta^4 = Inf used to manufacture NaN even though the
# empirical polynomial was only linear or quadratic.
ar_deg <- clusterIV:::.cjar_point_eval(
  c(1, 2, 0), c(1, 0, 1, 0, 0), 4, 1e200, "normal", "zero"
)
check(relnear(ar_deg$statistic, -1, 2e-12) && is.finite(ar_deg$p.value),
      "degree-degenerate CJAR tuple remains evaluable beyond beta^4 overflow")
sc_deg <- clusterIV:::.cjs_point_eval(c(1, 2), c(1, 0, 0), 1e200)
check(is.infinite(sc_deg$statistic) && sc_deg$p.value == 0,
      "degree-degenerate CJS tuple returns its infinite LM limit, not NaN")

if (length(.failures)) {
  stop(length(.failures), " inversion-safety check(s) failed: ",
       paste(.failures, collapse = "; "), call. = FALSE)
}
cat("\nAll focused inversion-safety tests passed.\n")
