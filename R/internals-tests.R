# Internal helpers: the cluster jackknife test kernels -- per-cluster
# sums, CJAR and CJS polynomial coefficients, point evaluation,
# confidence-set inversion, and the cross-fit variance estimator
# (Ligtenberg 2025).

# Whitened per-cluster instrument sums and their O(Gk) contractions, shared
# by the CJAR and CJS coefficient kernels: one C-level rowsum pass per
# variable, then whitening by a triangular backsolve, so columns of Yt, Xt
# are R^-T Z_[g]' y_[g] etc. (k x G). R is the upper-triangular Cholesky
# factor of Z'Z (shared with the CJIVE path); when NULL it is factorised
# here. Each caller forms only the k x k Gram contractions it needs.
.cluster_sums <- function(y, x, Z, cluster, R = NULL) {
  Z <- as.matrix(Z)
  k <- ncol(Z)
  cl <- droplevels(as.factor(cluster))
  G <- nlevels(cl)

  if (is.null(R)) R <- .gram_chol(Z)

  Yt <- backsolve(R, t(rowsum(Z * y, cl)), transpose = TRUE)
  Xt <- backsolve(R, t(rowsum(Z * x, cl)), transpose = TRUE)

  # The coefficient kernels need only the two whitened cluster matrices and
  # their dimensions.  Their raw/normalized fallbacks form the required
  # row/column contractions in one place; caching duplicates here was unused.
  list(Yt = Yt, Xt = Xt, k = k, G = G)
}

# CJAR polynomial coefficients (Ligtenberg 2025, Section 3 and p. 8). On the
# partialled data, Q(b) = eps(b)' Pddot eps(b) = n0 - n1*b + n2*b^2 and
# Vhat(b) = (2/k) sum_{g != h} c_gh(b)^2 = w0 + w1*b + ... + w4*b^4, where
# c_gh(b) = eps_[g](b)' P_Z[g,h] eps_[h](b). Every pairwise scalar is an inner
# product of the whitened per-cluster instrument sums from .cluster_sums, so
# the g != h sums collapse to k x k trace contractions: O(nk + Gk^2) time,
# O(Gk + k^2) memory, no G-fold R loop and no explicit (Z'Z)^-1.
# Takes the .cluster_sums() object itself, so one per-cluster-sums pass can
# feed both this kernel and .cjs_coefs (the iv_infer() path).
# With a non-NULL `cf` (a .crossfit_var() result) the cross-fit variance
# polynomial replaces the plain one; the statistic coefficients n0..n2 are
# untouched -- cross-fitting changes only the variance (Ligtenberg 2025,
# Section 5.4).
.cjar_coefs <- function(cs, cf = NULL) {
  k <- cs$k
  if (!is.null(cf) && is.null(cf$w)) {
    stop("internal error: the supplied cross-fit object carries no AR ",
         "variance polynomial (built with which = \"score\").",
         call. = FALSE)
  }

  # One pass of the full closed-form kernel for a given (Yt, Xt) pair; the
  # O(Gk) vectors are recomputed so that the same code serves the raw-scale
  # call (bit-identical to the historical path) and the normalized fallback.
  kernel <- function(Yt, Xt) {
    sy <- rowSums(Yt); sx <- rowSums(Xt)
    agg <- colSums(Yt^2); dgg <- colSums(Xt^2); mgg <- colSums(Xt * Yt)

    # Numerator coefficients in closed form, O(Gk).
    n0 <- sum(sy^2) - sum(agg)
    n1 <- 2 * (sum(sx * sy) - sum(mgg))
    n2 <- sum(sx^2) - sum(dgg)

    # k x k Gram contractions give the full g,h double sums; the diagonal
    # g = h terms are subtracted with the O(Gk) vectors above.  The names
    # record the g,h scalars each contraction sums, via the trace
    # identities on the k x k products (NOT the elementwise squares one
    # might guess): sum(AXY^2) = tr(Yt Xt' Xt Yt') = sum_{g,h} a_gh d_gh,
    # and sum(AXX * AYY) = tr((Xt Xt')(Yt Yt')) = ||Yt'Xt||_F^2
    # = sum_{g,h} m_gh^2.  The two enter w2 with a common factor; the names
    # below follow these trace identities.
    AYY <- tcrossprod(Yt)
    AXX <- tcrossprod(Xt)
    AXY <- tcrossprod(Xt, Yt)
    S_aa <- sum(AYY^2)
    S_dd <- sum(AXX^2)
    S_ad <- sum(AXY^2)            # sum over g,h of a_gh * d_gh
    S_mm <- sum(AXX * AYY)        # sum over g,h of m_gh^2
    S_mM <- sum(AXY * t(AXY))     # sum over g,h of m_gh * m_hg
    S_am <- sum(AXY * AYY)        # sum over g,h of a_gh * m_gh
    S_dm <- sum(AXX * AXY)        # sum over g,h of d_gh * m_gh

    w0 <- (2 / k) * (S_aa - sum(agg^2))
    w1 <- (-8 / k) * (S_am - sum(agg * mgg))
    w2 <- (2 / k) * (2 * S_mm + 2 * S_mM + 2 * S_ad -
                       sum(4 * mgg^2 + 2 * agg * dgg))
    w3 <- (-8 / k) * (S_dm - sum(dgg * mgg))
    w4 <- (2 / k) * (S_dd - sum(dgg^2))
    list(n = c(n0, n1, n2), w = c(w0, w1, w2, w3, w4),
         wref = (2 / k) * (S_aa + S_dd) + 1e-300)
  }

  # Sum-of-squares sign guard and clamp for a plain-path kernel result.
  # Runs only on finite values (the caller checks finiteness first), so the
  # comparison can no longer see NaN.
  guard <- function(z) {
    if (z$w[1L] < -1e-8 * z$wref || z$w[5L] < -1e-8 * z$wref) {
      stop("internal error: a sum-of-squares variance coefficient is negative.",
           call. = FALSE)
    }
    z$w[1L] <- max(z$w[1L], 0)
    z$w[5L] <- max(z$w[5L], 0)
    z
  }

  raw <- kernel(cs$Yt, cs$Xt)

  # Cross-fit path: the summands of w0^CF and w4^CF are products of two
  # DIFFERENT scalars, not squares, so no coefficient carries a sign
  # guarantee -- the sum-of-squares clamp and internal-error branch must NOT
  # run here (Ligtenberg 2025, Section 5.4).  cf$w is
  # finiteness-guarded at its source (.crossfit_var).
  if (!is.null(cf)) {
    if (all(is.finite(raw$n))) {
      return(list(n = raw$n, w = cf$w, k = k, G = cs$G))
    }
  } else if (all(is.finite(c(raw$n, raw$w))) &&
             !(all(raw$w == 0) && any(raw$n != 0))) {
    # Ordinary designs keep the historical floating-point path bit for bit.
    # The second condition detects quartic underflow: in exact arithmetic
    # Vhat is identically zero only when Q is too (every cross-cluster
    # scalar vanishes), so an all-zero w next to a non-zero n means the
    # variance coefficients (quartic in the data scale) underflowed while
    # the quadratic numerator survived -- previously a silent
    # accept-everything/statistic contradiction.
    raw <- guard(raw)
    return(list(n = raw$n, w = raw$w, k = k, G = cs$G))
  }

  # Overflow/underflow fallback: redo the contractions on unit-normalized
  # sums (all intermediates O(1)), then map back to the original beta scale
  # by the exact monomial powers of the two norms.  Where the mapped-back
  # coefficients are themselves not representable in double precision the
  # test family is still equivariant to a change of units, but the documented
  # self-contained coefficient contract is not -- fail closed with the
  # remedy rather than crash on NaN or silently return a degenerate zero
  # polynomial.
  a <- .stable_norm(as.numeric(cs$Yt))
  d <- .stable_norm(as.numeric(cs$Xt))
  if (!(a > 0) || !(d > 0) || !is.finite(a) || !is.finite(d)) {
    stop("internal error: degenerate per-cluster instrument sums in the ",
         "CJAR coefficient kernel.", call. = FALSE)
  }
  nz <- kernel(cs$Yt / a, cs$Xt / d)
  if (is.null(cf)) nz <- guard(nz)
  n <- nz$n * c(a * a, a * d, d * d)
  w <- nz$w * c(a^4, a^3 * d, a^2 * d^2, a * d^3, d^4)
  bad <- any(!is.finite(n)) || any(n == 0 & nz$n != 0) ||
    (is.null(cf) && (any(!is.finite(w)) || any(w == 0 & nz$w != 0)))
  if (bad) {
    stop("the CJAR/CJS coefficient polynomials are not representable in ",
         "double precision at this data scale (the variance coefficients ",
         "are quartic in the scale of the residualized y and x). The ",
         "statistics are scale-invariant under common rescaling. Rescale ",
         "`y` and `x` by a common factor to preserve beta, or if scaling ",
         "them separately transform `beta0` and the reported endpoints by ",
         "the outcome-scale/regressor-scale ratio, then refit.", call. = FALSE)
  }
  if (!is.null(cf)) {
    return(list(n = n, w = cf$w, k = k, G = cs$G))
  }
  list(n = n, w = w, k = k, G = cs$G)
}

# One-sided CJAR critical value: shifted-and-scaled chi-square (the paper's
# p. 8 recommendation, matching chi-square for small k and z_{1-a} as k grows)
# or the plain normal quantile.
.cjar_crit <- function(level, k, calibration) {
  switch(calibration,
    chisq  = (stats::qchisq(level, df = k) - k) / sqrt(2 * k),
    normal = stats::qnorm(level),
    stop("unknown `calibration`: ", calibration, call. = FALSE)
  )
}

# Stable point evaluation for the CJAR statistic and p-value.  The ordinary
# finite path deliberately retains the historical arithmetic.  When a large
# finite beta makes Q(beta), Vhat(beta), or k*Vhat(beta) overflow (or a small
# value underflow), the same scale-free ratio is evaluated from signed log
# magnitudes.  This is also the single implementation used by plot.cjar, so a
# fitted point null and the plotted p-value curve cannot diverge numerically.
.cjar_point_eval <- function(nc, wc, k, beta, calibration,
                             nonpos_var = c("zero", "accept")) {
  nonpos_var <- match.arg(nonpos_var)
  qco <- c(nc[1L], -nc[2L], nc[3L])

  pfun <- function(statistic) switch(calibration,
    chisq = stats::pchisq(k + sqrt(2 * k) * statistic, df = k,
                          lower.tail = FALSE),
    normal = stats::pnorm(statistic, lower.tail = FALSE),
    stop("unknown `calibration`: ", calibration, call. = FALSE))

  one <- function(b) {
    # Keep these two expressions in exactly the same operation order used by
    # cjar() before the stable fallback was added.  Besides avoiding an
    # unnecessary rounding change for ordinary inputs, this preserves the
    # ordinary finite-input arithmetic.  Horner/signed-log evaluation is
    # needed only after this direct expression ceases to be representable.
    qd <- nc[1L] - nc[2L] * b + nc[3L] * b^2
    vd <- sum(wc * b^(seq_along(wc) - 1L))
    qs <- .poly_signed_log(qco, b)
    vs <- .poly_signed_log(wc, b)

    # A finite non-zero direct value has the most faithful sign.  At zero or
    # non-finite values the signed-log path distinguishes genuine zeros from
    # underflow and retains the polynomial's leading-tail sign.
    vsgn <- if (is.finite(vd) && vd != 0) sign(vd) else unname(vs[1L])
    qval <- if (is.finite(qd)) qd else .signed_log_value(qs)
    vval <- if (is.finite(vd)) vd else .signed_log_value(vs)

    if (vsgn <= 0) {
      if (nonpos_var == "accept") {
        return(c(statistic = NA_real_, p.value = 1,
                 numerator = qval, variance = vval))
      }
      # Plain Vhat is non-negative in exact arithmetic.  Vhat = 0 forces
      # Q = 0, and the package convention defines the resulting 0/0 as T = 0.
      return(c(statistic = 0, p.value = pfun(0),
               numerator = qval, variance = max(vval, 0)))
    }

    kv <- k * vd
    statistic <- if (is.finite(qd) &&
                      (qd != 0 || unname(qs[1L]) == 0) &&
                      is.finite(vd) && vd > 0 &&
                      is.finite(kv) && kv > 0) {
      qd / sqrt(kv)
    } else if (unname(qs[1L]) == 0) {
      0
    } else {
      .signed_log_value(c(sign = unname(qs[1L]),
                          logabs = unname(qs[2L]) -
                            0.5 * (log(k) + unname(vs[2L]))))
    }
    c(statistic = statistic, p.value = pfun(statistic),
      numerator = qval, variance = vval)
  }

  ans <- vapply(as.numeric(beta), one, numeric(4))
  if (is.null(dim(ans))) ans <- matrix(ans, nrow = 4L)
  list(statistic = as.numeric(ans["statistic", ]),
       p.value = as.numeric(ans["p.value", ]),
       numerator = as.numeric(ans["numerator", ]),
       variance = as.numeric(ans["variance", ]))
}

# Stable point evaluation for the CJS LM statistic.  The paper's sqrt(n) and
# n scalings cancel in LM; they enter only the separately reported score and
# variance.  As above, the historical direct path is kept whenever every
# intermediate is representable, with a signed-log ratio as the fallback.
.cjs_point_eval <- function(sc, vc, beta, n = 1) {
  sco <- c(sc[1L], -sc[2L])

  one <- function(b) {
    # Match cjscore()'s historical operation order on its ordinary finite
    # path; the signed-log path below is solely an overflow/underflow fallback.
    sd <- sc[1L] - sc[2L] * b
    vd <- vc[1L] + vc[2L] * b + vc[3L] * b^2
    ss <- .poly_signed_log(sco, b)
    vs <- .poly_signed_log(vc, b)
    vsgn <- if (is.finite(vd) && vd != 0) sign(vd) else unname(vs[1L])

    score <- sd / sqrt(n)
    if (!is.finite(score)) score <- .signed_log_value(ss, -0.5 * log(n))
    variance <- vd / n
    if (!is.finite(variance)) variance <- .signed_log_value(vs, -log(n))

    if (vsgn <= 0) {
      return(c(statistic = NA_real_, p.value = 1,
               score = score, variance = variance))
    }

    score_sq <- score^2
    direct_ok <- is.finite(sd) && (sd != 0 || unname(ss[1L]) == 0) &&
      is.finite(vd) && vd > 0 && is.finite(score) &&
      is.finite(variance) && variance > 0 && is.finite(score_sq)
    statistic <- if (direct_ok) score_sq / variance else NA_real_
    if (!is.finite(statistic)) {
      statistic <- if (unname(ss[1L]) == 0) 0 else
        .signed_log_value(c(sign = 1,
                            logabs = 2 * unname(ss[2L]) -
                              unname(vs[2L])))
    }
    c(statistic = statistic,
      p.value = stats::pchisq(statistic, df = 1, lower.tail = FALSE),
      score = score, variance = variance)
  }

  ans <- vapply(as.numeric(beta), one, numeric(4))
  if (is.null(dim(ans))) ans <- matrix(ans, nrow = 4L)
  list(statistic = as.numeric(ans["statistic", ]),
       p.value = as.numeric(ans["p.value", ]),
       score = as.numeric(ans["score", ]),
       variance = as.numeric(ans["variance", ]))
}

.cjar_invert <- function(nc, wc, k, c_a, nonpos_var = c("zero", "accept")) {
  nonpos_var <- match.arg(nonpos_var)
  F_CJ <- if (wc[5] > 0) nc[3] / sqrt(k * wc[5]) else NA_real_

  qco <- c(nc[1L], -nc[2L], nc[3L])
  mult <- c_a^2 * k
  centre <- .poly_remote_center(qco, wc, mult)
  qw <- if (centre == 0) qco else .poly_translate(qco, centre)
  vw <- if (centre == 0) wc else .poly_translate(wc, centre)
  env <- .poly_square_envelope(qw, vw, mult)
  zg0 <- .poly_log_balance(env)
  gt0 <- .poly_square_minus_scaled(qw, vw, mult, zg0)
  gb <- .poly_balanced_form(gt0)
  gt <- gb$co
  zg <- zg0 + gb$z
  rg <- vapply(.poly_real_roots_low(gt), .scale_root2, numeric(1), z = zg)
  qf <- .poly_balanced_form(qw)
  vf <- .poly_balanced_form(vw)

  # Direct unsquared comparison by signs plus the conditioned squared
  # polynomial.  Once Q and c_a*sqrt(kV) have the same non-zero sign, their
  # ordering is exactly the ordering of their squares (reversed when both are
  # negative).  This retains differences smaller than the rounding error of
  # separately taking logs or square roots at a very remote boundary.
  cmp <- function(b) {
    vs <- .poly_balanced_sign(vf$co, vf$z, b)
    if (vs <= 0) {
      if (nonpos_var == "accept") {
        return(list(accept = TRUE, variance_nonpos = TRUE, residual = 0))
      }
      # Plain path: Vhat >= 0 pointwise in exact arithmetic, so a
      # numerically non-positive value is the degenerate tangency case
      # (Vhat = 0 forces Q = 0 there).  Mirror the point-statistic
      # convention exactly (.cjar_build: V0 <= 0 => T := 0): accepted iff
      # 0 <= c_a, which reproduces the statistic's decision under both
      # calibrations.  Without this branch, roundoff-negative micro-windows
      # between the split halves of a variance tangency root classified by
      # sign(Q) and could reject measure-zero slivers inside an accepted
      # region (visible at G = 2, where the boundary lives on {Vhat = 0}).
      return(list(accept = c_a >= 0, variance_nonpos = TRUE, residual = 0))
    }
    qs <- .poly_balanced_sign(qf$co, qf$z, b)
    rs <- if (vs > 0 && c_a != 0) sign(c_a) else 0
    gr <- .poly_balanced_relative(gt, zg, b)
    gs <- .poly_balanced_sign(gt, zg, b)
    acc <- if (qs < rs) TRUE
           else if (qs > rs) FALSE
           else if (qs == 0) TRUE
           else if (qs > 0) gs <= 0
           else gs >= 0
    residual <- if (c_a == 0)
      .poly_balanced_relative(qf$co, qf$z, b) else gr
    list(accept = acc, variance_nonpos = vs <= 0, residual = residual,
         same_sign = qs == rs)
  }
  accept <- function(b) cmp(b)$accept

  # Filter the squared equation against the unsquared boundary.  Variance
  # roots are pooled as breakpoint candidates on BOTH paths.  On the plain
  # path Vhat >= 0 pointwise, so a real root of Vhat is a tangency at which
  # every cross-cluster scalar -- hence also Q -- vanishes; such a point is
  # accepted (the T := 0 convention) and is a legitimate breakpoint.  The
  # acceptance boundary can lie exactly on {Vhat = 0} (with positive
  # probability only at G = 2, where Vhat is a single square and
  # T(b) = sign(c_12(b))): the squared-equation filter above rightly drops
  # those candidates as variance-degenerate, so without pooling them the
  # interior acceptance interval was lost entirely (an empty set returned
  # against a non-empty pointwise region for every level with c_a < 1).
  # Pooling is harmless elsewhere: spurious near-tangency roots only
  # subdivide an interval that the midpoint classification then re-merges.
  if (length(rg)) {
    keep <- vapply(rg, function(b) {
      z <- cmp(b)
      !z$variance_nonpos && z$same_sign && abs(z$residual) <= 1e-7
    }, logical(1))
    rg <- rg[keep]
  }
  rv <- vapply(.poly_real_roots_low(vf$co), .scale_root2, numeric(1),
               z = vf$z)
  r <- sort(c(rg, rv))
  if (length(r) > 1L) {
    tol <- 128 * .Machine$double.eps *
      pmax(abs(r[-length(r)]), abs(r[-1L]))
    r <- r[c(TRUE, diff(r) > tol)]
  }

  # Because rg and rv contain every real zero of the only two expressions
  # that can change acceptance, one direct probe beyond each outermost zero
  # gives the exact signed tail classification (including negative leading
  # cross-fit variance and every degree-degenerate case).
  probes <- .outside_probes(r)
  out <- .assemble_set(r, accept, accept(probes[1L]), accept(probes[2L]),
                       root_accept = accept, shift = centre)
  list(conf_set = out$conf_set, shape = out$shape, F_CJ = F_CJ)
}

# CJS polynomial coefficients (Ligtenberg 2025, Section 4 and eq. 6).  On the partialled
# data, sqrt(n) S(b) = x' Pddot eps(b) = s0 - s1*b is linear and
# n Vhat(b) = v(b) = v0 + v1*b + v2*b^2 is quadratic, with
#   v(b) = sum_g u_g(b)^2 + sum_{g != h} c_gh(b) c_hg(b),
#   u_g(b) = p_g - b q_g,  c_gh(b) = m_gh - b d_gh,
# reduced to the same whitened per-cluster instrument sums (.cluster_sums)
# as the CJAR path: O(nk + Gk^2) time, O(Gk + k^2) memory, no G-fold R loop,
# no explicit (Z'Z)^-1. Only two of the three k x k Grams are needed (no
# A_YY, deliberately not computed on this path). v2 is a sum of squares; v0
# and the value v(b) are NOT (the cross term can be negative) -- see the
# conservative non-positive-variance convention in cjscore()'s documentation.
# Takes the .cluster_sums() object itself, like .cjar_coefs.
# With a non-NULL `cf` (a .crossfit_var() result) the cross-fit variance
# polynomial replaces the plain one; s0, s1 are untouched (Ligtenberg 2025,
# Section 5.4).
.cjs_coefs <- function(cs, cf = NULL) {
  if (!is.null(cf) && is.null(cf$v)) {
    stop("internal error: the supplied cross-fit object carries no score ",
         "variance polynomial (built with which = \"ar\").",
         call. = FALSE)
  }
  # One pass of the closed-form kernel for a given (Yt, Xt) pair, shared by
  # the raw-scale call (bit-identical to the historical path) and the
  # normalized overflow/underflow fallback below.
  kernel <- function(Yt, Xt) {
    sy <- rowSums(Yt); sx <- rowSums(Xt)
    dgg <- colSums(Xt^2); mgg <- colSums(Xt * Yt)

    s0 <- sum(sx * sy) - sum(mgg)
    s1 <- sum(sx^2) - sum(dgg)

    # column sums of the cross-cluster scalars, O(Gk)
    p <- drop(crossprod(Yt, sx)) - mgg    # p_g = sum_{h != g} m_hg
    q <- drop(crossprod(Xt, sx)) - dgg    # q_g = sum_{h != g} d_gh

    AXX <- tcrossprod(Xt)
    AXY <- tcrossprod(Xt, Yt)
    S_mM <- sum(AXY * t(AXY))     # sum over g,h of m_gh * m_hg
    S_dm <- sum(AXX * AXY)        # sum over g,h of d_gh * m_gh
    S_dd <- sum(AXX^2)            # sum over g,h of d_gh^2

    v0 <- sum(p^2) + (S_mM - sum(mgg^2))
    v1 <- -2 * (sum(p * q) + (S_dm - sum(dgg * mgg)))
    v2 <- sum(q^2) + (S_dd - sum(dgg^2))
    list(s = c(s0, s1), v = c(v0, v1, v2),
         vref = sum(q^2) + S_dd + 1e-300)
  }

  # v2 is a sum of squares; anything beyond roundoff-negative is an internal
  # error, tiny negatives are clamped to zero. (v0 and v1 carry no sign
  # guarantee -- deliberately not clamped.)  Runs only on finite values, so
  # the comparison can no longer see NaN.
  guard <- function(z) {
    if (z$v[3L] < -1e-8 * z$vref) {
      stop("internal error: a sum-of-squares variance coefficient is negative.",
           call. = FALSE)
    }
    z$v[3L] <- max(z$v[3L], 0)
    z
  }

  raw <- kernel(cs$Yt, cs$Xt)

  # Cross-fit path: v2^CF carries no sum-of-squares sign guarantee, so the
  # clamp and internal-error branch must NOT run here
  # (Ligtenberg 2025, Section 5.4). The inverter and the conservative
  # non-positive-variance convention already work with signed
  # coefficients; cf$v is finiteness-guarded at its source (.crossfit_var).
  if (!is.null(cf)) {
    if (all(is.finite(raw$s))) {
      return(list(s = raw$s, v = cf$v, k = cs$k, G = cs$G))
    }
  } else if (all(is.finite(c(raw$s, raw$v))) &&
             !(all(raw$v == 0) && any(raw$s != 0))) {
    # Ordinary designs keep the historical floating-point path bit for bit;
    # the second condition detects quartic-in-scale underflow of the
    # variance coefficients next to a surviving score (see .cjar_coefs).
    raw <- guard(raw)
    return(list(s = raw$s, v = raw$v, k = cs$k, G = cs$G))
  }

  # Overflow/underflow fallback: unit-normalized recompute plus the exact
  # monomial map back; fail closed where the mapped coefficients are not
  # representable (see .cjar_coefs for the rationale and remedy message).
  a <- .stable_norm(as.numeric(cs$Yt))
  d <- .stable_norm(as.numeric(cs$Xt))
  if (!(a > 0) || !(d > 0) || !is.finite(a) || !is.finite(d)) {
    stop("internal error: degenerate per-cluster instrument sums in the ",
         "CJS coefficient kernel.", call. = FALSE)
  }
  nz <- kernel(cs$Yt / a, cs$Xt / d)
  if (is.null(cf)) nz <- guard(nz)
  s <- nz$s * c(a * d, d * d)
  v <- nz$v * c(a^2 * d^2, a * d^3, d^4)
  bad <- any(!is.finite(s)) || any(s == 0 & nz$s != 0) ||
    (is.null(cf) && (any(!is.finite(v)) || any(v == 0 & nz$v != 0)))
  if (bad) {
    stop("the CJAR/CJS coefficient polynomials are not representable in ",
         "double precision at this data scale (the variance coefficients ",
         "are quartic in the scale of the residualized y and x). The ",
         "statistics are scale-invariant under common rescaling. Rescale ",
         "`y` and `x` by a common factor to preserve beta, or if scaling ",
         "them separately transform `beta0` and the reported endpoints by ",
         "the outcome-scale/regressor-scale ratio, then refit.", call. = FALSE)
  }
  if (!is.null(cf)) {
    return(list(s = s, v = cf$v, k = cs$k, G = cs$G))
  }
  list(s = s, v = v, k = cs$k, G = cs$G)
}

# Analytic inversion of the CJS acceptance region under the conservative
# non-positive-variance convention:
#   reject b  iff  v(b) > 0 and A(b) > 0,
# with A(b) = (s0 - s1*b)^2 - crit * v(b); every identity used here is
# spot-checked against a dense computation.
# Where v > 0 this is the usual LM test {A <= 0}; where v <= 0 the statistic
# is undefined and b is accepted -- we decline to reject. The acceptance set
# is the complement of the intersection of two quadratic positivity regions,
# so it is still exactly invertible, no grid: the breakpoints are the pooled
# real roots of A and of v. Unlike the CJAR quartic no squaring step is
# involved, so every real root is a genuine breakpoint -- no spurious-root
# filter is needed; the same trimming, polish and midpoint/tail sign
# discipline as .cjar_invert applies, per polynomial. Every breakpoint is
# accepted by construction (at a root neither polynomial is > 0), so the
# acceptance set is closed; a breakpoint not adjacent to an accepted interval
# joins the set as a degenerate row. Because points with v(b) <= 0 are
# accepted, the set can pick up a closed interval where the variance estimate
# is non-positive, possibly disjoint from the main region (see ?cjscore).
.cjs_invert <- function(sc, vc, crit) {
  F_CJS <- if (vc[3] > 0) sc[2] / sqrt(vc[3]) else NA_real_

  sco <- c(sc[1L], -sc[2L])
  centre <- .poly_remote_center(sco, vc, crit)
  sw <- if (centre == 0) sco else .poly_translate(sco, centre)
  vw <- if (centre == 0) vc else .poly_translate(vc, centre)
  env <- .poly_square_envelope(sw, vw, crit)
  za0 <- .poly_log_balance(env)
  At0 <- .poly_square_minus_scaled(sw, vw, crit, za0)
  ab <- .poly_balanced_form(At0)
  At <- ab$co
  za <- za0 + ab$z
  ra <- vapply(.poly_real_roots_low(At), .scale_root2, numeric(1), z = za)
  vf <- .poly_balanced_form(vw)
  rv <- vapply(.poly_real_roots_low(vf$co), .scale_root2, numeric(1),
               z = vf$z)
  r <- sort(c(ra, rv))
  if (length(r) > 1L) {
    tol <- 128 * .Machine$double.eps *
      pmax(abs(r[-length(r)]), abs(r[-1L]))
    r <- r[c(TRUE, diff(r) > tol)]
  }

  accept <- function(b) {
    .poly_balanced_sign(vf$co, vf$z, b) <= 0 ||
      .poly_balanced_sign(At, za, b) <= 0
  }
  probes <- .outside_probes(r)
  out <- .assemble_set(r, accept, accept(probes[1L]), accept(probes[2L]),
                       root_accept = accept, shift = centre)
  list(conf_set = out$conf_set, shape = out$shape, F_CJS = F_CJS)
}

# Cross-fit variance polynomial coefficients for the CJAR and CJS tests
# (Ligtenberg 2025, Section 5.4, eq. 7), in whitened coordinates. The
# leave-out fitted values are affine in b, so Vhat_AR_CF stays quartic and
# n * Vhat_S_CF quadratic. Each leave-two-clusters-out projection is a
# downdate of the whitened Gram, I_k - S_g - S_h with S_g = Ztil_[g]' Ztil_[g],
# and every per-pair scalar reduces to k-vector contractions of the per-cluster
# whitened sums. The per-pair solve is dispatched on m = n_g + n_h like
# .cluster_block: for m < k the Woodbury identity on the m x m capacitance,
# for m >= k the k x k Cholesky downdate. S_g Grams are cached lazily, only for
# clusters entering an m >= k solve. Callers that ran .first_stage pass Ztil to
# avoid re-whitening; direct calls may omit it.
#
# `which` restricts work to the consumed polynomial: "ar" skips the score
# single-cluster loop, "score" skips the AR accumulation, "both" shares one
# pair sweep. Unrequested slots are NULL; .cjar_coefs/.cjs_coefs fail closed on
# misuse. `w` are the coefficients of Vhat_AR_CF(b) (2/k included) and `v` of
# n * Vhat_S_CF(b); no sign clamp, since cross-fit summands are products of two
# different scalars with no sum-of-squares guarantee.
#
# Guard: each leave-out Gram must be numerically usable, not merely PD in exact
# arithmetic. The eigenvalues of I_k - S_g - S_h are basis-invariant (the full
# whitened Gram is I_k); require reciprocal condition number > sqrt(eps),
# retaining ~half the double-precision digits per solve. Cheap certificates
# skip the eigendecomposition on safe sets (the O(1) Weyl bound from the
# per-cluster spectral norms, then Gershgorin on the k x k side); uncertified
# sets get the exact symmetric spectral check (on the m x m capacitance when
# m < k, whose spectrum is the downdate's plus k - m ones). Certificates only
# skip work -- the accept/stop classification stays spectral and basis-
# invariant. Fail closed, naming the cluster(s).
.crossfit_var <- function(y, x, Z, cluster, R = NULL, Ztil = NULL,
                          which = c("both", "ar", "score")) {
  which <- match.arg(which)
  do_ar <- which != "score"
  do_score <- which != "ar"
  Z <- as.matrix(Z)
  k <- ncol(Z)
  n <- length(y)
  cl <- droplevels(as.factor(cluster))
  G <- nlevels(cl)
  if (G < 3L) {
    stop("the cross-fit variance estimator requires at least 3 clusters: ",
         "its leave-two-clusters-out projections must retain data.",
         call. = FALSE)
  }
  if (is.null(Ztil)) {
    if (is.null(R)) R <- .gram_chol(Z)
    Ztil <- .whiten(Z, R)
  } else {
    Ztil <- as.matrix(Ztil)
    if (!is.numeric(Ztil) || nrow(Ztil) != n || ncol(Ztil) != k ||
        anyNA(Ztil) || any(!is.finite(Ztil))) {
      stop("internal error: supplied whitened instruments have incompatible dimensions or non-finite values.",
           call. = FALSE)
    }
  }
  groups <- split(seq_len(n), cl)
  labs <- levels(cl)

  ty <- drop(crossprod(Ztil, y))
  tx <- drop(crossprod(Ztil, x))
  Ay <- t(rowsum(Ztil * y, cl))          # column g = Ztil_[g]' y_[g], k x G
  Ax <- t(rowsum(Ztil * x, cl))
  # Per-cluster whitened row blocks, sizes, and spectral norms
  # lambda_max(S_g) = ||P_Z[g,g]||_2 via the shared size-dispatched kernel
  # (.cluster_block with u = NULL).  The norms feed the Weyl certificate in
  # lo_solve(); they are recomputed here rather than threaded from the
  # callers' leverage pass, which retains only the maximum.
  Zrows <- lapply(groups, function(i) Ztil[i, , drop = FALSE])
  ng <- lengths(groups)
  lev <- vapply(Zrows, function(Zg) .cluster_block(Zg, k)$lev, numeric(1))

  # The pair products below are quartic in the scale of the whitened cluster
  # sums.  Fail closed with the scale remedy before the loop rather than
  # propagate Inf/NaN or silently underflow the whole variance polynomial to
  # zero (which the conservative non-positive-variance convention would then
  # silently convert into accept-everything).
  ny <- .stable_norm(as.numeric(Ay))
  nx <- .stable_norm(as.numeric(Ax))
  quartic_bad <- function(s) {
    s2 <- s * s
    !is.finite(s2 * s2) || (s > 0 && s2 * s2 == 0)
  }
  if (quartic_bad(ny) || quartic_bad(nx)) {
    stop("the cross-fit variance polynomial is not representable in double ",
         "precision at this data scale (its coefficients are quartic in the ",
         "scale of the residualized y and x). The statistics are ",
         "scale-invariant under common rescaling. Rescale `y` and `x` by a ",
         "common factor to preserve beta, or if scaling them separately ",
         "transform `beta0` and the reported endpoints by the outcome-scale/",
         "regressor-scale ratio, then refit.", call. = FALSE)
  }

  # One k x k Cholesky per leave-out set; triangular backsolves, never an
  # explicit inverse.  sqrt(eps) is deliberately much larger than the LAPACK
  # positivity threshold: at rcond_2 <= sqrt(eps), a formally invertible
  # downdate can retain fewer than about eight decimal digits, and the AR/LM
  # variance coefficients multiply quantities obtained from these solves.
  rcond_tol <- sqrt(.Machine$double.eps)
  spectrally_usable <- function(M) {
    # Gershgorin lower/upper bounds are only a safe fast certificate.  If the
    # bounds cannot prove rcond_2 > rcond_tol, use eigenvalues; therefore a
    # change of instrument basis can affect cost, but not classification.
    am <- abs(M)
    rad <- rowSums(am) - diag(am)
    lower <- min(diag(M) - rad)
    upper <- max(diag(M) + rad)
    if (is.finite(lower) && is.finite(upper) && upper > 0 &&
        lower > rcond_tol * upper) {
      return(TRUE)
    }
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    lmax <- max(ev)
    lmin <- min(ev)
    is.finite(lmax) && is.finite(lmin) && lmax > 0 &&
      lmin > rcond_tol * lmax
  }

  usable_stop <- function(labels) {
    stop("cross-fit variance is undefined: the leave-out instrument Gram ",
         "is singular or numerically ill-conditioned when cluster(s) ",
         paste(labels, collapse = ", "), " are left out (spectral ",
         "reciprocal condition number <= sqrt(.Machine$double.eps)). ",
         "This is the leave-out analogue of the within-cluster leverage ",
         "advisory.", call. = FALSE)
  }

  # Exact spectral check on the m x m capacitance C = I_m - U U' of a
  # stacked row block U with m < k: the k x k downdate I_k - U'U has
  # eigenvalues {eigen(C)} plus k - m exact ones, so its reciprocal
  # condition number is decided on the small side without forming it.
  cap_usable <- function(cap) {
    ev <- eigen(cap, symmetric = TRUE, only.values = TRUE)$values
    lmax <- max(ev, 1)
    lmin <- min(ev, 1)
    is.finite(lmax) && is.finite(lmin) && lmax > 0 && lmin > rcond_tol * lmax
  }

  # k x k Grams are needed only by m >= k solves; cache them lazily so the
  # small-cluster regime never allocates the former G-element k x k list.
  Slist <- vector("list", G)
  getS <- function(g) {
    Sgg <- Slist[[g]]
    if (is.null(Sgg)) {
      Sgg <- crossprod(Zrows[[g]])
      Slist[[g]] <<- Sgg
    }
    Sgg
  }
  Ik <- NULL

  # Leave-out solve (I_k - sum_{g in gs} S_g)^{-1} rhs for one or two
  # clusters, dispatched on the stacked row count (see the header comment).
  # The Weyl certificate is checked first; an uncertified set falls through
  # to the exact spectral check on whichever side the branch already holds.
  lo_solve <- function(gs, rhs) {
    m <- sum(ng[gs])
    certified <- 1 - sum(lev[gs]) > rcond_tol
    if (m < k) {
      U <- if (length(gs) == 1L) Zrows[[gs]] else
        rbind(Zrows[[gs[1L]]], Zrows[[gs[2L]]])
      cap <- diag(m) - tcrossprod(U)
      if (!certified && !cap_usable(cap)) usable_stop(labs[gs])
      # The spectral margin makes a Cholesky failure unexpected; do not mask
      # a programming/LAPACK error behind the documented availability
      # condition.
      ch <- chol(cap)
      rhs + crossprod(U, backsolve(ch, backsolve(ch, U %*% rhs,
                                                 transpose = TRUE)))
    } else {
      if (is.null(Ik)) Ik <<- diag(k)
      Smat <- getS(gs[1L])
      if (length(gs) > 1L) Smat <- Smat + getS(gs[2L])
      M <- Ik - Smat
      if (!certified && !spectrally_usable(M)) usable_stop(labs[gs])
      ch <- chol(M)
      backsolve(ch, backsolve(ch, rhs, transpose = TRUE))
    }
  }

  w <- if (do_ar) numeric(5) else NULL   # Vhat_AR_CF quartic (2/k at the end)
  t2 <- if (do_score) numeric(3) else NULL   # second Vhat_S_CF term
  for (g in seq_len(G - 1L)) {
    ayg <- Ay[, g]; axg <- Ax[, g]
    for (h in (g + 1L):G) {
      ayh <- Ay[, h]; axh <- Ax[, h]
      small <- ng[g] + ng[h] < k
      sol <- lo_solve(c(g, h), cbind(ty - ayg - ayh, tx - axg - axh))
      # S_g sol without the k x k Gram on the small side; the m >= k side
      # keeps the cached-Gram product of the historical path.
      if (small) {
        SgS <- crossprod(Zrows[[g]], Zrows[[g]] %*% sol)
        ShS <- crossprod(Zrows[[h]], Zrows[[h]] %*% sol)
      } else {
        SgS <- getS(g) %*% sol
        ShS <- getS(h) %*% sol
      }
      # Ztil_[g]' u_g^{gh} etc.: whitened sums of the pair-specific
      # cross-fit residual vectors u = y - ytil(g,h), v = x - xtil(g,h).
      zug <- ayg - SgS[, 1L]; zvg <- axg - SgS[, 2L]
      zuh <- ayh - ShS[, 1L]; zvh <- axh - ShS[, 2L]

      if (do_ar) {
        # AR: a_gh(b) = A0 + A1 b + A2 b^2 and its (h, g) mirror; the ordered
        # sum counts each unordered pair twice with an identical product.
        A0 <- sum(zug * ayh)
        A1 <- -(sum(zug * axh) + sum(zvg * ayh))
        A2 <- sum(zvg * axh)
        B0 <- sum(zuh * ayg)
        B1 <- -(sum(zuh * axg) + sum(zvh * ayg))
        B2 <- sum(zvh * axg)
        w <- w + 2 * c(A0 * B0,
                       A0 * B1 + A1 * B0,
                       A0 * B2 + A1 * B1 + A2 * B0,
                       A1 * B2 + A2 * B1,
                       A2 * B2)
      }

      if (do_score) {
        # Score second term: affine [x'_g P_gh eps_h] times affine
        # [(u - v b)'_g P_gh x_h], both orders (the term is not symmetric in
        # g, h; the x'_g P_gh x_h factor F1 is).
        F0 <- sum(axg * ayh); F1 <- -sum(axg * axh)
        G0 <- sum(zug * axh); G1 <- -sum(zvg * axh)
        F0r <- sum(axh * ayg)
        G0r <- sum(zuh * axg); G1r <- -sum(zvh * axg)
        t2 <- t2 + c(F0 * G0 + F0r * G0r,
                     F0 * G1 + F1 * G0 + F0r * G1r + F1 * G0r,
                     F1 * G1 + F1 * G1r)
      }
    }
  }

  # Score first term: single-cluster leave-out (N-3 reading), factorised
  # over h via q_h = sum_{g != h} P_Z[h, g] x_[g] = Ztil_[h] (t_x - a_x_h);
  # every inner product reduces to the k-vector qx.  Skipped entirely when
  # the caller consumes only the AR polynomial.
  t1 <- if (do_score) numeric(3) else NULL
  if (do_score) {
    for (h in seq_len(G)) {
      ayh <- Ay[, h]; axh <- Ax[, h]
      qx <- tx - axh
      sol1 <- lo_solve(h, cbind(ty - ayh, tx - axh))
      ShS1 <- if (ng[h] < k) {
        crossprod(Zrows[[h]], Zrows[[h]] %*% sol1)
      } else {
        getS(h) %*% sol1
      }
      zu1 <- ayh - ShS1[, 1L]
      zv1 <- axh - ShS1[, 2L]
      L0 <- sum(zu1 * qx); L1 <- -sum(zv1 * qx)
      R0 <- sum(ayh * qx); R1 <- -sum(axh * qx)
      t1 <- t1 + c(L0 * R0, L0 * R1 + L1 * R0, L1 * R1)
    }
  }

  out <- list(w = if (do_ar) (2 / k) * w else NULL,
              v = if (do_score) t1 + t2 else NULL, k = k, G = G)
  if ((do_ar && any(!is.finite(out$w))) ||
      (do_score && any(!is.finite(out$v)))) {
    stop("internal error: non-finite cross-fit variance coefficients ",
         "survived the scale guard.", call. = FALSE)
  }
  out
}
