# Internal helpers: polynomial and root machinery for the analytic
# confidence-set inversions -- balanced/signed-log evaluation, real
# roots of low-degree polynomials, and acceptance-set assembly.  No
# test semantics live here.

# Log2 magnitude, with -Inf at exact zero.  The low-degree inverters do all
# beta and coefficient balancing in this representation, so multiplying a
# coefficient by s^j never creates an intermediate overflow or underflow.
.log2abs <- function(x) {
  out <- rep.int(-Inf, length(x))
  nz <- x != 0
  out[nz] <- log2(abs(x[nz]))
  out
}

# Choose beta = 2^z t by equilibrating the first and last non-zero empirical
# degrees of a magnitude envelope.  Unlike the old constant/leading-only
# rules this works when either endpoint coefficient is zero and it is blind
# to coefficient signs.  Rounding z to an integer loses at most half a binary
# exponent of equilibration, while making every coordinate scaling exact; an
# integer-power units change x -> 2^a x therefore subtracts a from z exactly.
.poly_log_balance <- function(logmag) {
  ii <- which(is.finite(logmag))
  if (length(ii) < 2L) return(0)
  i <- ii[1L] - 1L
  j <- ii[length(ii)] - 1L
  round((logmag[ii[1L]] - logmag[ii[length(ii)]]) / (j - i))
}

# Multiply by an integral power of two in bounded chunks.  Unlike a
# log2()/exp2() round trip this preserves every representable coefficient bit;
# chunking also handles exponents whose power of two cannot itself be stored.
.mul_pow2 <- function(x, z) {
  if (!is.finite(z) || z != round(z))
    stop("internal error: a non-integral binary scale was requested.",
         call. = FALSE)
  out <- x
  while (z > 512) {
    out <- out * 2^512
    z <- z - 512
  }
  while (z < -512) {
    out <- out * 2^-512
    z <- z + 512
  }
  out * 2^z
}

# Map a root in the balanced coordinate back by an exact binary scaling.
# Overflow/underflow occurs only when the mapped root is outside the finite
# double range.
.scale_root2 <- function(x, z) .mul_pow2(x, z)

# Horner evaluation and its absolute roundoff scale.
.poly_horner <- function(co, x) {
  v <- co[length(co)]
  if (length(co) > 1L) {
    for (j in (length(co) - 1L):1L) v <- v * x + co[j]
  }
  v
}

.poly_abs_horner <- function(co, x) .poly_horner(abs(co), abs(x))

# Coefficients after the exact change beta = centre + delta, lowest degree
# first.  `sum()` uses R's extended accumulator where available.  This small
# degree (at most four) composition is used only when a large translation,
# rather than a units change, is what makes the input coordinate ill-scaled.
.poly_translate <- function(co, centre) {
  d <- length(co) - 1L
  vapply(0:d, function(l) {
    j <- l:d
    sum(co[j + 1L] * choose(j, l) * centre^(j - l))
  }, numeric(1))
}

# Detect a remote polynomial location independently of a pure units change.
# The primary polynomial supplies its root centroid (or its sole linear root),
# while the jointly translated primary/secondary envelope supplies the local
# width.  A centre is used only when it is at least 64 local widths from zero;
# ordinary paths therefore keep their established floating-point route.
.poly_remote_center <- function(primary, secondary, mult) {
  ii <- which(primary != 0)
  if (length(ii) < 2L) return(0)
  d <- ii[length(ii)] - 1L
  centre <- if (d >= 2L && primary[d] != 0) {
    -primary[d] / (d * primary[d + 1L])
  } else if (primary[2L] != 0) {
    -primary[1L] / primary[2L]
  } else 0
  if (!is.finite(centre) || abs(centre) <= 64) return(0)
  pc <- .poly_translate(primary, centre)
  sc <- .poly_translate(secondary, centre)
  if (any(!is.finite(c(pc, sc)))) return(0)
  z <- .poly_log_balance(.poly_square_envelope(pc, sc, mult))
  width <- if (z > log2(.Machine$double.xmax)) Inf else 2^z
  if (is.finite(width) && abs(centre) > 64 * max(1, width)) centre else 0
}

# Complete real-root isolation for the degree-at-most-four polynomials used by
# CJAR/CJS.  Derivative roots partition the Cauchy-bounded real line into
# monotone intervals; sign-changing roots are bracketed by uniroot and a
# derivative root whose residual is at roundoff is retained as a repeated
# root.  This is not a beta grid.
.poly_real_roots_low <- function(co) {
  while (length(co) > 1L && co[length(co)] == 0) co <- co[-length(co)]
  d <- length(co) - 1L
  if (d < 1L) return(numeric(0))
  if (co[1L] == 0) {
    return(sort(c(0, .poly_real_roots_low(co[-1L]))))
  }
  if (d == 1L) return(-co[1L] / co[2L])

  dr <- .poly_real_roots_low((seq_len(d)) * co[-1L])
  lead <- abs(co[d + 1L])
  R <- 1 + max(abs(co[seq_len(d)]) / lead)
  if (!is.finite(R)) R <- .Machine$double.xmax / 4
  dr <- sort(dr[is.finite(dr) & dr > -R & dr < R])
  pts <- c(-R, dr, R)
  vals <- vapply(pts, .poly_horner, numeric(1), co = co)
  ztol <- function(x) 128 * (d + 1) * .Machine$double.eps *
    max(.poly_abs_horner(co, x), .Machine$double.xmin)

  out <- numeric(0)
  if (length(dr)) {
    keep <- vapply(seq_along(dr), function(i)
      abs(vals[i + 1L]) <= ztol(dr[i]), logical(1))
    out <- dr[keep]
  }
  for (i in seq_len(length(pts) - 1L)) {
    a <- pts[i]; b <- pts[i + 1L]
    fa <- vals[i]; fb <- vals[i + 1L]
    if (!is.finite(fa) || !is.finite(fb) || fa == 0 || fb == 0 ||
        sign(fa) == sign(fb)) next
    # zeroin's own stopping rule is 2 * eps * |x| + tol / 2, i.e. relative to
    # the current iterate.  A tolerance tied to the bracket ends instead
    # (the former 4 * eps * min(|a|, |b|)) stopped a root near zero early when
    # a roundoff-sized leading coefficient pushed a derivative root, and so a
    # bracket end, out to ~1e12.
    rr <- stats::uniroot(function(x) .poly_horner(co, x), c(a, b),
                         tol = .Machine$double.xmin, maxiter = 5000L)$root
    out <- c(out, rr)
  }
  out <- sort(out[is.finite(out)])
  if (length(out) > 1L) {
    # The recursive construction returns a repeated root only once.  This
    # final merge is at an ulp-scale tolerance, not the former fixed 1e-8
    # relative threshold that could collapse distinct nearby boundaries.
    tol <- 128 * .Machine$double.eps *
      pmax(abs(out[-length(out)]), abs(out[-1L]))
    out <- out[c(TRUE, diff(out) > tol)]
  }
  out
}

# Signed log-magnitude of a polynomial at beta.  Terms are summed only after
# their largest log magnitude is removed; signs (including beta^j for a
# negative beta) are retained.  This supplies stable direct acceptance
# predicates at extreme finite scales without changing either statistic.
.poly_signed_log <- function(co, beta) {
  if (beta == 0) {
    if (co[1L] == 0) return(c(sign = 0, logabs = -Inf))
    return(c(sign = sign(co[1L]), logabs = log(abs(co[1L]))))
  }
  j <- seq_along(co) - 1L
  nz <- co != 0
  if (!any(nz)) return(c(sign = 0, logabs = -Inf))
  lm <- log(abs(co[nz])) + j[nz] * log(abs(beta))
  sg <- sign(co[nz]) * if (beta < 0) (-1)^j[nz] else 1
  top <- max(lm)
  sm <- sum(sg * exp(lm - top))
  if (sm == 0) return(c(sign = 0, logabs = -Inf))
  c(sign = sign(sm), logabs = top + log(abs(sm)))
}

# Convert a signed-log value back to double precision after an optional
# multiplicative scale.  Overflow and underflow of the VALUE are legitimate
# here (Inf/0); the point-test ratios below remain evaluable because they are
# formed by subtracting log magnitudes before this conversion.
.signed_log_value <- function(x, log_scale = 0) {
  s <- unname(x[1L])
  if (s == 0) return(0)
  lv <- unname(x[2L]) + log_scale
  if (lv > log(.Machine$double.xmax)) return(s * Inf)
  s * exp(lv)
}

# A midpoint that cannot overflow when two finite roots have the same sign.
.safe_midpoint <- function(a, b) a / 2 + b / 2

# Shared assembly of an inverted confidence set.  `root_accept` is evaluated
# separately from the open root-induced intervals: only an actually accepted
# boundary can become an endpoint or an isolated [r,r] component.  No
# tolerance-based merging is performed here; rejected gaps, however narrow,
# remain gaps.
.assemble_set <- function(roots, accept, acc_left, acc_right, scale = 1,
                          root_accept = accept, shift = 0) {
  r <- sort(roots[is.finite(roots)])
  if (length(r) > 1L) r <- r[c(TRUE, diff(r) != 0)]
  m <- length(r)
  lo <- c(-Inf, r)
  hi <- c(r, Inf)
  acc <- logical(m + 1L)
  for (i in seq_len(m + 1L)) {
    acc[i] <- if (i == 1L) acc_left
              else if (i == m + 1L) acc_right
              else accept(.safe_midpoint(lo[i], hi[i]))
  }

  ivs <- list()
  i <- 1L
  while (i <= m + 1L) {
    if (acc[i]) {
      j <- i
      while (j < m + 1L && acc[j + 1L]) j <- j + 1L
      ivs[[length(ivs) + 1L]] <- c(lo[i], hi[j])
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  cs <- if (length(ivs)) do.call(rbind, ivs) else matrix(numeric(0), 0L, 2L)

  if (m > 0L) {
    covered <- vapply(seq_len(m), function(i)
      nrow(cs) > 0L && any(r[i] >= cs[, 1L] & r[i] <= cs[, 2L]),
      logical(1))
    add <- !covered & vapply(r, root_accept, logical(1))
    if (any(add)) {
      cs <- rbind(cs, cbind(r[add], r[add]))
      cs <- cs[order(cs[, 1L]), , drop = FALSE]
    }
  }
  cs <- cs * scale + shift                  # back to the input beta scale
  colnames(cs) <- c("lower", "upper")

  # Shape from the endpoint pattern alone; see the callers for how each
  # degenerate route (single ray, two rays) arises.
  list(conf_set = cs, shape = .shape_from_set(cs))
}

# Analytic inversion of the CJAR acceptance region
#   {b : Q(b) <= c_a * sqrt(k * Vhat(b))},
# the sublevel set of h(b) = Q(b) - c_a * sqrt(k * Vhat(b)). Boundary points
# are real roots of the quartic g(b) = Q(b)^2 - c_a^2 k Vhat(b); direct
# evaluation of h at each candidate removes the spurious roots of squaring
# (where Q = -c_a sqrt(k Vhat)) and is robust to c_a <= 0 and n2 <= 0 alike.
# nc = c(n0, n1, n2) with Q(b) = n0 - n1*b + n2*b^2; wc = c(w0, ..., w4).
# `nonpos_var` selects the convention where Vhat(b) <= 0. On the plain path
# ("zero", the default) Vhat >= 0 pointwise in exact arithmetic, so clamping
# at 0 is only roundoff hygiene and the acceptance rule is h(b) <= 0
# throughout. On the cross-fit path ("accept") Vhat carries no sign
# guarantee, and the CJS conservative convention applies verbatim
# (Ligtenberg 2025, Section 5.4): b is rejected iff
# Vhat(b) > 0 AND Q(b) > c_a * sqrt(k * Vhat(b)); where the variance
# estimate is non-positive the statistic is undefined and b is accepted.
# The boundary then also includes the real roots of Vhat itself, and tail
# classification uses the SIGNED leading variance coefficient, never
# max(., 0).
.poly_square_envelope <- function(qco, vco, mult) {
  d <- max(2L * (length(qco) - 1L), length(vco) - 1L)
  out <- rep.int(-Inf, d + 1L)
  lq <- .log2abs(qco)
  for (i in seq_along(qco)) for (j in seq_along(qco)) {
    if (is.finite(lq[i]) && is.finite(lq[j])) {
      h <- i + j - 1L
      out[h] <- max(out[h], lq[i] + lq[j])
    }
  }
  if (mult > 0) {
    lv <- .log2abs(vco) + log2(mult)
    ii <- which(is.finite(lv))
    out[ii] <- pmax(out[ii], lv[ii])
  }
  out
}

# Coefficients of q(2^z t)^2 - mult*v(2^z t), divided by one common positive
# power of two.  The common scale is load-bearing: q^2 and mult*v are never
# normalised independently.
.poly_square_minus_scaled <- function(qco, vco, mult, z) {
  jq <- seq_along(qco) - 1L
  jv <- seq_along(vco) - 1L
  lq <- .log2abs(qco) + jq * z
  lv <- if (mult > 0) .log2abs(vco) + log2(mult) + jv * z
        else rep.int(-Inf, length(vco))
  e <- ceiling(max(c(lq, lv / 2), na.rm = TRUE))
  if (!is.finite(e)) return(0)
  qs <- vapply(seq_along(qco), function(i)
    .mul_pow2(qco[i], jq[i] * z - e), numeric(1))
  if (mult > 0) {
    em <- floor(log2(mult))
    mm <- .mul_pow2(mult, -em)
    vs <- vapply(seq_along(vco), function(i)
      .mul_pow2(vco[i], jv[i] * z + em - 2 * e), numeric(1)) * mm
  } else {
    vs <- numeric(length(vco))
  }
  out <- numeric(2L * length(qco) - 1L)
  for (i in seq_along(qs)) {
    jj <- i:(i + length(qs) - 1L)
    out[jj] <- out[jj] + qs[i] * qs
  }
  length(out) <- max(length(out), length(vs))
  out[seq_along(vs)] <- out[seq_along(vs)] - vs
  out
}

# Balance and solve one low-degree polynomial, returning roots on the input
# beta scale.  Exact zeros, rather than a fixed relative leading trim, define
# the empirical degree after scale conditioning.
.poly_balanced_form <- function(co) {
  lm <- .log2abs(co)
  z <- .poly_log_balance(lm)
  j <- seq_along(co) - 1L
  lt <- lm + j * z
  top <- ceiling(max(lt))
  if (!is.finite(top)) return(list(co = rep.int(0, length(co)), z = 0))
  ct <- vapply(seq_along(co), function(i)
    .mul_pow2(co[i], j[i] * z - top), numeric(1))
  list(co = ct, z = z)
}

.poly_balanced_relative <- function(co, z, beta) {
  ii <- which(co != 0)
  if (!length(ii)) return(0)
  t <- .scale_root2(beta, -z)
  if (is.infinite(t)) {
    d <- max(ii) - 1L
    return(sign(co[d + 1L]) * if (t < 0) (-1)^d else 1)
  }
  v <- .poly_horner(co, t)
  den <- max(.poly_abs_horner(co, t), .Machine$double.xmin)
  if (!is.finite(v) || !is.finite(den)) {
    d <- max(ii) - 1L
    return(sign(co[d + 1L]) * if (t < 0) (-1)^d else 1)
  }
  v / den
}

.poly_balanced_sign <- function(co, z, beta) {
  v <- .poly_balanced_relative(co, z, beta)
  tol <- 512 * length(co) * .Machine$double.eps
  if (abs(v) <= tol) 0 else sign(v)
}

.outside_probes <- function(roots) {
  if (!length(roots)) return(c(-1, 1))
  left <- roots[1L]
  right <- roots[length(roots)]
  dl <- max(1, abs(left))
  dr <- max(1, abs(right))
  pl <- left - dl
  pr <- right + dr
  if (!is.finite(pl)) pl <- -.Machine$double.xmax
  if (!is.finite(pr)) pr <- .Machine$double.xmax
  c(pl, pr)
}
