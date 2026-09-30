# Internal helpers: input validation, data preparation and the shared
# estimation kernel -- FWL partialling, the instrument Gram Cholesky,
# whitening, the leave-cluster-out kernel and cluster-robust IV
# inference used by every public interface.

# FWL residualisation. Weights enter as a sqrt(w) transform; lm.fit is pivoted,
# so rank-deficient controls are fine. High-dimensional fixed effects are
# absorbed first by .demean(); absorption spans the intercept within groups,
# so the explicit intercept column is dropped on that path and only the
# remaining dense controls go through lm.fit.
.partial_out <- function(y, x, Z, controls = NULL, weights = NULL,
                         intercept = TRUE, fixed_effects = NULL) {
  n <- length(y)
  Z <- as.matrix(Z)
  rw <- if (is.null(weights)) rep.int(1, n) else sqrt(weights)

  yw <- y * rw
  xw <- x * rw
  Zw <- Z * rw

  if (!is.null(fixed_effects)) {
    Cw <- if (is.null(controls)) NULL else as.matrix(controls) * rw
    Cw0 <- Cw
    M <- cbind(yw, xw, Zw, Cw)
    M <- .demean(M, fixed_effects, weights = weights)
    yw <- M[, 1L]
    xw <- M[, 2L]
    Zw <- M[, 2L + seq_len(ncol(Z)), drop = FALSE]
    if (is.null(Cw)) {
      return(list(y = yw, x = as.numeric(xw), Z = Zw, rw = rw))
    }
    Cw <- M[, 2L + ncol(Z) + seq_len(ncol(Cw)), drop = FALSE]

    # A control in the joint FE span should be a zero residual.  The certified
    # matrix-free projection leaves at most relative numerical noise, however,
    # and a one-column QR would otherwise normalize that noise and project the
    # model variables onto it.  Treat columns inside the public 1e-10 HDFE
    # contract as absorbed before the dense QR.  This is scale relative and
    # therefore unchanged by the units of a control.
    absorb_tol <- max(1e-10, .gamma_tol(max(n, ncol(Cw))))
    absorbed <- vapply(seq_len(ncol(Cw)), function(j) {
      .relative_vector_zero(Cw[, j], Cw0[, j], absorb_tol)
    }, logical(1))
    Cw <- Cw[, !absorbed, drop = FALSE]
    if (!ncol(Cw)) {
      return(list(y = yw, x = as.numeric(xw), Z = Zw, rw = rw))
    }
  } else {
    C <- if (intercept) matrix(1, n, 1L) else NULL
    if (!is.null(controls)) {
      controls <- as.matrix(controls)
      C <- if (is.null(C)) controls else cbind(C, controls)
    }
    if (is.null(C)) {
      return(list(y = yw, x = xw, Z = Zw, rw = rw))
    }
    Cw <- C * rw
  }

  # Normalize every surviving control independently before QR.  The span is
  # unchanged, while numerical-rank decisions no longer depend on arbitrary
  # control units.  Zero columns on the no-FE path remain zero and are handled
  # by the pivoted fit.
  cn <- colnames(Cw)
  Cw <- vapply(seq_len(ncol(Cw)), function(j) .unit_norm(Cw[, j])$unit,
               numeric(n))
  dimnames(Cw) <- list(NULL, cn)

  # One matrix-response fit reuses a single pivoted QR factorisation of Cw.
  # Calling lm.fit() separately for y, x and Z performs the same factorisation
  # three times and is especially wasteful when the dense-control block is
  # wide.  C_Cdqrls applies the shared QR to every response column.  The
  # dimension-aware rank frontier matches the policy used for the instrument
  # Gram factorisation below.
  M <- cbind(yw, xw, Zw)
  rank_tol <- sqrt(max(nrow(Cw), ncol(Cw)) * .Machine$double.eps)
  Mr <- M - stats::lm.fit(Cw, M, tol = rank_tol)$fitted.values
  list(y = Mr[, 1L], x = as.numeric(Mr[, 2L]),
       Z = as.matrix(Mr[, 2L + seq_len(ncol(Z)), drop = FALSE]), rw = rw)
}

# Cholesky factor of the instrument Gram matrix Z'Z, with the package's
# singularity diagnosis. The single factorisation shared by every path
# (CJIVE first stage, CJAR/CJS coefficient kernels, the leverage diagnostic).
.stable_norm <- function(x) {
  if (!length(x)) return(0)
  ax <- abs(as.numeric(x))
  m <- max(ax)
  if (!is.finite(m) || m == 0) return(m)
  m * sqrt(sum((ax / m)^2))
}

# Unit vector plus a log-scale representation of its Euclidean norm.  Scaling
# by the largest entry first means the normalization itself cannot overflow,
# even when every input entry is finite but the raw Euclidean norm is not
# representable.  The log scale is used only when an ordinary norm ratio would
# overflow or underflow.
.unit_norm <- function(x) {
  x <- as.numeric(x)
  if (anyNA(x) || any(!is.finite(x))) {
    stop("internal inference input contains non-finite values.", call. = FALSE)
  }
  m <- if (length(x)) max(abs(x)) else 0
  if (!(m > 0)) {
    return(list(unit = rep.int(0, length(x)), max = 0,
                scaled_norm = 0, log2_norm = -Inf))
  }
  u <- x / m
  sn <- sqrt(sum(u * u))
  list(unit = u / sn, max = m, scaled_norm = sn,
       log2_norm = log2(m) + log2(sn))
}

# value * (||numerator|| / ||denominator||)^power without requiring either raw
# norm, or their ratio, to be representable.  Ordinary-scale arithmetic is
# retained whenever it is safe; the binary-log fallback is used only at the
# floating-point exponent boundary.
.rescale_by_norms <- function(value, numerator, denominator, power = 1) {
  if (value == 0) return(value)
  if (!(numerator$max > 0) || !(denominator$max > 0)) {
    stop("internal error: cannot rescale by a zero vector norm.", call. = FALSE)
  }
  ratio <- (numerator$max / denominator$max) *
    (numerator$scaled_norm / denominator$scaled_norm)
  if (is.finite(ratio) && ratio > 0) {
    out <- value * ratio^power
    if (is.finite(out) && out != 0) return(out)
  }

  logmag <- log2(abs(value)) +
    power * (numerator$log2_norm - denominator$log2_norm)
  if (logmag > log2(.Machine$double.xmax)) return(sign(value) * Inf)
  min_log <- log2(.Machine$double.xmin) + log2(.Machine$double.eps)
  if (logmag < min_log) return(sign(value) * 0)
  exponent <- floor(logmag)
  mantissa <- sign(value) * 2^(logmag - exponent)
  .mul_pow2(mantissa, exponent)
}

.col_norms <- function(M) {
  M <- as.matrix(M)
  vapply(seq_len(ncol(M)), function(j) .stable_norm(M[, j]), numeric(1))
}

.gamma_tol <- function(d) {
  de <- max(1, as.double(d)) * .Machine$double.eps
  if (de >= 1) Inf else de / (1 - de)
}

# Scale-free comparison of two vector norms.  Comparing their binary log norms
# avoids Inf/Inf and 0/0 when the raw norms lie outside the normal exponent
# range.  `after` is numerically zero relative to `before` at tolerance `tol`.
.relative_vector_zero <- function(after, before, tol) {
  an <- .unit_norm(after)
  bn <- .unit_norm(before)
  if (!(bn$max > 0)) return(!(an$max > 0))
  if (!(an$max > 0)) return(TRUE)
  an$log2_norm - bn$log2_norm <= log2(tol)
}

.gram_error <- function() {
  stop("instrument matrix is singular or numerically rank deficient after ",
       "partialling out covariates; its columns are collinear in the ",
       "residualized design.",
       call. = FALSE)
}

.gram_chol <- function(Z) {
  Z <- as.matrix(Z)
  scales <- .col_norms(Z)
  if (any(!is.finite(scales))) {
    stop("an instrument column norm is outside the representable ",
         "double-precision range; rescale the instrument columns and refit ",
         "(nonsingular column rescaling does not change their span).",
         call. = FALSE)
  }
  if (any(scales == 0)) .gram_error()

  gram <- crossprod(Z)
  ch <- if (all(is.finite(gram))) {
    tryCatch(chol(gram), error = function(e) NULL)
  } else {
    NULL
  }

  # Keep the original factorization on ordinary designs.  If the raw Gram
  # matrix over/underflows, factor the column-equilibrated design instead and
  # map the factor back to the original units.
  if (is.null(ch)) {
    Zunit <- sweep(Z, 2L, scales, "/")
    ch_unit <- tryCatch(chol(crossprod(Zunit)), error = function(e) NULL)
    if (is.null(ch_unit)) .gram_error()
    ch <- sweep(ch_unit, 2L, scales, "*")
  } else {
    ch_unit <- sweep(ch, 2L, scales, "/")
  }

  # ch_unit' ch_unit is the correlation Gram matrix.  Since the estimator
  # uses normal equations, rcond(ch_unit)^2 below d * eps is the numerical
  # rank frontier.  This check is invariant to the units of every column.
  rank_tol <- sqrt(max(nrow(Z), ncol(Z)) * .Machine$double.eps)
  rc <- tryCatch(rcond(ch_unit), error = function(e) NA_real_)
  if (!is.finite(rc) || rc <= rank_tol) .gram_error()
  ch
}

# Whitened instruments Ztil = Z R^-1 (so Ztil'Ztil = I_k). The k x k
# upper-triangular factor is inverted once (a triangular backsolve against
# the identity, O(k^3)) and applied with one GEMM: the same O(nk^2)
# arithmetic as the former BLAS-3 triangular solve of t(Z), but one n x k
# allocation instead of three (t(Z), the solve, the transpose back) --
# the difference between fitting and swapping at HDFE scale.  R is
# conditioning-gated by .gram_chol before it reaches this point, which is
# what makes the explicit triangular inverse as accurate here as the
# row-by-row solve (both are bounded by the same cond(R)).  In these
# coordinates every projection block is a plain Gram product,
# P_Z[g,h] = Ztil_g Ztil_h', which is what makes the per-cluster kernel
# below cheap on both sides of its size dispatch.
.whiten <- function(Z, R) Z %*% backsolve(R, diag(ncol(as.matrix(Z))))

# Full-sample first stage in whitened coordinates, computed once and reused.
# With t = Ztil'D: Xhat = Ztil t is the first-stage fit, h_i = ||Ztil_i||^2
# = z_i'(Z'Z)^-1 z_i the hat-matrix diagonal (used by JIVE), and R the
# Cholesky factor of Z'Z shared with the CJAR/CJS paths. No second inverse
# and no ZQinv: Ztil carries everything.
.first_stage <- function(D, Z, R = NULL) {
  if (is.null(R)) R <- .gram_chol(Z)
  Ztil <- .whiten(Z, R)
  t <- crossprod(Ztil, D)
  Xhat <- as.numeric(Ztil %*% t)
  e <- D - Xhat
  # Hat diagonal h_i = ||Ztil_i||^2, accumulated over column blocks: the
  # same O(nk) arithmetic as rowSums(Ztil^2) without materialising a second
  # full n x k matrix next to Ztil (the peak allocation that matters at
  # HDFE scale).
  h <- numeric(nrow(Ztil))
  k <- ncol(Ztil)
  for (jj in seq.int(1L, k, by = 64L)) {
    cols <- jj:min(jj + 63L, k)
    h <- h + rowSums(Ztil[, cols, drop = FALSE]^2)
  }
  list(R = R, Ztil = Ztil, t = t, Xhat = Xhat, e = e, h = h)
}

.leaveout_singular <- function(ng) {
  stop("leave-cluster-out fit is undefined for a cluster (I - H_g is ",
       "singular or numerically ill-conditioned): an instrument has no or ",
       "too little usable variation outside that cluster. ",
       sprintf("The offending cluster has %d observation(s).", ng),
       call. = FALSE)
}

# Solve (I - H)v = rhs for a symmetric projection block H.  Every caller uses
# the same Cholesky-based conditioning gate, so cjive() and iv_compare() cannot
# disagree at a numerical frontier.  The triangular reciprocal-condition
# estimate is compared on its square-root scale.  The spectral identity
# cond_2(I - H) = cond_2(chol(I - H))^2 motivates the eps^(1/4) cutoff;
# rcond() supplies a triangular norm estimate rather than an exact spectral
# condition number.  Below that frontier a formally invertible downdate cannot
# reliably retain the numerical accuracy required for agreement with the
# dense leave-cluster-out definition.  When requested, maxlev is then computed from H;
# iv_compare's comparison-only path skips that otherwise unused eigen pass.
# The solve uses the positive-definite structure rather than a general LU.
.projection_complement_solve <- function(H, rhs, ng, leverage) {
  m <- nrow(H)
  A <- diag(m) - H
  rcond_tol <- sqrt(.Machine$double.eps)
  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(ch)) .leaveout_singular(ng)
  rc_ch <- tryCatch(rcond(ch), error = function(e) NA_real_)
  if (!is.finite(rc_ch) || rc_ch <= sqrt(rcond_tol)) {
    .leaveout_singular(ng)
  }

  lev <- if (leverage) {
    max(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
  } else NA_real_
  sol <- backsolve(ch, backsolve(ch, rhs, transpose = TRUE))
  list(solution = as.numeric(sol), lev = lev)
}

# Per-cluster kernel on the whitened rows Zg = Ztil_g (n_g x k), shared by
# .leaveout_fit and .cluster_leverage so the size dispatch and the scalar
# shortcuts are written exactly once. G_g = Zg'Zg and H_g = Zg Zg' = P_Z[g,g]
# have the same non-zero eigenvalues, so the leverage lambda_max(P_Z[g,g])
# is read off whichever matrix the branch already formed:
#   n_g >  k: solve the k x k system, phat_g = Zg (I_k - G_g)^-1 u,
#             O(n_g k^2 + k^3);
#   n_g <= k: Woodbury on the small side, phat_g = a + H_g (I - H_g)^-1 a
#             with a = Zg u, O(n_g^2 k + n_g^3) -- the old block update,
#             surviving as the small-cluster branch.
# I_k - G_g is singular exactly when lambda_max = 1, the same condition under
# which the old I - H_g solve failed, so the error semantics carry over. The
# scalar cases (k = 1, n_g = 1) skip the matrix factorisation.  `leverage =
# FALSE` is used only where the caller discards maxlev, avoiding a cubic eigen
# pass without changing the fitted values.
.cluster_block <- function(Zg, k, u = NULL, leverage = TRUE) {
  ng <- nrow(Zg)
  if (ng > k) {
    if (k == 1L) {
      Gg <- sum(Zg^2)
      if (is.null(u)) return(list(lev = Gg))
      m <- 1 - Gg
      if (!is.finite(m) || m <= sqrt(.Machine$double.eps))
        .leaveout_singular(ng)
      return(list(lev = if (leverage) Gg else NA_real_,
                  phat = as.numeric(Zg) * (as.numeric(u) / m)))
    }
    Gg <- crossprod(Zg)
    if (is.null(u)) {
      return(list(lev = max(eigen(Gg, symmetric = TRUE,
                                  only.values = TRUE)$values)))
    }
    fit <- .projection_complement_solve(Gg, u, ng, leverage)
    list(lev = fit$lev, phat = as.numeric(Zg %*% fit$solution))
  } else {
    if (ng == 1L) {
      Hg <- sum(Zg^2)
      if (is.null(u)) return(list(lev = Hg))
      m <- 1 - Hg
      if (!is.finite(m) || m <= sqrt(.Machine$double.eps))
        .leaveout_singular(ng)
      return(list(lev = if (leverage) Hg else NA_real_,
                  phat = sum(as.numeric(Zg) * as.numeric(u)) / m))
    }
    Hg <- tcrossprod(Zg)
    if (is.null(u)) {
      return(list(lev = max(eigen(Hg, symmetric = TRUE,
                                  only.values = TRUE)$values)))
    }
    a <- as.numeric(Zg %*% u)
    fit <- .projection_complement_solve(Hg, a, ng, leverage)
    # (I - H)^-1 a = a + H(I - H)^-1 a, so the solution is already the
    # Woodbury fitted value; do not perform the redundant final matrix product.
    list(lev = fit$lev, phat = fit$solution)
  }
}

# Leave-cluster-out fitted values for a partition `groups` (list of row-index
# sets), in whitened coordinates. With u = t - t_g (t_g = Ztil_g' D_g) the R
# factors of pi_(-g) = (Z'Z - Z_g'Z_g)^-1 (Z'D - Z_g'D_g) cancel exactly:
#   phat_[g] = Ztil_g (I_k - G_g)^-1 u,
# dispatched per cluster to whichever side is smaller (see .cluster_block);
# no n_g x n_g object is formed when n_g > k. Exact against the brute-force
# definition. For singleton clusters this collapses to the improved JIVE
# (IJIVE; Ackerberg-Devereux 2009), the leave-one-out fit; Z is already
# residualised on the covariates (see .partial_out), which is exactly what
# distinguishes IJIVE from the original JIVE (Angrist-Imbens-Krueger 1999).
# maxlev = max_g lambda_max(P_Z[g,g]) comes free from the same pass.
.leaveout_fit <- function(D, Ztil, t, groups, k, leverage = TRUE) {
  phat <- numeric(length(D))
  maxlev <- if (leverage) 0 else NA_real_
  for (idx in groups) {
    Zg <- Ztil[idx, , drop = FALSE]
    blk <- .cluster_block(Zg, k, t - crossprod(Zg, D[idx]), leverage)
    phat[idx] <- blk$phat
    if (leverage && blk$lev > maxlev) maxlev <- blk$lev
  }
  list(phat = phat, maxlev = maxlev)
}

# Assumption A2(ii) diagnostic: max_g ||P_Z[g,g]||_2, the largest spectral
# norm of a diagonal cluster block of the instrument projection. The paper's
# A2(ii) bounds the SQUARE of this norm, ||P_Z[g,g]||_2^2 <= C < 1; the
# eigenvalues of a projection block lie in [0, 1], so bounding the norm and
# bounding its square are equivalent -- what is reported is the norm itself.
# Values near 1 mean one cluster nearly spans the instrument space (an
# instrument has almost no variation outside that cluster) and the CLT behind
# Theorem 1 is strained. The leverage-only pass over the same per-cluster
# branch as .leaveout_fit -- one shared size dispatch, not two copies.
.cluster_leverage <- function(Ztil, groups, k) {
  maxlev <- 0
  for (idx in groups) {
    lev <- .cluster_block(Ztil[idx, , drop = FALSE], k)$lev
    if (lev > maxlev) maxlev <- lev
  }
  maxlev
}

# Constructed instruments. Only p-hat changes between estimators; inference does
# not. .phat_jive is the leave-one-out fit; because fs is built from the
# residualised Z (.first_stage on partialled-out data), it is the improved JIVE
# (IJIVE; Ackerberg-Devereux 2009), not the original JIVE (Angrist-Imbens-
# Krueger 1999) that jackknifes the covariates with the instruments.
.phat_2sls <- function(fs) fs$Xhat

.phat_jive <- function(fs, D) {
  if (max(fs$h) >= 1 - 1e-12) {
    stop("JIVE is undefined: an observation's leverage is numerically 1 (an ",
         "instrument has no variation outside that observation). ",
         sprintf("max(h) = %.15g at row %d.", max(fs$h), which.max(fs$h)),
         call. = FALSE)
  }
  (fs$Xhat - fs$h * D) / (1 - fs$h)
}

.phat_cjive <- function(fs, D, groups)
  .leaveout_fit(D, fs$Ztil, fs$t, groups, ncol(fs$Ztil),
                leverage = FALSE)$phat

# Just-identified cluster-robust IV sandwich on the residualised data:
#   beta = (p'Y)/(p'D),  SE = sqrt(sum_g S_g^2 * G/(G-1)) / |p'D|,  S_g = sum_g p_i e_i.
# The critical-value/p-value step is dispatched on `inference`: "asymptotic"
# uses the standard normal, "t" uses t(G-1) (standard few-cluster practice).
# A future inference method is a new branch here, not a refactor.
.iv_inference <- function(phat, x, y, cluster, level = 0.95,
                          inference = "asymptotic") {
  cl <- droplevels(as.factor(cluster))
  G <- nlevels(cl)
  pu <- .unit_norm(phat)
  xu <- .unit_norm(x)
  yu <- .unit_norm(y)
  den_unit <- if (pu$max > 0 && xu$max > 0) {
    sum(pu$unit * xu$unit)
  } else 0
  if (!is.finite(den_unit) ||
      abs(den_unit) <= .gamma_tol(length(x))) {
    stop("IV estimator denominator is numerically zero; the residualized ",
         "first stage does not identify `x`.", call. = FALSE)
  }

  # Work entirely in unit-norm coordinates.  This is algebraically the same
  # covariance ratio and cluster sandwich, while remaining invariant to the
  # units of y and x and to a common rescaling of precision weights.  Taking a
  # stable norm of the cluster scores avoids the underflow-prone sum(Sg^2).
  theta <- sum(pu$unit * yu$unit) / den_unit
  score <- pu$unit * (yu$unit - theta * xu$unit)
  Sg <- as.numeric(rowsum(score, cl, reorder = FALSE))
  se_theta <- .stable_norm(Sg) * sqrt(G / (G - 1)) / abs(den_unit)
  beta <- .rescale_by_norms(theta, yu, xu)
  se <- .rescale_by_norms(se_theta, yu, xu)
  if (!is.finite(beta) || !is.finite(se)) {
    stop("IV coefficient or standard error is outside the representable ",
         "double-precision range; rescale `y` and `x` by a common factor to ",
         "preserve coefficient units, or track the unit conversion if they ",
         "are scaled separately, then refit.",
         call. = FALSE)
  }

  # A zero estimated variance and a nonzero coefficient imply an infinite
  # zero-null statistic.  At coefficient == SE == 0 the ratio is genuinely
  # undefined; return NA explicitly rather than leaking a NaN from 0/0.
  stat <- if (se_theta > 0) theta / se_theta
          else if (theta > 0) Inf
          else if (theta < 0) -Inf
          else NA_real_
  a <- 1 - (1 - level) / 2
  switch(inference,
    asymptotic = {
      zc <- stats::qnorm(a)
      pval <- if (is.na(stat)) NA_real_ else 2 * stats::pnorm(-abs(stat))
    },
    t = {
      zc <- stats::qt(a, df = G - 1)
      pval <- if (is.na(stat)) NA_real_
              else 2 * stats::pt(-abs(stat), df = G - 1)
    },
    stop("unknown `inference` method: ", inference, call. = FALSE)
  )
  list(coefficient = beta, se = se, statistic = stat,
       p.value = pval,
       conf.low = beta - zc * se, conf.high = beta + zc * se, G = G)
}

# A grouping (factor/character) z becomes a dummy design; numeric z is used as
# is. Here `intercept` is the tested internal coding flag: when true, one
# reference level is dropped; otherwise every level is retained, including the
# sole level of a one-level factor.
.expand_z <- function(z, intercept) {
  if (is.factor(z) || is.character(z)) {
    f <- as.factor(z)
    Z <- matrix(0, length(f), nlevels(f),
                dimnames = list(NULL, paste0("f", levels(f))))
    if (length(f)) Z[cbind(seq_along(f), as.integer(f))] <- 1
    if (intercept) Z <- Z[, -1L, drop = FALSE]
    return(list(Z = Z, grouping = TRUE, group = f))
  }
  if (!is.numeric(z))
    stop("`z` must be numeric, or a factor/character grouping variable.",
         call. = FALSE)
  Z <- as.matrix(z)
  storage.mode(Z) <- "double"
  list(Z = Z, grouping = FALSE, group = NULL)
}

.assert_finite <- function(v, name) {
  if (anyNA(v) || any(!is.finite(v)))
    stop(sprintf("`%s` contains missing or non-finite values; remove or impute them first.",
                 name), call. = FALSE)
}

.require_numeric <- function(v, name) {
  if (!is.numeric(v) || is.factor(v))
    stop(sprintf("`%s` must be numeric; factor and character values are not allowed.",
                 name), call. = FALSE)
  invisible(NULL)
}

.check_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop(sprintf("`%s` must be a single non-missing logical value.", name),
         call. = FALSE)
  invisible(NULL)
}

.check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  invisible(NULL)
}

# Stable label for the single endogenous regressor.  Bare names retain their
# spelling (including non-syntactic names); `$` and subsetting expressions use
# the underlying column/object name.  Computed expressions deliberately fall
# back to "x" rather than leaking a long deparse into tables.
.term_label <- function(expr, fallback = "x") {
  if (is.name(expr)) {
    out <- as.character(expr)
    if (length(out) == 1L && nzchar(out)) return(out)
  }
  if (is.call(expr) && length(expr) >= 2L) {
    op <- as.character(expr[[1L]])
    if (op %in% c("[", "[[")) {
      # A literal character index names the selected endogenous column, as in
      # d[["treatment"]] or d[, "treatment"]. Dynamic/computed indices keep
      # the underlying object label rather than guessing at run time.
      idx <- as.list(expr)[-(1:2)]
      for (j in rev(seq_along(idx))) {
        if (is.character(idx[[j]]) && length(idx[[j]]) == 1L &&
            nzchar(idx[[j]])) return(idx[[j]])
      }
      return(.term_label(expr[[2L]], fallback))
    }
    if (op == "$" && length(expr) == 3L) {
      rhs <- expr[[3L]]
      if (is.name(rhs)) return(as.character(rhs))
      if (is.character(rhs) && length(rhs) == 1L && nzchar(rhs)) return(rhs)
    }
  }
  fallback
}

# Complete topology of a confidence set.  `shape` is retained for backwards
# compatibility as a tail-classification label; these fields are the complete
# description when bounded middle components coexist with one or both tails.
.set_topology <- function(conf_set) {
  cs <- as.matrix(conf_set)
  if (ncol(cs) != 2L || anyNA(cs) ||
      (nrow(cs) && any(cs[, 1L] > cs[, 2L]))) {
    stop("internal error: malformed confidence-set endpoint matrix.",
         call. = FALSE)
  }
  list(
    n_components = nrow(cs),
    unbounded_left = nrow(cs) > 0L && is.infinite(cs[1L, 1L]) &&
      cs[1L, 1L] < 0,
    unbounded_right = nrow(cs) > 0L &&
      is.infinite(cs[nrow(cs), 2L]) && cs[nrow(cs), 2L] > 0
  )
}

# The backwards-compatible tail-classification label, derived from the
# topology above.  The single shape classifier of the package: the
# inverters and the tidiers both call this.
.shape_from_set <- function(conf_set) {
  tp <- .set_topology(conf_set)
  if (tp$n_components == 0L) "empty"
  else if (tp$n_components == 1L && tp$unbounded_left &&
           tp$unbounded_right) "whole_line"
  else if (xor(tp$unbounded_left, tp$unbounded_right)) "ray"
  else if (tp$unbounded_left && tp$unbounded_right) "two_rays"
  else "bounded"
}

.resolve_iv_tests <- function(tests, omitted = FALSE) {
  allowed <- c("cjar", "cjscore")
  if (isTRUE(omitted)) return(c("cjar", "cjscore"))
  if (is.null(tests) || length(tests) == 0L) return(character(0))
  if (!is.character(tests)) {
    stop("`tests` must be NULL or a character subset of ",
         "c(\"cjar\", \"cjscore\").",
         call. = FALSE)
  }
  if (anyNA(tests) || any(!nzchar(tests))) {
    stop("`tests` must not contain NA or empty strings.", call. = FALSE)
  }
  resolved <- vapply(tests, function(value) {
    if (value %in% allowed) return(value)
    hit <- allowed[startsWith(allowed, value)]
    if (length(hit) == 1L) return(hit)
    if (length(hit) > 1L) {
      stop("ambiguous partial `tests` value ", sQuote(value), ": matches ",
           paste(hit, collapse = ", "), ".", call. = FALSE)
    }
    stop("unknown `tests` value ", sQuote(value), ".", call. = FALSE)
  }, character(1))
  allowed[allowed %in% unique(resolved)]
}

.check_dots <- function(dots, fn) {
  if (is.null(dots) || length(dots) == 0L) return(invisible(NULL))
  nms <- names(dots)
  if (is.null(nms)) nms <- rep.int("", length(dots))
  labels <- ifelse(nzchar(nms), nms, paste0("..", seq_along(dots)))
  stop("unused argument(s) in `", fn, "(... )`: ",
       paste(labels, collapse = ", "), ".", call. = FALSE)
}

# UseMethod() partially matches names against a method's formals before the
# method body can inspect `...`.  Check named arguments at the generic first;
# unnamed arguments are left for ordinary positional matching and the method's
# stricter .check_dots() guard.
.check_named_dots <- function(dots, allowed, fn) {
  if (is.null(dots) || length(dots) == 0L) return(invisible(NULL))
  nms <- names(dots)
  if (is.null(nms)) return(invisible(NULL))
  bad <- unique(nms[nzchar(nms) & !(nms %in% allowed)])
  if (!length(bad)) return(invisible(NULL))
  stop("unused argument(s) in `", fn, "(...)`: ",
       paste(bad, collapse = ", "), ".", call. = FALSE)
}

# Validate names only when UseMethod() will enter one of this package's two
# fitting methods.  A third-party S3 method may intentionally define different
# arguments and must remain free to receive them.
.check_dispatch_dots <- function(y, dots, allowed, fn, caller) {
  generic <- get(fn, mode = "function")
  table <- get(".__S3MethodsTable__.", envir = environment(generic),
               inherits = FALSE)
  lookup <- function(method) {
    out <- get0(method, envir = caller, mode = "function", inherits = TRUE)
    if (is.null(out)) get0(method, envir = table, mode = "function",
                           inherits = FALSE) else out
  }

  method <- NULL
  for (cl in class(y)) {
    method <- lookup(paste(fn, cl, sep = "."))
    if (!is.null(method)) break
  }
  if (is.null(method)) method <- lookup(paste0(fn, ".default"))

  builtins <- list(
    get0(paste0(fn, ".default"), envir = table, inherits = FALSE),
    get0(paste0(fn, ".formula"), envir = table, inherits = FALSE)
  )
  if (!any(vapply(builtins, identical, logical(1), y = method)))
    return(invisible(NULL))
  .check_named_dots(dots, allowed, fn)
}

# Coerce `fixed_effects` (a factor, or a list/data.frame of factors) to a list
# of factors of length n.
.coerce_fe <- function(fixed_effects, n) {
  if (is.null(fixed_effects)) return(NULL)
  fe <- if (is.list(fixed_effects)) as.list(fixed_effects)
        else list(fixed_effects)
  if (length(fe) == 0L) return(NULL)
  for (k in seq_along(fe)) {
    f <- fe[[k]]
    if (length(f) != n)
      stop("each fixed-effect factor must have length n.", call. = FALSE)
    if (anyNA(f))
      stop("`fixed_effects` contains missing values.", call. = FALSE)
    fe[[k]] <- droplevels(as.factor(f))
  }
  fe
}

# Validate inputs, expand z, partial out, build the cluster partition. Shared by
# cjive() and iv_compare().
.prep_data <- function(y, x, z, cluster, controls, weights, intercept,
                       fixed_effects = NULL) {
  .check_flag(intercept, "intercept")
  .require_numeric(y, "y")
  .require_numeric(x, "x")
  y <- as.numeric(y)
  x <- as.numeric(x)
  n <- length(y)

  if (length(x) != n) stop("`x` and `y` must have the same length.", call. = FALSE)
  if (length(cluster) != n) stop("`cluster` must have length n.", call. = FALSE)
  .assert_finite(y, "y")
  .assert_finite(x, "x")
  if (anyNA(cluster)) stop("`cluster` contains missing values.", call. = FALSE)
  if (!is.null(weights)) {
    .require_numeric(weights, "weights")
    weights <- as.numeric(weights)
    if (length(weights) != n) stop("`weights` must have length n.", call. = FALSE)
    if (any(!is.finite(weights)) || any(weights <= 0))
      stop("`weights` must be finite and strictly positive.", call. = FALSE)
  }

  fe <- .coerce_fe(fixed_effects, n)

  if (anyNA(z)) stop("`z` contains missing values.", call. = FALSE)
  # Every non-empty FE dummy span contains the global intercept.  Full-level
  # coding of a grouping instrument would therefore become singular after FE
  # absorption even when the user requested intercept = FALSE.
  zinfo <- .expand_z(z, intercept = intercept || !is.null(fe))
  if (nrow(zinfo$Z) != n) stop("`z` must have n rows.", call. = FALSE)
  if (ncol(zinfo$Z) < 1L)
    stop("`z` supplies no instrument columns (a grouping factor needs at least 2 levels).",
         call. = FALSE)
  .assert_finite(zinfo$Z, "z")

  if (!is.null(controls)) {
    if (NROW(controls) != n) stop("`controls` must have n rows.", call. = FALSE)
    if (anyNA(controls)) stop("`controls` contains missing values.", call. = FALSE)
    controls <- if (is.data.frame(controls)) {
      if (ncol(controls) == 0L) {
        matrix(numeric(0), nrow = n, ncol = 0L)
      } else {
        # Absorbed FE already span the intercept, so reference coding is the
        # non-redundant factor-control basis even when intercept = FALSE.
        cf <- if (intercept || !is.null(fe)) ~ . else ~ 0 + .
        mm <- stats::model.matrix(cf, data = controls)
        mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
      }
    } else as.matrix(controls)
    if (!is.numeric(controls))
      stop("`controls` must expand to a numeric matrix.", call. = FALSE)
    .assert_finite(controls, "controls")
    if (ncol(controls) == 0L) controls <- NULL
  }

  cl <- droplevels(as.factor(cluster))
  G <- nlevels(cl)
  if (G < 2L) stop("at least 2 clusters are required.", call. = FALSE)

  # Do not impose a raw group-by-cluster spread rule on the dense paths.  They
  # operate on the encoded, globally residualised instrument columns, and raw
  # support is not equivalent to rank of that design: after intercept FWL, a
  # dummy whose level lies in one cluster is a nonzero constant outside that
  # cluster and its leave-out Gram can still be full rank.  The common Gram
  # checks below and the leave-out kernel therefore decide support from the
  # matrix actually used.  The explicit leaveout_mean path has a different
  # group-mean requirement and checks positive outside-group mass itself.

  # Unlike an instrument group, a fixed-effect level lying entirely inside one
  # cluster is legitimate: FE are absorbed globally by .partial_out before the
  # leave-cluster-out step, consistent with the dense-controls convention (each
  # observation's own cluster enters the partialling out; see the many-controls
  # caveat in ?cjive).
  rw <- if (is.null(weights)) rep.int(1, n) else sqrt(weights)
  yw0 <- y * rw
  xw0 <- x * rw
  Zw0 <- zinfo$Z * rw
  .assert_finite(yw0, "weighted y")
  .assert_finite(xw0, "weighted x")
  .assert_finite(Zw0, "weighted z")

  po <- .partial_out(y, x, zinfo$Z, controls, weights, intercept,
                     fixed_effects = fe)
  .assert_finite(po$y, "residualized y")
  .assert_finite(po$x, "residualized x")
  .assert_finite(po$Z, "residualized z")

  dense_count <- if (is.null(controls)) 0L else ncol(as.matrix(controls))
  fe_count <- if (is.null(fe)) {
    if (intercept) 1L else 0L
  } else {
    # Reference coding spans the same joint FE space with one intercept and
    # L_j - 1 columns per supplied dimension.  This remains a raw upper bound:
    # nesting, duplicate dimensions and dense-control collinearity reduce rank.
    1L + sum(vapply(fe, nlevels, 0L) - 1L)
  }
  k_controls <- dense_count + fe_count
  proj_tol <- .gamma_tol(max(n, ncol(zinfo$Z), k_controls))
  if (!is.null(fe)) proj_tol <- max(proj_tol, 1e-10)

  if (.relative_vector_zero(po$x, xw0, proj_tol))
    stop("`x` has no usable variation after partialling out covariates.",
         call. = FALSE)
  if (.relative_vector_zero(po$y, yw0, proj_tol))
    stop("`y` has no usable variation after partialling out covariates.",
         call. = FALSE)

  bad_z <- which(vapply(seq_len(ncol(po$Z)), function(j) {
    .relative_vector_zero(po$Z[, j], Zw0[, j], proj_tol)
  }, logical(1)))
  if (length(bad_z)) {
    zn <- colnames(zinfo$Z)
    if (is.null(zn)) zn <- paste0("z", seq_len(ncol(zinfo$Z)))
    zn[!nzchar(zn)] <- paste0("z", which(!nzchar(zn)))
    stop("instrument column(s) have no usable variation after partialling ",
         "out covariates: ", paste(zn[bad_z], collapse = ", "), ".",
         call. = FALSE)
  }

  R <- .gram_chol(po$Z)
  list(y = po$y, x = po$x, Z = po$Z, cluster = cl,
       groups = split(seq_len(n), cl), n = n, G = G, k = ncol(zinfo$Z),
       grouping = zinfo$grouping, group = zinfo$group, rw = po$rw,
       weights = weights, R = R,
       fe_dims = if (is.null(fe)) 0L else length(fe),
       fe_levels = if (is.null(fe)) 0L else sum(vapply(fe, nlevels, 0L)),
       k_controls = k_controls)
}

# Leave-cluster-out group mean, FLM's closed form for the pure judge design:
# p_i = (S_j - S_{j,g}) / (N_j - N_{j,g}), the weighted mean of x over group j
# outside i's cluster. Differs from the dense route by an intercept term of
# order n_g/n (~ 1/G balanced): the leave-out sample re-estimates the
# intercept direction, and the gap does not shrink in n_g alone.
.leaveout_mean <- function(x, group, cluster, weights) {
  g <- as.factor(group)
  cl <- as.factor(cluster)
  w <- if (is.null(weights)) rep.int(1, length(x)) else weights

  Sj <- tapply(w * x, g, sum)[g]
  Nj <- tapply(w, g, sum)[g]
  gc <- interaction(g, cl, drop = TRUE)
  Sjg <- tapply(w * x, gc, sum)[gc]
  Njg <- tapply(w, gc, sum)[gc]

  denom <- Nj - Njg
  if (any(denom <= 0)) {
    stop("leave-cluster-out group mean is undefined: an instrument group lies ",
         "entirely in one cluster.", call. = FALSE)
  }
  as.numeric((Sj - Sjg) / denom)
}

# The three test-path advisories (small G, dominating cluster, near-unit
# within-cluster leverage), shared by cjar(), cjscore() and iv_infer() so the
# thresholds and wording live in exactly one place. `label` names the
# calibration in the message ("CJAR", "CJS", or "CJAR/CJS" from iv_infer(),
# which fires each advisory once for the whole workflow).
.advise <- function(G, groups, n, maxlev, label) {
  if (G < 20L)
    warning(sprintf("only %d clusters: the %s calibration is asymptotic in the number of clusters.",
                    G, label), call. = FALSE)
  ngmax <- max(lengths(groups))
  if (ngmax / n > 0.2)
    warning(sprintf("the largest cluster holds %.0f%% of the sample; the %s asymptotics require no dominating cluster.",
                    100 * ngmax / n, label), call. = FALSE)
  if (maxlev > 0.99)
    warning(sprintf("max within-cluster leverage = %.4f: a cluster nearly spans the instrument space (Assumption A2(ii) of Ligtenberg 2025 bounds its square below 1, equivalently the norm itself).",
                    maxlev), call. = FALSE)
  invisible(NULL)
}

# The design/diagnostic fields common to every fitted object (cjive, cjar,
# cjscore). Spliced into each *_build constructor so the three classes carry
# an identical block; the class-specific fields (coefficients, sets, k/p, path)
# stay in the callers. Field-name access downstream is order-independent.
.fit_common <- function(d, eff, maxlev) {
  list(n = d$n, G = d$G, maxlev = maxlev,
       F_eff = eff$F_eff, K_eff = eff$K_eff, F_eff_crit = eff$F_eff_crit,
       ng_max = max(lengths(d$groups)),
       fe_dims = d$fe_dims, fe_levels = d$fe_levels,
       k_controls = d$k_controls, n_dropped = 0L)
}

# Scalar-null and level validation shared by the test paths (cjar, cjscore,
# iv_infer): beta0 a single finite number, level strictly inside (0, 1).
.check_beta0_level <- function(beta0, level) {
  if (!is.numeric(beta0) || length(beta0) != 1L || !is.finite(beta0))
    stop("`beta0` must be a single finite number.", call. = FALSE)
  .check_level(level)
  invisible(NULL)
}
