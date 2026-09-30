# Internal helpers: high-dimensional fixed-effect absorption.  The
# matrix-free alternating-projection solver (rowsum demeaning, CG
# refinement, certificates) behind `fixed_effects =`, ending in .demean.

# Weighted projector data for the fixed-effect solver.  M is already in the
# sqrt(w) coordinates used by .partial_out.  For one FE dimension j, P_j X is
# the projection onto its transformed dummy span.  The groupwise-normalised
# representation below is algebraically identical to
#
#   sqrt(w) * [rowsum(sqrt(w) * X, f) / rowsum(w, f)][f, ],
#
# but it is unchanged by a common rescaling of the weights and avoids forming
# raw weight sums that can overflow or underflow.  Every application is one
# C-level rowsum call; there are no per-level R loops or dummy matrices.
.hdfe_problem <- function(M, fe_list, weights, tol, maxit) {
  if (!is.numeric(tol) || is.factor(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a single finite positive number.", call. = FALSE)
  }
  if (!is.numeric(maxit) || is.factor(maxit) || length(maxit) != 1L ||
      !is.finite(maxit) || maxit <= 0 || maxit != floor(maxit) ||
      maxit > .Machine$integer.max) {
    stop("`maxit` must be a single positive finite integer-like value.",
         call. = FALSE)
  }
  if (!is.numeric(M) || is.factor(M))
    stop("`M` must be a numeric matrix.", call. = FALSE)
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  if (anyNA(M) || any(!is.finite(M)))
    stop("`M` must contain only finite values.", call. = FALSE)
  n <- nrow(M)
  if (n < 1L)
    stop("`M` must have at least one row.", call. = FALSE)
  if (ncol(M) < 1L)
    stop("`M` must have at least one column.", call. = FALSE)
  if (!is.list(fe_list) || length(fe_list) < 1L)
    stop("`fe_list` must be a non-empty list of fixed-effect vectors.",
         call. = FALSE)

  if (is.null(weights)) {
    w <- rep.int(1, n)
  } else {
    if (!is.numeric(weights) || is.factor(weights) || length(weights) != n ||
        anyNA(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
      stop("`weights` must have length n and contain finite, strictly positive numeric values.",
           call. = FALSE)
    }
    w <- as.numeric(weights)
  }

  fe <- vector("list", length(fe_list))
  for (j in seq_along(fe_list)) {
    f <- fe_list[[j]]
    if (length(f) != n)
      stop("each fixed-effect vector must have length n.", call. = FALSE)
    if (anyNA(f))
      stop("fixed-effect vectors must not contain missing values.",
           call. = FALSE)
    fe[[j]] <- droplevels(as.factor(f))
  }

  # A coarser partition's dummy span is contained in a finer partition's span.
  # Remove exact containment (including relabelled duplicates) before choosing
  # the solver: it changes neither the joint projection nor the omitted FE's
  # orthogonality condition, and prevents an algebraically zero reduced
  # two-factor operator from being mistaken for numerical ill-conditioning.
  # The comparisons loop only over FE dimensions, never over levels.
  if (length(fe) > 1L) {
    codes <- lapply(fe, function(f) match(f, unique(f)))
    is_subset <- function(a, b) {
      # Va is contained in Vb iff every level of b maps to one level of a.
      first <- a[match(seq_len(max(b)), b)]
      all(a == first[b])
    }
    keep <- rep(TRUE, length(fe))
    for (j in seq_along(fe)) {
      for (k in seq_along(fe)) {
        if (j == k) next
        jk <- is_subset(codes[[j]], codes[[k]])
        if (!jk) next
        equivalent <- is_subset(codes[[k]], codes[[j]])
        if (!equivalent || k < j) {
          keep[j] <- FALSE
          break
        }
      }
    }
    fe <- fe[keep]
    codes <- codes[keep]

    # The reduced two-space solve is algebraically symmetric but its finite-
    # precision Krylov path is not.  Put retained partitions in one canonical
    # order so reversing fixed_effects cannot change a return into a maxit
    # stop (or conversely).  Fewer-level partitions come first; ties are
    # broken lexicographically on first-occurrence canonical codes, which is
    # invariant to factor labels.  This loops only over FE dimensions.
    if (length(fe) > 1L) {
      before <- function(a, b) {
        na <- max(a)
        nb <- max(b)
        if (na != nb) return(na < nb)
        different <- a != b
        if (!any(different)) return(FALSE)
        first <- which.max(different)
        a[first] < b[first]
      }
      ord <- seq_along(fe)
      for (j in seq.int(2L, length(ord))) {
        pos <- j
        while (pos > 1L &&
               before(codes[[ord[pos]]], codes[[ord[pos - 1L]]])) {
          tmp <- ord[pos - 1L]
          ord[pos - 1L] <- ord[pos]
          ord[pos] <- tmp
          pos <- pos - 1L
        }
      }
      fe <- fe[ord]
    }
  }

  # One common square-root-weight scale is harmless.  A second, groupwise
  # scale keeps every group's denominator representable even when weights
  # span nearly the full floating-point exponent range.
  rw <- sqrt(w)
  rw <- rw / max(rw)
  specs <- lapply(fe, function(f) {
    g <- as.integer(f)
    group_scale <- as.numeric(rowsum(rw, g, reorder = TRUE))
    a <- rw / group_scale[g]
    denom <- as.numeric(rowsum(a * a, g, reorder = TRUE))
    if (any(!is.finite(a)) || any(!is.finite(denom)) || any(denom <= 0)) {
      stop("fixed-effect projection weights are non-finite or numerically degenerate.",
           call. = FALSE)
    }
    list(g = g, a = a, denom = denom)
  })

  # A pair of distinct weighted FE dummies can be indistinguishable at double
  # precision even when neither factor has a globally extreme weight range.
  # .hdfe_pair_gap_witness() computes a stable local principal-angle witness
  # for every unordered FE pair.  Taking the minimum is invariant to FE-list
  # order.  Relabelling can choose another cell when weights tie, but every
  # tied-cell gap is far above the numerical stop threshold, so the safety
  # decision is label-invariant.  A cheap exact prefilter avoids the pair-cell
  # aggregation on ordinary geometries: if q_min is the smallest
  # observation loading in any unit transformed dummy, every distinct cell
  # witness has lambda >= q_min^2 / 2.  Hence the scan cannot fire when that
  # lower bound already exceeds the same numerical safety floor.
  work_tol <- min(tol, 1e-12)
  geometry_floor <- max(2 * work_tol^2,
                        128 * .Machine$double.eps * length(specs))
  loading_floor <- max(work_tol,
                       128 * .Machine$double.eps * length(specs))
  minimum_loading <- min(vapply(specs, function(spec) {
    min(spec$a / sqrt(spec$denom[spec$g]))
  }, numeric(1)))
  if (!is.finite(minimum_loading) || minimum_loading < 0) {
    stop("fixed-effect projection loadings are non-finite or numerically degenerate.",
         call. = FALSE)
  }
  # When the loading already sits at or below its own safety floor, .demean()
  # stops on the loading condition before the pair gap is ever consulted, so
  # the O(J^2 n log n) pair scan cannot change the decision and is skipped.
  pair_gap <- Inf
  if (length(specs) > 1L &&
      minimum_loading > loading_floor &&
      minimum_loading <= sqrt(2 * geometry_floor)) {
    for (j in seq_len(length(specs) - 1L)) {
      for (k in seq.int(j + 1L, length(specs))) {
        pair_gap <- min(pair_gap,
                        .hdfe_pair_gap_witness(specs[[j]], specs[[k]]))
      }
    }
  }
  if (!is.finite(pair_gap) && !is.infinite(pair_gap)) {
    stop("fixed-effect geometry diagnostic produced a non-finite value.",
         call. = FALSE)
  }

  list(M = M, specs = specs, n = n, p = ncol(M),
       tol = tol, maxit = as.integer(maxit),
       minimum_loading = minimum_loading,
       minimum_pair_gap = pair_gap)
}

# Stable local principal-angle witness for two weighted FE spaces.  For unit
# transformed dummies u_g and v_h sharing cell c,
#
#   rho_c^2 = (W_c / W_g) (W_c / W_h),
#   lambda_c = 1 - rho_c.
#
# If c is the largest-weight cell for both levels, compute the two outside
# shares directly and form
#
#   delta_c = 1 - rho_c^2 = o_g + o_h - o_g o_h,
#   lambda_c = delta_c / (1 + sqrt(1 - delta_c)).
#
# Direct outside sums retain a weak incidence edge that `1 - share` would
# round to zero.  A small lambda is an explicit near-common-direction witness:
# (u_g - v_h) / sqrt(2 lambda) has unit forward norm but individual projection
# corrections sqrt(lambda / 2).  Exact equal supports are harmless duplicate
# directions and are excluded.  All aggregation uses rowsum; the only loops
# in the caller are over FE dimensions, never levels.
.hdfe_pair_gap_witness <- function(left, right) {
  gl <- left$g
  gr <- right$g
  # Code only observed (g_l, g_r) pairs.  interaction(..., drop = TRUE) is not
  # HDFE-safe because it constructs all Cartesian level labels before dropping
  # empty cells.  Sorting is O(n) memory regardless of the number of possible
  # level combinations.
  ord <- order(gl, gr)
  gl_ord <- gl[ord]
  gr_ord <- gr[ord]
  boundary <- c(TRUE, gl_ord[-1L] != gl_ord[-length(gl_ord)] |
                       gr_ord[-1L] != gr_ord[-length(gr_ord)])
  cell <- integer(length(gl))
  cell[ord] <- cumsum(boundary)
  nc <- sum(boundary)
  first <- ord[boundary]
  gl_cell <- gl[first]
  gr_cell <- gr[first]

  mass_left <- left$a * left$a / left$denom[gl]
  mass_right <- right$a * right$a / right$denom[gr]
  share_left <- as.numeric(rowsum(mass_left, cell, reorder = TRUE))
  share_right <- as.numeric(rowsum(mass_right, cell, reorder = TRUE))

  dominant <- function(level, share) {
    ord <- order(level, -share, seq_along(level))
    out <- logical(length(level))
    out[ord[!duplicated(level[ord])]] <- TRUE
    out
  }
  top_left <- dominant(gl_cell, share_left)
  top_right <- dominant(gr_cell, share_right)

  outside_left <- as.numeric(rowsum(
    replace(share_left, top_left, 0), gl_cell, reorder = TRUE))
  outside_right <- as.numeric(rowsum(
    replace(share_right, top_right, 0), gr_cell, reorder = TRUE))
  o_left <- outside_left[gl_cell]
  o_right <- outside_right[gr_cell]

  cell_n <- tabulate(cell, nbins = nc)
  left_n <- tabulate(gl)[gl_cell]
  right_n <- tabulate(gr)[gr_cell]
  distinct_support <- cell_n != left_n | cell_n != right_n
  use <- top_left & top_right & distinct_support
  if (!any(use)) return(Inf)

  delta <- o_left[use] + o_right[use] - o_left[use] * o_right[use]
  if (any(!is.finite(delta)) || any(delta < 0)) return(0)
  delta <- pmin(delta, 1)
  lambda <- delta / (1 + sqrt(pmax(0, 1 - delta)))
  min(lambda)
}

# Use a fast C-level sum of squares when it is safely in the normal floating-
# point range, but fall back to the package's scale-first norm if squaring
# underflows or overflows.  Near-null FE directions can produce an O(1e-200)
# operator image from a unit input; treating its squared norm as exact zero
# would bypass the conditioning guard and falsely certify a forward error.
.hdfe_norms <- function(X) {
  X <- as.matrix(X)
  ss <- colSums(X * X)
  out <- sqrt(ss)
  fallback <- !is.finite(ss) | ss < 1024 * .Machine$double.xmin
  if (any(fallback)) {
    jj <- which(fallback)
    out[jj] <- vapply(jj, function(j) .stable_norm(X[, j]), numeric(1))
  }
  if (any(!is.finite(out))) {
    stop("non-finite norm produced inside the fixed-effect solver.",
         call. = FALSE)
  }
  out
}

.hdfe_col_scale <- function(X, scale) {
  X * rep(as.numeric(scale), each = nrow(X))
}

.hdfe_project_one <- function(X, spec) {
  X <- as.matrix(X)
  sx <- rowsum(spec$a * X, spec$g, reorder = TRUE)
  mu <- sx / spec$denom
  out <- spec$a * mu[spec$g, , drop = FALSE]
  out
}

# S = sum_j P_j is self-adjoint positive semidefinite and
# <v, S v> = sum_j ||P_j v||^2.  Thus ker(S) is the joint FE-orthogonal
# space and range(S) is the sum of the FE dummy spaces.
.hdfe_project_sum <- function(X, pb) {
  out <- matrix(0, nrow(X), ncol(X))
  for (spec in pb$specs) out <- out + .hdfe_project_one(X, spec)
  if (any(!is.finite(out)))
    stop("non-finite value produced by the joint fixed-effect projector.",
         call. = FALSE)
  out
}

# For two spaces, V1 + V2 is the orthogonal direct sum
#
#   V1 + (I - P1)V2.
#
# After P1 has been removed once, A = (I-P1) P2 (I-P1) is self-adjoint PSD,
# has range (I-P1)V2, and has the same harmless coefficient non-identification
# as the full dummy system.  Solving with A therefore obtains the exact joint
# projection with three projector applications per Krylov step instead of
# repeatedly applying both full spaces until their slow angle contracts.  The
# conditioning guard and final per-FE certificate remain mandatory.
.hdfe_project_two_reduced <- function(X, pb) {
  R1X <- X - .hdfe_project_one(X, pb$specs[[1L]])
  P2R1X <- .hdfe_project_one(R1X, pb$specs[[2L]])
  out <- P2R1X - .hdfe_project_one(P2R1X, pb$specs[[1L]])
  if (any(!is.finite(out))) {
    stop("non-finite value produced by the reduced two-factor projector.",
         call. = FALSE)
  }
  out
}

.hdfe_cyclic_sweep <- function(X, pb) {
  out <- X
  for (spec in pb$specs) out <- out - .hdfe_project_one(out, spec)
  if (any(!is.finite(out)))
    stop("non-finite value produced by a fixed-effect stability sweep.",
         call. = FALSE)
  out
}

# Independent optimality diagnostics.  The group score for level g is
# rowsum(a * E)_g / sqrt(denom_g); its Euclidean norm equals the norm of the
# correction P_j E in exact arithmetic.  These are first-order/orthogonality
# residuals.  Without a certified lower spectral bound for S they are a
# backward certificate, not a universal forward-error bound; the joint CG
# solve below targets the projection itself, and an uncertified result stops.
.hdfe_metrics <- function(E, pb) {
  en <- .hdfe_norms(E)
  corr <- score <- matrix(0, length(pb$specs), ncol(E))
  for (j in seq_along(pb$specs)) {
    spec <- pb$specs[[j]]
    sums <- rowsum(spec$a * E, spec$g, reorder = TRUE)
    scores <- sums / sqrt(spec$denom)
    C <- spec$a * (sums / spec$denom)[spec$g, , drop = FALSE]
    if (any(!is.finite(scores)) || any(!is.finite(C)))
      stop("non-finite value produced while certifying fixed-effect orthogonality.",
           call. = FALSE)
    corr[j, ] <- .hdfe_norms(C)
    score[j, ] <- .hdfe_norms(scores)
  }
  absolute <- pmax(apply(corr, 2L, max), apply(score, 2L, max))
  relative <- absolute / pmax(en, .Machine$double.eps)
  list(norm = en, correction = corr, group_score = score,
       absolute = absolute, relative = relative)
}

# A candidate is accepted only after every individual FE correction and every
# weighted within-level score is small, one additional complete cyclic sweep
# is stable, and the same checks pass again after that sweep.
.hdfe_certificate <- function(E, pb, work_tol) {
  before <- .hdfe_metrics(E, pb)
  E2 <- .hdfe_cyclic_sweep(E, pb)
  change <- .hdfe_norms(E2 - E)
  after <- .hdfe_metrics(E2, pb)
  # E's input columns were independently scaled to unit norm.  `absolute` and
  # `change` are therefore already scale-free relative diagnostics.  Dividing
  # again by ||E|| would be inappropriate near an exactly absorbed column,
  # where E consists only of floating-point projection noise.
  column_ok <- before$absolute <= work_tol & change <= work_tol &
    after$absolute <= work_tol
  list(ok = all(column_ok), column_ok = column_ok, residuals = E2,
       change = max(change),
       orthogonality = max(c(before$absolute, after$absolute)),
       before = before, after = after)
}

# Matrix-free CG correction for A Q = A E, where A is either S = sum_j P_j or
# the reduced self-adjoint two-FE operator above.  Each nonzero data column and
# each nonzero right-hand side is scaled independently, so units of
# y/x/z/control columns cannot change either the path or the accept/error
# decision.  Starting at zero keeps every iterate in range(A), avoiding
# unidentified FE coefficients in disconnected or rank-deficient designs.
.hdfe_cg_correction <- function(E, pb, work_tol, budget,
                                project = .hdfe_project_sum) {
  escale <- .hdfe_norms(E)
  U <- E
  nz <- escale > 0
  if (any(nz)) U[, nz] <- .hdfe_col_scale(U[, nz, drop = FALSE],
                                          1 / escale[nz])

  B <- project(U, pb)
  bscale <- .hdfe_norms(B)
  active <- bscale > 0
  if (!any(active)) {
    return(list(correction = matrix(0, nrow(E), ncol(E)), iterations = 0L,
                converged = TRUE, breakdown = FALSE,
                linear_residual = 0, min_rayleigh = Inf))
  }

  rhs <- matrix(0, nrow(E), ncol(E))
  rhs[, active] <- .hdfe_col_scale(B[, active, drop = FALSE],
                                   1 / bscale[active])
  H <- matrix(0, nrow(E), ncol(E))
  R <- rhs
  P <- rhs
  rr <- .hdfe_norms(R)^2
  # The independent projection certificate is the acceptance gate.  Solving
  # the auxiliary system an extra decimal place beyond that gate only adds
  # rowsum passes; if this matched tolerance is insufficient, the outer
  # refinement loop solves the remaining FE component and certifies again.
  linear_tol <- max(32 * .Machine$double.eps, work_tol)
  # A merely positive Rayleigh quotient proves algebraic progress, not enough
  # forward accuracy in floating point.  For an operator with norm at most J,
  # a backward perturbation of order 128 * eps * J can be amplified by the
  # inverse Rayleigh quotient.  Requiring that amplification to stay below the
  # caller's requested tolerance gives this fail-closed floor.  It is
  # intentionally conservative: a design below it stops instead of relying on
  # a tiny normal-equation residual as a forward-error certificate.
  condition_floor <- 128 * .Machine$double.eps * length(pb$specs) /
    min(pb$tol, 1)
  min_rayleigh <- Inf
  iterations <- 0L
  breakdown <- FALSE

  while (any(active) && iterations < budget) {
    jj <- which(active)
    AP <- project(P[, jj, drop = FALSE], pb)
    iterations <- iterations + 1L
    denom <- colSums(P[, jj, drop = FALSE] * AP)
    pnorm2 <- .hdfe_norms(P[, jj, drop = FALSE])^2
    rayleigh <- denom / pnorm2
    if (any(!is.finite(denom)) || any(!is.finite(rayleigh)) ||
        any(denom <= 0) || any(rayleigh <= condition_floor)) {
      breakdown <- TRUE
      break
    }
    min_rayleigh <- min(min_rayleigh, rayleigh)

    alpha <- rr[jj] / denom
    if (any(!is.finite(alpha))) {
      breakdown <- TRUE
      break
    }
    H[, jj] <- H[, jj, drop = FALSE] +
      .hdfe_col_scale(P[, jj, drop = FALSE], alpha)
    R[, jj] <- R[, jj, drop = FALSE] - .hdfe_col_scale(AP, alpha)

    # Recurrences are cheap; periodically, and whenever they appear done,
    # replace them with the explicitly recomputed linear-system residual.
    rnew <- .hdfe_norms(R[, jj, drop = FALSE])
    apparent <- rnew <= linear_tol
    if (iterations %% 8L == 0L || any(apparent)) {
      exact <- rhs[, jj, drop = FALSE] -
        project(H[, jj, drop = FALSE], pb)
      R[, jj] <- exact
      rnew <- .hdfe_norms(exact)
    }
    rr_new <- rnew^2
    done <- rnew <= linear_tol
    if (any(done)) active[jj[done]] <- FALSE

    keep <- !done
    if (any(keep)) {
      beta <- rr_new[keep] / rr[jj[keep]]
      if (any(!is.finite(beta))) {
        breakdown <- TRUE
        break
      }
      P[, jj[keep]] <- R[, jj[keep], drop = FALSE] +
        .hdfe_col_scale(P[, jj[keep], drop = FALSE], beta)
      rr[jj[keep]] <- rr_new[keep]
    }
    if (any(done)) {
      P[, jj[done]] <- 0
      rr[jj[done]] <- rr_new[done]
    }
  }

  mult <- bscale * escale
  correction <- .hdfe_col_scale(H, mult)
  if (any(!is.finite(correction))) breakdown <- TRUE
  lin <- if (any(active)) max(sqrt(rr[active])) else 0
  list(correction = correction, iterations = iterations,
       converged = !any(active) && !breakdown, breakdown = breakdown,
       linear_residual = lin, min_rayleigh = min_rayleigh)
}

.hdfe_failure <- function(reason, iterations, maxit, change, orthogonality,
                          tol, work_tol, linear_residual = NA_real_) {
  stop(sprintf(paste0(
    "fixed-effect demeaning failed to converge/certify after %d solver ",
    "iteration(s) (maximum %d): %s; final relative full-sweep change = %.6g; ",
    "final joint orthogonality residual = %.6g; final linear-system ",
    "residual = %.6g; requested tolerance = %.6g; binding internal ",
    "certificate tolerance = %.6g. The fixed-effect ",
    "system may be ill-conditioned or near-collinear."),
    iterations, maxit, reason, change, orthogonality, linear_residual, tol,
    work_tol),
    call. = FALSE)
}

# Numerically certified joint FE projection.  The production path stays
# rowsum-only and never expands dummies.
.demean <- function(M, fe_list, weights = NULL, tol = 1e-10,
                    maxit = 10000L) {
  pb <- .hdfe_problem(M, fe_list, weights, tol, maxit)
  M <- pb$M
  scales <- .col_norms(M)
  X <- M
  nz <- scales > 0
  if (any(nz)) X[, nz] <- .hdfe_col_scale(X[, nz, drop = FALSE],
                                          1 / scales[nz])

  # The default asks for a substantially tighter internal certificate than
  # the public 1e-10 projection contract.  A smaller caller-supplied tolerance
  # remains binding.
  work_tol <- min(tol, 1e-12)

  # A witnessed principal-angle gap below either the backward-certificate
  # resolution (2 * work_tol^2) or the CG spectral floor is not numerically
  # distinguishable from an exact common FE direction.  Stop before applying
  # projectors: otherwise cancellation can erase that direction and make both
  # a normalized refinement and the orthogonality certificate falsely pass.
  # This is a one-sided conditioning witness, not a claimed global lower bound
  # for every possible multi-space dependence.
  geometry_floor <- max(2 * work_tol^2,
                        128 * .Machine$double.eps * length(pb$specs))
  loading_floor <- max(work_tol,
                       128 * .Machine$double.eps * length(pb$specs))
  if (length(pb$specs) > 1L &&
      pb$minimum_loading <= loading_floor) {
    stop(sprintf(paste0(
      "fixed-effect demeaning failed to certify before iteration: the ",
      "weighted FE geometry contains an unresolved weak incidence link ",
      "(minimum normalized within-level dummy loading = %.6g; safety ",
      "threshold = %.6g; requested tolerance = %.6g). Extreme within-level ",
      "weight imbalance can hide a material joint fixed-effect direction."),
      pb$minimum_loading, loading_floor, tol), call. = FALSE)
  }
  if (length(pb$specs) > 1L &&
      pb$minimum_pair_gap <= geometry_floor) {
    stop(sprintf(paste0(
      "fixed-effect demeaning failed to certify before iteration: the ",
      "weighted FE geometry is numerically unresolved (witnessed ",
      "principal-angle gap = %.6g; safety threshold = %.6g; ",
      "requested tolerance = %.6g). Extreme within-level weight imbalance ",
      "can hide a material joint fixed-effect direction."),
      pb$minimum_pair_gap, geometry_floor, tol), call. = FALSE)
  }
  solver_iterations <- 0L
  linear_residual <- 0
  min_rayleigh <- Inf
  remaining_correction <- Inf

  if (length(pb$specs) == 1L) {
    # The mathematical projection itself is exactly one FE sweep.  The
    # following pass is only the mandatory stability/certificate check.
    E <- X - .hdfe_project_one(X, pb$specs[[1L]])
    cert <- .hdfe_certificate(E, pb, work_tol)
    if (!cert$ok) {
      .hdfe_failure("the exact one-factor projection did not pass its numerical certificate",
                    1L, maxit, cert$change, cert$orthogonality, tol,
                    work_tol, 0)
    }
    E <- cert$residuals
    solver_iterations <- 1L
    remaining_correction <- 0
  } else {
    if (length(pb$specs) == 2L) {
      E <- X
      joint_projector <- .hdfe_project_two_reduced
      reduced_two <- TRUE
    } else {
      E <- X
      joint_projector <- .hdfe_project_sum
      reduced_two <- FALSE
    }
    repeat {
      if (reduced_two) {
        # The reduced operator acts on V1's orthogonal complement.  Restore
        # that invariant before every refinement; the final cyclic stability
        # sweep ends in the second FE space and need not preserve it exactly.
        E <- E - .hdfe_project_one(E, pb$specs[[1L]])
      }
      # Always enter the normalized joint solve, even when the raw FE scores
      # are already small.  A unit-sized FE-span component can have an
      # arbitrarily small score in a near-null direction of S; accepting it
      # from the backward certificate alone would be a false projection.
      cg <- .hdfe_cg_correction(E, pb, work_tol,
                                maxit - solver_iterations,
                                project = joint_projector)
      solver_iterations <- solver_iterations + cg$iterations
      linear_residual <- cg$linear_residual
      min_rayleigh <- min(min_rayleigh, cg$min_rayleigh)
      remaining_correction <- max(.hdfe_norms(cg$correction))
      candidate <- E - cg$correction
      if (any(!is.finite(candidate))) {
        .hdfe_failure("a non-finite joint-projection candidate was produced",
                      solver_iterations, maxit, Inf, Inf, tol,
                      work_tol, linear_residual)
      }
      cert <- .hdfe_certificate(candidate, pb, work_tol)
      E <- cert$residuals
      # A successful first solve can remove a well-conditioned FE component
      # while a near-null component is hidden below its relative RHS tolerance.
      # Do not return until a freshly normalised solve of the candidate finds
      # that the remaining *projection correction itself* is below tolerance.
      # This iterative probe exposes mixed fast/slow spectra to the Rayleigh
      # conditioning guard instead of relying on the backward certificate.
      if (cert$ok && cg$converged &&
          remaining_correction <= work_tol) break
      if (!cert$ok && cg$converged && reduced_two) {
        # The reduced solve is exact in real arithmetic, but cancellation in
        # an extreme orientation can leave a tiny component above the tighter
        # internal certificate.  Refine that component with the symmetric
        # full-space operator; it is never accepted merely from the reduced
        # solve's small change.
        joint_projector <- .hdfe_project_sum
        reduced_two <- FALSE
      }
      if (cg$breakdown && reduced_two &&
          is.finite(cg$min_rayleigh)) {
        # A later loss of conjugacy after at least one well-resolved reduced
        # direction is retried with S.  A first-direction conditioning failure
        # has min_rayleigh = Inf and still stops immediately; the full-space
        # retry also retains its own Rayleigh guard.
        joint_projector <- .hdfe_project_sum
        reduced_two <- FALSE
        next
      }
      if (cg$breakdown) {
        .hdfe_failure("the matrix-free CG solve encountered a non-positive or numerically singular direction",
                      solver_iterations, maxit, cert$change,
                      cert$orthogonality, tol, work_tol, linear_residual)
      }
      if (solver_iterations >= maxit || cg$iterations == 0L) {
        .hdfe_failure("the joint projection did not meet the final certificate",
                      solver_iterations, maxit, cert$change,
                      cert$orthogonality, tol, work_tol, linear_residual)
      }
      # Joint iterative refinement: solve for the FE component remaining in
      # the candidate, with that residual column-normalised afresh.  A small
      # S(E) is never by itself treated as convergence.
    }
  }

  out <- .hdfe_col_scale(E, scales)
  dimnames(out) <- dimnames(M)
  if (any(!is.finite(out))) {
    .hdfe_failure("rescaling the certified projection produced non-finite values",
                  solver_iterations, maxit, cert$change,
                  cert$orthogonality, tol, work_tol, linear_residual)
  }
  out
}
