# Internal machinery for the clustered-jackknife family.
#
# The package is organised around the seam described in the design notes: a
# shared front end (FWL + weighting), a leave-out kernel that takes the
# leave-out assignment as an argument, a set of thin instrument constructors,
# and a single inference back end.  CJIVE is the only routine wired up in v1;
# the factoring is what lets siblings (a cluster-jackknife AR test, a multi-way
# jackknife) slot in without reworking the core.

# ---------------------------------------------------------------------------
# Front end: Frisch-Waugh-Lovell residualisation in the weighted geometry.
#
# Optional precision weights w enter as a sqrt(w) left transform on every
# variable (outcome, endogenous regressor, instruments, covariates); the whole
# estimator then operates on one consistent weighted geometry and centring in
# the covariance form is the weighted mean.  An intercept is partialled out
# unless `intercept = FALSE`; the covariate design may be rank deficient
# (FE/dummy controls are allowed), so we residualise with a pivoted least
# squares fit rather than forming a normal-equations inverse.
# ---------------------------------------------------------------------------
.partial_out <- function(y, x, Z, controls = NULL, weights = NULL,
                         intercept = TRUE) {
  n <- length(y)
  Z <- as.matrix(Z)
  rw <- if (is.null(weights)) rep.int(1, n) else sqrt(weights)

  # Covariate design to partial out: intercept (optional) plus any controls.
  C <- if (intercept) matrix(1, n, 1L) else NULL
  if (!is.null(controls)) {
    controls <- as.matrix(controls)
    C <- if (is.null(C)) controls else cbind(C, controls)
  }

  # Weighted geometry.
  yw <- y * rw
  xw <- x * rw
  Zw <- Z * rw

  if (is.null(C)) {
    return(list(y = yw, x = xw, Z = Zw, rw = rw))
  }

  Cw <- C * rw
  # Pivoted LS tolerates rank deficiency; fitted values are the projection
  # regardless of which aliased columns are dropped.
  yr <- yw - stats::lm.fit(Cw, yw)$fitted.values
  xr <- xw - stats::lm.fit(Cw, xw)$fitted.values
  Zr <- Zw - stats::lm.fit(Cw, Zw)$fitted.values
  list(y = yr, x = as.numeric(xr), Z = as.matrix(Zr), rw = rw)
}

# ---------------------------------------------------------------------------
# Full-sample first stage, computed once and reused by every constructor.
#   Q  = Z'Z,  Qinv via Cholesky;  pihat = Qinv Z'D;  Xhat = Z pihat;
#   e  = D - Xhat;  ZQinv = Z Qinv (sliced per cluster later);
#   h  = diag(Z Qinv Z') the observation leverages (for JIVE).
# ---------------------------------------------------------------------------
.first_stage <- function(D, Z) {
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
       e = D - Xhat, ZQinv = ZQinv, h = rowSums(ZQinv * Z))
}

# ---------------------------------------------------------------------------
# Leave-out kernel (the required optimisation).
#
# For a leave-out assignment given as a list of row-index sets, return the
# leave-set-out first-stage fitted values.  For CJIVE the assignment is the
# cluster partition, so the sets are disjoint and the Woodbury block identity
#
#     p_g = Xhat_g - H_g (I - H_g)^{-1} e_g,     H_g = Z_g Q^{-1} Z_g',
#
# is exact: it reproduces, for every cluster, the fit that re-estimates
# (sum_{l not in g} Z_l Z_l')^{-1} sum_{l not in g} Z_l D_l and applies it to
# Z_g -- but at the cost of one Cholesky of Q plus a small n_g x n_g solve per
# cluster, never G refactorisations of Q.  For singleton clusters the block
# collapses to the scalar JIVE update automatically.
#
# The kernel takes the assignment as an argument and does not assume "exactly
# one partition" beyond using the disjoint-block update, so a covering can be
# supplied by a future multi-way sibling without changing this interface.
# `maxlev` = max_g lambda_max(H_g) is a conditioning diagnostic.
# ---------------------------------------------------------------------------
.leaveout_fit <- function(D, Z, groups, fs = NULL) {
  if (is.null(fs)) fs <- .first_stage(D, Z)
  phat <- fs$Xhat            # overwritten on the rows of each leave-out set
  maxlev <- 0

  for (idx in groups) {
    ng <- length(idx)
    Zg <- Z[idx, , drop = FALSE]
    Hg <- fs$ZQinv[idx, , drop = FALSE] %*% t(Zg)   # n_g x n_g leverage block
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

# ---------------------------------------------------------------------------
# Instrument constructors: each returns a constructed instrument p-hat on the
# (already residualised) data.  Only the constructed instrument distinguishes
# the estimators; inference is identical.
# ---------------------------------------------------------------------------
.phat_2sls <- function(fs) fs$Xhat
.phat_jive <- function(fs, D) (fs$Xhat - fs$h * D) / (1 - fs$h)
.phat_cjive <- function(fs, D, Z, groups) .leaveout_fit(D, Z, groups, fs)$phat

# ---------------------------------------------------------------------------
# Inference back end: the just-identified cluster-robust IV sandwich.
#
#   s_i = p_i (Y_i - beta D_i),  S_g = sum_{i in g} s_i,  den = sum_i p_i D_i,
#   SE  = sqrt( sum_g S_g^2 * G/(G-1) ) / |den|.
#
# Operates on the residualised quantities; the point estimate is the
# covariance-ratio (p'Y)/(p'D), automatic once covariates (incl. intercept)
# are partialled out.  Kept separate from the constructors so a test /
# confidence-set back end (CJAR) can later reuse the same residualised inputs.
# ---------------------------------------------------------------------------
.iv_inference <- function(phat, x, y, cluster, level = 0.95) {
  cl <- droplevels(as.factor(cluster))
  G <- nlevels(cl)
  den <- sum(phat * x)
  beta <- sum(phat * y) / den

  s <- phat * (y - beta * x)
  Sg <- tapply(s, cl, sum)
  Sg[is.na(Sg)] <- 0
  se <- sqrt(sum(Sg^2) * G / (G - 1)) / abs(den)

  stat <- beta / se
  zc <- stats::qnorm(1 - (1 - level) / 2)
  list(coefficient = beta, se = se, statistic = stat,
       p.value = 2 * stats::pnorm(-abs(stat)),
       conf.low = beta - zc * se, conf.high = beta + zc * se, G = G)
}

# ---------------------------------------------------------------------------
# Expand the instrument argument.  A factor/character grouping vector (a judge
# design) is expanded to a dummy design with one reference level dropped, the
# intercept supplying the rest; a numeric vector/matrix is taken as is.
# ---------------------------------------------------------------------------
.expand_z <- function(z) {
  if (is.factor(z) || is.character(z)) {
    f <- as.factor(z)
    Z <- stats::model.matrix(~f)[, -1L, drop = FALSE]
    return(list(Z = Z, grouping = TRUE, group = f))
  }
  Z <- as.matrix(z)
  storage.mode(Z) <- "double"
  list(Z = Z, grouping = FALSE, group = NULL)
}

# Reject missing or non-finite values up front, with a message that names the
# offending argument -- otherwise NAs surface as a cryptic LS failure, or (for
# the cluster id) are silently dropped from the partition while still counted
# in n, returning a quietly wrong estimate.
.assert_finite <- function(v, name) {
  if (anyNA(v) || any(!is.finite(v)))
    stop(sprintf("`%s` contains missing or non-finite values; remove or impute them first.",
                 name), call. = FALSE)
}

# Shared preparation: validate, expand z, partial out, build the cluster
# partition.  Used by both cjive() and iv_compare() so they cannot diverge.
.prep_data <- function(y, x, z, cluster, controls, weights, intercept) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  n <- length(y)

  if (length(x) != n) stop("`x` and `y` must have the same length.", call. = FALSE)
  if (length(cluster) != n) stop("`cluster` must have length n.", call. = FALSE)
  .assert_finite(y, "y")
  .assert_finite(x, "x")
  if (anyNA(cluster)) stop("`cluster` contains missing values.", call. = FALSE)
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != n) stop("`weights` must have length n.", call. = FALSE)
    if (any(!is.finite(weights)) || any(weights <= 0))
      stop("`weights` must be finite and strictly positive.", call. = FALSE)
  }

  if (anyNA(z)) stop("`z` contains missing values.", call. = FALSE)
  zinfo <- .expand_z(z)
  if (nrow(zinfo$Z) != n) stop("`z` must have n rows.", call. = FALSE)
  .assert_finite(zinfo$Z, "z")

  # Covariates: a data frame (possibly with factor columns) is expanded to a
  # numeric design via model.matrix; a matrix is taken as is.  Rank deficiency
  # is tolerated downstream by the pivoted residualisation.
  if (!is.null(controls)) {
    if (NROW(controls) != n) stop("`controls` must have n rows.", call. = FALSE)
    if (anyNA(controls)) stop("`controls` contains missing values.", call. = FALSE)
    controls <- if (is.data.frame(controls))
      stats::model.matrix(~., data = controls)[, -1L, drop = FALSE]
    else as.matrix(controls)
    .assert_finite(controls, "controls")
  }

  cl <- droplevels(as.factor(cluster))
  G <- nlevels(cl)
  if (G < 2L) stop("at least 2 clusters are required.", call. = FALSE)

  # FLM require several clusters per instrument group; for a grouping z we can
  # check this directly (a judge confined to one cluster has no leave-out fit).
  if (zinfo$grouping) {
    spread <- rowSums(table(zinfo$group, cl) > 0L)
    if (any(spread < 2L)) {
      bad <- names(spread)[spread < 2L]
      shown <- bad[seq_len(min(5L, length(bad)))]
      stop("each instrument group must appear in at least 2 clusters; ",
           "offending group(s): ", paste(shown, collapse = ", "),
           if (length(bad) > 5L) ", ..." else "", ".", call. = FALSE)
    }
  }

  po <- .partial_out(y, x, zinfo$Z, controls, weights, intercept)
  list(y = po$y, x = po$x, Z = po$Z, cluster = cl,
       groups = split(seq_len(n), cl), n = n, G = G, p = ncol(zinfo$Z),
       grouping = zinfo$grouping, group = zinfo$group, rw = po$rw,
       weights = weights)
}

# ---------------------------------------------------------------------------
# Leave-cluster-out group mean: FLM's closed form for the pure judge design.
# For observation i in instrument group j and cluster g,
#   p_i = (S_j - S_{j,g}) / (N_j - N_{j,g}),
# the (weighted) mean of x over members of j outside i's cluster.  This removes
# the cluster from the group *level*, so it differs from the dense/FWL route by
# the small O(1/n_g) intercept term.  O(n), never selected automatically.
# ---------------------------------------------------------------------------
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
