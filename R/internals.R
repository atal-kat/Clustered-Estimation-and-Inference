# Internal helpers: a shared front end, the leave-out kernel, the instrument
# constructors and the inference back end. Only CJIVE is exposed in this version.

# FWL residualisation. Weights enter as a sqrt(w) transform; lm.fit is pivoted,
# so rank-deficient controls (e.g. fixed effects) are fine.
.partial_out <- function(y, x, Z, controls = NULL, weights = NULL,
                         intercept = TRUE) {
  n <- length(y)
  Z <- as.matrix(Z)
  rw <- if (is.null(weights)) rep.int(1, n) else sqrt(weights)

  C <- if (intercept) matrix(1, n, 1L) else NULL
  if (!is.null(controls)) {
    controls <- as.matrix(controls)
    C <- if (is.null(C)) controls else cbind(C, controls)
  }

  yw <- y * rw
  xw <- x * rw
  Zw <- Z * rw

  if (is.null(C)) {
    return(list(y = yw, x = xw, Z = Zw, rw = rw))
  }

  Cw <- C * rw
  yr <- yw - stats::lm.fit(Cw, yw)$fitted.values
  xr <- xw - stats::lm.fit(Cw, xw)$fitted.values
  Zr <- Zw - stats::lm.fit(Cw, Zw)$fitted.values
  list(y = yr, x = as.numeric(xr), Z = as.matrix(Zr), rw = rw)
}

# Full-sample first stage, computed once and reused. h is the diagonal of the
# hat matrix Z (Z'Z)^-1 Z' (used by JIVE).
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

# Leave-out fitted values for a partition `groups` (list of row-index sets).
# Woodbury block identity, exact and one Cholesky cheaper than refitting G times:
#   p_g = Xhat_g - H_g (I - H_g)^-1 e_g,   H_g = Z_g (Z'Z)^-1 Z_g'.
# Collapses to JIVE for singleton clusters. maxlev = max_g lambda_max(H_g).
.leaveout_fit <- function(D, Z, groups, fs = NULL) {
  if (is.null(fs)) fs <- .first_stage(D, Z)
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

# Constructed instruments. Only p-hat changes between estimators; inference does
# not.
.phat_2sls <- function(fs) fs$Xhat
.phat_jive <- function(fs, D) (fs$Xhat - fs$h * D) / (1 - fs$h)
.phat_cjive <- function(fs, D, Z, groups) .leaveout_fit(D, Z, groups, fs)$phat

# Just-identified cluster-robust IV sandwich on the residualised data:
#   beta = (p'Y)/(p'D),  SE = sqrt(sum_g S_g^2 * G/(G-1)) / |p'D|,  S_g = sum_g p_i e_i.
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

# A grouping (factor/character) z becomes a dummy design, one level dropped;
# numeric z is used as is.
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

.assert_finite <- function(v, name) {
  if (anyNA(v) || any(!is.finite(v)))
    stop(sprintf("`%s` contains missing or non-finite values; remove or impute them first.",
                 name), call. = FALSE)
}

# Validate inputs, expand z, partial out, build the cluster partition. Shared by
# cjive() and iv_compare().
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

  # A grouping z needs every group spread over at least two clusters.
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

# Leave-cluster-out group mean, FLM's closed form for the pure judge design:
# p_i = (S_j - S_{j,g}) / (N_j - N_{j,g}), the weighted mean of x over group j
# outside i's cluster. Differs from the dense route by an O(1/n_g) term.
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
