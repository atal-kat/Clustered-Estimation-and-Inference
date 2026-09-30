# Effective first-stage F of Montiel Olea & Pflueger (2013, JBES 31(3)),
# clustered variant, with the simplified TSLS critical value.  The formulas
# below follow equations 3, 16, 32, and 33 of the paper.
#
# The paper's Step 1 requires the exogenous regressors partialled out of
# every variable (already done by .partial_out) and the instruments
# orthonormalised, Z*'Z*/S = I_K. The first-stage whitened instruments satisfy
# Ztil'Ztil = I_k, so Z* = sqrt(n) * Ztil and both sqrt(n) factors cancel in
# the statistic: everything reduces to objects the package already holds.
# Instrument-basis invariance depends on this orthonormalisation.

# Patnaik critical value (eq. 32): the upper-alpha quantile of
# chisq_{K_eff}(x * K_eff) / K_eff, with x = 1/tau. K_eff may be fractional;
# base R's qchisq accepts non-integer df and an ncp. Reject the null of weak
# instruments (Nagar bias > tau of the worst-case benchmark) when
# F_eff > crit.
.eff_f_crit <- function(K_eff, x = 10, alpha = 0.05) {
  stats::qchisq(1 - alpha, df = K_eff, ncp = x * K_eff) / K_eff
}

# The statistic and its simplified-TSLS critical value (Def. 2, eq. 33), on
# the whitened path: t_x = Ztil'x (the first-stage fit coefficients in
# whitened coordinates, `t` in .first_stage) and vhat = x - Ztil t_x (the
# first-stage residuals, `e` in .first_stage). Then
#   x' P_Z x   = sum(t_x^2),
#   What_2     = crossprod(rowsum(Ztil * vhat, cluster)),  # clustered meat
#   F_eff      = sum(t_x^2) / tr(What_2),
# at O(nk + Gk * min(G, k)) plus an eigenproblem on the smaller of the G x G
# and k x k Gram matrices, no O(n_g^2) object and no second pass over Z. What_2 carries
# no small-sample correction (faithful to the paper; Stata's weakivtest
# multiplies in the usual cluster finite-sample factor, so its statistic
# differs by that known scalar). K_eff is scale-invariant in What_2, so a
# correction there could not move the critical value anyway.
#
# tau = 0.10 and alpha = 0.05 are fixed (the paper's baseline and
# weakivtest's default); the simplified procedure sets x = 1/tau = 10.
# A nonzero numerator and an exactly zero clustered-moment denominator imply
# F_eff = Inf. K_eff and its critical value are then undefined; a zero numerator
# and denominator makes F_eff itself undefined. No tolerance-classified
# "exact fit" shortcut is used: a small but nonzero denominator produces the
# corresponding large finite statistic.
.eff_f <- function(Ztil, t_x, vhat, cluster) {
  Sg <- rowsum(Ztil * as.numeric(vhat), cluster)     # G x k
  su <- .unit_norm(as.numeric(Sg))
  tu <- .unit_norm(as.numeric(t_x))
  if (!(su$max > 0)) {
    return(list(F_eff = if (tu$max > 0) Inf else NA_real_,
                K_eff = NA_real_, F_eff_crit = NA_real_))
  }

  # Normalize before every square: this removes both the quadratic underflow
  # in tr(W2) and the quartic overflow in tr(W2^2). The nonzero eigenvalues of
  # Sg'Sg and Sg Sg' coincide, so form whichever Gram matrix is smaller.
  Sgn <- matrix(su$unit, nrow = nrow(Sg), ncol = ncol(Sg))
  W2n <- if (ncol(Sgn) <= nrow(Sgn)) crossprod(Sgn) else tcrossprod(Sgn)
  lmaxn <- max(eigen(W2n, symmetric = TRUE, only.values = TRUE)$values)
  trW2sq <- sum(W2n * W2n)
  x <- 10                                            # 1 / tau, tau = 0.10
  K_eff <- (1 + 2 * x) / (trW2sq + 2 * x * lmaxn)
  F_eff <- if (tu$max > 0) {
    .rescale_by_norms(1, tu, su, power = 2)
  } else 0
  list(F_eff = F_eff, K_eff = K_eff,
       F_eff_crit = .eff_f_crit(K_eff, x = x, alpha = 0.05))
}

# The NA triple for paths that carry no whitened first stage (the
# leaveout_mean path of cjive).
.eff_f_na <- function() {
  list(F_eff = NA_real_, K_eff = NA_real_, F_eff_crit = NA_real_)
}
