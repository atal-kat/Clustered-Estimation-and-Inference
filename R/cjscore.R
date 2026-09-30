#' Cluster jackknife score test (CJS)
#'
#' Tests point hypotheses on the structural coefficient of a single endogenous
#' regressor using the cluster jackknife score (LM) statistic of Ligtenberg
#' (2025, Section 4), with the confidence set obtained by analytic,
#' grid-free test inversion by sign-classifying the quadratic rejection
#' boundary and quadratic variance polynomial; roots and endpoints are
#' floating-point quantities, and coverage is the test's asymptotic coverage.
#' Like \code{\link{cjar}} the
#' test is cluster-robust and remains valid under weak and under many
#' instruments.  The paper finds complementary power: CJS can perform better
#' near the tested value and under high heteroskedasticity, while CJAR can
#' perform better with weak instruments.  For reporting confidence sets,
#' prefer \code{\link{cjar}}: see
#' \sQuote{Score confidence sets and spurious regions} below.
#'
#' @param y Outcome (numeric vector).  Factor and character outcomes are
#'   rejected.  The generic dispatches on this argument; a formula as the
#'   first argument selects the formula method.
#' @param x Single endogenous regressor (numeric vector).  Factor and
#'   character regressors are rejected; for the print method, a fitted
#'   \code{"cjscore"} object.
#' @param z Instruments: a numeric vector/matrix, or a factor/character
#'   grouping vector for a judge design.  A grouping instrument uses reference
#'   coding whenever an intercept or any fixed effect is partialled out.  It
#'   uses one column for every level only when \code{intercept = FALSE} and no
#'   fixed effects are supplied.
#'   Identical to the \code{z} argument of \code{\link{cjar}} and
#'   \code{\link{cjive}}.
#' @param cluster Cluster identifiers (length n).  In the formula method all
#'   four forms work: a bare column name, a one-sided formula (\code{~ g}), a
#'   column name as a string, or a full vector (see \code{\link{cjive}}).  The
#'   test allows arbitrary dependence within clusters; validity rests on
#'   independence across them.  A cluster formula must contain exactly one
#'   bare variable; multiway formulas and arithmetic combinations of cluster
#'   codes are not supported.
#' @param formula A model formula in either restricted layout documented in
#'   \code{\link{cjive}}: \code{y ~ x | z} (optionally
#'   \code{y ~ x | z | fe}) or \code{y ~ exog | fe | endo ~ inst}.  Formula
#'   sections are additive lists of bare variable names, with exactly one
#'   outcome and one endogenous regressor; arithmetic, transformed terms and
#'   subtraction are rejected except for the \code{0}/\code{-1} intercept
#'   markers.
#' @param data A data frame in which to evaluate the formula.
#' @param subset Optional expression selecting the rows to use, evaluated in
#'   \code{data} (formula method only), as in \code{\link[stats]{lm}}.
#' @param na.action How to treat missing values (formula method only); see
#'   \code{\link{cjive}}.  The default \code{na.omit} drops incomplete rows
#'   in a single complete-cases filter and records the count as
#'   \code{n_dropped}; the vector interface requires clean input.
#' @param controls Optional exogenous covariates: a matrix or data frame, or a
#'   one-sided formula in the formula method.  Partialled out of \code{y},
#'   \code{x} and every instrument by Frisch-Waugh-Lovell before the test, the
#'   package-wide convention; the covariate regimes and caveats documented in
#'   \code{\link{cjar}} (\sQuote{Controls and fixed effects}) apply unchanged.
#' @param fixed_effects Optional high-dimensional fixed effects to absorb: a
#'   factor, or a list/data frame of factors (one per dimension).  In the
#'   formula method a third \code{|}-separated part
#'   (\code{y ~ x | z | fe1 + fe2}) or a one-sided formula is also accepted.
#'   Absorbed by the same guarded numerical projection as in
#'   \code{\link{cjar}} and \code{\link{cjive}}.
#' @param weights Optional finite, strictly positive numeric precision weights;
#'   factor and character weights are rejected.  \code{y},
#'   \code{x} and the instruments are transformed by \code{sqrt(weights)}
#'   before partialling out, so the result is the CJS test of the transformed
#'   (weighted) model.  In the formula method the same four forms as
#'   \code{cluster} work, including a bare data-frame column name.
#' @param beta0 One finite numeric null value of the structural coefficient at
#'   which the statistic and p-value are reported (default 0, the no-effect
#'   null).  The confidence set is always computed regardless of \code{beta0}.
#'   Testing a
#'   specific point null is this function's primary use case.
#' @param level One finite numeric confidence level strictly between 0 and 1
#'   (default 0.95), as in \code{\link{cjar}}.
#' @param intercept One non-missing logical value; partial out an intercept
#'   (default \code{TRUE}).  Disable only when the inputs are already
#'   residualised.
#' @param variance Variance estimator studentising the score.  \code{"plain"}
#'   (the default) is the feasible estimator of Ligtenberg (2025, eq. 6),
#'   the unchanged default path; \code{"crossfit"} is the cross-fit
#'   estimator of Section 5.4 (eq. 7), unbiased under the paper's assumptions
#'   at every \eqn{\beta} rather than only at the truth.  The score
#'   coefficients are untouched; only the
#'   variance polynomial changes, and the conservative non-positive-variance
#'   convention below applies unchanged (the cross-fit estimate, like the
#'   plain one, carries no sign guarantee).  Costs one leave-out solve per
#'   cluster pair, dispatched on the stacked size of the two left-out
#'   clusters (\eqn{O(m^2 k + m^3)} for \eqn{m = n_g + n_h < k},
#'   \eqn{O(k^3)} otherwise); see \sQuote{Cross-fit variance} in
#'   \code{\link{cjar}} for the shared construction.
#' @param ... Must be empty.  Unknown, misspelled and unnamed arguments are
#'   errors.
#'
#' @return An object of class \code{"cjscore"}: a list with
#'   \describe{
#'     \item{\code{statistic}, \code{p.value}}{the score statistic
#'       \eqn{LM(\beta_0) = S(\beta_0)^2 / \hat V^S(\beta_0)} and its p-value
#'       under the asymptotic \eqn{\chi^2_1} calibration.  When
#'       \eqn{\hat V^S(\beta_0) \le 0} the statistic is not defined: it is
#'       reported as \code{NA} and the p-value is 1, so \eqn{\beta_0} is
#'       accepted by the conservative convention (see \sQuote{Non-positive
#'       variance estimates}).}
#'     \item{\code{score}, \code{variance}}{the paper-scaled ingredients
#'       \eqn{S(\beta_0) = X'\ddot P_Z\,\varepsilon(\beta_0)/\sqrt{n}} and
#'       \eqn{\hat V^S(\beta_0)}; the variance may be non-positive in finite
#'       samples (its cross-cluster term is not a sum of squares).}
#'     \item{\code{term}, \code{beta0}, \code{level}, \code{variance_estimator},
#'       \code{crit}}{the endogenous-regressor label, null value, confidence
#'       level, variance estimator (\code{"plain"} or
#'       \code{"crossfit"}), and the critical
#'       value \eqn{\chi^2_{1,1-\alpha}} = \code{qchisq(level, 1)}.}
#'     \item{\code{conf_set}}{an m x 2 matrix of interval endpoints (0 rows
#'       when empty; \code{-Inf}/\code{Inf} allowed in the first/last row; a
#'       row with equal endpoints is an isolated accepted point).  Because
#'       points where the variance estimate is non-positive are accepted, the
#'       set can contain more than one component: the usual LM region plus
#'       regions on which \eqn{\hat V^S \le 0} (see \sQuote{Non-positive
#'       variance estimates}).  Under the plain variance such a region is
#'       bounded; under cross-fit it may extend into a tail.}
#'     \item{\code{shape}}{one of \code{"bounded"}, \code{"ray"},
#'       \code{"two_rays"}, \code{"whole_line"}, \code{"empty"}, sharing
#'       the \code{\link{cjar}} tail-classification vocabulary.  It is not a
#'       complete topology when bounded middle components coexist with tails.}
#'     \item{\code{n_components}, \code{unbounded_left},
#'       \code{unbounded_right}}{the complete number-of-components and tail
#'       topology read directly from \code{conf_set}.}
#'     \item{\code{bounded}}{logical; \code{TRUE} when every endpoint is
#'       finite (including the empty set).}
#'     \item{\code{F_CJS}}{the score-form first-stage strength statistic
#'       (\code{NA} if its variance term is zero); see \sQuote{Instrument
#'       strength}.}
#'     \item{\code{n}, \code{G}, \code{k}}{observations, clusters, and number
#'       of instrument columns after partialling.}
#'     \item{\code{maxlev}}{\eqn{\max_g \|P_{Z,[g,g]}\|_2}, the same
#'       cluster-leverage diagnostic reported by \code{\link{cjar}}.}
#'     \item{\code{F_eff}, \code{K_eff}, \code{F_eff_crit},
#'       \code{ng_max}}{the Montiel Olea-Pflueger effective first-stage F
#'       and its simplified-TSLS critical value, and the largest cluster
#'       size, exactly as in \code{\link{cjar}} (whose documentation carries
#'       the details and the \sQuote{Where is the first-stage F?} FAQ).  If
#'       the estimated denominator is exactly zero with nonzero signal,
#'       \code{F_eff = Inf} and both \code{K_eff} and \code{F_eff_crit} are
#'       \code{NA}.}
#'     \item{\code{fe_dims}, \code{fe_levels}, \code{k_controls},
#'       \code{n_dropped}}{as in \code{\link{cjive}}.}
#'     \item{\code{coef_score}, \code{coef_var}}{the coefficients
#'       \eqn{(s_0, s_1)} of \eqn{\sqrt{n}\,S(\beta) = s_0 - s_1\beta} and
#'       \eqn{(v_0, v_1, v_2)} of \eqn{n\hat V^S(\beta) = v_0 + v_1\beta +
#'       v_2\beta^2}, on the original beta scale; they make the object
#'       self-contained (\code{confint} re-inverts at a new level without
#'       refitting).}
#'     \item{\code{call}}{the matched call.}
#'   }
#'
#' @details
#' \subsection{Validated design}{
#' Weighting and dense or fixed-effect residualisation happen before the
#' identifying checks.  The transformed data must be finite; \code{y},
#' \code{x} and every instrument column must retain usable variation; and the
#' residualised instrument matrix must have full numerical column rank.
#' Scale- and dimension-aware checks report absorbed or collinear instruments
#' rather than silently dropping them.  See \code{\link{cjive}} for the shared
#' input and restricted-formula contract.
#' }
#'
#' \subsection{The test}{
#' With \eqn{\varepsilon(\beta) = y - x\beta} on the partialled data, the
#' statistic is \eqn{S(\beta) = X'\ddot P_Z\,\varepsilon(\beta)/\sqrt{n}}
#' (Ligtenberg 2025, Section 4), where \eqn{\ddot P_Z} is the
#' instrument projection with its diagonal cluster blocks removed -- the same
#' cluster jackknife as \code{\link{cjar}}, applied to the score
#' \eqn{X'P_Z\varepsilon} instead of the quadratic form
#' \eqn{\varepsilon'P_Z\varepsilon}.  It is studentised by the feasible
#' variance estimator of eq. (6) (conditionally unbiased and consistent by
#' the paper's Theorem 3, under its assumptions) and compared, two-sided,
#' against the asymptotic normal limit given
#' by the paper's Theorem 2: \eqn{H_0\colon\beta = \beta_0} is rejected when
#' \eqn{LM(\beta_0) > \chi^2_{1,1-\alpha}}.  The chi-square(1) reference is
#' asymptotic under the paper's assumptions, not an exact finite-sample
#' law.  The paper's shifted-and-scaled
#' \eqn{\chi^2_k} critical value is specific to the AR statistic (whose
#' reference distribution depends on \eqn{k}) and is deliberately not carried
#' over; a \code{calibration} argument is not offered because the two
#' candidate calibrations -- two-sided normal on \eqn{S/\sqrt{\hat V}} and
#' \eqn{\chi^2_1} on \eqn{LM} -- are the same test.
#' }
#'
#' \subsection{When to reach for the score test}{
#' The score test offers a different power profile while retaining weak-ID
#' robustness (Ligtenberg 2025, Section 4, after Kleibergen 2002).  In the
#' paper's simulations its advantage concentrates near the tested value with
#' a moderate number of strong instruments and under strong
#' heteroskedasticity; under weak instruments the cluster jackknife AR has
#' better power far from the null.  Its
#' appendix also finds the cluster jackknife score notably robust to a
#' dominating cluster, where the cluster jackknife AR over-rejects.
#' }
#'
#' \subsection{Score confidence sets and spurious regions}{
#' Inverting a score test can produce confidence sets containing regions far
#' from any plausible parameter value: a score statistic vanishes wherever
#' the underlying objective is stationary, not only near the truth -- here,
#' \eqn{S(\beta)} is linear in \eqn{\beta} and crosses zero at
#' \eqn{\beta = s_0/s_1}, a point the test always accepts, however the data
#' were generated (where the variance estimate is positive the statistic is
#' zero there; where it is non-positive the point is accepted by the
#' conservative convention).  In addition, the conservative convention can
#' add a closed interval on which the variance estimate is non-positive,
#' possibly disjoint from the main region (see \sQuote{Non-positive variance
#' estimates}).  The set reported by \code{cjscore()} has asymptotic
#' coverage under the paper's assumptions (no unconditional finite-sample
#' coverage guarantee is claimed), but for interval reporting prefer
#' \code{\link{cjar}}, whose statistic diverges along any fixed alternative
#' and whose plain variance estimator is a sum of squares (its cross-fit
#' variance, like the cross-fit CJS variance, has no pointwise sign guarantee);
#' use \code{cjscore()} primarily to test pre-specified point nulls, where its
#' complementary power properties are relevant.
#' }
#'
#' \subsection{Non-positive variance estimates}{
#' Under \code{variance = "plain"}, unlike the plain AR variance (a sum of
#' squares), the score variance estimator's
#' cross-cluster term \eqn{\sum_{g\neq h} c_{gh} c_{hg}} can be negative, so
#' \eqn{\hat V^S(\beta) \le 0} is possible in finite samples (it requires at
#' least three clusters; the estimator is correct on average and consistent
#' at the true \eqn{\beta_0} by the paper's Theorem 3, so a non-positive
#' value at a plausible \eqn{\beta} is itself a red flag).  Automatic
#' acceptance at such points is the \emph{package's} conservative
#' convention -- the cited papers do not prove or prescribe this rule.  The
#' convention, applied identically in the statistic and the inversion, is:
#' \eqn{\beta} is rejected if and only if \eqn{\hat V^S(\beta) > 0}
#' \emph{and} \eqn{S(\beta)^2 > \chi^2_{1,1-\alpha} \hat V^S(\beta)}.  Where
#' the variance estimate is positive this is the usual LM test; where it is
#' non-positive the statistic is not defined and the point is accepted --
#' the test declines to reject.  For the plain estimator, the leading variance
#' coefficient satisfies \eqn{v_2 \ge 0}, so negativity is confined to a
#' bounded beta-window and does not change the tails.  The cross-fit variance
#' has no corresponding sign guarantee: a non-positive region can extend into
#' a tail.  The same conservative acceptance convention is applied, but the
#' plain-estimator tail statement and boundedness formula below do not carry
#' over without qualification.
#' }
#'
#' \subsection{Instrument strength: F_CJS}{
#' For \code{variance = "plain"}, in exact arithmetic the confidence set is
#' bounded if and only if
#' \eqn{v_2 > 0} and
#' \code{F_CJS^2 > crit} (both tails rejected), with
#' \eqn{F_{CJS} = s_1/\sqrt{v_2}} -- the algebraic leading-coefficient
#' boundedness criterion.  Algebraically, \eqn{F_{CJS}^2} is the CJS
#' statistic applied to the first stage (testing instrument irrelevance,
#' \eqn{\Pi = 0}, with \eqn{x} in the role of the outcome) -- the score
#' analogue of \code{\link{cjar}}'s \code{F_CJ}, with the same numerator
#' \eqn{x'\ddot P_Z x} and the score-form variance in the denominator.  Under
#' cross-fit, \code{F_CJS} is recomputed from the selected variance leading
#' coefficient.  It is \code{NA} when that coefficient is non-positive; when
#' finite, its comparison with \code{crit} describes the selected procedure's
#' two tails rather than a variance-invariant first-stage quantity.  Both
#' forms are free by-products of the stored coefficients.
#' }
#'
#' \subsection{Relation to cjar and cjive}{
#' \code{cjscore()} shares \code{\link{cjar}}'s entire infrastructure: the
#' partialling conventions, the instrument Gram factorisation, the
#' \code{maxlev} diagnostic and the advisory thresholds, so results from the
#' two tests and from \code{\link{cjive}} refer to the same design.  The
#' recommended workflow remains: \code{\link{cjive}} for the point estimate,
#' \code{\link{cjar}} for the reported confidence set, \code{cjscore()} for
#' complementary tests of specific point nulls (and as a cross-check: under
#' strong identification all three agree closely).
#' }
#'
#' \subsection{FAQ}{
#' \describe{
#'   \item{My data aren't clustered -- can I still use this?}{Yes: set
#'     \code{cluster = seq_len(n)}, every cluster a singleton.  \eqn{\ddot
#'     P_Z} is then the projection matrix with a zeroed diagonal and
#'     \code{cjscore()} computes the jackknife Lagrange multiplier test of
#'     Matsushita and Otsu (2024): statistic, variance estimator and
#'     chi-square(1) calibration coincide algebraically (their eqs. 3--4).
#'     The print
#'     method labels this case \dQuote{independent data; Matsushita-Otsu
#'     2024}.  The identity holds for the default \code{variance = "plain"}
#'     only -- \code{variance = "crossfit"} substitutes Ligtenberg's (2025,
#'     Section 5.4) leave-out variance construction, which is not the
#'     Matsushita-Otsu estimator, and the print label says so -- and it is
#'     an identity for the no-control (or previously partialled) design.
#'     Two package conventions sit on top of their test: the
#'     conservative acceptance at non-positive variance estimates
#'     (\sQuote{Non-positive variance estimates} above), which their paper
#'     does not discuss, and the ex-ante global Frisch-Waugh-Lovell
#'     partialling of supplied controls, which is not claimed to equal the
#'     included-regressor correction analysed by Matsushita and Otsu -- with
#'     many controls see the caveat in \code{\link{cjar}} (\sQuote{Controls
#'     and fixed effects}).  The validity
#'     citation for clustered data remains Ligtenberg (2025, Theorem 2),
#'     which nests this test as the singleton case.}
#' }
#' }
#'
#' \subsection{Computation}{
#' The coefficients of the linear statistic and quadratic variance are
#' computed once from whitened per-cluster instrument sums in
#' \eqn{O(nk + Gk^2)} time and \eqn{O(Gk + k^2)} memory -- the same kernel
#' shape as \code{\link{cjar}}, one Gram matrix cheaper.  The confidence set
#' is the sign classification of two quadratics (the acceptance polynomial
#' and the variance polynomial); no grid search anywhere.  The \code{maxlev}
#' diagnostic runs in the whitened leave-cluster-out kernel shared with
#' \code{\link{cjive}} and \code{\link{cjar}}, at
#' \eqn{O(nk^2 + \sum_g \min(n_g, k)^3)}.
#' The implementation is checked against a direct dense calculation with the
#' diagonal cluster blocks of the instrument projector explicitly removed.
#' }
#'
#' @references
#' Ligtenberg, J. W. (2025). Inference in clustered IV models with many and
#' weak instruments. arXiv:2306.08559v3.  Section 4 for the statistic,
#' Assumption 3 and Theorem 2; eq. (6) for the feasible variance
#' estimator and Theorem 3 for its unbiasedness and consistency; Section 5.4
#' (eq. 7) for the cross-fit variance estimator behind
#' \code{variance = "crossfit"}; Section 6 and Appendix F for the power and
#' dominant-cluster evidence cited above.
#'
#' Kleibergen, F. (2002). Pivotal statistics for testing structural
#' parameters in instrumental variables regression. \emph{Econometrica},
#' 70(5), 1781--1803.  The score test's ancestor for fixed k.
#'
#' Matsushita, Y. and Otsu, T. (2024). A jackknife Lagrange multiplier test
#' with many weak instruments. \emph{Econometric Theory}, 40(2), 447--470.
#' The independent-data jackknife LM test (eqs. 3--4 for the statistic and
#' its variance estimator), compared against \eqn{\chi^2} critical values --
#' the calibration convention adopted here.  At singleton clusters with
#' \code{variance = "plain"} and no supplied controls, \code{cjscore()}
#' reduces to this test algebraically (see the FAQ).
#'
#' Dufour, J.-M. (1997). Some impossibility theorems in econometrics with
#' applications to structural and dynamic models. \emph{Econometrica}, 65(6),
#' 1365--1387.  Why weak-instrument-valid confidence sets must sometimes be
#' unbounded.
#'
#' @seealso \code{\link{iv_infer}} for the recommended one-call workflow
#'   (CJIVE estimate + CJAR confidence set + CJS point-null test with
#'   shared preprocessing); \code{\link{cjar}} for the companion AR test and the recommended
#'   confidence sets; \code{\link{cjive}} for the point estimate;
#'   \code{\link{iv_compare}} for estimator comparison.  Two precomputed
#'   worked examples on real data:
#'   \code{vignette("queens-workflow", package = "clusterIV")} and
#'   \code{vignette("miami-bail", package = "clusterIV")}.
#'
#' @examples
#' ## Simulated judge design, as in ?cjar
#' set.seed(42)
#' G  <- 30; ng <- 8; n <- G * ng
#' cl <- rep(seq_len(G), each = ng)            # cluster identifier
#' judge <- factor(rep(1:6, length.out = n))   # judge identity (grouping z)
#' u  <- rnorm(G)[cl]                          # cluster-level shock
#' x  <- 0.6 * as.numeric(judge) + u + rnorm(n)
#' y  <- 0.5 * x + u + rnorm(n)
#'
#' ## Test a point null with the score test
#' sc <- cjscore(y, x, judge,
#'               cluster = cl,     # arbitrary dependence within clusters
#'               beta0   = 0,      # the point null of interest (default)
#'               level   = 0.95)   # level of the accompanying set (default)
#' sc
#'
#' ## Recommended workflow: CJIVE estimate, CJAR set, CJS point-null test
#' est <- cjive(y, x, judge, cluster = cl)
#' ar  <- cjar(y, x, judge, cluster = cl)
#' coef(est)
#' confint(ar)                                 # report this set
#' cjscore(y, x, judge, cluster = cl, beta0 = 1)$p.value  # test beta = 1
#'
#' @export
cjscore <- function(y, ...) {
  .check_dispatch_dots(
    y,
    match.call(expand.dots = FALSE)$...,
    c("x", "z", "cluster", "controls", "fixed_effects", "weights",
      "intercept", "level", "beta0", "variance", "data", "subset",
      "na.action"),
    "cjscore", parent.frame()
  )
  UseMethod("cjscore")
}

# Assemble a "cjscore" object from the prepared design, the shared
# per-cluster instrument sums and the leverage diagnostic; shared by
# cjscore.default and iv_infer() so the two produce identical fields.
.cjs_build <- function(d, csums, maxlev, eff, beta0, level, call,
                       variance = "plain", cf = NULL, term = "x") {
  cs <- .cjs_coefs(csums, cf = cf)
  k <- cs$k

  crit <- stats::qchisq(level, df = 1)
  inv <- .cjs_invert(cs$s, cs$v, crit)

  # Paper-scaled score and variance at beta0, with the LM ratio evaluated in
  # signed-log form when a large finite null would overflow the raw quadratic.
  # The conservative non-positive-variance rule is inside the shared helper,
  # which is also used by the plotted p-value curve.
  ev <- .cjs_point_eval(cs$s, cs$v, beta0, n = d$n)
  statistic <- ev$statistic
  p.value <- ev$p.value
  S0 <- ev$score
  V0 <- ev$variance

  topology <- .set_topology(inv$conf_set)
  structure(c(list(
    statistic = statistic, p.value = p.value,
    score = S0, variance = V0,
    term = term, beta0 = beta0, level = level,
    variance_estimator = variance, crit = crit,
    conf_set = inv$conf_set, shape = inv$shape,
    n_components = topology$n_components,
    unbounded_left = topology$unbounded_left,
    unbounded_right = topology$unbounded_right,
    bounded = nrow(inv$conf_set) == 0L || all(is.finite(inv$conf_set)),
    F_CJS = inv$F_CJS, k = k,
    coef_score = cs$s, coef_var = cs$v, call = call),
    .fit_common(d, eff, maxlev)), class = "cjscore")
}

#' @rdname cjscore
#' @export
cjscore.default <- function(y, x, z, cluster, controls = NULL,
                            fixed_effects = NULL, weights = NULL,
                            intercept = TRUE, level = 0.95, beta0 = 0,
                            variance = c("plain", "crossfit"), ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "cjscore")
  cl <- match.call()
  .check_flag(intercept, "intercept")
  variance <- match.arg(variance)
  term <- .term_label(substitute(x))
  .check_beta0_level(beta0, level)

  d <- .prep_data(y, x, z, cluster, controls, weights, intercept,
                  fixed_effects = fixed_effects)

  # The instrument Gram factorisation shared with the CJIVE/CJAR paths; no
  # first-stage fit is needed on the test path.
  R <- d$R

  csums <- .cluster_sums(d$y, d$x, d$Z, d$cluster, R = R)
  # As in cjar(): one whitened first-stage pass supplies the leverage
  # diagnostic and the effective-F ingredients.
  fs <- .first_stage(d$x, d$Z, R = R)
  maxlev <- .cluster_leverage(fs$Ztil, d$groups, d$k)
  eff <- .eff_f(fs$Ztil, fs$t, fs$e, d$cluster)
  .advise(d$G, d$groups, d$n, maxlev, "CJS")

  # The cross-fit variance replaces only the variance polynomial; the score
  # coefficients and every diagnostic above are shared.  Only the score
  # polynomial is consumed here (which = "score").
  cf <- if (variance == "crossfit") {
    .crossfit_var(d$y, d$x, d$Z, d$cluster, R = R, Ztil = fs$Ztil,
                  which = "score")
  }

  .cjs_build(d, csums, maxlev, eff, beta0, level, cl,
             variance = variance, cf = cf, term = term)
}

#' @rdname cjscore
#' @export
cjscore.formula <- function(formula, data, cluster, controls = NULL,
                            fixed_effects = NULL, weights = NULL,
                            subset, na.action = stats::na.omit,
                            intercept = TRUE, level = 0.95, beta0 = 0,
                            variance = c("plain", "crossfit"), ...) {
  .check_dots(match.call(expand.dots = FALSE)$..., "cjscore")
  .check_flag(intercept, "intercept")
  .check_beta0_level(beta0, level)
  .formula_method(formula, data, missing(data), substitute(cluster),
                  substitute(weights),
                  if (missing(subset)) NULL else substitute(subset),
                  controls, fixed_effects, na.action, intercept,
                  parent.frame(), match.call(),
                  fit = function(prep, intercept)
                    cjscore.default(prep$y, prep$x, prep$z,
                                    cluster = prep$cluster,
                                    controls = prep$controls,
                                    fixed_effects = prep$fixed_effects,
                                    weights = prep$weights,
                                    intercept = intercept, level = level,
                                    beta0 = beta0, variance = variance))
}
