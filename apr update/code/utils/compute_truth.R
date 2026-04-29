## utils/compute_truth.R
## Analytic Weibull g-formula truth for the marginal survival functions
## S^a(t) = E_L [ exp(-Lambda_a(t, L)) ] under the DGM in utils/sim_data.R.
##
## We approximate the expectation over L by Monte Carlo on a very large
## sample drawn under the same DGM; this matches the approach used in
## gen_truth.R / aiptw_discrete.R and avoids the need for analytic
## marginalisation over the multivariate-normal covariate distribution.

compute_truth <- function(t_eval = c(1, 5, 10),
                           dgm,
                           N_truth     = NULL,
                           truth_seed  = 1234567L) {

  # Honour dgm$N_truth if the caller passed it inside the DGM list (e.g.
  # from a YAML config). Explicit N_truth argument wins; otherwise fall
  # back to a 500000-row default.
  if (is.null(N_truth)) {
    N_truth <- if (!is.null(dgm$N_truth)) as.integer(dgm$N_truth) else 500000L
  }

  ## dgm must contain at least: N.Lcovs.linear, N.Lcovs.sq, mu, sigma,
  ##   alpha.L, alpha.W, coeff.A, coeff.L, coeff.Lsq, coeff.O, coeff.W,
  ##   gamma.tte, lambda.tte, lambda.cens, admin.cens

  big <- sim_surv_data(
    seed         = truth_seed,
    N            = N_truth,
    Lcovs.linear = dgm$N.Lcovs.linear,
    Lcovs.sq     = dgm$N.Lcovs.sq,
    mu           = dgm$mu,
    sigma        = dgm$sigma,
    alpha.L      = dgm$alpha.L,
    alpha.W      = dgm$alpha.W,
    coeff.A      = dgm$coeff.A,
    coeff.L      = dgm$coeff.L,
    coeff.Lsq    = dgm$coeff.Lsq,
    coeff.O      = dgm$coeff.O,
    coeff.W      = dgm$coeff.W,
    gamma.tte    = dgm$gamma.tte,
    lambda.tte   = dgm$lambda.tte,
    lambda.cens  = dgm$lambda.cens,
    admin.cens   = dgm$admin.cens,
    gen.truth    = 1
  )

  cov.mat   <- big$cov.mat
  cov.mat.L <- cov.mat[, seq_len(dgm$N.Lcovs.linear + dgm$N.Lcovs.sq)]

  # Linear predictor under control (excludes coeff.A * A); the treated arm
  # adds coeff.A in line.
  lp_base <- as.vector(
    cov.mat.L %*% dgm$coeff.L +
    cov.mat.L[, 1:2]^2 %*% dgm$coeff.Lsq[1:2] +
    dgm$coeff.O * cov.mat[, "O"] +
    dgm$coeff.W * cov.mat[, "W"]
  )

  S0 <- vapply(t_eval, function(t) {
    mean(exp(-dgm$lambda.tte * t^dgm$gamma.tte * exp(lp_base)))
  }, numeric(1))

  S1 <- vapply(t_eval, function(t) {
    mean(exp(-dgm$lambda.tte * t^dgm$gamma.tte *
             exp(dgm$coeff.A + lp_base)))
  }, numeric(1))

  data.frame(t = t_eval, S0 = S0, S1 = S1, RD = S1 - S0)
}
