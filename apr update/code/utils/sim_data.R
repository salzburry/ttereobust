sim_surv_data <- function(seed = 2026,
                          N = 100,
                          Lcovs.linear = 4,
                          Lcovs.sq = 2,
                          mu = rep(0,8), # Lcovs.linear + Lcovs.sq + 2 (W and O)
                          sigma = diag(8),
                          alpha.L = c(0.5, 0.4, -0.6, 0.4, 0.5, 0.3), # Apr 28 doc Section 6.2 default; per-DGM files override
                          alpha.W = 1,
                          # coefficient / log-HRs for outcome model
                          coeff.A = log(0.9),
                          coeff.L = log(c(0.95, 1.1, 0.95, 0.9, 1.05, 1.01)),
                          coeff.Lsq = log(c(0.945, 1.05)),
                          coeff.O = log(0.9), # this will only appear in outcome model
                          coeff.W = log(0.9), # this will only appear in exposure model
                          gamma.tte = 1,
                          lambda.tte = 0.3,
                          lambda.cens = 0.2,
                          admin.cens = 10,
                          gen.truth = NA
                          ) {

  # Reproducibility: if a non-NA seed is supplied, set RNG state so that the
  # simulated dataset is fully determined by (seed, parameters). Pass seed = NA
  # (or NULL) to inherit the ambient RNG state — useful when an outer driver
  # is managing seeds across replications.
  if (!is.null(seed) && !is.na(seed)) {
    set.seed(seed)
  }

  ## --- Simulate covs
  tot.Lcovs <- Lcovs.linear + Lcovs.sq
  Ncovs <- tot.Lcovs + 2

  cov.mat <- mvrnorm(n = N, mu = mu, Sigma = sigma)
  L.names <- paste0("L", seq_len(tot.Lcovs))   # base R, no purrr dependency
  dimnames(cov.mat) <- list(NULL, c(L.names, "W", "O"))
  cov.mat.L <- cov.mat[,1:(Lcovs.linear+Lcovs.sq)]

  ## --- Set model params

  model.probA.L <- alpha.L
  model.probA.W <- alpha.W # this is param that only goes in exposure model

  gamma <- gamma.tte
  lambda <- lambda.tte

  cens_lambda <- lambda.cens # tweak this to get certain % of censoring

  maxT <- admin.cens

  ## --- Define exposure model
  expit <- function(x) { exp(x) / (1 + exp(x)) }

  # Exposure model — matches Apr 28 methods document Section 6.2:
  #   logit P(A = 1) = -2 + 0.5*L1 + 0.4*L2 + (-0.6)*L3 + 0.4*L4 + 0.5*L5 + 0.3*L6 + W
  # Coefficients are supplied via alpha.L (length 6, for L1..L6) and alpha.W.
  # L1, L2 enter linearly here; the quadratic L1, L2 terms appear only in the
  # outcome hazard (per the doc).
  if(!is.na(gen.truth)) {
    A <- gen.truth
  } else {
    A <- rbinom(N, 1, expit(-2 +
                              model.probA.L[1]*cov.mat[,"L1"] +
                              model.probA.L[2]*cov.mat[,"L2"] +
                              model.probA.L[3]*cov.mat[,"L3"] +
                              model.probA.L[4]*cov.mat[,"L4"] +
                              model.probA.L[5]*cov.mat[,"L5"] +
                              model.probA.L[6]*cov.mat[,"L6"] +
                              model.probA.W*cov.mat[,"W"]
    ) )
  }


  # Generate event times X, censoring times C and event indicators D.obs
  X <- rep(NA,N)

  u.t <- runif(N,0,1)

  weihaz.denom <- lambda * exp(coeff.A*A +
                                 cov.mat.L %*% coeff.L +
                                 cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
                                 coeff.W*cov.mat[,"W"] + coeff.O*cov.mat[,"O"]
  )


  # simulate time from a weibull distribution.
  # Tiny jitter is applied to the raw event/censoring times BEFORE the pmin with
  # admin.cens so the administrative boundary at maxT is preserved exactly. If
  # we jitter T.obs after pmin, admin-censored individuals can be pushed past
  # admin.cens and then mis-counted as still-at-risk by the >= IPCW outcome.
  X <- (-log(u.t)/weihaz.denom)^(1/gamma) + runif(N, 0, 1) / 1000

  C <- rexp(N, cens_lambda) + runif(N, 0, 1) / 1000

  T.obs <- pmin(X, C, maxT)

  D.obs <- ifelse(T.obs == X, 1, 0)   # event flag preserved post-jitter

  covs.df <- as.data.frame(cov.mat) %>% mutate(id = row_number()) %>%
    mutate(L1sq = L1^2, L2sq = L2^2) # add non-linear effects

  # The final simulated dataset
  surv.df <- data.frame(id=1:N, eventtime = T.obs, event = D.obs, A = A) %>%
    left_join(., covs.df, by = c("id"))

  return(list(data = surv.df, cov.mat = cov.mat))
  }
