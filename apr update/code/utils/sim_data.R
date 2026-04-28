sim_surv_data <- function(seed = 2026,
                          N = 100,
                          Lcovs.linear = 4,
                          Lcovs.sq = 2,
                          mu = rep(0,8), # Lcovs.linear + Lcovs.sq + 2 (W and O)
                          sigma = diag(8),
                          alpha.L = c(0.005, 0.001, -0.3, 0.1, 0.2, 0.15), # params for exposure model c(L_1, ..., L_{Lcovs.linear + Lcovs.sq})
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

  ## --- Simulate covs
  tot.Lcovs <- Lcovs.linear + Lcovs.sq
  Ncovs <- tot.Lcovs + 2

  cov.mat <- mvrnorm(n = N, mu = mu, Sigma = sigma)
  L.names <- lapply(c(1:tot.Lcovs), function(l){ paste0("L",l) }) %>% list_c
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

  if(!is.na(gen.truth)) {
    A <- gen.truth
  } else {
    A <- rbinom(N, 1, expit(-1 +
                              model.probA.L[1]*cov.mat[,"L1"]^2 +
                              model.probA.L[2]*cov.mat[,"L2"]^2 +
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


  # simulate time from a weibull distribution
  X <- (-log(u.t)/weihaz.denom)^(1/gamma)

  C <- rexp(N, cens_lambda)

  T.obs <- pmin(X, C, maxT)

  D.obs <- ifelse(pmin(X, C, maxT) == X,
                  1,
                  0
  )

  covs.df <- as.data.frame(cov.mat) %>% mutate(id = row_number()) %>%
    mutate(L1sq = L1^2, L2sq = L2^2) # add non-linear effects

  # The final simulate dataset
  surv.df <- data.frame(id=1:N, eventtime = T.obs, event = D.obs, A = A) %>%
    left_join(., covs.df, by = c("id")) %>%
    mutate(eventtime = eventtime + runif(n(), 0, 1) / 1000) # ensure no ties, so rescale eventtime

  return(list(data = surv.df, cov.mat = cov.mat))
  }
