## aiptw.R
##
## AIPTW-Cox: continuous-time variant of the AIPTW estimator that uses an
## unweighted Cox PH outcome model combined with the explicit AIPTW
## augmentation term and IPCW correction (marginal KM under independent
## censoring).
##
## Adds value over aiptw_discrete.R / utils/aiptw_estimator.R, which use
## a discrete-time logit-PLR outcome model (the canonical link for the
## AIPTW algebraic-equivalence result, Gabriel et al. 2024). This script
## is the "Cox outcome" comparison arm.
##
## Formula:
##   S^a_AIPTW(t) = (1/n) sum_i {
##       Q_i(a, t)  +  [I(A_i = a) / pi_i(a)]  *  [Y~_i(t) - Q_i(a, t)]
##   }
## with
##   Q_i(a, t) = exp(-H0(t) * exp(LP_i(a)))
##   pi_i(a)   = P(A_i = a | L_i)                trimmed propensity
##   Y~_i(t)   = I(T~_i >= t) / G(t-)            IPCW-corrected outcome
##
## For protocol-grade replicated runs use simulate_aiptw.R (logit-PLR
## variant). The Cox variant has not been wired into the harness; this
## script is a single-dataset interactive demo.

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(ggplot2); library(MASS)
})

source("utils/sim_data.R")
source("utils/compute_truth.R")
source("utils/scenarios.R")


## ---- AIPTW-Cox estimator --------------------------------------------------
##
## Pure function: takes one dataset and the spec strings, returns a data
## frame with t / S0 / S1 / RD on the requested grid.

run_aiptw_cox <- function(surv.df, ps_spec, out_spec,
                           t_eval, admin.cens,
                           ps_trim = c(0.01, 0.99)) {

  # Propensity score (trimmed)
  ps.mod <- glm(as.formula(paste("A ~", ps_spec)),
                family = "binomial", data = surv.df)
  pi_1 <- pmin(pmax(predict(ps.mod, type = "response"),
                    ps_trim[1]), ps_trim[2])
  ipw_ind1 <- ifelse(surv.df$A == 1L, 1 / pi_1,       0)
  ipw_ind0 <- ifelse(surv.df$A == 0L, 1 / (1 - pi_1), 0)

  # Cox outcome model + baseline cumulative hazard (with H0(0) = 0
  # prepended so survival starts at 1 exactly).
  cox.mod <- coxph(as.formula(paste("Surv(eventtime, event) ~", out_spec)),
                   data = surv.df, ties = "breslow", x = TRUE)
  H0 <- rbind(data.frame(time = 0, hazard = 0),
              basehaz(cox.mod, centered = FALSE))
  lp0 <- predict(cox.mod, newdata = mutate(surv.df, A = 0), type = "lp")
  lp1 <- predict(cox.mod, newdata = mutate(surv.df, A = 1), type = "lp")
  h0_all <- approx(H0$time, H0$hazard, xout = t_eval,
                   method = "constant", f = 0, rule = 2)$y

  # IPCW G(t-) (left limit; avoids divide-by-zero at admin.cens)
  cens.km <- survfit(Surv(eventtime, 1L - event) ~ 1, data = surv.df)
  G_tminus <- pmax(
    summary(cens.km, times = pmax(t_eval - 1e-8, 0), extend = TRUE)$surv,
    1e-6
  )

  obs_time <- surv.df$eventtime
  do.call(rbind, lapply(seq_along(t_eval), function(k) {
    t_k <- t_eval[k]
    if (t_k == 0) return(data.frame(t = 0, S0 = 1, S1 = 1, RD = 0))
    Q0 <- exp(-h0_all[k] * exp(lp0))
    Q1 <- exp(-h0_all[k] * exp(lp1))
    Y_ipcw <- as.numeric(obs_time >= t_k) / G_tminus[k]
    S0 <- mean(Q0 + ipw_ind0 * (Y_ipcw - Q0))
    S1 <- mean(Q1 + ipw_ind1 * (Y_ipcw - Q1))
    data.frame(t = t_k, S0 = S0, S1 = S1, RD = S1 - S0)
  }))
}


## ---- Configuration --------------------------------------------------------

simN   <- 2500                          # protocol Section 6.2
t.eval <- seq(0, 10, 0.1)               # plot grid
t.perf <- c(1, 5, 10)                   # protocol Section 6.3

# Active scenario (any pair of names from utils/scenarios.R)
ps.idx  <- "correct"
out.idx <- "correct"

# Uncomment desired DGM:
# source("get_params_waning.R")
# source("get_params_delayed.R")
source("get_params_ph.R")


## ---- Truth & dataset ------------------------------------------------------

dgm_params <- list(
  N.Lcovs.linear = N.Lcovs.linear, N.Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma,
  alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens
)

message("Computing analytic Weibull truth...")
truth.fine <- compute_truth(t_eval = t.eval, dgm = dgm_params)
truth.perf <- compute_truth(t_eval = t.perf, dgm = dgm_params)

message("Simulating analysis dataset (N = ", simN, ")...")
sim <- sim_surv_data(
  seed = this.seed, N = simN,
  Lcovs.linear = N.Lcovs.linear, Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma,
  alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens
)
surv.df <- sim$data


## ---- Run AIPTW-Cox --------------------------------------------------------

est <- run_aiptw_cox(
  surv.df,
  ps_spec    = ps.specs[[ps.idx]],
  out_spec   = out.specs[[out.idx]],
  t_eval     = t.eval,
  admin.cens = admin.cens
)

est.perf <- run_aiptw_cox(
  surv.df,
  ps_spec    = ps.specs[[ps.idx]],
  out_spec   = out.specs[[out.idx]],
  t_eval     = t.perf,
  admin.cens = admin.cens
)


## ---- Performance ---------------------------------------------------------

cat(sprintf("\nAIPTW-Cox  (PS = %s | Outcome = %s)\n", ps.idx, out.idx))
cat(strrep("-", 40), "\n", sep = "")
for (k in seq_len(nrow(truth.perf))) {
  e <- est.perf[which.min(abs(est.perf$t - truth.perf$t[k])), ]
  cat(sprintf("  t=%2g  S0=%.4f (true %.4f)  S1=%.4f (true %.4f)  RD=%+.4f (true %+.4f)\n",
              truth.perf$t[k],
              e$S0, truth.perf$S0[k],
              e$S1, truth.perf$S1[k],
              e$RD, truth.perf$RD[k]))
}


## ---- Plot ----------------------------------------------------------------

plot.df <- bind_rows(
  est        %>% transmute(time = t, A = 0L, surv = S0, method = "AIPTW-Cox"),
  est        %>% transmute(time = t, A = 1L, surv = S1, method = "AIPTW-Cox"),
  truth.fine %>% transmute(time = t, A = 0L, surv = S0, method = "True Weibull"),
  truth.fine %>% transmute(time = t, A = 1L, surv = S1, method = "True Weibull")
) %>%
  mutate(
    # Display-only clipping; raw estimates are on est / est.perf above.
    surv = pmin(pmax(surv, 0), 1),
    A    = factor(A, levels = c(0L, 1L),
                  labels = c("Control (A=0)", "Treatment (A=1)"))
  )

ggplot(plot.df, aes(x = time, y = surv, colour = method)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~A) +
  ylim(0, 1) + xlim(0, admin.cens) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  xlab("Time (years)") + ylab("Marginal Survival Probability") +
  ggtitle(sprintf("AIPTW-Cox  (PS = %s | Outcome = %s)",
                  ps.idx, out.idx))
