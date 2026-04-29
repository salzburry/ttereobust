## aiptw.R
##
## Interactive single-dataset demo of the AIPTW-Cox arm. This script
## now uses utils/aiptw_cox_estimator.R (the same estimator the
## protocol harness runs when --method aiptw_cox or --method both is
## passed to simulate_aiptw.R), so the demo and the production harness
## are guaranteed to compute identical numbers on the same dataset.
##
## For protocol-grade replicated runs (R x scenario grid x bootstrap
## CIs) use simulate_aiptw.R --method aiptw_cox or run_aiptw_cox_analysis.R.

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(ggplot2); library(MASS)
})

source("utils/sim_data.R")
source("utils/compute_truth.R")
source("utils/scenarios.R")


## ---- AIPTW-Cox estimator --------------------------------------------------
## Sourced from utils/ so this demo and the production harness use the
## SAME estimator (riskRegression::ate-based AIPTW with Cox outcome and
## a marginal censoring model). The local wrapper below adapts the
## utility's wide-format return to a slim (t, S0, S1, RD) data frame for
## plotting.

source("utils/aiptw_cox_estimator.R")

run_aiptw_cox <- function(surv.df, ps_spec, out_spec,
                           t_eval, admin.cens,
                           ps_trim = c(0.01, 0.99)) {
  est <- aiptw_cox_estimate(surv.df, ps_spec, out_spec,
                             t_eval = t_eval,
                             admin.cens = admin.cens,
                             ps_trim = ps_trim)
  est[, c("t", "S0", "S1", "RD")]
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
## Saved to figs/ so it is visible from Rscript as well as RStudio.

dir.create("figs", showWarnings = FALSE)

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

p <- ggplot(plot.df, aes(x = time, y = surv, colour = method)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~A) +
  ylim(0, 1) + xlim(0, admin.cens) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  xlab("Time (years)") + ylab("Marginal Survival Probability") +
  ggtitle(sprintf("AIPTW-Cox  (PS = %s | Outcome = %s)",
                  ps.idx, out.idx))

print(p)
out_path <- file.path("figs", sprintf("aiptw_cox_%s_%s.png", ps.idx, out.idx))
ggsave(out_path, p, width = 9, height = 5, dpi = 120)
message("[plot] wrote ", out_path)
