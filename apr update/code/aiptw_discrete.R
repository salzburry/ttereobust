## aiptw_discrete.R
##
## Interactive single-DGM demo of AIPTW-PLR (logit, discrete-time) across
## the five main misspecification scenarios from the Apr 28 methods
## document Section 6.4. One simulated dataset, no bootstrap, plot of
## survival curves vs analytic Weibull truth.
##
## For protocol-grade replicated runs (1900 replicates with bootstrap
## CIs and the full scenario grid) use simulate_aiptw.R +
## summarise_aiptw.R; the canonical estimator lives in
## utils/aiptw_estimator.R and the spec strings in utils/scenarios.R.

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(tidyr)
  library(ggplot2);  library(MASS)
})

source("utils/sim_data.R")
source("utils/aiptw_estimator.R")
source("utils/compute_truth.R")
source("utils/scenarios.R")


## ---- Configuration ---------------------------------------------------------

simN   <- 2500                           # protocol Section 6.2
t.eval <- seq(0, 10, 0.25)               # discrete-grid for plotting
t.perf <- c(1, 5, 10)                    # protocol Section 6.3

# Uncomment desired DGM scenario:
# source("get_params_waning.R")
# source("get_params_delayed.R")
source("get_params_ph.R")

dgm_params <- list(
  N.Lcovs.linear = N.Lcovs.linear, N.Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma,
  alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens
)


## ---- Truth ----------------------------------------------------------------

message("Computing analytic Weibull truth...")
truth.fine <- compute_truth(t_eval = seq(0, admin.cens, 0.1),
                             dgm = dgm_params)
truth.perf <- compute_truth(t_eval = t.perf, dgm = dgm_params)


## ---- Analysis dataset -----------------------------------------------------

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


## ---- Run AIPTW for each main misspec scenario -----------------------------

all.est <- bind_rows(lapply(names(misspec_patterns), function(sc_name) {
  pat <- misspec_patterns[[sc_name]]
  message("  [", sc_name, "]  PS = ", pat$ps, "  |  Outcome = ", pat$out)
  est <- aiptw_estimate(
    surv.df,
    ps_spec    = ps.specs[[pat$ps]],
    out_spec   = out.specs[[pat$out]],
    t_eval     = t.eval,
    admin.cens = admin.cens
  )
  est$scenario <- pat$label
  est
}))


## ---- Performance at t = 1, 5, 10 ------------------------------------------

est.perf <- bind_rows(lapply(names(misspec_patterns), function(sc_name) {
  pat <- misspec_patterns[[sc_name]]
  est <- aiptw_estimate(
    surv.df,
    ps_spec    = ps.specs[[pat$ps]],
    out_spec   = out.specs[[pat$out]],
    t_eval     = t.perf,
    admin.cens = admin.cens
  )
  est$scenario <- pat$label
  est
}))

for (sc in unique(est.perf$scenario)) {
  cat("\n", sc, "\n", strrep("-", nchar(sc)), "\n", sep = "")
  sub <- dplyr::filter(est.perf, scenario == sc)
  for (k in seq_len(nrow(sub))) {
    tr <- truth.perf[truth.perf$t == sub$t[k], ]
    cat(sprintf("  t=%2g  S0=%.4f (true %.4f)  S1=%.4f (true %.4f)  RD=%+.4f (true %+.4f)\n",
                sub$t[k],
                sub$S0[k], tr$S0,
                sub$S1[k], tr$S1,
                sub$RD[k], tr$RD))
  }
}


## ---- Plot -----------------------------------------------------------------
## Saved to figs/ so it is visible whether the script is run via Rscript
## (non-interactive, no graphics device) or sourced inside RStudio.

dir.create("figs", showWarnings = FALSE)

plot.df <- bind_rows(
  all.est %>% transmute(time = t, A = 0L, surv = S0, scenario),
  all.est %>% transmute(time = t, A = 1L, surv = S1, scenario),
  truth.fine %>% transmute(time = t, A = 0L, surv = S0,
                            scenario = "True Weibull"),
  truth.fine %>% transmute(time = t, A = 1L, surv = S1,
                            scenario = "True Weibull")
) %>%
  mutate(
    # Display-only clipping; raw estimates are preserved on est.perf above.
    surv = pmin(pmax(surv, 0), 1),
    A    = factor(A, levels = c(0L, 1L),
                  labels = c("Control (A=0)", "Treatment (A=1)"))
  )

p <- ggplot(plot.df, aes(x = time, y = surv, colour = scenario)) +
  geom_step(linewidth = 0.6) +
  facet_wrap(~A) +
  ylim(0, 1) + xlim(0, admin.cens) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 3)) +
  xlab("Time (years)") + ylab("Marginal Survival Probability") +
  ggtitle("AIPTW-PLR: Marginal Survival by Misspecification Scenario",
          subtitle = sprintf("N = %d  |  one dataset  |  no bootstrap",
                              simN))

# Print so it shows up in interactive mode, AND save so it shows up in
# non-interactive mode.
print(p)
out_path <- file.path("figs", "aiptw_discrete_misspec_curves.png")
ggsave(out_path, p, width = 9, height = 5, dpi = 120)
message("[plot] wrote ", out_path)
