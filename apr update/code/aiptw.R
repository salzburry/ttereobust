## =============================================================================
## aiptw.R
## Augmented Inverse Probability of Treatment Weighting (AIPTW)
## for Survival Outcomes
##
## Two estimators are implemented:
##
##  1. AIPTW-PLR
##     Weighted pooled logistic regression (PLR) with the canonical logit link
##     followed by regression standardisation. By the result of Blanche et al.
##     (2023) and Gabriel et al. (2024, Stat Med) this is algebraically
##     equivalent to the classical AIPTW / AIPW estimator. The logit link is
##     the canonical link for the Bernoulli likelihood; using a non-canonical
##     link such as cloglog (as in the DR-Standardisation in analysis.R) breaks
##     the double-robustness property.
##
##  2. AIPTW-Cox
##     Unweighted Cox PH outcome model combined with the explicit AIPTW
##     augmentation term. Independent censoring is handled via IPCW using the
##     Kaplan-Meier estimator of the censoring distribution.
##
##     Formula (Bang & Robins 2005; Kurz 2022):
##       S^a_AIPTW(t) = (1/n) sum_i {
##           Q_i(a,t)  +  [I(A_i=a)/pi_i(a)]  *  [Y~_i(t) - Q_i(a,t)]
##       }
##     where
##       Q_i(a,t)  = exp(-H0(t) * exp(LP_i(a)))     Cox predicted survival
##       pi_i(a)   = P(A_i = a | L_i)                propensity score
##       Y~_i(t)   = I(T~_i > t) / G(t)              IPCW-corrected outcome
##       G(t)      = KM censoring survival at t
##
## Double robustness: consistent if EITHER the propensity score model OR the
## outcome model is correctly specified.
##
## Misspecification scenarios (per protocol Section 6.4):
##   Exposure:  correct / missing W / wrong functional form / heavy
##   Outcome:   correct / missing O,W / wrong functional form / heavy
##
## References:
##   Bang & Robins (2005) Biometrics
##   Funk et al. (2011) Am J Epidemiol
##   Kurz (2022) Med Decis Making
##   Blanche et al. (2023) Lifetime Data Anal
##   Gabriel et al. (2024) Stat Med
## =============================================================================


# ---- Libraries ---------------------------------------------------------------

library(survival)
library(dplyr)
library(ggplot2)
library(MASS)
library(purrr)

source("utils/sim_data.R")


# ---- Simulation parameters ---------------------------------------------------

simN <- 10000

# Uncomment the desired data-generating scenario:
# source("get_params_waning.R")
# source("get_params_delayed.R")
source("get_params_ph.R")

rescale_time <- 1/4                           # weekly discrete-time intervals
cutpoints    <- seq(0, admin.cens, rescale_time)
t.eval       <- seq(0, admin.cens, 0.1)       # grid for evaluation / plotting


# ---- Model specification scenarios (per protocol Section 6.4) ----------------

exposure.covs <- c(
  "L1sq + L2sq + L3 + L4 + L5 + L6 + W",          # 1: correct
  "L1sq + L2sq + L3 + L4 + L5 + L6",               # 2: missing W
  "L1 + L2 + L3 + L4 + L5 + L6 + W",               # 3: wrong functional form
  "L3 + L4 + L5"                                     # 4: heavy misspecification
)

outcome.covs <- c(
  "A + L1sq + L2sq + L3 + L4 + L5 + L6 + O + W",  # 1: correct
  "A + L1sq + L2sq + L3 + L4 + L5 + L6",           # 2: missing O and W
  "A + L1 + L2 + L3 + L4 + L5 + L6 + O + W",      # 3: wrong functional form
  "A + L3 + L4 + L5"                                # 4: heavy misspecification
)

# Active scenario indices (1 = correctly specified)
exp.idx <- 1
out.idx <- 1


# ---- Generate true marginal survival curves ----------------------------------
# Use a very large sample to approximate the population-level truth via
# the Weibull g-computation integral (same method as analysis.R / gen_truth.R).

message("Generating true survival curves (N = 500,000)...")

true.df1 <- sim_surv_data(
  seed          = this.seed,  N             = 500000,
  Lcovs.linear  = N.Lcovs.linear,  Lcovs.sq = N.Lcovs.sq,
  mu            = mu,               sigma    = sigma,
  alpha.L       = alpha.L,          alpha.W  = alpha.W,
  coeff.A       = coeff.A,          coeff.L  = coeff.L,
  coeff.Lsq     = coeff.Lsq,        coeff.O  = coeff.O,
  coeff.W       = coeff.W,
  gamma.tte     = gamma.tte,  lambda.tte  = lambda.tte,
  lambda.cens   = lambda.cens, admin.cens = admin.cens,
  gen.truth     = 1
)

cov.mat   <- true.df1$cov.mat
cov.mat.L <- cov.mat[, 1:(N.Lcovs.linear + N.Lcovs.sq)]

# Pre-compute individual linear predictors (avoids recomputation inside sapply)
lin_pred_1 <- as.vector(
  coeff.A +
  cov.mat.L %*% coeff.L +
  cov.mat.L[, 1:2]^2 %*% coeff.Lsq[1:2] +
  coeff.O * cov.mat[, "O"] +
  coeff.W * cov.mat[, "W"]
)
lin_pred_0 <- lin_pred_1 - coeff.A   # set A = 0

surv1_wei <- sapply(t.eval, function(t)
  mean(exp(-lambda.tte * t^gamma.tte * exp(lin_pred_1))))
surv0_wei <- sapply(t.eval, function(t)
  mean(exp(-lambda.tte * t^gamma.tte * exp(lin_pred_0))))

true.surv.df <- data.frame(
  time   = rep(t.eval, 2),
  A      = rep(c(0L, 1L), each = length(t.eval)),
  surv   = c(surv0_wei, surv1_wei),
  method = "True Weibull"
)


# ---- Generate analysis dataset -----------------------------------------------

message("Simulating analysis dataset (N = ", simN, ")...")

sim.df <- sim_surv_data(
  seed         = this.seed,  N            = simN,
  Lcovs.linear = N.Lcovs.linear,  Lcovs.sq = N.Lcovs.sq,
  mu           = mu,               sigma    = sigma,
  alpha.L      = alpha.L,          alpha.W  = alpha.W,
  coeff.A      = coeff.A,          coeff.L  = coeff.L,
  coeff.Lsq    = coeff.Lsq,        coeff.O  = coeff.O,
  coeff.W      = coeff.W,
  gamma.tte    = gamma.tte,  lambda.tte  = lambda.tte,
  lambda.cens  = lambda.cens, admin.cens = admin.cens,
  gen.truth    = NA
)

surv.df <- sim.df$data

# Long format (person-period) for discrete-time models
surv.long.df <- survSplit(
  Surv(eventtime, event) ~ .,
  data    = surv.df,
  cut     = cutpoints,
  episode = "time_period"
)

# Time-period to calendar-time mapping
interval_mapping <- data.frame(
  time_period = seq_along(cutpoints),
  time        = cutpoints
)

# Baseline covariate lookup (used when building prediction grids)
baseline.covs <- dplyr::select(
  surv.df, id, L1, L1sq, L2, L2sq, L3, L4, L5, L6, O, W
)


# ---- Censoring survival for IPCW (AIPTW-Cox only) ---------------------------
# Reverse the event indicator: censored observations become "events" for KM.
# This gives G(t) = P(C > t), the marginal censoring survival under independent
# censoring.  Pre-compute at every evaluation point for efficiency.

cens.km <- survfit(Surv(eventtime, 1L - event) ~ 1, data = surv.df)
G_t_vec <- summary(cens.km, times = t.eval, extend = TRUE)$surv
G_t_vec <- pmax(G_t_vec, 1e-6)   # guard against zero division at late times


# ---- IPTW weights (reused from colleague's IPTW code in analysis.R) ----------

compute_iptw <- function(exposure.formula.str, surv.df, surv.long.df) {

  # Denominator: P(A | L)
  ps.mod  <- glm(as.formula(paste("A ~", exposure.formula.str)),
                 family = "binomial", data = surv.df)
  pi_1    <- predict(ps.mod, type = "response")   # P(A=1|L), length = N

  # Marginal P(A=1) for stabilisation (scalar)
  pi_marg <- mean(surv.df$A == 1L)

  # Stabilised weight: P(A_i) / P(A_i | L_i)
  wt_denom <- ifelse(surv.df$A == 1L, pi_1,       1 - pi_1)
  wt_num   <- ifelse(surv.df$A == 1L, pi_marg, 1 - pi_marg)
  ipw_s    <- wt_num / wt_denom

  # Merge weights into long-format data (treatment is time-fixed, so the
  # same weight applies to every time-row of a given individual)
  wt.df <- data.frame(id = surv.df$id, ipw_s = ipw_s)
  survwt.long.df <- left_join(surv.long.df, wt.df, by = "id")

  list(
    ps.mod         = ps.mod,
    pi_1           = pi_1,           # P(A=1|L) for every individual
    ipw_s          = ipw_s,          # stabilised weights (wide, one per person)
    survwt.long.df = survwt.long.df  # long-format with ipw_s column
  )
}


# ---- AIPTW-PLR ---------------------------------------------------------------
#
# Step 1: Fit PLR on long-format data weighted by IPTW using the CANONICAL
#         logit link.  The logit link is the natural exponential family link for
#         the Bernoulli distribution.  Gabriel et al. (2024, Stat Med) prove
#         that IPTW-weighted GLM with a canonical link, followed by regression
#         standardisation, is algebraically equivalent to the AIPW estimator.
#         Using a non-canonical link (cloglog) breaks this equivalence.
#
# Step 2: Build a full counterfactual prediction grid (every individual x every
#         time-period x each treatment arm) and obtain discrete hazard
#         predictions from the fitted model.
#
# Step 3: Compute individual cumulative survival: cumprod(1 - hazard).
#
# Step 4: Standardise (average) individual survival curves over the empirical
#         covariate distribution to obtain the marginal survival curve.

run_aiptw_plr <- function(outcome.covs.str, survwt.long.df,
                           baseline.covs, interval_mapping, method.name) {

  plr.formula <- as.formula(
    paste0("event ~ ", outcome.covs.str, " + as.factor(time_period)")
  )

  # Weighted PLR — quasibinomial to accommodate non-integer weights
  # CRITICAL: logit link (canonical), NOT cloglog
  plr.mod <- glm(
    plr.formula,
    data    = survwt.long.df,
    family  = quasibinomial(link = "logit"),
    weights = ipw_s
  )

  # Full counterfactual prediction grid
  time_periods <- sort(unique(survwt.long.df$time_period))

  pred.df <- expand.grid(
    id          = baseline.covs$id,
    time_period = time_periods,
    A           = c(0L, 1L),
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    left_join(baseline.covs, by = "id") %>%
    mutate(event = 0L)

  pred.df$haz <- predict(plr.mod, newdata = pred.df, type = "response")

  # Individual cumulative survival
  pred.df <- pred.df %>%
    arrange(id, A, time_period) %>%
    group_by(id, A) %>%
    mutate(csurv = cumprod(1 - haz)) %>%
    ungroup()

  # Marginalise over individuals, map to calendar time, prepend t = 0
  pred.df %>%
    group_by(time_period, A) %>%
    summarise(surv = mean(csurv), .groups = "drop") %>%
    left_join(dplyr::select(interval_mapping, time_period, time),
              by = "time_period") %>%
    bind_rows(
      data.frame(time_period = NA_integer_, A = c(0L, 1L), surv = 1, time = 0),
      .
    ) %>%
    arrange(A, time) %>%
    mutate(method = method.name) %>%
    dplyr::select(time, A, surv, method)
}


# ---- AIPTW-Cox ---------------------------------------------------------------
#
# Fits an unweighted Cox PH outcome model and applies the explicit AIPTW
# augmentation term with IPCW correction for independent right censoring:
#
#   S^a_AIPTW(t) = (1/n) sum_i {
#       Q_i(a,t)  +  [I(A_i=a) / pi_i(a)]  *  [Y~_i(t)  -  Q_i(a,t)]
#   }
#
# Q_i(a,t) is obtained from the Cox baseline cumulative hazard and the
# individual linear predictors, avoiding expensive per-individual survfit calls:
#   Q_i(a,t) = exp( -H0(t) * exp(LP_i(a)) )
#
# Y~_i(t) = I(T~_i > t) / G(t) is the IPCW-corrected observed outcome.
# Under independent censoring the marginal KM estimate G(t) suffices.
#
# NOTE: the Cox model is unweighted — the IPW enters only through the
# explicit augmentation term, not through the outcome model itself.

run_aiptw_cox <- function(outcome.covs.str, surv.df, pi_1, G_t_vec,
                           t.eval, method.name) {

  cox.formula <- as.formula(
    paste0("Surv(eventtime, event) ~ ", outcome.covs.str)
  )

  # Unweighted Cox outcome model
  cox.mod <- coxph(cox.formula, data = surv.df, ties = "breslow", x = TRUE)

  # Baseline cumulative hazard H0(t), un-centred
  H0 <- basehaz(cox.mod, centered = FALSE)

  # Linear predictors under each counterfactual treatment arm
  lp0 <- predict(cox.mod, newdata = mutate(surv.df, A = 0), type = "lp")
  lp1 <- predict(cox.mod, newdata = mutate(surv.df, A = 1), type = "lp")

  # Unstabilised IPW indicator weights: I(A_i = a) / P(A_i = a | L_i)
  # (unstabilised is correct for the AIPTW augmentation term)
  ipw_ind0 <- ifelse(surv.df$A == 0L, 1 / (1 - pi_1), 0)
  ipw_ind1 <- ifelse(surv.df$A == 1L, 1 / pi_1,       0)

  obs_time <- surv.df$eventtime

  # Pre-compute baseline cumulative hazard at all evaluation points
  h0_all <- approx(H0$time, H0$hazard, xout = t.eval, rule = 2)$y

  # Vectorised loop over time points (lapply + rbind faster than map_dfr here)
  out.list <- lapply(seq_along(t.eval), function(k) {

    h0t <- h0_all[k]
    G_t <- G_t_vec[k]

    # Individual Cox survival predictions: S_i(t|a) = exp(-H0(t)*exp(LP_i(a)))
    Q0 <- exp(-h0t * exp(lp0))
    Q1 <- exp(-h0t * exp(lp1))

    # IPCW-corrected observed outcome
    Y_ipcw <- as.numeric(obs_time > t.eval[k]) / G_t

    # AIPTW: g-computation + IPW-weighted residual correction
    surv0 <- mean(Q0 + ipw_ind0 * (Y_ipcw - Q0))
    surv1 <- mean(Q1 + ipw_ind1 * (Y_ipcw - Q1))

    # Clip to [0,1]: can deviate slightly due to sampling noise
    data.frame(
      time   = t.eval[k],
      A      = c(0L, 1L),
      surv   = pmin(pmax(c(surv0, surv1), 0), 1),
      method = method.name
    )
  })

  do.call(rbind, out.list)
}


# ---- Performance evaluation helper -------------------------------------------

eval_performance <- function(res.df, true_S0, true_S1,
                              t.perf = c(1, 5, 10)) {

  label <- unique(res.df$method)
  cat("\n", label, "\n", strrep("-", nchar(label)), "\n", sep = "")

  true_list <- list("0" = true_S0, "1" = true_S1)

  for (a in c(0L, 1L)) {
    sub    <- dplyr::filter(res.df, A == a)
    true_S <- true_list[[as.character(a)]]

    for (t_pt in t.perf) {
      idx_true <- which.min(abs(t.eval - t_pt))
      idx_est  <- which.min(abs(sub$time - t_pt))
      truth    <- true_S[idx_true]
      est      <- sub$surv[idx_est]
      rel_bias <- 100 * (est - truth) / truth

      cat(sprintf("  A=%d  t=%-3g  Est=%.4f  True=%.4f  RelBias=%+.2f%%\n",
                  a, t_pt, est, truth, rel_bias))
    }
  }
}


# ==============================================================================
# Main analysis (default: both models correctly specified, exp.idx=1, out.idx=1)
# ==============================================================================

message("\n--- IPTW weights (exposure model ", exp.idx, ") ---")
iptw.res <- compute_iptw(exposure.covs[exp.idx], surv.df, surv.long.df)

message("--- AIPTW-PLR (outcome model ", out.idx, ") ---")
aiptw.plr.df <- run_aiptw_plr(
  outcome.covs.str = outcome.covs[out.idx],
  survwt.long.df   = iptw.res$survwt.long.df,
  baseline.covs    = baseline.covs,
  interval_mapping = interval_mapping,
  method.name      = sprintf("AIPTW-PLR (exp=%d, out=%d)", exp.idx, out.idx)
)

message("--- AIPTW-Cox (outcome model ", out.idx, ") ---")
aiptw.cox.df <- run_aiptw_cox(
  outcome.covs.str = outcome.covs[out.idx],
  surv.df          = surv.df,
  pi_1             = iptw.res$pi_1,
  G_t_vec          = G_t_vec,
  t.eval           = t.eval,
  method.name      = sprintf("AIPTW-Cox (exp=%d, out=%d)", exp.idx, out.idx)
)

# Performance at t = 1, 5, 10
message("\n===== Performance measures =====")
eval_performance(aiptw.plr.df, surv0_wei, surv1_wei)
eval_performance(aiptw.cox.df, surv0_wei, surv1_wei)


# ==============================================================================
# Optional: Loop over all 4 x 4 misspecification scenarios
# (demonstrates double robustness: estimate stays consistent when at least one
#  model is correctly specified)
# Uncomment to run.
# ==============================================================================

# all.aiptw.plr <- list()
# all.aiptw.cox <- list()
#
# for (ei in 1:4) {
#   for (oi in 1:4) {
#     key <- paste0("e", ei, "_o", oi)
#     message("Running scenario: exp=", ei, " out=", oi)
#
#     iptw.i <- compute_iptw(exposure.covs[ei], surv.df, surv.long.df)
#
#     all.aiptw.plr[[key]] <- run_aiptw_plr(
#       outcome.covs[oi], iptw.i$survwt.long.df, baseline.covs,
#       interval_mapping, sprintf("AIPTW-PLR (e%d,o%d)", ei, oi)
#     )
#     all.aiptw.cox[[key]] <- run_aiptw_cox(
#       outcome.covs[oi], surv.df, iptw.i$pi_1, G_t_vec, t.eval,
#       sprintf("AIPTW-Cox (e%d,o%d)", ei, oi)
#     )
#   }
# }


# ==============================================================================
# Visualisation
# ==============================================================================

plot.df <- bind_rows(true.surv.df, aiptw.plr.df, aiptw.cox.df) %>%
  mutate(
    A           = factor(A, levels = c(0L, 1L),
                         labels = c("Control (A=0)", "Treatment (A=1)")),
    method_grp  = sub(" \\(.*", "", method)   # strip scenario label for colour
  )

method.colors <- c(
  "True Weibull" = "black",
  "AIPTW-PLR"    = "#E41A1C",
  "AIPTW-Cox"    = "#377EB8"
)

ggplot(plot.df,
       aes(x = time, y = surv,
           colour   = method_grp,
           linetype = method)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~A) +
  scale_colour_manual(values = method.colors, name = "Estimator") +
  scale_linetype_discrete(name = "Model spec.") +
  ylim(0, 1) + xlim(0, admin.cens) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.box      = "vertical",
        legend.text     = element_text(size = 8)) +
  xlab("Time (years)") +
  ylab("Marginal Survival Probability") +
  ggtitle(
    sprintf("AIPTW: Marginal Survival Curves  [exp=%d, out=%d]",
            exp.idx, out.idx),
    subtitle = paste0(
      "Exposure model: ", exposure.covs[exp.idx], "\n",
      "Outcome model:  ", outcome.covs[out.idx]
    )
  )
