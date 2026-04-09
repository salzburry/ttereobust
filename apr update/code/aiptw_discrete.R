## =============================================================================
## aiptw_discrete.R
## Augmented Inverse Probability of Treatment Weighting (AIPTW)
## Discrete-time implementation for survival outcomes
##
## Estimand (protocol Section 6.3):
##   S^a(t) = P(T^a > t)   for a in {0, 1}
##   RD(t)  = S^1(t) - S^0(t)   at t = 1, 5, 10
##
## AIPTW estimator (Bang & Robins 2005; Funk et al. 2011; Kurz 2022):
##
##   S^a_AIPTW(t) = (1/n) * sum_i {
##       Q_i(a, t)
##     + [I(A_i = a) / pi_i(a)] * [Ytilde_i(t) - Q_i(a, t)]
##   }
##
##   Q_i(a, t)    : individual survival P(T>t | A=a, L_i) from UNWEIGHTED PLR
##   pi_i(a)      : P(A_i=a | L_i) from the propensity score model
##   Ytilde_i(t)  : I(Ttilde_i > t) / G(t)  — IPCW-corrected observed outcome
##   G(t)         : KM estimate of censoring survival (independent censoring)
##
## Double robustness: consistent if EITHER the propensity score model OR the
## outcome model is correctly specified (but not necessarily both).
##
## This is NOT the same as DR Standardisation in analysis.R:
##   - DR Stand uses a WEIGHTED outcome model (IPTW applied to the GLM)
##     with a cloglog link, which is non-canonical and is NOT doubly robust
##     (Gabriel et al. 2024, Stat Med; Denz et al. 2023, Stat Med).
##   - AIPTW uses an UNWEIGHTED outcome model; the IPW enters only through
##     the explicit augmentation term. With a canonical logit link the two
##     are algebraically equivalent (Gabriel et al. 2024), but the explicit
##     formula here makes the double-robustness mechanism transparent.
##
## Misspecification scenarios follow protocol Section 6.4:
##   Each scenario isolates misspecification in ONE nuisance model only.
##
## DGM note on W and O:
##   W appears in both the exposure and outcome models in the DGM as coded
##   (coeff.W != 0 in the survival linear predictor). The get_params comment
##   "only appear in exposure model" is inconsistent with the truth formula;
##   this inconsistency is noted but not changed here to preserve comparability
##   with analysis.R. The key distinction between W and O is that O is
##   UNMEASURED (so it cannot appear in the exposure model), while W is
##   MEASURED (and used in both models).
##
## References:
##   Bang & Robins (2005) Biometrics 61:962-973
##   Funk et al. (2011) Am J Epidemiol 173:761-767
##   Kurz (2022) Med Decis Making 42:156-167
##   Gabriel et al. (2024) Stat Med 43:534-547
## =============================================================================


# ---- Libraries ---------------------------------------------------------------
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)

source("utils/sim_data.R")


# ---- Parameters --------------------------------------------------------------
simN <- 2500   # N per protocol Section 6.2 (NOT 10,000)

# Uncomment desired DGM scenario:
# source("get_params_waning.R")
# source("get_params_delayed.R")
source("get_params_ph.R")

rescale_time <- 1/4                          # weekly discrete-time intervals
cutpoints    <- seq(0, admin.cens, rescale_time)
t.eval       <- seq(0, admin.cens, 0.1)      # fine grid for plotting
t.perf       <- c(1, 5, 10)                  # protocol performance time points


# ---- Propensity score model specifications -----------------------------------
# Note: separate from outcome model — do NOT share a single index.
# Scenario 2 changes only the PS model; scenario 3 changes only the outcome.

ps.specs <- c(
  correct  = "L1sq + L2sq + L3 + L4 + L5 + L6 + W",   # 1
  no_W     = "L1sq + L2sq + L3 + L4 + L5 + L6",         # 2: missing W
  wrong_ff = "L1 + L2 + L3 + L4 + L5 + L6 + W",         # 3: wrong functional form
  heavy    = "L3 + L4 + L5"                               # 4: heavy misspec
)

# ---- Outcome model specifications --------------------------------------------
# Correct model includes O (the measured unmeasured confounder available in
# the outcome model) and W (which also affects outcome per the DGM truth formula).
# Misspecification scenarios remove one variable at a time or change functional form.

out.specs <- c(
  correct  = "A + L1sq + L2sq + L3 + L4 + L5 + L6 + O + W",   # 1
  no_O     = "A + L1sq + L2sq + L3 + L4 + L5 + L6 + W",         # 2: missing O
  wrong_ff = "A + L1 + L2 + L3 + L4 + L5 + L6 + O + W",        # 3: wrong functional form
  heavy    = "A + L3 + L4 + L5"                                   # 4: heavy misspec
)

# ---- Protocol misspecification scenarios (Section 6.4) ----------------------
# Each row isolates ONE misspecified model — this is the key design the
# shared cov.index in analysis.R cannot represent.

scenarios <- list(
  both_correct    = list(ps = "correct",  out = "correct",
                         label = "Both correct"),
  miss_W_ps_only  = list(ps = "no_W",     out = "correct",
                         label = "PS: missing W"),
  miss_O_out_only = list(ps = "correct",  out = "no_O",
                         label = "Outcome: missing O"),
  wrong_ff_out    = list(ps = "correct",  out = "wrong_ff",
                         label = "Outcome: wrong form"),
  wrong_ff_ps     = list(ps = "wrong_ff", out = "correct",
                         label = "PS: wrong form")
)


# ---- Generate true marginal survival curves ----------------------------------
# Analytical Weibull g-computation on a large population (same as analysis.R).
# Provides the benchmark S^0(t), S^1(t), RD(t) against which AIPTW is assessed.

message("Generating true survival curves (N = 500,000)...")

true.df1 <- sim_surv_data(
  seed = this.seed, N = 500000,
  Lcovs.linear = N.Lcovs.linear, Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma, alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens,
  gen.truth = 1
)

cov.mat   <- true.df1$cov.mat
cov.mat.L <- cov.mat[, 1:(N.Lcovs.linear + N.Lcovs.sq)]

# Individual log-linear predictors: shared across S^0 and S^1
lin_pred_base <- as.vector(
  cov.mat.L %*% coeff.L +
  cov.mat.L[, 1:2]^2 %*% coeff.Lsq[1:2] +
  coeff.O * cov.mat[, "O"] +
  coeff.W * cov.mat[, "W"]
)

surv1_wei <- sapply(t.eval, function(t)
  mean(exp(-lambda.tte * t^gamma.tte * exp(coeff.A + lin_pred_base))))
surv0_wei <- sapply(t.eval, function(t)
  mean(exp(-lambda.tte * t^gamma.tte * exp(lin_pred_base))))
rd_wei    <- surv1_wei - surv0_wei

true.surv.df <- data.frame(
  time   = rep(t.eval, 2),
  A      = rep(c(0L, 1L), each = length(t.eval)),
  surv   = c(surv0_wei, surv1_wei),
  method = "True Weibull"
)


# ---- Generate analysis dataset -----------------------------------------------

message("Simulating analysis dataset (N = ", simN, ")...")

sim.df <- sim_surv_data(
  seed = this.seed, N = simN,
  Lcovs.linear = N.Lcovs.linear, Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma, alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens,
  gen.truth = NA
)

surv.df <- sim.df$data    # one row per individual

# Long (person-period) format for discrete-time models
surv.long.df <- survSplit(
  Surv(eventtime, event) ~ .,
  data    = surv.df,
  cut     = cutpoints,
  episode = "time_period"
)

# Time-period → calendar-time mapping
interval_mapping <- data.frame(
  time_period = seq_along(cutpoints),
  time        = cutpoints
)

# Baseline covariate lookup: used to build counterfactual prediction grids
baseline.covs <- dplyr::select(surv.df, id, L1, L1sq, L2, L2sq, L3, L4, L5, L6, O, W)


# ---- Censoring KM for IPCW ---------------------------------------------------
# Reverse the event indicator so that censored observations become "events".
# G(t) = P(C > t) under independent censoring.
# Pre-computed once and reused across all scenarios.

cens.km <- survfit(Surv(eventtime, 1L - event) ~ 1, data = surv.df)
G_t_vec <- summary(cens.km, times = t.eval, extend = TRUE)$surv
G_t_vec <- pmax(G_t_vec, 1e-6)   # numerical safety for late time points


# ---- Core AIPTW function -----------------------------------------------------
#
# Takes one PS specification and one outcome specification.
# Returns marginal survival curves S^0(t), S^1(t) and risk differences RD(t).
#
# Algorithm:
#   1. Fit propensity score model  →  pi_1 = P(A=1|L)
#   2. Fit UNWEIGHTED PLR (logit)  →  conditional hazard h_k(a, L_i)
#   3. Build counterfactual grid   →  Q_i(a, t) = cumprod(1 - h_k(a, L_i))
#   4. Compute IPCW outcome        →  Ytilde_i(t) = I(Ttilde_i > t) / G(t)
#   5. Combine via AIPTW formula   →  S^a(t) + augmentation term

run_aiptw_discrete <- function(ps.spec.str, out.spec.str,
                                surv.df, surv.long.df,
                                baseline.covs, interval_mapping,
                                G_t_vec, t.eval,
                                method.label) {

  n <- nrow(surv.df)

  # ------ Step 1: Propensity score model ------
  ps.mod <- glm(as.formula(paste("A ~", ps.spec.str)),
                family = "binomial", data = surv.df)
  pi_1   <- predict(ps.mod, type = "response")   # P(A=1|L_i), length n

  # Unstabilised IPW indicators I(A_i=a) / P(A_i=a|L_i) for each arm.
  # These are used ONLY in the augmentation term, not in the outcome model.
  ipw_ind1 <- ifelse(surv.df$A == 1L, 1 / pi_1,       0)   # arm A=1
  ipw_ind0 <- ifelse(surv.df$A == 0L, 1 / (1 - pi_1), 0)   # arm A=0

  # ------ Step 2: Unweighted outcome model (PLR, canonical logit link) ------
  # Model is fitted on the observed data WITHOUT IPTW weighting.
  # Covariate adjustment alone provides the g-computation component Q_i(a,t).
  # Note: using a weighted PLR with logit link + standardisation would also be
  # a valid AIPTW implementation (algebraically equivalent under canonical link,
  # per Gabriel et al. 2024), but the explicit unweighted form here is clearer.
  plr.formula <- as.formula(
    paste0("event ~ ", out.spec.str, " + as.factor(time_period)")
  )
  plr.mod <- glm(plr.formula, data = surv.long.df,
                 family = binomial(link = "logit"))

  # ------ Step 3: Counterfactual prediction grid ------
  # All individuals x all time periods x each treatment arm.
  # Setting A = 0 or 1 for everyone (regardless of observed treatment) gives
  # counterfactual hazard predictions.
  time_periods <- sort(unique(surv.long.df$time_period))

  pred.df <- expand.grid(
    id          = baseline.covs$id,
    time_period = time_periods,
    A           = c(0L, 1L),
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    left_join(baseline.covs, by = "id") %>%
    mutate(event = 0L)

  pred.df$haz <- predict(plr.mod, newdata = pred.df, type = "response")

  # Individual survival Q_i(a, t_k) = product_{j <= k} (1 - h_j(a, L_i))
  pred.df <- pred.df %>%
    arrange(id, A, time_period) %>%
    group_by(id, A) %>%
    mutate(csurv = cumprod(1 - haz)) %>%
    ungroup()

  # ------ Step 4 & 5: AIPTW at each evaluation time point ------
  # For each t, extract Q_i(a, t) from the prediction grid, form the
  # IPCW outcome, and apply the augmentation formula.

  obs_time <- surv.df$eventtime   # T_tilde_i

  out.list <- lapply(seq_along(t.eval), function(k) {

    t   <- t.eval[k]
    G_t <- G_t_vec[k]

    if (t == 0) {
      return(data.frame(time = 0, A = c(0L, 1L), surv = 1,
                        method = method.label))
    }

    # Time period corresponding to calendar time t
    # (the last completed discrete interval at or before t)
    tp <- interval_mapping$time_period[max(which(interval_mapping$time <= t))]

    # Q_i(a, t): extract and align to surv.df row order via id join
    Q_at_t <- pred.df %>%
      dplyr::filter(time_period == tp) %>%
      dplyr::select(id, A, csurv)

    Q0 <- left_join(dplyr::select(surv.df, id),
                    dplyr::filter(Q_at_t, A == 0L),
                    by = "id")$csurv
    Q1 <- left_join(dplyr::select(surv.df, id),
                    dplyr::filter(Q_at_t, A == 1L),
                    by = "id")$csurv

    # IPCW-corrected observed outcome:  Ytilde_i(t) = I(T_tilde_i > t) / G(t)
    Y_ipcw <- as.numeric(obs_time > t) / G_t

    # AIPTW formula:
    #   S^a(t) = mean_i { Q_i(a,t) + [I(A_i=a)/pi_i(a)] * [Ytilde_i(t) - Q_i(a,t)] }
    #          = g-computation term + IPW-weighted residual correction
    surv0 <- mean(Q0 + ipw_ind0 * (Y_ipcw - Q0))
    surv1 <- mean(Q1 + ipw_ind1 * (Y_ipcw - Q1))

    # Clip to [0,1]: small deviations can occur with finite samples
    data.frame(
      time   = t,
      A      = c(0L, 1L),
      surv   = pmin(pmax(c(surv0, surv1), 0), 1),
      method = method.label
    )
  })

  surv.result <- do.call(rbind, out.list)

  # Risk difference RD(t) = S^1(t) - S^0(t)
  rd.result <- surv.result %>%
    pivot_wider(id_cols = time, names_from = A, values_from = surv,
                names_prefix = "S") %>%
    rename(S0 = S0, S1 = S1) %>%
    mutate(RD = S1 - S0, method = method.label)

  list(surv = surv.result, rd = rd.result)
}


# ---- Performance reporting helper --------------------------------------------

report_performance <- function(result, true_S0, true_S1, rd_wei,
                                t.eval, t.perf = c(1, 5, 10)) {

  lbl <- unique(result$surv$method)
  cat("\n", lbl, "\n", strrep("-", nchar(lbl)), "\n", sep = "")
  cat(sprintf("  %-4s  %-3s  %-8s  %-8s  %s\n",
              "t", "A", "Est", "True", "RelBias%"))

  for (t_pt in t.perf) {
    for (a in c(0L, 1L)) {
      sub      <- dplyr::filter(result$surv, A == a)
      idx_est  <- which.min(abs(sub$time - t_pt))
      idx_true <- which.min(abs(t.eval - t_pt))
      est      <- sub$surv[idx_est]
      truth    <- if (a == 0L) true_S0[idx_true] else true_S1[idx_true]
      rb       <- 100 * (est - truth) / truth
      cat(sprintf("  %-4g  A=%d  %-8.4f  %-8.4f  %+.2f%%\n",
                  t_pt, a, est, truth, rb))
    }

    rd_sub   <- dplyr::filter(result$rd, abs(time - t_pt) == min(abs(result$rd$time - t_pt)))
    rd_true  <- rd_wei[which.min(abs(t.eval - t_pt))]
    cat(sprintf("        RD    %-8.4f  %-8.4f\n",
                rd_sub$RD[1], rd_true))
  }
}


# ==============================================================================
# Run AIPTW across all protocol misspecification scenarios
# ==============================================================================

message("\n=== Running AIPTW for all misspecification scenarios ===\n")

all.results <- lapply(names(scenarios), function(sc_name) {
  sc <- scenarios[[sc_name]]
  message("  [", sc_name, "]  PS: ", sc$ps, "  |  Outcome: ", sc$out)

  run_aiptw_discrete(
    ps.spec.str  = ps.specs[sc$ps],
    out.spec.str = out.specs[sc$out],
    surv.df          = surv.df,
    surv.long.df     = surv.long.df,
    baseline.covs    = baseline.covs,
    interval_mapping = interval_mapping,
    G_t_vec          = G_t_vec,
    t.eval           = t.eval,
    method.label     = paste0("AIPTW (", sc$label, ")")
  )
})
names(all.results) <- names(scenarios)

# Print performance for every scenario
invisible(lapply(all.results, report_performance,
                 true_S0 = surv0_wei, true_S1 = surv1_wei,
                 rd_wei  = rd_wei,    t.eval  = t.eval))


# ==============================================================================
# Visualisation 1: Marginal survival curves by scenario
# ==============================================================================

plot.surv.df <- bind_rows(
  true.surv.df,
  lapply(all.results, function(r) r$surv)
) %>%
  mutate(A = factor(A, levels = c(0L, 1L),
                    labels = c("Control (A=0)", "Treatment (A=1)")))

ggplot(plot.surv.df,
       aes(x = time, y = surv, colour = method, linetype = method)) +
  geom_step(linewidth = 0.6) +
  facet_wrap(~A) +
  scale_colour_manual(
    values = c(
      "True Weibull"                = "black",
      "AIPTW (Both correct)"        = "#1A9850",
      "AIPTW (PS: missing W)"       = "#D73027",
      "AIPTW (Outcome: missing O)"  = "#FC8D59",
      "AIPTW (Outcome: wrong form)" = "#91BFDB",
      "AIPTW (PS: wrong form)"      = "#4575B4"
    ),
    name = "Method / Scenario"
  ) +
  scale_linetype_manual(
    values = c(
      "True Weibull"                = "solid",
      "AIPTW (Both correct)"        = "solid",
      "AIPTW (PS: missing W)"       = "dashed",
      "AIPTW (Outcome: missing O)"  = "dashed",
      "AIPTW (Outcome: wrong form)" = "dotdash",
      "AIPTW (PS: wrong form)"      = "dotdash"
    ),
    name = "Method / Scenario"
  ) +
  ylim(0, 1) + xlim(0, admin.cens) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 8),
        legend.key.width = unit(1.5, "cm")) +
  guides(colour   = guide_legend(nrow = 3),
         linetype = guide_legend(nrow = 3)) +
  xlab("Time (years)") + ylab("Marginal Survival Probability") +
  ggtitle("AIPTW: Marginal Survival Curves by Misspecification Scenario",
          subtitle = sprintf("N=%d, DGM: PH scenario", simN))


# ==============================================================================
# Visualisation 2: Risk difference RD(t) at t = 1, 5, 10
# ==============================================================================

rd.perf.df <- bind_rows(lapply(names(all.results), function(sc_name) {
  r   <- all.results[[sc_name]]
  sub <- r$rd %>% dplyr::filter(time %in% t.perf)
  sub$scenario <- scenarios[[sc_name]]$label
  sub
})) %>%
  mutate(t_label = factor(paste("t =", time), levels = paste("t =", t.perf)))

true.rd.df <- data.frame(
  time    = t.perf,
  true_RD = rd_wei[sapply(t.perf, function(t) which.min(abs(t.eval - t)))],
  t_label = factor(paste("t =", t.perf), levels = paste("t =", t.perf))
)

ggplot() +
  geom_col(data = rd.perf.df,
           aes(x = scenario, y = RD, fill = scenario),
           position = "dodge", alpha = 0.8) +
  geom_hline(data = true.rd.df,
             aes(yintercept = true_RD), lty = 2, colour = "black", linewidth = 0.8) +
  facet_wrap(~t_label) +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 40, hjust = 1, size = 8),
        legend.position = "none") +
  xlab("Misspecification scenario") +
  ylab(expression(hat(RD)(t) == hat(S)^1(t) - hat(S)^0(t))) +
  ggtitle("AIPTW Risk Differences by Scenario  (dashed = truth)")
