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
##   pi_i(a)      : P(A_i=a | L_i) from the propensity score model (trimmed)
##   Ytilde_i(t)  : I(Ttilde_i >= t) / G(t-)  — IPCW-corrected observed outcome
##   G(t-)        : left-limit of KM censoring survival (handles admin. censoring)
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
## Outcome model family note:
##   This implementation uses a logit-link PLR as the working outcome model.
##   The underlying DGM is a continuous-time Weibull process for which a
##   cloglog discrete-time hazard model is the exact interval hazard.
##   Logit is chosen here because the canonical link is required for the
##   AIPTW algebraic equivalence result (Gabriel et al. 2024); it is a
##   working model rather than the true conditional hazard model.
##   The "covariate sets correct" scenarios below refer to correct covariate
##   inclusion, not correct link-function family.
##
## Misspecification scenarios follow protocol Section 6.4:
##   Each scenario isolates misspecification in ONE nuisance model only.
##
## DGM note on W and O:
##   To preserve comparability with analysis.R / gen_truth.R, the truth
##   formula here follows those files exactly: both W and O enter the
##   outcome hazard (coeff.W != 0 and coeff.O != 0 in the linear predictor).
##   The distinction used in the misspecification scenarios is that O cannot
##   appear in the PS model (it is not used for treatment assignment) while
##   W can appear in both, consistent with the truth formula.
##   The get_params comment "coeff.W will only appear in exposure model" is
##   inconsistent with the truth code and is flagged for future resolution.
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
library(grid)     # required for unit() in legend.key.width
library(MASS)

source("utils/sim_data.R")


# ---- Parameters --------------------------------------------------------------
simN <- 2500   # N per protocol Section 6.2

# Uncomment desired DGM scenario:
# source("get_params_waning.R")
# source("get_params_delayed.R")
source("get_params_ph.R")

rescale_time <- 1/4                      # quarter-year discrete intervals
t.eval       <- seq(0, admin.cens, 0.1)  # fine evaluation grid
t.perf       <- c(1, 5, 10)             # protocol performance time points

# Discrete-time interval setup.
# interval_ends: right endpoints of each interval  → 0.25, 0.50, ..., 10
# split_cuts:    interior cut times passed to survSplit (excludes admin.cens
#                endpoint so survSplit does not create an empty final interval)
interval_ends <- seq(rescale_time, admin.cens, by = rescale_time)
split_cuts    <- interval_ends[-length(interval_ends)]

interval_mapping <- data.frame(
  time_period = seq_along(interval_ends),
  time        = interval_ends
)


# ---- Propensity score model specifications -----------------------------------
# Kept separate from outcome model specs.  Sharing a single index (as in the
# original analysis.R cov.index) prevents representing "misspecify one model
# only" scenarios.

ps.specs <- c(
  correct  = "L1sq + L2sq + L3 + L4 + L5 + L6 + W",   # correct covariate set
  no_W     = "L1sq + L2sq + L3 + L4 + L5 + L6",         # missing W
  wrong_ff = "L1 + L2 + L3 + L4 + L5 + L6 + W",         # wrong functional form (L1,L2 linear)
  heavy    = "L3 + L4 + L5"                               # heavy misspecification
)


# ---- Outcome model specifications --------------------------------------------
# Correct covariate set includes O and W consistent with the truth formula.
# See DGM note in file header.

out.specs <- c(
  correct  = "A + L1sq + L2sq + L3 + L4 + L5 + L6 + O + W",
  no_O     = "A + L1sq + L2sq + L3 + L4 + L5 + L6 + W",          # missing O
  wrong_ff = "A + L1 + L2 + L3 + L4 + L5 + L6 + O + W",          # wrong functional form
  heavy    = "A + L3 + L4 + L5"                                     # heavy misspecification
)


# ---- Protocol misspecification scenarios (Section 6.4) ----------------------
# Every scenario changes exactly ONE nuisance model.

scenarios <- list(
  both_cvs_correct = list(ps = "correct",  out = "correct",
                           label = "Both covariate sets correct"),
  miss_W_ps_only   = list(ps = "no_W",     out = "correct",
                           label = "PS: missing W"),
  miss_O_out_only  = list(ps = "correct",  out = "no_O",
                           label = "Outcome: missing O"),
  wrong_ff_out     = list(ps = "correct",  out = "wrong_ff",
                           label = "Outcome: wrong form"),
  wrong_ff_ps      = list(ps = "wrong_ff", out = "correct",
                           label = "PS: wrong form"),
  heavy_ps         = list(ps = "heavy",    out = "correct",
                           label = "PS: heavy misspec"),
  heavy_out        = list(ps = "correct",  out = "heavy",
                           label = "Outcome: heavy misspec")
)


# ---- True marginal survival curves -------------------------------------------
# Analytical Weibull g-computation over a large population.
# Identical to analysis.R / gen_truth.R so benchmarks are comparable.

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


# ---- Analysis dataset --------------------------------------------------------

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

surv.df <- sim.df$data   # one row per individual

# Long (person-period) format using corrected interior cut points.
# split_cuts excludes 0 and admin.cens, so every interval is non-empty and
# the time_period labels align cleanly with interval_mapping$time.
surv.long.df <- survSplit(
  Surv(eventtime, event) ~ .,
  data    = surv.df,
  cut     = split_cuts,
  episode = "time_period"
)

baseline.covs <- dplyr::select(surv.df, id, L1, L1sq, L2, L2sq, L3, L4, L5, L6, O, W)


# ---- Censoring KM for IPCW ---------------------------------------------------
# Reverse the event indicator: censored observations become "events" for this KM.
# G(t) = P(C > t) under independent censoring.
#
# IPCW uses the left-limit G(t-) evaluated just before t.  This matters at the
# administrative censoring boundary (t = admin.cens = 10): at exactly t = 10,
# G(10) can be 0 or near-0 due to the mass of admin-censoring, making G(t)
# an unstable denominator.  Using G(t-) = G(10 - eps) avoids this.
# The outcome indicator is also changed to >= (individual observed to survive
# through the start of the last interval is counted as surviving past t).

cens.km      <- survfit(Surv(eventtime, 1L - event) ~ 1, data = surv.df)
G_tminus_vec <- summary(cens.km,
                         times  = pmax(t.eval - 1e-8, 0),
                         extend = TRUE)$surv
G_tminus_vec <- pmax(G_tminus_vec, 1e-6)   # guard against exact zero


# ---- Core AIPTW function -----------------------------------------------------
#
# Inputs:  one PS spec string, one outcome spec string, pre-built data objects
# Output:  list(surv = survival-curve data frame, rd = risk-difference data frame)
#
# Steps:
#   1. Fit PS model → trimmed pi_1 = P(A=1|L)
#   2. Fit UNWEIGHTED PLR (logit) → conditional discrete hazard h_k(a, L_i)
#   3. Build full counterfactual grid → Q_i(a,t) = cumprod(1 - h_k(a,L_i))
#   4. At each t: form IPCW outcome Ytilde_i(t) = I(T_i >= t) / G(t-)
#   5. Apply AIPTW: S^a(t) = mean_i{Q_i + [I(A=a)/pi_i] * [Ytilde - Q_i]}

run_aiptw_discrete <- function(ps.spec.str, out.spec.str,
                                surv.df, surv.long.df,
                                baseline.covs, interval_mapping,
                                G_tminus_vec, t.eval,
                                method.label) {

  # --- Step 1: Propensity score model ---
  ps.mod <- glm(as.formula(paste("A ~", ps.spec.str)),
                family = "binomial", data = surv.df)

  # Trim predicted probabilities away from 0/1 to prevent extreme weights
  pi_1 <- pmin(pmax(predict(ps.mod, type = "response"), 0.01), 0.99)

  # Unstabilised IPW indicators  I(A_i=a) / P(A_i=a|L_i)
  # Used only in the augmentation term; the outcome model is unweighted.
  ipw_ind1 <- ifelse(surv.df$A == 1L, 1 / pi_1,       0)
  ipw_ind0 <- ifelse(surv.df$A == 0L, 1 / (1 - pi_1), 0)

  # --- Step 2: Unweighted PLR outcome model (canonical logit link) ---
  # IPW does NOT enter here — it enters only through the augmentation term.
  # logit is the canonical link for the Bernoulli, which is the link required
  # for the AIPTW algebraic equivalence result (Gabriel et al. 2024).
  plr.formula <- as.formula(
    paste0("event ~ ", out.spec.str, " + as.factor(time_period)")
  )
  plr.mod <- glm(plr.formula, data = surv.long.df,
                 family = binomial(link = "logit"))

  # --- Step 3: Counterfactual prediction grid ---
  # Every individual x every time period x A in {0, 1}.
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

  # Individual cumulative survival Q_i(a, t_k) = prod_{j<=k} (1 - h_j(a,L_i))
  pred.df <- pred.df %>%
    arrange(id, A, time_period) %>%
    group_by(id, A) %>%
    mutate(csurv = cumprod(1 - haz)) %>%
    ungroup()

  obs_time <- surv.df$eventtime

  # --- Steps 4 & 5: AIPTW at each evaluation time point ---
  out.list <- lapply(seq_along(t.eval), function(k) {

    t        <- t.eval[k]
    G_tminus <- G_tminus_vec[k]

    if (t == 0) {
      return(data.frame(time = 0, A = c(0L, 1L), surv = 1,
                        method = method.label))
    }

    # Last completed interval at or before t
    tp <- interval_mapping$time_period[max(which(interval_mapping$time <= t))]

    # Q_i(a, t): align to surv.df row order via id join
    Q_at_t <- pred.df %>%
      dplyr::filter(time_period == tp) %>%
      dplyr::select(id, A, csurv)

    Q0 <- left_join(dplyr::select(surv.df, id),
                    dplyr::filter(Q_at_t, A == 0L), by = "id")$csurv
    Q1 <- left_join(dplyr::select(surv.df, id),
                    dplyr::filter(Q_at_t, A == 1L), by = "id")$csurv

    # IPCW outcome using >= and G(t-):
    #   I(T_tilde_i >= t) / G(t-)
    # Using >= (not >) ensures individuals admin-censored exactly at t = 10
    # are correctly treated as "observed to survive through t".
    # G(t-) avoids division by near-zero at the admin-censoring boundary.
    Y_ipcw <- as.numeric(obs_time >= t) / G_tminus

    # AIPTW formula
    surv0 <- mean(Q0 + ipw_ind0 * (Y_ipcw - Q0))
    surv1 <- mean(Q1 + ipw_ind1 * (Y_ipcw - Q1))

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
    mutate(RD = S1 - S0, method = method.label)

  list(surv = surv.result, rd = rd.result)
}


# ---- Performance reporting ---------------------------------------------------

report_performance <- function(result, true_S0, true_S1, rd_wei,
                                t.eval, t.perf = c(1, 5, 10)) {

  lbl <- unique(result$surv$method)
  cat("\n", lbl, "\n", strrep("-", nchar(lbl)), "\n", sep = "")
  cat(sprintf("  %-4s  %-3s  %-8s  %-8s  %s\n", "t", "A", "Est", "True", "RelBias%"))

  for (t_pt in t.perf) {
    for (a in c(0L, 1L)) {
      sub      <- dplyr::filter(result$surv, A == a)
      idx_est  <- which.min(abs(sub$time - t_pt))
      idx_true <- which.min(abs(t.eval - t_pt))
      est      <- sub$surv[idx_est]
      truth    <- if (a == 0L) true_S0[idx_true] else true_S1[idx_true]
      cat(sprintf("  %-4g  A=%d  %-8.4f  %-8.4f  %+.2f%%\n",
                  t_pt, a, est, truth, 100 * (est - truth) / truth))
    }
    rd_sub  <- result$rd[which.min(abs(result$rd$time - t_pt)), ]
    rd_true <- rd_wei[which.min(abs(t.eval - t_pt))]
    cat(sprintf("        RD    %-8.4f  %-8.4f\n", rd_sub$RD, rd_true))
  }
}


# ==============================================================================
# Run AIPTW across all misspecification scenarios
# ==============================================================================

message("\n=== Running AIPTW for all misspecification scenarios ===\n")

all.results <- lapply(names(scenarios), function(sc_name) {
  sc <- scenarios[[sc_name]]
  message("  [", sc_name, "]  PS: ", sc$ps, "  |  Outcome: ", sc$out)
  run_aiptw_discrete(
    ps.spec.str      = ps.specs[sc$ps],
    out.spec.str     = out.specs[sc$out],
    surv.df          = surv.df,
    surv.long.df     = surv.long.df,
    baseline.covs    = baseline.covs,
    interval_mapping = interval_mapping,
    G_tminus_vec     = G_tminus_vec,
    t.eval           = t.eval,
    method.label     = paste0("AIPTW (", sc$label, ")")
  )
})
names(all.results) <- names(scenarios)

invisible(lapply(all.results, report_performance,
                 true_S0 = surv0_wei, true_S1 = surv1_wei,
                 rd_wei  = rd_wei,    t.eval  = t.eval))


# ==============================================================================
# Visualisation 1: Marginal survival curves
# ==============================================================================

sc_colours <- c(
  "True Weibull"                         = "black",
  "AIPTW (Both covariate sets correct)"  = "#1A9850",
  "AIPTW (PS: missing W)"                = "#D73027",
  "AIPTW (Outcome: missing O)"           = "#FC8D59",
  "AIPTW (Outcome: wrong form)"          = "#91BFDB",
  "AIPTW (PS: wrong form)"               = "#4575B4",
  "AIPTW (PS: heavy misspec)"            = "#762A83",
  "AIPTW (Outcome: heavy misspec)"       = "#E7298A"
)

sc_linetypes <- c(
  "True Weibull"                         = "solid",
  "AIPTW (Both covariate sets correct)"  = "solid",
  "AIPTW (PS: missing W)"                = "dashed",
  "AIPTW (Outcome: missing O)"           = "dashed",
  "AIPTW (Outcome: wrong form)"          = "dotdash",
  "AIPTW (PS: wrong form)"               = "dotdash",
  "AIPTW (PS: heavy misspec)"            = "dotted",
  "AIPTW (Outcome: heavy misspec)"       = "dotted"
)

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
  scale_colour_manual(values = sc_colours,   name = "Method / Scenario") +
  scale_linetype_manual(values = sc_linetypes, name = "Method / Scenario") +
  ylim(0, 1) + xlim(0, admin.cens) +
  theme_classic(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 8),
        legend.key.width = unit(1.5, "cm")) +
  guides(colour   = guide_legend(nrow = 4),
         linetype = guide_legend(nrow = 4)) +
  xlab("Time (years)") + ylab("Marginal Survival Probability") +
  ggtitle("AIPTW: Marginal Survival Curves by Misspecification Scenario",
          subtitle = sprintf("N = %d  |  DGM: PH  |  PLR working model (logit link)", simN))


# ==============================================================================
# Visualisation 2: Risk differences RD(t) at t = 1, 5, 10
# ==============================================================================

rd.perf.df <- bind_rows(lapply(names(all.results), function(sc_name) {
  r        <- all.results[[sc_name]]
  sub      <- r$rd[r$rd$time %in% t.perf, ]
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
           alpha = 0.8) +
  geom_hline(data = true.rd.df,
             aes(yintercept = true_RD),
             lty = 2, colour = "black", linewidth = 0.8) +
  facet_wrap(~t_label) +
  theme_classic(base_size = 11) +
  theme(axis.text.x    = element_text(angle = 40, hjust = 1, size = 8),
        legend.position = "none") +
  xlab("Misspecification scenario") +
  ylab(expression(hat(RD)(t) == hat(S)^1(t) - hat(S)^0(t))) +
  ggtitle("AIPTW Risk Differences by Scenario  (dashed = truth)")
