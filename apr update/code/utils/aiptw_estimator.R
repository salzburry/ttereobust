## utils/aiptw_estimator.R
## Pure-function AIPTW survival estimator for the discrete-time / logit-PLR
## formulation (Bang & Robins 2005; Funk et al. 2011; Kurz 2022; Gabriel et
## al. 2024). Designed to be called many times from a replication driver and
## from a bootstrap loop without any top-level side effects.
##
## Inputs:
##   surv.df      one row per individual; must contain columns
##                  id, eventtime, event, A, L1..L6, L1sq, L2sq, W, O
##   ps_spec      RHS string for the propensity model, e.g. "L1 + L2 + ... + W"
##   out_spec     RHS string for the discrete-time hazard model, e.g.
##                  "A + L1 + L1sq + L2 + L2sq + L3 + L4 + L5 + L6 + O + W"
##   t_eval       vector of evaluation time points (default 1, 5, 10)
##   admin.cens   administrative censoring boundary (default 10)
##   rescale_time discrete-time interval width (default 1/4)
##   ps_trim      P(A|L) clipping range (default c(0.01, 0.99))
##
## Returns: data.frame with columns t, S0, S1, RD on RAW scale (no clipping).
##
## Notes:
##   - IPW does NOT enter the outcome PLR; it enters only via the augmentation.
##   - Logit link is required for the AIPTW algebraic equivalence (Gabriel
##     et al. 2024); cloglog would break double-robustness.
##   - IPCW uses the left-limit G(t-) and >= so that admin-censored survivors
##     at t = admin.cens are correctly retained.
##   - Returns NA on any internal error, so a bootstrap loop can keep going.

aiptw_estimate <- function(surv.df,
                            ps_spec,
                            out_spec,
                            t_eval        = c(1, 5, 10),
                            admin.cens    = 10,
                            rescale_time  = 1/4,
                            ps_trim       = c(0.01, 0.99)) {

  tryCatch({

    interval_ends <- seq(rescale_time, admin.cens, by = rescale_time)
    split_cuts    <- interval_ends[-length(interval_ends)]
    interval_mapping <- data.frame(
      time_period = seq_along(interval_ends),
      time        = interval_ends
    )

    # Long-format (person-period)
    surv.long.df <- survival::survSplit(
      survival::Surv(eventtime, event) ~ .,
      data    = surv.df,
      cut     = split_cuts,
      episode = "time_period"
    )

    baseline.covs <- dplyr::select(
      surv.df, id, L1, L1sq, L2, L2sq, L3, L4, L5, L6, O, W
    )

    # Censoring KM at left limit G(t-)
    cens.km <- survival::survfit(
      survival::Surv(eventtime, 1L - event) ~ 1, data = surv.df
    )
    G_tminus <- summary(
      cens.km, times = pmax(t_eval - 1e-8, 0), extend = TRUE
    )$surv
    G_tminus <- pmax(G_tminus, 1e-6)

    # Propensity score model
    ps.mod <- stats::glm(
      stats::as.formula(paste("A ~", ps_spec)),
      family = "binomial", data = surv.df
    )
    pi_1 <- pmin(pmax(stats::predict(ps.mod, type = "response"),
                      ps_trim[1]), ps_trim[2])

    # Unstabilised indicator weights for the augmentation term
    ipw_ind1 <- ifelse(surv.df$A == 1L, 1 / pi_1,       0)
    ipw_ind0 <- ifelse(surv.df$A == 0L, 1 / (1 - pi_1), 0)

    # Outcome model: unweighted PLR with canonical logit link
    plr.formula <- stats::as.formula(
      paste0("event ~ ", out_spec, " + as.factor(time_period)")
    )
    plr.mod <- stats::glm(
      plr.formula, data = surv.long.df,
      family = stats::binomial(link = "logit")
    )

    # Counterfactual prediction grid
    time_periods <- sort(unique(surv.long.df$time_period))
    pred.df <- expand.grid(
      id          = baseline.covs$id,
      time_period = time_periods,
      A           = c(0L, 1L),
      KEEP.OUT.ATTRS = FALSE
    )
    pred.df <- dplyr::left_join(pred.df, baseline.covs, by = "id")
    pred.df$event <- 0L
    pred.df$haz   <- stats::predict(plr.mod, newdata = pred.df,
                                    type = "response")

    pred.df <- pred.df %>%
      dplyr::arrange(id, A, time_period) %>%
      dplyr::group_by(id, A) %>%
      dplyr::mutate(csurv = cumprod(1 - haz)) %>%
      dplyr::ungroup()

    obs_time <- surv.df$eventtime

    out <- lapply(seq_along(t_eval), function(k) {
      t <- t_eval[k]
      if (t == 0) return(data.frame(t = 0, S0 = 1, S1 = 1, RD = 0))

      tp <- interval_mapping$time_period[
        max(which(interval_mapping$time <= t))
      ]

      Q_at_t <- pred.df %>%
        dplyr::filter(time_period == tp) %>%
        dplyr::select(id, A, csurv)

      Q0 <- dplyr::left_join(
        dplyr::select(surv.df, id),
        dplyr::filter(Q_at_t, A == 0L), by = "id"
      )$csurv
      Q1 <- dplyr::left_join(
        dplyr::select(surv.df, id),
        dplyr::filter(Q_at_t, A == 1L), by = "id"
      )$csurv

      Y_ipcw <- as.numeric(obs_time >= t) / G_tminus[k]

      S0 <- mean(Q0 + ipw_ind0 * (Y_ipcw - Q0))
      S1 <- mean(Q1 + ipw_ind1 * (Y_ipcw - Q1))
      data.frame(t = t, S0 = S0, S1 = S1, RD = S1 - S0)
    })

    out_df <- do.call(rbind, out)
    out_df$status    <- "ok"
    out_df$error_msg <- NA_character_
    out_df

  }, error = function(e) {
    # Surface the failure to the caller instead of silently returning NA.
    # The replication driver and aggregator both inspect $status so a
    # scenario's success rate can be reported separately from n_reps.
    data.frame(
      t         = t_eval,
      S0        = NA_real_,
      S1        = NA_real_,
      RD        = NA_real_,
      status    = "error",
      error_msg = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}


## Bootstrap wrapper -- non-parametric individual-level resampling.
##
## Returns a data.frame with one row per (t, target) and columns
##   t, target, se, ci_lo, ci_hi   (target in {"S0", "S1", "RD"})

aiptw_bootstrap <- function(surv.df, ps_spec, out_spec,
                             t_eval       = c(1, 5, 10),
                             admin.cens   = 10,
                             rescale_time = 1/4,
                             ps_trim      = c(0.01, 0.99),
                             B            = 200,
                             ci_level     = 0.95) {

  N <- nrow(surv.df)
  alpha <- (1 - ci_level) / 2

  # Each bootstrap returns a (length(t_eval) x 3) data frame; stack into
  # a long table keyed by (b, t, target).
  boot_long <- vector("list", B)
  for (b in seq_len(B)) {
    idx     <- sample.int(N, N, replace = TRUE)
    boot.df <- surv.df[idx, ]
    boot.df$id <- seq_len(N)            # re-id so survSplit etc. behave
    est <- aiptw_estimate(boot.df, ps_spec, out_spec, t_eval,
                          admin.cens, rescale_time, ps_trim)
    boot_long[[b]] <- data.frame(
      b      = b,
      t      = rep(est$t, 3),
      target = rep(c("S0", "S1", "RD"), each = nrow(est)),
      val    = c(est$S0, est$S1, est$RD)
    )
  }
  boot_long_df <- do.call(rbind, boot_long)

  boot_long_df %>%
    dplyr::group_by(t, target) %>%
    dplyr::summarise(
      n_boot_ok = sum(!is.na(val)),
      se        = stats::sd(val, na.rm = TRUE),
      ci_lo     = stats::quantile(val, alpha,     na.rm = TRUE, names = FALSE),
      ci_hi     = stats::quantile(val, 1 - alpha, na.rm = TRUE, names = FALSE),
      .groups   = "drop"
    )
}
