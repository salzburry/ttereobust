## utils/aiptw_cox_estimator.R
##
## AIPTW with a Cox PH outcome model, wrapping riskRegression::ate(). Sister
## to utils/aiptw_estimator.R (discrete-time logit-PLR variant).
##
## Why this exists alongside aiptw_estimate():
##   - The methods document specifies discrete-time canonical-link PLR for
##     the AIPTW double-robustness algebraic-equivalence result (Gabriel et
##     al. 2024). That is what aiptw_estimate() does.
##   - In practice many applied analysts reach for a Cox outcome model. This
##     file provides AIPTW with that choice so the project can compare both
##     formulations on the same protocol grid.
##   - The Cox-AIPTW variant is implemented via riskRegression::ate(), the
##     well-maintained CRAN implementation (also wrapped by
##     adjustedCurves::surv_aiptw). We call ate() directly to avoid an extra
##     layer of indirection.
##
## Function signatures match aiptw_estimate() / aiptw_bootstrap() so the
## driver can dispatch on --method without estimator-specific glue.
##
## Bootstrap: produces the same percentile-based SE / CI columns as the PLR
## arm. Additionally, the point-estimate function returns ate()'s
## influence-function-based SEs and analytical CIs in *_if_se / *_if_ci_*
## columns so the report can compare bootstrap vs IF inference on the
## Cox arm.

suppressPackageStartupMessages({
  if (!requireNamespace("riskRegression", quietly = TRUE)) {
    stop("Package 'riskRegression' is required for the AIPTW-Cox arm. ",
         "Install with install.packages('riskRegression').")
  }
})


## Internal point-estimate worker. with_if_se = FALSE skips the
## influence-function variance computation, which dominates ate()'s wall
## time; bootstrap calls set this to FALSE since bootstrap uses percentile
## intervals from refit estimates rather than IF-based SEs.
aiptw_cox_estimate_inner <- function(surv.df, ps_spec, out_spec,
                                       t_eval, admin.cens,
                                       ps_trim    = c(0.01, 0.99),
                                       with_if_se = TRUE) {

  # Work on a copy so we can add A_fct without mutating the caller's df.
  df          <- surv.df
  df$A_fct    <- factor(df$A, levels = c(0L, 1L))
  out_spec_cox <- gsub("\\bA\\b", "A_fct", out_spec)

  # Treatment model — logistic regression. x=TRUE/y=TRUE keep design /
  # response on the fitted object so ate() can recompute internals.
  ps_fit <- glm(stats::as.formula(paste("A_fct ~", ps_spec)),
                family = stats::binomial(link = "logit"),
                data = df, x = TRUE, y = TRUE)

  # Outcome model — Cox PH on A_fct + covariates.
  cox_formula <- stats::as.formula(
    paste("survival::Surv(eventtime, event) ~", out_spec_cox)
  )
  cox_fit <- survival::coxph(cox_formula, data = df,
                              x = TRUE, y = TRUE, ties = "breslow")

  # ate() requires strictly positive evaluation times.
  t_eval_pos <- t_eval[t_eval > 0]
  if (length(t_eval_pos) == 0L) {
    stop("aiptw_cox_estimate(): t_eval must contain at least one t > 0")
  }

  ate_fit <- riskRegression::ate(
    event     = cox_fit,
    treatment = ps_fit,
    data      = df,
    times     = t_eval_pos,
    cause     = 1,
    estimator = "AIPTW",
    se        = with_if_se,
    band      = FALSE,
    verbose   = FALSE
  )

  # Convert risk -> survival (S = 1 - risk). RD on survival = -RD on risk.
  mr <- as.data.frame(ate_fit$meanRisk)
  dr <- as.data.frame(ate_fit$diffRisk)

  out <- lapply(t_eval, function(t) {
    if (t == 0) {
      return(data.frame(
        t           = 0,
        S0          = 1, S1 = 1, RD = 0,
        S0_if_se    = 0, S1_if_se = 0, RD_if_se = 0,
        S0_if_ci_lo = 1, S0_if_ci_hi = 1,
        S1_if_ci_lo = 1, S1_if_ci_hi = 1,
        RD_if_ci_lo = 0, RD_if_ci_hi = 0
      ))
    }
    r0 <- mr[mr$time == t & mr$treatment == 0, ]
    r1 <- mr[mr$time == t & mr$treatment == 1, ]
    drow <- dr[dr$time == t, ]
    if (nrow(r0) == 0L || nrow(r1) == 0L || nrow(drow) == 0L) {
      stop("ate() did not return rows for t = ", t)
    }
    se_S0 <- if (with_if_se && !is.null(r0$se))    r0$se    else NA_real_
    se_S1 <- if (with_if_se && !is.null(r1$se))    r1$se    else NA_real_
    se_RD <- if (with_if_se && !is.null(drow$se))  drow$se  else NA_real_
    lo_S0 <- if (with_if_se && !is.null(r0$lower)) 1 - r0$upper else NA_real_
    hi_S0 <- if (with_if_se && !is.null(r0$upper)) 1 - r0$lower else NA_real_
    lo_S1 <- if (with_if_se && !is.null(r1$lower)) 1 - r1$upper else NA_real_
    hi_S1 <- if (with_if_se && !is.null(r1$upper)) 1 - r1$lower else NA_real_
    lo_RD <- if (with_if_se && !is.null(drow$lower)) -drow$upper else NA_real_
    hi_RD <- if (with_if_se && !is.null(drow$upper)) -drow$lower else NA_real_
    data.frame(
      t           = t,
      S0          = 1 - r0$estimate,
      S1          = 1 - r1$estimate,
      RD          = -drow$estimate,
      S0_if_se    = se_S0,
      S1_if_se    = se_S1,
      RD_if_se    = se_RD,
      S0_if_ci_lo = lo_S0, S0_if_ci_hi = hi_S0,
      S1_if_ci_lo = lo_S1, S1_if_ci_hi = hi_S1,
      RD_if_ci_lo = lo_RD, RD_if_ci_hi = hi_RD
    )
  })
  do.call(rbind, out)
}


## Public point-estimate API. tryCatch wrapper mirrors aiptw_estimate()
## so the replication driver sees consistent status / error_msg columns
## across both arms.
aiptw_cox_estimate <- function(surv.df, ps_spec, out_spec,
                                t_eval       = c(1, 5, 10),
                                admin.cens   = 10,
                                rescale_time = NULL,    # accepted for API parity, unused
                                ps_trim      = c(0.01, 0.99)) {
  tryCatch({
    out_df <- aiptw_cox_estimate_inner(surv.df, ps_spec, out_spec,
                                        t_eval, admin.cens,
                                        ps_trim, with_if_se = TRUE)
    out_df$status    <- "ok"
    out_df$error_msg <- NA_character_
    out_df
  }, error = function(e) {
    nas <- rep(NA_real_, length(t_eval))
    data.frame(
      t           = t_eval,
      S0          = nas, S1 = nas, RD = nas,
      S0_if_se    = nas, S1_if_se = nas, RD_if_se = nas,
      S0_if_ci_lo = nas, S0_if_ci_hi = nas,
      S1_if_ci_lo = nas, S1_if_ci_hi = nas,
      RD_if_ci_lo = nas, RD_if_ci_hi = nas,
      status      = "error",
      error_msg   = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}


## Bootstrap wrapper — non-parametric individual-level resampling, same
## design as aiptw_bootstrap(). Returns one row per (t, target) with the
## same n_boot_total / n_boot_ok / se / ci_lo / ci_hi / boot_errors
## schema, so the rep CSV layer is identical between arms.
aiptw_cox_bootstrap <- function(surv.df, ps_spec, out_spec,
                                 t_eval       = c(1, 5, 10),
                                 admin.cens   = 10,
                                 rescale_time = NULL,
                                 ps_trim      = c(0.01, 0.99),
                                 B            = 200,
                                 ci_level     = 0.95) {

  stopifnot(is.numeric(B), length(B) == 1L, B >= 0L)
  N     <- nrow(surv.df)
  alpha <- (1 - ci_level) / 2

  if (B == 0L) {
    return(data.frame(
      t             = rep(t_eval, 3),
      target        = rep(c("S0", "S1", "RD"), each = length(t_eval)),
      n_boot_total  = 0L, n_boot_ok = 0L, n_boot_failed = 0L,
      se            = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      boot_errors   = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  boot_long <- vector("list", B)
  boot_errs <- character(0)
  for (b in seq_len(B)) {
    idx      <- sample.int(N, N, replace = TRUE)
    boot.df  <- surv.df[idx, ]
    boot.df$id <- seq_len(N)
    est <- tryCatch(
      aiptw_cox_estimate_inner(boot.df, ps_spec, out_spec,
                                t_eval, admin.cens, ps_trim,
                                with_if_se = FALSE),
      error = function(e) {
        msg <- conditionMessage(e)
        if (!msg %in% boot_errs) boot_errs <<- c(boot_errs, msg)
        nas <- rep(NA_real_, length(t_eval))
        data.frame(
          t = t_eval, S0 = nas, S1 = nas, RD = nas,
          S0_if_se = nas, S1_if_se = nas, RD_if_se = nas,
          S0_if_ci_lo = nas, S0_if_ci_hi = nas,
          S1_if_ci_lo = nas, S1_if_ci_hi = nas,
          RD_if_ci_lo = nas, RD_if_ci_hi = nas
        )
      })
    boot_long[[b]] <- data.frame(
      b      = b,
      t      = rep(est$t, 3),
      target = rep(c("S0", "S1", "RD"), each = nrow(est)),
      val    = c(est$S0, est$S1, est$RD)
    )
  }
  boot_long_df <- do.call(rbind, boot_long)
  boot_err_str <- if (length(boot_errs) == 0L) NA_character_ else
                   paste(boot_errs, collapse = " | ")

  safe_sd <- function(v) {
    if (sum(!is.na(v)) < 2L) NA_real_ else stats::sd(v, na.rm = TRUE)
  }
  safe_q <- function(v, p) {
    if (sum(!is.na(v)) < 1L) NA_real_ else
      suppressWarnings(stats::quantile(v, p, na.rm = TRUE, names = FALSE))
  }

  res <- aggregate(boot_long_df$val,
                    by  = list(t = boot_long_df$t,
                               target = boot_long_df$target),
                    FUN = function(v) {
                      c(n_boot_ok = sum(!is.na(v)),
                        se        = safe_sd(v),
                        ci_lo     = safe_q(v, alpha),
                        ci_hi     = safe_q(v, 1 - alpha))
                    })
  out <- data.frame(
    t             = res$t,
    target        = res$target,
    n_boot_total  = B,
    n_boot_ok     = as.integer(res$x[, "n_boot_ok"]),
    n_boot_failed = B - as.integer(res$x[, "n_boot_ok"]),
    se            = res$x[, "se"],
    ci_lo         = res$x[, "ci_lo"],
    ci_hi         = res$x[, "ci_hi"],
    boot_errors   = boot_err_str,
    stringsAsFactors = FALSE
  )
  out
}
