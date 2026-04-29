# get working dir
wd <- getwd()

# load libraries
library(survival)
library(dplyr)
library(ggplot2)
library(flexsurv)      # provides flexsurvspline() and standsurv()
library(adjustedCurves)
library(MASS)
library(purrr)
library(broom)         # tidy()

source("utils/sim_data.R")
source("utils/scenarios.R")   # ps.specs / out.specs (single source of truth)

# Local positional aliases for backwards compatibility with this script's
# index-based access pattern (exposure.covs[ps.index] etc.). Names are
# kept so downstream paste0() output uses the same scenario labels as
# the AIPTW harness.
exposure.covs        <- ps.specs
exposure.covs.sceN   <- length(exposure.covs)
outcome.covs         <- out.specs
outcome.covs.sceN    <- length(outcome.covs)
 
# Simulate data - set params for analysis data
simN = 2500   # N per protocol Section 6.2
 
#source("get_params_waning.R")
#source("get_params_delayed.R")
source("get_params_ph.R")
 
## ---- Generate "true" survival curves:
 
true.df0 <- sim_surv_data(seed = this.seed,
                         N = 500000,
                         Lcovs.linear = N.Lcovs.linear,
                         Lcovs.sq = N.Lcovs.sq,
                         mu = mu, # Lcovs.linear + Lcovs.sq + 2 (W and O)
                         sigma = sigma,
                         alpha.L = alpha.L, # params for exposure model c(L_1, ..., L_{Lcovs.linear + Lcovs.sq})
                         alpha.W = alpha.W,
                         # coefficient / log-HRs for outcome model
                         coeff.A = coeff.A,
                         coeff.L = coeff.L,
                         coeff.Lsq = coeff.Lsq,
                         coeff.O = coeff.O, # outcome-only confounder
                         coeff.W = coeff.W, # confounder of BOTH exposure and outcome (matches sim_data.R DGM)
                         gamma.tte = gamma.tte,
                         lambda.tte = lambda.tte,
                         lambda.cens = lambda.cens,
                         admin.cens = admin.cens,
                         gen.truth = 0 # set to 1 to gen truth for trt arm
)
true.df1 <- sim_surv_data(seed = this.seed,
                         N = 500000,
                         Lcovs.linear = N.Lcovs.linear,
                         Lcovs.sq = N.Lcovs.sq,
                         mu = mu, # Lcovs.linear + Lcovs.sq + 2 (W and O)
                         sigma = sigma,
                         alpha.L = alpha.L, # params for exposure model c(L_1, ..., L_{Lcovs.linear + Lcovs.sq})
                         alpha.W = alpha.W,
                         # coefficient / log-HRs for outcome model
                         coeff.A = coeff.A,
                         coeff.L = coeff.L,
                         coeff.Lsq = coeff.Lsq,
                         coeff.O = coeff.O, # outcome-only confounder
                         coeff.W = coeff.W, # confounder of BOTH exposure and outcome (matches sim_data.R DGM)
                         gamma.tte = gamma.tte,
                         lambda.tte = lambda.tte,
                         lambda.cens = lambda.cens,
                         admin.cens = admin.cens,
                         gen.truth = 1 # set to 0 to gen truth for control
)
 
true.surv.df <- rbind(true.df0$data, true.df1$data)
cov.mat <- true.df1$cov.mat # just need to extract the covartiate matrix from either simulation
cov.mat.L <- cov.mat[,1:(N.Lcovs.linear+N.Lcovs.sq)]
 
t.list <- seq(0,admin.cens, 0.1)
surv1_wei <- sapply(t.list, function(t) {
  mean(exp(-lambda.tte*(t^gamma.tte)*exp(
    coeff.A +
      cov.mat.L %*% coeff.L +
      cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
      coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
  )))
})
surv0_wei <- sapply(t.list, function(t) {
  mean(exp(-lambda.tte*(t^gamma.tte)*exp(
    cov.mat.L %*% coeff.L +
      cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
      coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
  )))
})
 
weib.surv.df <- data.frame(
  time = c(t.list, t.list),
  A = c(rep(0,length(t.list)), rep(1,length(t.list))),
  surv = c(surv0_wei, surv1_wei),
  method = "True Weibull Survival"
)
 
## ---- Generate analysis dataset:
 
sim.df <- sim_surv_data(seed = this.seed,
                         N = simN,
                         Lcovs.linear = N.Lcovs.linear,
                         Lcovs.sq = N.Lcovs.sq,
                         mu = mu, # Lcovs.linear + Lcovs.sq + 2 (W and O)
                         sigma = sigma,
                         alpha.L = alpha.L, # params for exposure model c(L_1, ..., L_{Lcovs.linear + Lcovs.sq})
                         alpha.W = alpha.W,
                         # coefficient / log-HRs for outcome model
                         coeff.A = coeff.A,
                         coeff.L = coeff.L,
                         coeff.Lsq = coeff.Lsq,
                         coeff.O = coeff.O, # outcome-only confounder
                         coeff.W = coeff.W, # confounder of BOTH exposure and outcome (matches sim_data.R DGM)
                         gamma.tte = gamma.tte,
                         lambda.tte = lambda.tte,
                         lambda.cens = lambda.cens,
                         admin.cens = admin.cens,
                         gen.truth = NA
)
surv.df <- sim.df$data %>%
  mutate(A_fct = as.factor(A))
 
# Separate indices for propensity score and outcome model.
# Using one shared index (as previously) changes BOTH models together and does
# not represent the protocol's "misspecify one model only" scenarios.
ps.index  <- 1   # index into exposure.covs: 1=correct, 2=no W, 3=wrong form, 4=heavy
out.index <- 1   # index into outcome.covs:  1=correct, 2=no O,  3=wrong form, 4=heavy

iptw.formula <- paste("A", exposure.covs[ps.index], sep = " ~ ")
 
## ------ Singly Robust Methods:
 
## IPTW only analysis ##
iptw.model <- glm(as.formula(iptw.formula),
                  family="binomial",
                  data = surv.df)
pred.wt <- predict(iptw.model,
                   type = "response",
                   newdata = surv.df)
survwt.df <- surv.df
survwt.df$wt = ifelse(surv.df$A==1, pred.wt, 1-pred.wt)
survwt.df$wt.cum = ave(survwt.df$wt,survwt.df$id,FUN=cumprod)
 
iptw.model.num <- glm(A ~ 1, family="binomial", data = surv.df)
pred.wt.num <- predict(iptw.model.num,type = "response",newdata = surv.df)
 
survwt.df$wt.num = ifelse(survwt.df$A==1, pred.wt.num, 1-pred.wt.num)
survwt.df$wt.cum.num = ave(survwt.df$wt.num,survwt.df$id,FUN=cumprod)
 
#Stabilized weights
survwt.df$ipw.s <- survwt.df$wt.cum.num/survwt.df$wt.cum
 
# Fit the weighted cox model
iptw.cox <- coxph(Surv(eventtime, event) ~ A,
                  data = survwt.df,
                  weights = ipw.s,
                  ties = "breslow",
                  id = id)
 
# Build a treatment-specific marginal survival curve on the protocol time
# grid via survfit() rather than subject-level predict() at observed
# event times. The previous predict(..., type='survival', newdata=surv.df)
# returns one conditional survival per subject at that subject's observed
# event time, which is not a marginal S(t|A) curve and is not comparable
# with AIPTW or the protocol targets at t = 1, 5, 10.
cox.t.grid <- seq(0, admin.cens, 0.1)
# Two single-curve survfit calls instead of one multi-newdata call.
# summary.survfit() exposes $surv differently depending on the number of
# curves and the time grid; with two newdata rows it can be either a
# matrix or a vector pasted end-to-end. Splitting per arm avoids that
# ambiguity.
iptw.cox.fit0 <- survfit(iptw.cox, newdata = data.frame(A = 0))
iptw.cox.fit1 <- survfit(iptw.cox, newdata = data.frame(A = 1))
s0 <- summary(iptw.cox.fit0, times = cox.t.grid, extend = TRUE)$surv
s1 <- summary(iptw.cox.fit1, times = cox.t.grid, extend = TRUE)$surv

cox.surv.df <- data.frame(
  time   = rep(cox.t.grid, 2),
  A      = rep(c(0, 1), each = length(cox.t.grid)),
  surv   = c(s0, s1),
  method = paste0("IPTW Cox (", exposure.covs[ps.index], ")")
)
 
# Fit NPH FPM:
fpm.nph.fit <- flexsurvspline(Surv(eventtime, event) ~ A + gamma1(A),
                              data = survwt.df,
                              weights = ipw.s,
                              k = 2, scale = "hazard")
tidy(fpm.nph.fit)
 
fpm.surv.df <- summary(fpm.nph.fit, type="survival",
                    t = seq(0, admin.cens, 0.1),
                    newdata=data.frame(A = c(0,1)),
                    tidy = TRUE) %>%
  rename(surv = est) %>%
  mutate(method = paste0("IPTW FPM NPH (", exposure.covs[ps.index],")")) %>%
  dplyr::select(time, surv, A, method)
 
# Fit discrete-time pooled logistic model.
# Discrete-time interval setup mirrors utils/aiptw_estimator.R:
# interval_ends are the right endpoints (0.25, ..., admin.cens) and
# split_cuts are the interior cut times passed to survSplit. Excluding
# 0 and admin.cens prevents the first interval being mis-labelled at
# t = 0 and the last interval being empty.
rescale_time  <- 1/4
interval_ends <- seq(rescale_time, admin.cens, by = rescale_time)
split_cuts    <- interval_ends[-length(interval_ends)]
cutpoints     <- interval_ends   # legacy alias used by interval_mapping below
survwt.long.df <- survSplit(Surv(eventtime, event) ~ .,
                                       data = survwt.df,
                                       cut = split_cuts,
                                       episode = "time_period")
 
# Fit the pooled logistic regression
plr_model <- glm(event ~ A + as.factor(time_period),
                 data = survwt.long.df,
                 family = quasibinomial(link = "cloglog"), # need this due to weighting
                 weights = ipw.s)
 
# create prediction data
time_points <- levels(as.factor(survwt.long.df$time_period))
 
plr.surv.df <- data.frame(
  time_period = c(as.numeric(time_points),as.numeric(time_points)),
  A = c(rep(0,length(time_points)),
        rep(1,length(time_points)))
)
 
plr.surv.df$haz <- predict(plr_model, newdata = plr.surv.df, type = "response")
plr.surv.df$surv_interval <- 1 - plr.surv.df$haz # Calculate interval-specific survival
 
plr.surv.df$surv <- NA
plr.surv.df$surv[plr.surv.df$A==0] <- cumprod(plr.surv.df$surv_interval[plr.surv.df$A==0]) # Use cumprod to get the cumulative survival curve
plr.surv.df$surv[plr.surv.df$A==1] <- cumprod(plr.surv.df$surv_interval[plr.surv.df$A==1]) # Use cumprod to get the cumulative survival curve
 
interval_mapping <- data.frame(
  time_period = seq_along(interval_ends),
  start_time  = c(0, interval_ends[-length(interval_ends)]),
  time        = interval_ends
)
# Map period -> calendar time first, then prepend an explicit S(0) = 1
# row at time = 0. Adding the boundary row via time_period = 1 (as the
# previous code did) would map t = 0 onto interval_ends[1] = 0.25,
# duplicating the first interval endpoint and shifting the curve.
plr.surv.df <- plr.surv.df %>%
  left_join(interval_mapping %>% dplyr::select(time_period, time),
            by = "time_period") %>%
  bind_rows(data.frame(time_period = NA_integer_, surv = 1,
                        A = c(0L, 1L), time = 0)) %>%
  mutate(method = paste0("IPTW Discrete (", exposure.covs[ps.index], ")")) %>%
  dplyr::select(time, A, surv, method)
 
## Fit a super learner ??
 
 
ggplot() +
  geom_step(data=rbind(cox.surv.df, weib.surv.df, fpm.surv.df, plr.surv.df),
            aes(x = time, y = surv, color=method)) +
  facet_grid(~A) +
  ylim(0,1) + xlim(0,admin.cens) +
  theme_classic() +
  xlab("Time") + ylab("Marginal Survival Probability") +
  ggtitle("Compare Methods")
 
 
## Regression Standardisation / simple G-computation analysis ##
surv.formula <- "Surv(eventtime, event)"
surv.model.formula <- paste(surv.formula, outcome.covs[out.index], sep = " ~ ")
glm.model.formula <- paste0("event ~ ", outcome.covs[out.index], " + as.factor(time_period)")
fpm.model.formula <- paste(surv.model.formula, "gamma1(A)", sep = " + ") # add NPH for A
 
# For the cox model, make use of adjustedCurves(). Alternatively can use riskRegression::ate()
cond.cox <- coxph(as.formula(gsub("\\bA\\b", "A_fct", surv.model.formula)), # use the factor var for A
                  data = surv.df,
                  ties = "breslow",
                  id = id,
                  x=TRUE)
 
cox_adjsurv <- adjustedsurv(data=surv.df,
                        variable="A_fct",
                        ev_time="eventtime",
                        event="event",
                        method="direct",
                        outcome_model=cond.cox,
                        conf_int=FALSE)
 
cox.adjsurv.df <- cox_adjsurv$adj %>% rename(A = group) %>%
  mutate(method = paste0("Cox Reg Stand (", outcome.covs[out.index],")"))
 
# Fit NPH FPM:
cond.fpm.nph.fit <- flexsurvspline(as.formula(fpm.model.formula),
                              data = surv.df,
                              k = 2, scale = "hazard")
tidy(cond.fpm.nph.fit)
 
fpm.adjsurv.df <- standsurv(cond.fpm.nph.fit, type="survival",
                       t = seq(0, admin.cens, 0.1),
                       at=list(list(A=0),
                               list(A=1))) %>%
  tidyr::pivot_longer(., cols=c("at1","at2"), names_to="A", values_to = "surv") %>%
  mutate(A = ifelse(A == "at1", 0, 1)) %>%
  mutate(method = paste0("FPM NPH Reg Stand (", outcome.covs[out.index],")")) %>%
  dplyr::select(time, surv, A, method)
 
 
# Fit the pooled logistic regression for reg stand - need to do manual implementation
 
surv.long.df <- survSplit(Surv(eventtime, event) ~ .,
                            data = surv.df,
                            cut = split_cuts,           # interior cuts only
                            episode = "time_period")
 
cond.plr.mod <- glm(as.formula(glm.model.formula),
                    data = surv.long.df,
                    family = binomial(link = "cloglog")) # don't use quasi, as we don't have weighted df
 
baseline.covs <- surv.df %>%
  dplyr::select(id, L1, L1sq, L2, L2sq, L3, L4, L5, L6, O, W)
 
pred.df0 <- expand.grid(
  id          = surv.df$id,
  time_period = sort(unique(surv.long.df$time_period))
) %>%
  left_join(baseline.covs, by = "id") %>%
  mutate(
    A     = 0,
    event = 0L     # event column required by model
  )
 
pred.df1 <- expand.grid(
  id          = surv.df$id,
  time_period = sort(unique(surv.long.df$time_period))
) %>%
  left_join(baseline.covs, by = "id") %>%
  mutate(
    A     = 1,
    event = 0L     # event column required by model
  )
 
predgrid.df <- rbind(pred.df0, pred.df1)
 
# get discrete haz
predgrid.df$haz <- predict(cond.plr.mod, newdata = predgrid.df, type = "response")
 
# gte individual/conditional surv
plr.adjsurv.df_ <- predgrid.df %>%
  arrange(id, time_period) %>%
  group_by(id, A) %>%
  mutate(csurv = cumprod(1 - haz)) %>%
  ungroup()
 
# regression standaridsation step for adjusted surv
plr.adjsurv.df<- plr.adjsurv.df_ %>%
  group_by(time_period, A) %>%
  summarise(surv = mean(csurv), .groups = "drop") %>%
  left_join(
    interval_mapping %>% dplyr::select(time_period, time),
    by = "time_period"
  ) %>%
  mutate(method = paste0("Discrete Reg Stand (", outcome.covs[out.index],")")) %>%
  dplyr::select(-time_period)
 
 
ggplot() +
  geom_step(data=rbind(cox.surv.df, weib.surv.df, fpm.surv.df, plr.surv.df,
                       cox.adjsurv.df, fpm.adjsurv.df, plr.adjsurv.df),
            aes(x = time, y = surv, color=method)) +
  facet_wrap(~A) +
  ylim(0,1) + xlim(0,admin.cens) +
  theme_classic() +
  xlab("Time") + ylab("Marginal Survival Probability") +
  ggtitle("Compare Methods")
 
## ------ Weighted Standardisation (non-DR):
##
## NOTE: Weighted Cox / FPM / cloglog-PLR standardisation is NOT doubly
## robust (Gabriel et al. 2024; Apr 28 methods document Section 3.2).
## The canonical DR estimator uses a logit (canonical) link in the
## weighted PLR; that variant is implemented as AIPTW in
## utils/aiptw_estimator.R. The block below is kept as a comparison
## arm and is labelled accordingly.
 
# For the cox model, make use of adjustedCurves(). Alternatively can use riskRegression::ate()
cond.wtcox <- coxph(as.formula(gsub("\\bA\\b", "A_fct", surv.model.formula)), # use the factor var for A
                  data = survwt.df,
                  weights = ipw.s,
                  ties = "breslow",
                  id = id,
                  x=TRUE)
 
cox_drsurv <- adjustedsurv(data=surv.df,
                            variable="A_fct",
                            ev_time="eventtime",
                            event="event",
                            method="direct",
                            outcome_model=cond.wtcox,
                            conf_int=FALSE)
 
cox.drsurv.df <- cox_drsurv$adj %>% rename(A = group) %>%
  mutate(method = paste0("Cox Weighted Standardisation (non-DR) (", outcome.covs[out.index],")"))
 
# Fit FPM NPH
cond.wtfpm.nph.fit <- flexsurvspline(as.formula(fpm.model.formula),
                                   data = survwt.df,
                                   weights = ipw.s,
                                   k = 2, scale = "hazard")
fpm.drsurv.df <- standsurv(cond.wtfpm.nph.fit, type="survival",
                            t = seq(0, admin.cens, 0.1),
                            at=list(list(A=0),
                                    list(A=1))) %>%
  tidyr::pivot_longer(., cols=c("at1","at2"), names_to="A", values_to = "surv") %>%
  mutate(A = ifelse(A == "at1", 0, 1)) %>%
  mutate(method = paste0("FPM NPH Weighted Standardisation (non-DR) (", outcome.covs[out.index],")")) %>%
  dplyr::select(time, surv, A, method)
 
# Fit the pooled logistic regression for reg stand - need to do manual implementation
cond.wtplr.mod <- glm(as.formula(glm.model.formula),
                    data = survwt.long.df,
                    family = quasibinomial(link = "cloglog"),
                    weights = ipw.s)
 
dr.predgrid.df <- rbind(pred.df0, pred.df1)
 
# get discrete haz
dr.predgrid.df$haz <- predict(cond.wtplr.mod, newdata = dr.predgrid.df, type = "response")
 
# gte individual/conditional surv
plr.drsurv.df_ <- dr.predgrid.df %>%
  arrange(id, time_period) %>%
  group_by(id, A) %>%
  mutate(csurv = cumprod(1 - haz)) %>%
  ungroup()
 
plr.drsurv.df <- plr.drsurv.df_ %>%
  group_by(time_period, A) %>%
  summarise(surv = mean(csurv), .groups = "drop") %>%
  left_join(
    interval_mapping %>% dplyr::select(time_period, time),
    by = "time_period"
  ) %>%
  mutate(method = paste0("Discrete Weighted Standardisation (non-DR) (", outcome.covs[out.index],")")) %>%
  dplyr::select(-time_period)
 
ggplot() +
  geom_step(data=rbind(weib.surv.df,
                       cox.drsurv.df, plr.drsurv.df, fpm.drsurv.df),
            aes(x = time, y = surv, color=method)) +
  facet_wrap(~A) +
  ylim(0,1) + xlim(0,admin.cens) +
  theme_classic() +
  xlab("Time") + ylab("Marginal Survival Probability") +
  ggtitle("Compare Methods")
 
ggplot() +
  geom_step(data=rbind(cox.surv.df, weib.surv.df, fpm.surv.df, plr.surv.df,
                       cox.adjsurv.df, fpm.adjsurv.df, plr.adjsurv.df,
                       cox.drsurv.df, fpm.drsurv.df, plr.drsurv.df),
            aes(x = time, y = surv, color=method)) +
  facet_wrap(~A) +
  ylim(0,1) + xlim(0,admin.cens) +
  theme_classic() +
  xlab("Time") + ylab("Marginal Survival Probability") +
  ggtitle("Compare Methods")
