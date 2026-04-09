# get working dir
wd <- getwd()
 
# load libraries
library(survival)
library(dplyr)
library(ggplot2)
library(flexsurv)
library(adjustedCurves)
library(MASS)
library(purrr)
library(gfoRmula)
 
source("utils/sim_data.R")
 
# Specify covs and model formulae
exposure.covs <- c("L1sq + L2sq + L3 + L4 + L5 + L6 + W", # correct
                   "L1sq + L2sq + L3 + L4 + L5 + L6", # no W
                   "L1 + L2 + L3 + L4 + L5 + L6 + W", # wrong functional form
                   "L3 + L4 + L5") # heavy mis-specification
exposure.covs.sceN <- length(exposure.covs)
 
outcome.covs <- c("A + L1sq + L2sq + L3 + L4 + L5 + L6 + O + W",
                  "A + L1sq + L2sq + L3 + L4 + L5 + L6",
                  "A + L1 + L2 + L3 + L4 + L5 + L6 + O + W",
                  "A + L3 + L4 + L5")
outcome.covs.sceN <- length(outcome.covs)
 
# Simulate data - set params for analysis data
simN = 10000
 
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
                         coeff.O = coeff.O, # this will only appear in outcome model
                         coeff.W = coeff.W, # this will only appear in exposure model
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
                         coeff.O = coeff.O, # this will only appear in outcome model
                         coeff.W = coeff.W, # this will only appear in exposure model
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
                         coeff.O = coeff.O, # this will only appear in outcome model
                         coeff.W = coeff.W, # this will only appear in exposure model
                         gamma.tte = gamma.tte,
                         lambda.tte = lambda.tte,
                         lambda.cens = lambda.cens,
                         admin.cens = admin.cens,
                         gen.truth = NA
)
surv.df <- sim.df$data %>%
  mutate(A_fct = as.factor(A))
 
# get model formula
cov.index <- 1
 
iptw.formula <- paste("A", exposure.covs[cov.index], sep = " ~ ")
 
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
 
surv_iptwcox <- predict(iptw.cox, type = "survival", newdata = surv.df)
 
cox.surv.df <- data.frame(
  time = surv.df$eventtime,
  A = surv.df$A,
  surv = surv_iptwcox,
  method = paste0("IPTW Cox (", exposure.covs[cov.index],")")
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
  mutate(method = paste0("IPTW FPM NPH (", exposure.covs[cov.index],")")) %>%
  dplyr::select(time, surv, A, method)
 
# Fit discrete-time pooled logistic model
 
#cutpoints <- unique(survwt.df$eventtime[survwt.df$event == 1])
rescale_time <- 1/4
cutpoints <- seq(0,admin.cens, rescale_time) # weekly cutpoints
survwt.long.df <- survSplit(Surv(eventtime, event) ~ .,
                                       data = survwt.df,
                                       cut = cutpoints,
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
  time_period = 1:(length(cutpoints)),
  start_time = c(0, cutpoints[-length(cutpoints)]),
  time = cutpoints
)
plr.surv.df <- left_join(plr.surv.df %>%
                           add_row(data.frame(time_period = c(1,1), surv = c(1,1), A = c(0,1))),
                         interval_mapping %>% dplyr::select(time_period, time), by = c("time_period")) %>%
  mutate(method = paste0("IPTW Discrete (", exposure.covs[cov.index],")")) %>%
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
surv.model.formula <- paste(surv.formula, outcome.covs[cov.index], sep = " ~ ")
glm.model.formula <- paste0("event ~ ", outcome.covs[cov.index], " + as.factor(time_period)")
fpm.model.formula <- paste(surv.model.formula, "gamma1(A)", sep = " + ") # add NPH for A
 
# For the cox model, make use of adjustedCurves(). Alternatively can use riskRegression::ate()
cond.cox <- coxph(as.formula(gsub(\bA\b, "A_fct", surv.model.formula)), # use the factor var for A
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
  mutate(method = paste0("Cox Reg Stand (", outcome.covs[cov.index],")"))
 
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
  mutate(method = paste0("FPM NPH Reg Stand (", outcome.covs[cov.index],")")) %>%
  dplyr::select(time, surv, A, method)
 
 
# Fit the pooled logistic regression for reg stand - need to do manual implementation
 
surv.long.df <- survSplit(Surv(eventtime, event) ~ .,
                            data = surv.df,
                            cut = cutpoints,
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
  mutate(method = paste0("Discrete Reg Stand (", outcome.covs[cov.index],")")) %>%
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
 
## ------ Doubly Robust Standardisation:
 
# For the cox model, make use of adjustedCurves(). Alternatively can use riskRegression::ate()
cond.wtcox <- coxph(as.formula(gsub(\bA\b, "A_fct", surv.model.formula)), # use the factor var for A
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
  mutate(method = paste0("Cox Doubly Robust RS (", outcome.covs[cov.index],")"))
 
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
  mutate(method = paste0("FPM NPH Doubly Robust RS (", outcome.covs[cov.index],")")) %>%
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
  mutate(method = paste0("Discrete Doubly Robust RS (", outcome.covs[cov.index],")")) %>%
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
                       cox.drsurv.df, plr.drsurv.df),
            aes(x = time, y = surv, color=method)) +
  facet_wrap(~A) +
  ylim(0,1) + xlim(0,admin.cens) +
  theme_classic() +
  xlab("Time") + ylab("Marginal Survival Probability") +
  ggtitle("Compare Methods")
