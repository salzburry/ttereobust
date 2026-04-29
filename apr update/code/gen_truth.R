# get working dir
wd <- getwd()
 
# load libraries
library(survival)
library(dplyr)
library(purrr)
library(stringr)
library(MASS) # to simulate from a mvn
library(ggplot2)
library(flexsurv)
library(bshazard)
library(gridExtra)
library(grid)
 
source("utils/sim_data.R")
source("utils/compute_truth.R")

# Simulate data - set params
simN = 500000

#source("get_params_waning.R")
#source("get_params_delayed.R")
source("get_params_ph.R")

# DGM parameter bag passed to compute_truth() and reused for the analytic
# hazard / HR computation below.
dgm_params <- list(
  N.Lcovs.linear = N.Lcovs.linear, N.Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma,
  alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens
)

# Two A-forced samples used only for the KM Monte-Carlo truth check below.
# The analytic Weibull truth comes from compute_truth().
sim.df0 <- sim_surv_data(seed = this.seed, N = simN,
  Lcovs.linear = N.Lcovs.linear, Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma, alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens,
  gen.truth = 0
)
sim.df1 <- sim_surv_data(seed = this.seed, N = simN,
  Lcovs.linear = N.Lcovs.linear, Lcovs.sq = N.Lcovs.sq,
  mu = mu, sigma = sigma, alpha.L = alpha.L, alpha.W = alpha.W,
  coeff.A = coeff.A, coeff.L = coeff.L, coeff.Lsq = coeff.Lsq,
  coeff.O = coeff.O, coeff.W = coeff.W,
  gamma.tte = gamma.tte, lambda.tte = lambda.tte,
  lambda.cens = lambda.cens, admin.cens = admin.cens,
  gen.truth = 1
)

surv.df   <- rbind(sim.df0$data, sim.df1$data)
cov.mat   <- sim.df1$cov.mat
cov.mat.L <- cov.mat[, 1:(N.Lcovs.linear + N.Lcovs.sq)]

# Analytic Weibull truth. Replaces the two sapply blocks that recomputed
# what compute_truth() does. The hazard / HR diagnostics below still need
# cov.mat for the per-row hazard computation, so it stays available.
t.list   <- seq(0, admin.cens, 0.1)
truth_df <- compute_truth(t_eval = t.list, dgm = dgm_params)
surv0_wei <- truth_df$S0
surv1_wei <- truth_df$S1

weib.surv.df <- data.frame(
  time = c(t.list, t.list),
  A    = c(rep(0, length(t.list)), rep(1, length(t.list))),
  surv = c(surv0_wei, surv1_wei)
)
 
#Kaplan-Meier estimates of survival probabilities
surv.A <- survfit(Surv(eventtime,event) ~ A,data = surv.df)
surv.A.summ <- summary(surv.A, times=t.list)
 
#Obtain survival probabilities under the two treatment strategies: "always treated" (surv1), "never treated" (surv0)
surv0_km <- surv.A.summ$surv[surv.A.summ$strata=="A=0"]
surv1_km <- surv.A.summ$surv[surv.A.summ$strata=="A=1"]
 
#corresponding risk differences
survdiff_km <- surv1_km-surv0_km
 
# save data frame with truth
km.surv.df <- data.frame(
  time = c(surv.A.summ$time[surv.A.summ$strata=="A=0"], surv.A.summ$time[surv.A.summ$strata=="A=1"]),
  A = c(rep(0,length(surv.A.summ$time[surv.A.summ$strata=="A=0"])), rep(1,length(surv.A.summ$time[surv.A.summ$strata=="A=1"]))),
  surv = c(surv0_km, surv1_km)
)
 
surv.plot <- ggplot() +
  geom_step(data=km.surv.df, aes(x = time, y = surv, linetype=as.factor(A)), color="red") +
  geom_step(data=weib.surv.df, aes(x = time, y = surv, linetype=as.factor(A)), color="blue") +
  ylim(0,1) + xlim(0,admin.cens) +
  theme_classic() +
  xlab("Time") + ylab("Survival") +
  ggtitle("True Survival Curves")
#
# ggplot() +
#   geom_step(data=km.surv.df, aes(x = log(time), y = log(-log(surv)), linetype=as.factor(A)), color="red") +
#   geom_step(data=weib.surv.df, aes(x = log(time), y = log(-log(surv)), linetype=as.factor(A)), color="blue") +
#   theme_classic() +
#   xlab("log-Time") + ylab("Log-Cumulative Hazards") +
#   ggtitle("Check PH")
 
## get "true" HRs
haz1_wei <- sapply(t.list, function(t) {
 
  numer <- mean(gamma.tte*lambda.tte*(t^(gamma.tte-1))*exp(
    coeff.A +
      cov.mat.L %*% coeff.L +
      cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
      coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
  )*
    exp(-lambda.tte*(t^gamma.tte)*exp(
      coeff.A +
        cov.mat.L %*% coeff.L +
        cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
        coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
    ))
  )
 
  denom <- mean(exp(-lambda.tte*(t^gamma.tte)*exp(
    coeff.A +
      cov.mat.L %*% coeff.L +
      cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
      coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
  )))
 
  numer/denom
 
})
 
haz0_wei <- sapply(t.list, function(t) {
 
  numer <- mean(gamma.tte*lambda.tte*(t^(gamma.tte-1))*exp(
      cov.mat.L %*% coeff.L +
      cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
      coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
  )*
    exp(-lambda.tte*(t^gamma.tte)*exp(
        cov.mat.L %*% coeff.L +
        cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
        coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
    ))
  )
 
  denom <- mean(exp(-lambda.tte*(t^gamma.tte)*exp(
      cov.mat.L %*% coeff.L +
      cov.mat.L[,1:2]^2 %*% coeff.Lsq[1:2] +
      coeff.O*cov.mat[,"O"] + coeff.W*cov.mat[,"W"]
  )))
 
  numer/denom
 
})
 
weib.haz.df <- data.frame(
  time = c(t.list,t.list),
  haz = c(haz0_wei,haz1_wei),
  A = c(rep(0, length(haz0_wei)),rep(1, length(haz0_wei)))
)
 
weib.hr.df <- data.frame(
  time = t.list,
  haz0 = haz0_wei,
  haz1 = haz1_wei
)  %>% mutate(hr = haz1 / haz0)
 
hr.plot <- ggplot() +
  geom_line(data=weib.hr.df, aes(x = time, y = hr)) +
  geom_hline(yintercept = exp(coeff.A), lty=2, color="red") +
  theme_classic() +
  xlab("Time") + ylab("Marginal Hazard Ratio") +
  ggtitle("True weibull marginal HRs") +
  ylim(exp(coeff.A)-0.1,1)
 
haz.plot <- ggplot() +
  geom_line(data=weib.haz.df , aes(x = time, y = haz, color=as.factor(A))) +
  theme_classic() +
  xlab("Time") + ylab("Marginal Hazard Rates") +
  ggtitle("True weibull marginal hazards")
 
coef_table <- data.frame(
  Variable = c("A", "L", "L (non-linear)", "O", "W",
               "gamma", "lambda"),
  Coefficient = c(
    toString(round(coeff.A, 4)),
    toString(round(coeff.L, 4)),
    toString(round(coeff.Lsq, 4)),
    toString(round(coeff.O, 4)),
    toString(round(coeff.W, 4)),
    toString(round(gamma.tte, 4)),
    toString(round(lambda.tte, 4))
  )
)
 
table <- tableGrob(coef_table, rows = NULL)
 
grid.arrange(surv.plot, haz.plot, hr.plot, table, nrow = 2, ncol = 2)
 
## Get smoothed HRs
 
# fit0 <- bshazard(
#   Surv(eventtime,event) ~ 1,
#   data = subset(surv.df, A == 0),
#   nk = 31, nbin = 100
# )
# fit1 <- bshazard(
#   Surv(eventtime,event) ~ 1,
#   data = subset(surv.df, A == 1),
#   nk = 31, nbin = 100
# )
#
# haz.ratio.smooth <- bind_rows(
#   data.frame(id = 1:length(fit0$time), time = fit0$time, Hazard = fit0$hazard, A = 0),
#   data.frame(id = 1:length(fit0$time), time = fit1$time, Hazard = fit1$hazard, A = 1)
# )
#  
#  
# haz.ratio.smooth.df <- tidyr::pivot_wider(haz.ratio.smooth %>% mutate(time = round(time, 3)),
#                                           id_cols = c("time"),
#                                           names_from = "A",
#                                           values_from = "Hazard",
#                                           names_prefix = "haz") %>%
#   mutate(HR = haz1/haz0)
#  
# ggplot() +
#   geom_line(data=haz.ratio.smooth.df, aes(x = time, y = HR)) +
#   geom_hline(yintercept = 0.5, lty=2, color="red") +
#   theme_classic() +
#   xlab("Time") + ylab("Hazard Ratio") +
#   ggtitle("Smooth Hazard Ratios using bshazard()") +
#   ylim(0.5,1)
 
 
# Let us also fit a NPH unadjusted marginal model to visualise the time-varying HR
# fpm.nph.fit <- flexsurvspline(Surv(eventtime,event) ~ A + gamma1(A),
#                               data = surv.df,
#                               k = 2, scale = "hazard")
# tidy(fpm.nph.fit)
#
# haz.df <- summary(fpm.nph.fit, type="hazard", t = seq(0, admin.cens, 0.1), newdata=data.frame(A = c(0,1)), tidy = TRUE)
#
# hr.df <- tidyr::pivot_wider(haz.df,
#                             id_cols = c("time"),
#                             names_from = "A",
#                             values_from = "est",
#                             names_prefix = "haz") %>%
#   mutate(HR = haz1/haz0)
#
# ggplot() +
#   geom_line(data=hr.df, aes(x = time, y = HR)) +
#   geom_hline(yintercept = 0.5, lty=2, color="red") +
#   theme_classic() +
#   xlab("Time") + ylab("Hazard Ratio") +
#   ggtitle("NPH Flexible Parametric Model (unadjusted)") +
#   ylim(0.5,1)
 
 
# ## check for PHs - this looks roughly proportional
# cox.ph <- coxph(Surv(eventtime,event) ~ A, data = surv.df)
# summary(cox.ph)
# ph.test <- cox.zph(cox.ph)
