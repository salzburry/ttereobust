# See simulation plan located here: https://myteams.gsk.com/:w:/r/sites/SarwarMStatsTeams/SiteAssets/My%20Documents/Projects/Doubly%20Robust%20Methods%20for%20Survival%20Outcomes/Methods%20Summary.docx?d=we0eddb92cdba46bba99544d2c95a49ac&csf=1&web=1&e=4ohVZy
 
## --- Setup
 
# get working dir
wd <- getwd()
 
# load libraries
library(survival)
library(dplyr)
library(purrr)
library(stringr)
library(MASS) # to simulate from a mvn
library(ggplot2)
 
source("utils/sim_data.R")
 
# Simulate data - set params
simN = 100
 
source("get_params.R")
 
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
                         admin.cens = admin.cens
)
 
surv.df <- sim.df$data
cov.mat <- sim.df$cov.mat
cov.mat.L <- cov.mat[,1:(N.Lcovs.linear+N.Lcovs.sq)]
 
## --- Now let's just check out some of what we have generated
 
# No. treated vs untreated
table(surv.df$A)
 
# Check functional form of covs
 
Lcovs.sq <- c(1:N.Lcovs.sq)
HR.Lsq <- lapply(Lcovs.sq, function(l) {
  exp(cov.mat.L[,l]*coeff.L[l] + coeff.Lsq[l]*(cov.mat.L[,l])^2)
})
names(HR.Lsq) <- paste0("HR.L", Lcovs.sq)
 
Lcovs.linear <- c((N.Lcovs.sq+1):(N.Lcovs.sq+N.Lcovs.linear))
HR.L.linear <- lapply(Lcovs.linear, function(l) {
  exp(cov.mat.L[,l]*coeff.L[l])
})
names(HR.L.linear) <- paste0("HR.L", Lcovs.linear)
 
check.cov.funcs.HR <- cbind(as.data.frame(HR.L.linear),
      as.data.frame(HR.Lsq)) %>%
  mutate(id = row_number()) %>%
  left_join(., surv.df %>% mutate(val = "cov") %>% dplyr::select(id, matches("^L\\d+$")),
            by="id")
 
check.cov.funcs.HR1 <- check.cov.funcs.HR %>%
  dplyr::select(id, starts_with("HR.L")) %>%
  tidyr::pivot_longer(
    cols = -id,
    names_to = "key",
    names_prefix = "HR.L",
    values_to = "HR"
  )
 
check.cov.funcs.HR2 <- check.cov.funcs.HR %>%
  dplyr::select(id, matches("^L\\d+$")) %>%
  tidyr::pivot_longer(
    cols = -id,
    names_to = "key",
    names_prefix = "L",     
    values_to = "L"
  )
 
check.cov.funcs.HR.df <- full_join(check.cov.funcs.HR1, check.cov.funcs.HR2, by = c("id", "key"))
 
ggplot(data = check.cov.funcs.HR.df) +  
  geom_point(aes(x = L, y = HR)) +
  facet_wrap(~ key) +
  theme_classic() +
  ggtitle("HRs for L Covariates")
