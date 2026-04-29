# BELOW PARAMS GEN ROUGHLY Delayed Effects
 
### parameters for simulation:
this.seed <- 2026
# set covariate params
N.Lcovs.linear = 4
N.Lcovs.sq = 2
mu = rep(0,8) # Lcovs.linear + Lcovs.sq + 2 (W and O)
sigma = diag(8)
 
# exposure model params - obsolete for truth params
alpha.L = c(0.5, 0.4, -0.6, 0.4, 0.5, 0.3) # params for exposure model c(L_1, ..., L_{Lcovs.linear + Lcovs.sq})
alpha.W = 1 # coefficient / log-HRs for outcome model
 
# outcome model params - double these (apart from A) so that the ill die early leaving us with the healthiest
coeff.A = log(0.5)
coeff.L = log(c(0.75, 2, 0.5, 0.8, 5, 1.5))/3
coeff.Lsq = log(c(0.95, 1.5))/3
coeff.O = log(0.8)/3 # outcome-only confounder
coeff.W = log(1.5)/3 # confounder of BOTH exposure and outcome (matches sim_data.R DGM)
 
# params for generating event and censoring times
gamma.tte = 2
lambda.tte = 0.015
lambda.cens = 0.1
admin.cens = 10
