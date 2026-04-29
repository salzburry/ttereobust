# BELOW PARAMS GEN ROUGHLY PH: Note we need to have low heterogeneity of risk within the trt arm i.e. patients are very similar to each other w.r.t risk
 
### parameters for simulation:
this.seed <- 2026
# set covariate params
N.Lcovs.linear = 4
N.Lcovs.sq = 2
mu = rep(0,8) # Lcovs.linear + Lcovs.sq + 2 (W and O)
sigma = diag(8)
 
# exposure model params - obsolete for truth params
alpha.L = c(0.5, 0.4, -0.6, 0.4, 0.5, 0.3) # params for exposure model c(L_1, ..., L_{Lcovs.linear + Lcovs.sq})
alpha.W = 1 # coefficient on W in the exposure logit (sim_data.R)
 
# outcome model params - have low heterogeneity of risk to get approx proportional on conditional and marginal
coeff.A = log(0.5)
coeff.L = log(c(0.75, 2, 0.5, 0.8, 5, 1.5))/10
coeff.Lsq = log(c(0.95, 1.5))/10
coeff.O = log(0.8)/10 # outcome-only confounder
coeff.W = log(1.5)/10 # confounder of BOTH exposure and outcome (matches sim_data.R DGM)
 
# params for generating event and censoring times
gamma.tte = 1
lambda.tte = 0.2
lambda.cens = 0.1
admin.cens = 10
