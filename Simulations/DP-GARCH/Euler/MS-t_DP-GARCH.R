
################
## PLAYGROUND ##
################

source("Simulation Study.R")

# define ground truths

# simulation parameters
len = 1000
n_sim = 5

# set GARCH parameters (meaning depends on model type)
omega = 1
alpha = 0.3
beta  = 0.4
gamma = 0.3
nu    = 5


## set Value-at-Risk level of interest
VaR_q = c(0.01, 0.05)

# set process

truth = list(
  model  = "MS-GARCH",
  models = c("sGARCH", "gjrGARCH", "eGARCH"),
  innov  = c("std", "std", "std"),
  pars   = list("omega_1"=omega, "alpha_1"=alpha, "beta_1"=beta, "shape"=nu,
                "omega_2"=omega, "alpha_2"=alpha, "gamma_2"=gamma, "beta_2"=beta, "shape"=nu,
                "omega_3"=omega, "gamma_3"=gamma, "alpha_3"=alpha, "beta_3"=beta, "shape"=nu,
                "P_11"=0.8, "P_12"=0.1,
                "P_21"=0.1, "P_22"=0.8,
                "P_31"=0.1, "P_32"=0.1),   # order matters!!
  VaR_q = VaR_q
)

truth_name = "MS(t)"


## define model

model    = list(type="NPB-GARCH", law="DP", n_post=2000L, burnin=10000L, thin=5L, VaR_q=VaR_q)
model_name = "DP-GARCH"


## Go 

sim_study(truth=truth, truth_name=truth_name, model=model, model_name=model_name, len=len, n_sim=n_sim,
          save_results=TRUE, freeCores=1)

# plot results
#source("plots.R")
#res_plots(process="Exp-N", models=c("GARCH", "Bay-GARCH", "MS-GARCH", "B-MS-GARCH"))