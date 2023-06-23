
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
  model = "eGARCH",
  innov = "std",
  order = c(1,1),
  pars  = list(omega=omega, alpha1=alpha, beta1=beta, gamma1=gamma, shape=nu),
  VaR_q = VaR_q
)

truth_name = "Exp(t)"


## define model

model      = list(type="NPB-GARCH", law="DP", n_post=2000L, burnin=10000L, thin=5L, VaR_q=VaR_q)
model_name = "DP-GARCH"

sim_study(truth=truth, truth_name=truth_name, model=model, model_name=model_name, len=len, n_sim=n_sim,
          save_results=TRUE, freeCores=1)

# plot results
#source("plots.R")
#res_plots(process="Exp-N", models=c("GARCH", "Bay-GARCH", "MS-GARCH", "B-MS-GARCH"))