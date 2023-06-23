
source("Simulation Study.R")


## set length & no. of simulations
len   = 1000
n_sim = 40


## set Value-at-Risk level of interest
VaR_q = c(0.01, 0.05)

## define true process
omega = 1
alpha = 0.2
alpha_Exp = 0
beta  = 0.4
gamma = 0.5
nu    = 5

truth = list(
  model  = "MS-GARCH",
  models = c("sGARCH", "gjrGARCH", "eGARCH"),
  innov  = c("std", "std", "std"),
  pars   = list("omega_1"=omega, "alpha_1"=alpha, "beta_1"=beta, "shape"=nu,
                "omega_2"=omega, "alpha_2"=alpha, "gamma_2"=gamma, "beta_2"=beta, "shape"=nu,
                "omega_3"=omega, "gamma_3"=gamma, "alpha_3"=alpha_Exp, "beta_3"=beta, "shape"=nu,
                "P_11"=0.8, "P_12"=0.1,
                "P_21"=0.1, "P_22"=0.8,
                "P_31"=0.1, "P_32"=0.1),   # order matters!!
  VaR_q = VaR_q
)

truth_name = "MS(t)"


## define models
models = list(
  #"DP-GARCH"   = list(type="NPB-GARCH", law="DP", n_post=200L, burnin=200L, thin=1L, VaR_q=VaR_q),
  "BayMS-GARCH"  = list(type="MS-GARCH", est="MCMC", mix=FALSE, models=rep("sGARCH", 3L), innovs=rep("norm", 3L), n_post=10000L, burnin=10000L, thin=1L, VaR_q=VaR_q),
  "MS-GARCH"     = list(type="MS-GARCH", est="MLE",  mix=FALSE, models=rep("sGARCH", 3L), innovs=rep("norm", 3L), n_post=10000L, burnin=10000L, thin=1L, VaR_q=VaR_q),
  "BayMix-GARCH" = list(type="MS-GARCH", est="MCMC", mix=TRUE,  models=rep("sGARCH", 3L), innovs=rep("norm", 3L), n_post=10000L, burnin=10000L, thin=1L, VaR_q=VaR_q),
  "Mix-GARCH"    = list(type="MS-GARCH", est="MLE",  mix=TRUE,  models=rep("sGARCH", 3L), innovs=rep("norm", 3L), n_post=10000L, burnin=10000L, thin=1L, VaR_q=VaR_q),
  "Bay-GARCH"    = list(type="Bay-GARCH", n_post=10000L, burnin=10000L, thin=1L, innov="norm", VaR_q=VaR_q),
  "RF-GARCH"     = list(type="NP-GARCH", M=20L, K=5L, method="RF",  innov="norm", fin_smooth=FALSE, weigh=FALSE, VaR_q=VaR_q),
  "T-GARCH"      = list(type="GARCH", model="tGARCH", p=1L, q=1L, innov="norm", VaR_q=VaR_q),
  "GARCH(3,3)"   = list(type="GARCH", model="sGARCH", p=3L, q=3L, innov="norm", VaR_q=VaR_q),
  "GARCH"        = list(type="GARCH", model="sGARCH", p=1L, q=1L, innov="norm", VaR_q=VaR_q)
)


sim_study(truth=truth, truth_name=truth_name, models=models, len=len, n_sim=n_sim,
          save_results=TRUE, freeCores=1)

# plot results
#source("plots.R")
#res_plots(process="Exp-N", models=c("GARCH", "Bay-GARCH", "MS-GARCH", "B-MS-GARCH"))