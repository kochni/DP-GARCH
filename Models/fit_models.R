
#### helper functions
source("helpers.R")

## models
source("GARCH.R")
source("MS_GARCH.R")
source("NP_GARCH.R")
source("Bay_GARCH.R")
source("DP_GARCH.R")


GARCH_fit = function(X, args) {
  
  type = args$type
  
  if (type == "GARCH") {
    p = args$p
    q = args$q
    model = args$model
    innov = args$innov
    VaR_q = args$VaR_q
    fit = GARCH(X=X, model=model, p=p, q=q, innov=innov, VaR_q=VaR_q)
  
    } else if (type == "NP-GARCH") {
    M = args$M
    K = args$K
    method = args$method
    weigh = args$weigh
    fin_smooth = args$fin_smooth
    innov = args$innov
    VaR_q = args$VaR_q
    fit = NP_GARCH(X=X, M=M, K=K, method=method, weigh=weigh, fin_smooth=fin_smooth, innov=innov, VaR_q=VaR_q)
  
  } else if (type == "Bay-GARCH") {
    n_post = args$n_post
    burnin = args$burnin
    thin = args$thin
    innov = args$innov
    VaR_q = args$VaR_q
    fit = Bay_GARCH(X=X, innov=innov, n_post=n_post, burnin=burnin, thin=thin, VaR_q=VaR_q)
    
  } else if (type == "MS-GARCH") {
    models = args$models
    innovs = args$innovs
    est = args$est
    mix = args$mix
    n_post = args$n_post
    burnin = args$burnin
    thin = args$thin
    VaR_q = args$VaR_q
    fit = MS_GARCH(X=X, models=models, innovs=innovs, mix=mix, VaR_q=VaR_q, est=est, n_post=n_post, burnin=burnin, thin=thin)
  
  } else if (type == "NPB-GARCH") {
    law = args$law
    burnin = args$burnin
    n_post = args$n_post
    thin = args$thin
    VaR_q = args$VaR_q
    fit = DP_GARCH(X=X, law=law, n_post=n_post, burnin=burnin, thin=thin, VaR_q=VaR_q)
  
    } else {
    cat("GARCH type not implemented")
  }
  
  return(fit)
}
