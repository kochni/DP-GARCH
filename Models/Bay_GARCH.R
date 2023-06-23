

## Bayesian GARCH

Bay_GARCH = function(X, innov, VaR_q, n_post, burnin, thin, plot=FALSE) {
  ## X:      n-vector of observations
  ## n_post: no. of posterior samples
  
  cat("Fitting Bay-GARCH... ")
  
  length = n_post*thin + burnin  # elongate chain to preserve effective sample size
  n = length(X)
  
  if (innov == "norm") {
    lambda = 100
    delta = 500
  } else {
    lambda = NULL
    delta = NULL
  }
  
  control = list(l.chain=length, refresh=0)
  mcmc = bayesGARCH(X, control=control, lambda=lambda, delta=delta)
  post_samples = formSmpl(mcmc, l.bi=burnin, batch.size=thin)
  
  # average volatilities over posterior samples
  vol = matrix(rep(NA, n*n_post), nrow=n_post, ncol=n)
  X_sample = matrix(rep(NA, n*n_post), nrow=n_post, ncol=n)
  
  for (i in 1:n_post) {
    coef = post_samples[i,]
    omega = as.numeric(coef["alpha0"])
    alpha = as.numeric(coef["alpha1"])
    beta = as.numeric(coef["beta"])
    
    # get volatilities
    vol[i,] = get_vol(X, omega, alpha, beta)
    
    # sample 1 datapoint per time index for each posterior sample
    if (innov == "norm") {
      X_sample[i,] = rnorm(n=n, mean=0, sd=vol[i,])
      
    } else if (innov == "std") {
      df = as.numeric(coef["nu"])
      X_sample[i,] = r_std_t(n=n, mu=0, sd=vol[i,], df=df)
      
    } else {
      cat("Bay-GARCH: unspecified innovation distribution")
    }
  }
  
  vol_mean = colMeans(vol)
  vol_med = apply(vol, 2, median)
  vol = data.frame("mean"=vol_mean, "med"=vol_med)
  
  var = apply(X_sample, 2, quantile, probs=VaR_q)
  var = t(var)
  var = (-1)*var
  var = as.data.frame(var)
  colnames(var) = paste0(VaR_q*100, "%")
  rownames(var) = paste0("t=", 1:n)
  
  cat("done!\n")
  
  result = list("Vol_pred"=vol, "VaR_pred"=var)
  
  # clear memory
  rm(mcmc, post_samples, vol, vol_mean, vol_med, var, X_sample, Vol_pred)
  gc()
  
  return(result)
}
