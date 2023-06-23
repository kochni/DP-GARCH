

## Markov-Switching GARCH

MS_GARCH = function(X, models, innovs, mix, est, n_post, burnin, thin, VaR_q) {
  
  cat(paste0("Fitting MS-GARCH (", est, ")... "))
  
  n = length(X)
  
  spec = CreateSpec(
    variance.spec     = list(model = models),
    distribution.spec = list(distribution = innovs),
    switch.spec       = list(do.mix = mix)
  )
  ctr = list(n.burn=burnin, nmcmc=n_post*thin, nthin=thin,
             do.se=FALSE)
  
  fit = FitMCMC(data=X, spec=spec, ctr=ctr)
  
  ## extract volatility & VaR
  
  if (est == "MLE") { 
    # initialize MLE optimizer w MCMC mean
    # to increase numeric stability
    par0 = colMeans(fit$par)
    
    ctr = list(par0=par0, do.se=FALSE)
    fit = FitML(data=X, spec=spec, ctr=ctr)
  }
  
  vol = Volatility(fit)
  vol = as.numeric(vol)
  vol = data.frame("mean"=vol, "med"=vol)
  
  var = Risk(fit, alpha=VaR_q, do.its=TRUE, do.es=FALSE)$VaR
  var = (-1)*var
  var = as.data.frame(var)
  colnames(var) = paste0(VaR_q*100, "%")
  rownames(var) = paste0("t=", 1:n)
  
  cat("done!\n")
  
  result = list("Vol_pred"=vol, "VaR_pred"=var)
  
  # clear memory
  rm(fit, vol, var)
  gc()
  
  return(result)
}
