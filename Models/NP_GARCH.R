

## Nonparameric GARCH (Bühlmann & McNeil, 2002)


NP_GARCH = function(X, M, K, method, weigh=FALSE, fin_smooth, innov, VaR_q) {
  ## M:      no. of iterative updates of volatility estimates
  ## K:      no. of last few iterations over which volatility estimates are averaged
  ## method: (non-parametric) regression method for estimating f
  ## weigh: if TRUE, weights observations w/ inverse of previous variance estimates
  
  cat(paste0("Fitting ", method, "-GARCH... "))
  
  n = length(X)
  sigma = matrix(rep(NA, (M+1)*n), nrow=n, ncol=M+1)
  
  ## initiate volatility estimates with MLE of GARCH(1,1)
  spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                    distribution.model = innov)
  
  fit0 = ugarchfit(data=X, spec=spec,
                   fit.control = list(stationarity=0))
  
  sigma[,1] = sigma(fit0)
  sigma[1,] = sigma(fit0)[1]
  df_est = coef(fit0)["shape"]
  
  ## iteratively regress X^2 on backshifted X and volatility estimates
  ## and update volatility estimates
  
  target = X[2:n]^2
  pred1 = X[1:(n-1)]
  
  for (m in 2:(M+1)) {
    
    pred2 = sigma[1:(n-1), m-1]^2
    
    # regression weights
    weights = NULL
    if (weigh == TRUE) {
      weights = sigma[2:n, m-1]^(-2)
      weights = weights/mean(weights)
    }
    
    # regress target on predictors nonparametrically
    sigma[2:n, m] = np_regr(method, target, pred1, pred2, weights)$vol_pred
  }
  
  # avg. over last K volatility estimates
  if (K > 1) {
    vol_avg = rowMeans(sigma[, (M-K+2):(M+1)])
  } else {
    vol_avg = sigma[,(M+1)]
  }
  
  # final smoothing step
  if (fin_smooth == TRUE) {
    pred2 = vol_avg[1:(n-1)]^2
    vol_2_to_n = np_regr(method, target, pred1, pred2, weights)$vol_pred
    vol = c(vol_avg[1], vol_2_to_n)
    
  } else {
    vol = vol_avg
  }
  
  # estimate VaR as 1% quantile of innovation distribution
  # with standard deviation equal to estimated volatility
  # (and the estimated degrees of freedom if t)
  if (innov == "norm") {
    var = sapply(VaR_q, qnorm, mean=0, sd=vol)
    
  } else if (innov == "std") {
    var = sapply(VaR_q, q_std_t, mu=0, sd=vol, df=df_est)
    
  } else {
    cat("Innovation distribution not implemented in NP_GARCH")
  }
  
  vol = data.frame("mean"=vol, "med"=vol)
  
  var = (-1)*var
  var = as.data.frame(var)
  colnames(var) = paste0(VaR_q*100, "%")
  rownames(var) = paste0("t=", 1:n)
  
  cat("done!\n")
  
  result = list("Vol_pred"=vol, "VaR_pred"=var)
  
  # clear memory
  rm(sigma, vol, var, target, pred1, pred2)
  if (K>1) rm(vol_avg)
  if (fin_smooth) rm(vol_2_to_n)
  gc()
  
  return(result)
}