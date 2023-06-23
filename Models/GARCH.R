

## Standard GARCH model & extensions

GARCH = function(X, model, p, q, innov, VaR_q) {
  
  cat(paste0("Fitting ", model, "(", p, ",", q, ")... "))
  
  n = length(X)
  
  if (model == "tGARCH") {
    spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                      variance.model = list(model = "fGARCH", submodel="TGARCH", garchOrder = c(p,q)),
                      distribution.model = innov)
  
    } else {
    spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                      variance.model = list(model = model, garchOrder = c(p,q)),
                      distribution.model = innov)
  }
  
  cat(paste("Fitting GARCH w", innov, "errors"))
  
  fit = ugarchfit(data = X, spec = spec,
                  fit.control = list(stationarity=0),
                  solver.control = list(tol = 1e-10))
  
  # extract volatility from fit
  vol = sigma(fit)
  vol = as.numeric(vol)
  vol = data.frame("mean"=vol, "med"=vol)
  
  # extract VaR from fit
  var = lapply(VaR_q, FUN=quantile, x=fit)
  var = do.call(cbind.data.frame, var)
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
