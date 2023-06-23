

## packages

packages_modeling = c(
  "rugarch",       # GARCH & basic extensions
  "MSGARCH",       # Markov-switching GARCH
  "bayesGARCH",    # Bayesian GARCH
  #"stochvol",     # stochastic volatility model
  "randomForest",  # random forest regressor
  "mgcv",          # GAM regressor
  "earth",         # MARS regressor
  "metRology",     # scaled Student-t distribution
  "extraDistr",    # Laplace distribution
  "foreach",       # parallelization
  "doParallel"     # parallelization
  )

lapply(packages_modeling, library, character.only = TRUE)


## Standardized Student-t distribution
r_std_t = function(n, mu=0, sd, df) {
  X = rt.scaled(n=n, mean=mu, sd=sd*sqrt((df-2)/df), df=df)
  return(X)
}

d_std_t = function(x, mu=0, sd, df) {
  dens = dt.scaled(x, mean=mu, sd=sd*sqrt((df-2)/df), df=df)
  return(dens)
}

q_std_t = function(p, mu=0, sd, df) {
  quant = qt.scaled(p=p, mean=mu, sd=sd*sqrt((df-2)/df), df=df)
  return(quant)
}
# Note: x.scaled gives corresponds to t scaled by sd;
#       sd of scaled t is thus sd*sqrt(df/(df-2));
#       hence sd=sd/sqrt(df/(df-2)) yields "actual"
#       std. dev. of  1


## Standardized Laplace distribution
r_std_lap = function(n, mu=0, sd) {
  sample = rlaplace(n=n, mu=mu, sigma=sqrt(0.5*sd))
  return(sample)
}

d_std_lap = function(x, mu=0, sd) {
  dens = dlaplace(x=x, mu=mu, sigma=sqrt(0.5*sd))
  return(dens)
}

q_std_lap = function(p, mu=0, sd) {
  quant = qlaplace(p=p, mu=mu, sigma=sqrt(0.5*sd))
  return(quant)
}

## loss functions
loss = function(error, abs, root) {
  
  if (abs==TRUE) {
    loss = mean(abs(error))
  } else {
    loss = mean(error^2)
  }
  
  if (root==TRUE) {
    loss = sqrt(loss)
  }
  
  return(loss)
}

# non-parametric regressor
np_regr = function(method, target, pred1, pred2, weights) {
  
  if (method=="RF") {
    fit = randomForest(target ~ pred1 + pred2,
                       ntree = 2000,
                       weights = weights,
                       na.action = "na.omit")
    
  } else if (method=="MARS") {
    fit = earth(target ~ pred1 + pred2,
                #glm = list(family=Gamma),
                weights = weights,
                degree = 2)
    
  } else if (method=="GAM") {
    fit = gam(target ~ s(pred1) + s(pred2) + s(pred1, pred2),
              #family = gaussian(link="log"),
              weights = weights,
              control = list(irls.reg=0))
    
  } else if (method=="loess") {
    fit = loess(target ~ pred1 + pred2,
                weights = weights)
    
  } else {
    cat("NP-GARCH: regression method not implemented")
    stop()
  }
  
  # truncate predictions at 0
  var_pred = predict(fit)
  var_pred[var_pred<0] = 0
  vol_pred = sqrt(var_pred)
  
  return(list("vol_pred"=vol_pred))
}


## get volatility from GARCH(1,1) coeff's & obs

get_vol = function(X, omega, alpha, beta) {
  
  n = length(X)
  vol = rep(NA, n)
  vol[1] = omega
  
  for (t in 2:n) {
    vol[t] = sqrt(omega + alpha*X[t-1]^2 + beta*vol[t-1]^2)
  }
  
  return("Vol"=vol)
}