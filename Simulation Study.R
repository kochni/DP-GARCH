
source("fit_models.R")


sim_study = function(truth, truth_name, models, len, n_sim,
                     save_results, freeCores) {
  ## processes: ground truth processes to be simulated and estimated
  ## models:    models to be fitted to the simulated data
  ## len:       length of simulated processes
  ## n_sim:     no. of time each process is simulated
  
  start = Sys.time()
  
  n_models = length(models)
  names_models = names(models)
  
  cat("Simulation:\n")
  cat(paste("- Process:", truth_name, "\n"))
  cat(paste("- Length: ", len, "\n"))
  cat(paste("- Samples:", n_sim, "\n\n"))
  
  cat("Models:\n")
  for (j in 1:n_models) {
    cat(paste("-", names_models[j], "\n"))
  }
  cat("\n")
  
  # set up parallelization
  totalCores = parallelly::availableCores()
  usedCores = totalCores[1] - freeCores
  cluster <- makeCluster(usedCores)
  registerDoParallel(cluster)
  cat(paste0("Using ", usedCores, "/", totalCores[1], " CPU cores\n\n"))
  
  if (truth$model != "MS-GARCH" && truth$model != "SV") {
    truth_spec = ugarchspec(
      mean.model         = list(armaOrder = c(0,0), include.mean = FALSE),
      variance.model     = list(model=truth$model, garchOrder=truth$order),
      fixed.pars         = c(truth$pars, list("scale"=1, "mu"=0)),
      distribution.model = truth$innov
    )
  }
  
  
  ## simulate n_sim realizations from truth and estimate models
  
  result = foreach(sim = 1:n_sim, .combine="+") %dopar% {
    
    source("fit_models.R")
    
    ## Simulate time series
    
    if (truth$model == "MS-GARCH") {
      # due to bug in MSGARCH package
      # need to specify anew when inside 'dopar'
      truth_spec = CreateSpec(
        variance.spec     = list(model = truth$models),
        distribution.spec = list(distribution = truth$innov),
        switch.spec       = list(do.mix=FALSE)
      )
      
      simulation = simulate(object=truth_spec, nsim=1L, nahead=len, par=unlist(truth$pars), nburn=500L)
      
      X_sim = simulation$draw
      Vol_true = rep(NA, len)
      Vol_true[simulation$state==1] = simulation$CondVol[,,1][simulation$state==1]
      Vol_true[simulation$state==2] = simulation$CondVol[,,2][simulation$state==2]
      Vol_true[simulation$state==3] = simulation$CondVol[,,3][simulation$state==3]
      
    } else if (truth$model == "SV") {
      simulation = svsim(len = len, mu = truth$pars$mu,
                         phi = truth$pars$phi,
                         sigma = truth$pars$sigma,
                         nu = truth$pars$nu,
                         rho = truth$pars$rho)
      X_sim = simulation$y
      Vol_true = simulation$vol
      
    } else {
      simulation = ugarchpath(n.sim = len, n.start = 500L, spec = truth_spec)
      
      X_sim = fitted(simulation)
      Vol_true = sigma(simulation)
    }
    
    X_sim = as.numeric(X_sim)
    Vol_true = as.numeric(Vol_true)
    
    
    ## extract true VaR
    
    if (truth$innov[1] == "norm") {
      VaR_true = sapply(truth$VaR_q, qnorm, mean=0, sd=Vol_true)
      
    } else if (truth$innov[1] == "std") {
      if (truth$pars$shape > 2) {
        VaR_true = sapply(truth$VaR_q, q_std_t, mu=0, sd=Vol_true, df=truth$pars$shape)
      } else if (truth$pars$shape == 1) {
        VaR_true = sapply(truth$VaR_q, q_std_lap, mu=0, sd=Vol_true)
      }
      
    } else {
      cat("Innovation distribution not implemented")
      stop()
    }
    VaR_true = (-1)*VaR_true
    colnames(VaR_true) = paste0("VaR_", 100*truth$VaR_q, "%_true")
    rownames(VaR_true) = paste0("t=", 1:len)
    
    # ensure that valid values produced
    if (sum(is.na(c(X_sim, Vol_true, VaR_true))) > 0) {
      cat("Simulation produced NAs")
      stop()
    }
    
    # save truths of first simulation
    if (save_results==TRUE && sim==1) {
      sim_df = data.frame("X_sim"=X_sim, "Vol_true"=Vol_true)
      sim_df = cbind(sim_df, VaR_true)
      write.csv(sim_df, paste0("Simulations/", truth_name, "_truth.csv"))
    }
    
    # fit models to simulated data
    errors = foreach(j = 1:n_models, .combine="rbind") %do% {
      
      fit = GARCH_fit(X=X_sim, args=models[[j]])
      
      # extract estimated volatility & VaR
      Vol_pred      = fit$Vol_pred
      Vol_pred_mean = Vol_pred$mean
      Vol_pred_med  = Vol_pred$med
      VaR_pred      = fit$VaR_pred
      
      # compute errors
      Vol_mean_err = Vol_true - Vol_pred$mean
      Vol_med_err = Vol_true - Vol_pred$med
      VaR_err = VaR_true - VaR_pred
      
      # save predictions of first simulation
      if (save_results == TRUE && sim==1) {
        df_preds = cbind(Vol_pred, VaR_pred)
        write.csv(df_preds, paste0("Simulations/", truth_name, "_", names_models[j], "_pred.csv"))
      }
      
      vol_rmse = loss(Vol_mean_err, abs=FALSE, root=TRUE)
      vol_mad  = loss(Vol_med_err,  abs=TRUE,  root=FALSE)
      
      var_rmse = apply(VaR_err, 2, loss, abs=FALSE, root=TRUE)
      var_mad  = apply(VaR_err, 2, loss, abs=TRUE, root=FALSE)
      
      errors = c(vol_rmse, vol_mad, var_rmse, var_mad)
      names(errors) = c("Vol_RMSE", "Vol_MAD",
                        paste0("VaR_", truth$VaR_q, "_RMSE"), paste0("VaR_", truth$VaR_q, "_MAD"))
      
      # compute square of errors (allows to later compute std dev)
      errors_sq = errors^2
      names(errors_sq) = paste0(names(errors), "_Sq")
      
      errors_model = c(errors, errors_sq)
      
      rm(fit, Vol_pred, Vol_pred_mean, Vol_pred_med, VaR_pred, Vol_mean_err, Vol_med_err, VaR_err, df_preds, errors, errors_sq)
      gc()
      
      errors_model
    }
    errors
    
  }
  
  if (n_models == 1) {
    result = as.data.frame(t(result))
  } else {
    result = as.data.frame(result)
  }
  
  rownames(result) = paste(names_models)
  
  errors_sum = result[, 1:(ncol(result)/2)]    # sum of errors
  errors_sum_sq = result[, (ncol(result)/2+1):ncol(result)] # sum of squared errors
  
  # compute mean errors
  errors_mean = errors_sum/n_sim
  
  # compute standard deviations
  errors_sd = sqrt( (errors_sum_sq - 1/n_sim*errors_sum^2) / (n_sim-1) )
  colnames(errors_sd) = paste0(names(errors_mean), "_StD")
  
  errors = cbind(errors_mean, errors_sd)[order(c(seq_along(errors_mean), seq_along(errors_sd)))]
  
  if (save_results == TRUE) {
    write.csv(errors, paste0("Simulations/", truth_name, "_Results.csv"))
  }
  
  cat(paste(truth_name, "done!\n\n"))

  stopCluster(cluster)
  
  end = Sys.time()
  time = end-start
  time
  
}
  
