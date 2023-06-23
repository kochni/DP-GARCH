
## Nonparametric Bayesian GARCH Mixture
## with various Laws, incl. Dirichlet process


DP_GARCH = function(X, law, burnin, n_post, thin, VaR_q) {
  
  if (law == "DP") {
    Option <- 4
  } else if (law == "PD") {
    Option <- 3
  } else if (law == "NGG") {
    Option <- 2
  } else if (law == "None") {
    Option <- 1
  }
  
  Return.R = X
  
  ### Start: Code Lau & Cripps
  ### (be confident about changes!)
  
  if (Option==1) {
    filename.Rdata <- "SP500-NON.Rdata";
    filename.log   <- "SP500-NON.log";
    mix.ar <- FALSE; mix.garch <- FALSE; PD <- FALSE; AMDR <- FALSE; NON <- TRUE;
    set.PARA.alpha <- 0.5;
    set.PARA.theta <- 1.0;
    set.PARA.b     <- 0.25;
    set.PARA.q     <- 1.0;
  } else if (Option==2) {
    filename.Rdata <- "SP500-NGG.Rdata";
    filename.log   <- "SP500-NGG.log";
    mix.ar <- FALSE; mix.garch <- TRUE; PD <- FALSE; AMDR <- FALSE; NON <- FALSE;
    set.PARA.alpha <- 0.5;
    set.PARA.theta <- 1.0;
    set.PARA.b     <- 0.25;
    set.PARA.q     <- 1.0;
  } else if (Option==3) {
    filename.Rdata <- "SP500-XPD.Rdata";
    filename.log   <- "SP500-XPD.log";
    mix.ar <- FALSE; mix.garch <- TRUE; PD <- TRUE; AMDR <- FALSE; NON <- FALSE;
    set.PARA.alpha <- 0.5;
    set.PARA.theta <- 1.0;
    set.PARA.b     <- 1.0;
    set.PARA.q     <- 0.6769;
  } else if (Option==4) {
    filename.Rdata <- "SP500-XDP.Rdata";
    filename.log   <- "SP500-XDP.log";
    mix.ar <- FALSE; mix.garch <- TRUE; PD <- TRUE; AMDR <- FALSE; NON <- FALSE;
    set.PARA.alpha <- 0.0;
    set.PARA.theta <- 1.0;
    set.PARA.b     <- 1.0;
    set.PARA.q     <- 2.3538;
  } else {
    filename.Rdata <- "SP500-TEST.Rdata";
    filename.log   <- "SP500-TEST.log";
    mix.ar <- FALSE; mix.garch <- FALSE; PD <- TRUE; AMDR <- FALSE; NON <- FALSE;
    set.PARA.alpha <- 0.0;
    set.PARA.theta <- 1.0;
    set.PARA.b     <- 1.0;
    set.PARA.q     <- 2.3538;
  }
  
  
  ### Posterior Sampling ## Begin
  
  
  Y             <- Return.R; # **************************************************
  
  n.size        <- length(Y);
  ar.order      <- 1;
  garch.a.order <- 1;
  garch.b.order <- 1;
  max.order     <- max(ar.order,garch.a.order,garch.b.order);
  
  dim.ar1      <- 1 + ar.order;
  dim.ar2      <- dim.ar1 * dim.ar1;
  dim.garch1   <- 1 + garch.a.order + garch.b.order;
  dim.garch2   <- dim.garch1 * dim.garch1;
  dim.garch.a1 <- 1 + garch.a.order;
  dim.garch.a2 <- dim.garch.a1 * dim.garch.a1;
  dim.garch.b1 <- garch.b.order;
  dim.garch.b2 <- dim.garch.b1 * dim.garch.b1;
  
  Y.series               <- c(rep(0,max.order),Y);
  m.series               <- rep(0,n.size+max.order);
  e.series               <- rep(0,n.size+max.order);
  h2.series              <- rep(var(Y),n.size+max.order);
  proposal.m.series      <- m.series;
  proposal.e.series      <- e.series;
  proposal.h2.series     <- h2.series;
  
  ave.h2.series <- rep(0.0,n.size+max.order);
  ave.h1.series <- rep(0.0,n.size+max.order);
  ave.m1.series <- rep(0.0,n.size+max.order);
  
  ORDERs              <- c(max.order,ar.order,garch.a.order,garch.b.order);
  DIMs                <- c(dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2);
  PARA.ar             <- rep(c(0,rep(0.0,ar.order)),n.size);
  PARA.garch          <- rep(c(var(Y),rep(0.1,garch.a.order),rep(0.1,garch.b.order)),n.size);
  proposal.PARA.garch <- PARA.garch;
  proposal.PARA.ar    <- PARA.ar;
  
  mean.ar       <- 0.0;
  vari.ar       <- 10.0;
  mean.garch    <- 0.0;
  vari.garch    <- 10.0;
  
  PRIOR.ar      <- rep(c(mean.ar,vari.ar),dim.ar1);
  PRIOR.garch   <- rep(c(mean.garch,vari.garch),dim.garch1);
  
  AMDR.AVE.ARs        <- rep(rep(mean.ar,dim.ar1),n.size);
  AMDR.COV.ARs        <- rep(diag(vari.ar,dim.ar1),n.size);
  AMDR.AVE.GARCHs     <- rep(rep(mean.garch,dim.garch1),n.size);
  AMDR.COV.GARCHs     <- rep(diag(vari.garch,dim.garch1),n.size);
  
  # C.n      <- label;
  # label    <- apply(outer(C.n,1:max(C.n),"=="),2,sum);
  # np.n     <- max(C.n);
  
  PARA.alpha <- set.PARA.alpha;
  PARA.theta <- set.PARA.theta;
  PARA.b     <- set.PARA.b;
  PARA.q     <- set.PARA.q;
  PARAs      <- c(PARA.alpha,PARA.theta,PARA.b,PARA.q);
  
  C.n      <- 1:n.size;
  #C.n      <- c(rep(1, ceiling(n.size/2)), rep(2, floor(n.size/2)));
  e.n      <- c(apply(outer(C.n,1:max(C.n),"=="),2,sum),rep(0,n.size-max(C.n)));
  np.n     <- max(C.n);
  
  C.global      <- rep(1,n.size);
  e.global      <- c(apply(outer(C.global,1:max(C.n),"=="),2,sum),rep(0,n.size-max(C.global)));
  np.global     <- max(C.global);
  
  C.temp      <- C.global
  e.temp      <- e.global
  np.temp     <- np.global
  
  n.warmup    <- burnin ;       # **************************************************
  n.iter      <- n_post * thin; # **************************************************
  n.runs      <- n.iter + n.warmup
  
  all.garch      <- array(0,c(n.iter,n.size*dim.garch1));
  all.ar         <- array(0,c(n.iter,n.size*dim.ar1));
  all.np         <- array(0,c(n.iter,1));
  all.likelihood <- array(0,c(n.iter,1));
  all.C.n        <- array(0,c(n.iter,n.size));
  
  pred.X   <- seq(-5,5,0.1);
  pred.n   <- length(pred.X);
  pred.Y   <- rep(0,pred.n);
  pred.lnY <- rep(0,pred.n);
  
  Option    <- 1;
  max.R     <- 1000;
  
  approx.N  <- 100;
  
  if ((mix.ar==FALSE)&&(mix.garch==FALSE)) {
    C.n      <- rep(1,n.size);
    e.n      <- c(n.size,rep(0,n.size-1));
    np.n     <- 1;
  }
  
  if (PD) {
    sample.partition  <- "sample_PD_GIB_1STEP_MGARCH";
    pred.function     <- "logLikelihood_pred_PD";
  } else {
    sample.partition  <- "sample_NGG_GIB_2STEP_MGARCH";
    pred.function     <- "logLikelihood_pred_NGG";
  }
  
  if (NON) {
    sample.partition  <- "sample_NON_GIB_1STEP_MGARCH";
    pred.function     <- "logLikelihood_pred_NON";
  }
  
  
  if (AMDR) {
    sample.garch.para <- "sample_GARCH_AMDR_MCMC";
    sample.ar.para    <- "sample_AR_AMDR_MCMC";
  } else {
    sample.garch.para <- "sample_GARCH_pq_MCMC";
    sample.ar.para    <- "sample_AR_p_MCMC";
  }
  
  if(.Platform$OS.type == "unix") {
    dyn.text <- "GARCH.so"
  } else {
    dyn.text <- "GARCH.dll"
  }
  
  date1 <- Sys.time();
  date2 <- Sys.time();
  date3 <- Sys.time();
  
  # system("./lcompile")
  
  date0 <- Sys.time()
  
  
  dyn.load(dyn.text)
  system.time({
    
    for (i in 1:n.runs) {
      
      out <- .C(sample.partition,
                as.integer(C.n),as.integer(e.n),as.integer(n.size),
                as.integer(np.n),as.integer(ORDERs),as.integer(DIMs),
                as.double(PARAs),
                as.double(PARA.ar),as.double(proposal.PARA.ar),
                as.double(PARA.garch),as.double(proposal.PARA.garch),
                as.double(PRIOR.ar),as.double(PRIOR.garch),
                as.double(Y.series),
                as.double(m.series),as.double(proposal.m.series),
                as.double(e.series),as.double(proposal.e.series),
                as.double(h2.series),as.double(proposal.h2.series),
                as.integer(Option),as.integer(max.R),as.integer(approx.N),
                as.integer(mix.ar),as.integer(mix.garch)
      );
      C.n <- out[[1]]; e.n <- out[[2]]; np.n.past <- np.n; np.n <- out[[4]];
      PARA.ar <- out[[8]]; PARA.garch <- out[[10]];
      if (np.n.past>np.n) {
        AMDR.AVE.ARs[(np.n*dim.ar1+1):(np.n.past*dim.ar1)] <-
          rep(rep(mean.ar,dim.ar1),np.n.past-np.n);
        AMDR.COV.ARs[(np.n*dim.ar2+1):(np.n.past*dim.ar2)] <-
          rep(diag(vari.ar,dim.ar1),np.n.past-np.n);
        AMDR.AVE.GARCHs[(np.n*dim.garch1+1):(np.n.past*dim.garch1)] <-
          rep(rep(mean.garch,dim.garch1),np.n.past-np.n);
        AMDR.COV.GARCHs[(np.n*dim.garch2+1):(np.n.past*dim.garch2)] <-
          rep(diag(vari.garch,dim.garch1),np.n.past-np.n);
      }
      
      if (mix.garch) {
        C.temp <- C.n; e.temp <- e.n; np.temp <- np.n;
      } else {
        C.temp <- C.global; e.temp <- e.global; np.temp <- np.global;
      }
      out <- .C(sample.garch.para,
                as.integer(C.temp),as.integer(e.temp),as.integer(n.size),
                as.integer(np.temp),as.integer(ORDERs),as.integer(DIMs),
                as.double(PARA.ar),as.double(proposal.PARA.ar),
                as.double(PARA.garch),as.double(proposal.PARA.garch),
                as.double(PRIOR.ar),as.double(PRIOR.garch),
                as.double(Y.series),
                as.double(m.series),as.double(proposal.m.series),
                as.double(e.series),as.double(proposal.e.series),
                as.double(h2.series),as.double(proposal.h2.series),
                as.integer(Option),
                as.double(AMDR.AVE.ARs),as.double(AMDR.COV.ARs),
                as.double(AMDR.AVE.GARCHs),as.double(AMDR.COV.GARCHs),
                as.integer(i),as.integer(max.R)
      );
      PARA.garch <- out[[9]];
      if (mix.garch) {
        AMDR.AVE.GARCHs[1:(np.n*dim.garch1)] <- out[[23]][1:(np.n*dim.garch1)];
        AMDR.COV.GARCHs[1:(np.n*dim.garch2)] <- out[[24]][1:(np.n*dim.garch2)];
      } else {
        PARA.garch[1:(np.n*dim.garch1)] <- rep(PARA.garch[1:dim.garch1],np.n);
        AMDR.AVE.GARCHs[1:dim.garch1]   <- out[[23]][1:dim.garch1];
        AMDR.COV.GARCHs[1:dim.garch2]   <- out[[24]][1:dim.garch2];
      }
      
      # if (mix.ar) {
      # 	C.temp <- C.n; e.temp <- e.n; np.temp <- np.n;
      # 	} else {
      # 	C.temp <- C.global; e.temp <- e.global; np.temp <- np.global;
      # 	}
      # out <- .C(sample.ar.para,
      # 	as.integer(C.temp),as.integer(e.temp),as.integer(n.size),
      # 	as.integer(np.temp),as.integer(ORDERs),as.integer(DIMs),
      # 	as.double(PARA.ar),as.double(proposal.PARA.ar),
      # 	as.double(PARA.garch),as.double(proposal.PARA.garch),
      # 	as.double(PRIOR.ar),as.double(PRIOR.garch),
      # 	as.double(Y.series),
      # 	as.double(m.series),as.double(proposal.m.series),
      # 	as.double(e.series),as.double(proposal.e.series	),
      # 	as.double(h2.series),as.double(proposal.h2.series),
      # 	as.integer(Option),
      # 	as.double(AMDR.AVE.ARs),as.double(AMDR.COV.ARs),
      # 	as.double(AMDR.AVE.GARCHs),as.double(AMDR.COV.GARCHs),
      # 	as.integer(i),as.integer(max.R)
      # 	);
      # PARA.ar <- out[[7]];
      # if (mix.ar) {
      # 	AMDR.AVE.ARs[1:(np.n*dim.ar1)] <- out[[21]][1:(np.n*dim.ar1)];
      # 	AMDR.COV.ARs[1:(np.n*dim.ar2)] <- out[[22]][1:(np.n*dim.ar2)];
      # 	} else {
      #  	PARA.ar[1:(np.n*dim.ar1)] <- rep(PARA.ar[1:dim.ar1],np.n);
      # 	AMDR.AVE.ARs[1:dim.ar1] <- out[[21]][1:dim.ar1];
      # 	AMDR.COV.ARs[1:dim.ar2] <- out[[22]][1:dim.ar2];
      # 	}
      
      k <- i - n.warmup
      if (k>0 && (k-1) %% thin == 0) {
        out <- .C(pred.function,
                  as.integer(C.n),as.integer(e.n),as.integer(n.size),
                  as.integer(np.n),as.integer(ORDERs),as.integer(DIMs),
                  as.double(PARAs),
                  as.double(PARA.ar),	as.double(PARA.garch),
                  as.double(PRIOR.ar),as.double(PRIOR.garch),
                  as.double(Y.series),as.double(m.series),as.double(e.series),as.double(h2.series),
                  as.double(0),
                  as.double(pred.X),as.double(pred.lnY),as.integer(pred.n),
                  as.integer(mix.ar),as.integer(mix.garch)
        );
        pred.Y             <- pred.Y + exp(out[[18]]);
        # lines(pred.X,exp(out[[18]]),type="l",col=2,pch=16,lwd=2);
        ave.h2.series      <- ave.h2.series + out[[15]];
        ave.h1.series      <- ave.h1.series + sqrt(out[[15]]);
        ave.m1.series      <- ave.m1.series + out[[13]];
        all.garch[k,]      <- PARA.garch;
        all.ar[k,]         <- PARA.ar;
        all.np[k,]         <- np.n;
        all.likelihood[k,] <- out[[16]];
        all.C.n[k,]        <- C.n;
      }
      
      if (i%%100==0) {
        date3 <- Sys.time();
        cat("complete",i,"iterations...",format(date3-date2),"\n");
        print(t(rbind(
          table(C.n),
          array(PARA.ar[1:(np.n*(1+ar.order))],c(1+ar.order,np.n)),
          array(PARA.garch[1:(np.n*(1+garch.a.order+garch.b.order))],c(1+garch.a.order+garch.b.order,np.n))
        )));
        write(paste(
          filename.Rdata,"<-",
          paste("complete",i,"iterations...",format(date3-date2)),"Overall:",format(date3-date1)),
          filename.log,append=T);
        date2 <- Sys.time();
      }
      
    }
    
    ave.h2.series <- ave.h2.series / n.iter;
    ave.h1.series <- ave.h1.series / n.iter;
    ave.m1.series <- ave.m1.series / n.iter;
    pred.Y        <- pred.Y / n.iter;
    
    date3 <- Sys.time(); cat("Overall:",format(date3-date1),"\n");
    date1 <- Sys.time()
    
    write(paste(
      filename.Rdata,"<-",date1,
      paste("|",sep=""),
      format(date1-date0)),
      filename.log,append=T)
    
  })
  dyn.unload(dyn.text)
  
  ###---- Posterior Sampling End ----###
  ### End Code Lau & Cripps
  
  n = length(X)
  
  all.garch = as.data.frame(all.garch)
  rownames(all.garch) = paste0("i=", burnin+c(1:n_post))
  colnames(all.garch) = paste0(rep(c("omega_", "alpha_", "beta_"), n), rep(1:n, each=3))
  
  # save first (after burn-in) & last 50 posterior coefficient samples
  save_samples = rbind(head(all.garch, n=50), tail(all.garch, n=50))
  write.csv(save_samples, paste0("Simulations/", law, "PostSamples.csv"))

  # extract coefficients
  omega_ind = c(1:(3*n))[(c(1:(3*n))-1) %% 3 == 0]
  alpha_ind = c(1:(3*n))[(c(1:(3*n))-2) %% 3 == 0]
  beta_ind = c(1:(3*n))[(c(1:(3*n))-3) %% 3 == 0]
  
  post_samples_omega = all.garch[, omega_ind]
  post_samples_alpha = all.garch[, alpha_ind]
  post_samples_beta = all.garch[, beta_ind]
  
  ## compute volatilities & VaR
  vol = matrix(rep(NA, n*n_post), nrow=n_post, ncol=n)
  X_sample = matrix(rep(NA, n*n_post), nrow=n_post, ncol=n)
  
  for (i in 1:n_post) {
    
    omega = post_samples_omega[i,]
    alpha = post_samples_alpha[i,]
    beta = post_samples_beta[i,]
    
    vol[i,1] = omega[1]
    
    for (t in 2:n) {
      vol[i,t] = sqrt( omega[t] + alpha[t]*X[t-1]^2 + beta[t]*vol[i,t-1]^2 )
    }
    
    # sample 1 datapoint per time index for each posterior sample
    X_sample[i,] = rnorm(n=n, mean=0, sd=vol[i,])
  }
  
  vol_mean = colMeans(vol)
  vol_med = apply(vol, 2, median)
  Vol_pred = data.frame("mean"=vol_mean, "med"=vol_med)
  rownames(Vol_pred) = paste0("t=", 1:n)
  
  var = apply(X_sample, 2, quantile, probs=VaR_q)
  var = t(var)
  var = (-1)*var
  var = as.data.frame(var)
  colnames(var) = paste0(VaR_q*100, "%")
  rownames(var) = paste0("t=", 1:n)
  
  result = list("Vol_pred"=Vol_pred, "VaR_pred"=var)
  
  return(result)
}



