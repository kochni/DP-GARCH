
## packages

source("helpers.R")

packages_visual = c(
  "ggplot2",   # plotting
  "ggpubr",    # ggplot customization
  #"ggstance",  # horizontal barplot
  "scales"
)

lapply(packages_visual, library, character.only = TRUE)

## plot theme
#theme = theme_bw()
theme = theme_light(base_size=20)
#theme = theme_classic()
theme_set(theme)

## plot of true vs. predicted values

res_plots = function(process, models, quant, save=FALSE) {
  
  if (models[1]=="all") {
    models = c("GARCH", "Mix-GARCH",  "MS-GARCH",
               "Bay-GARCH", "BayMix-GARCH", "BayMS-GARCH",
               "GARCH(3,3)", "T-GARCH", "RF-GARCH")
  }
  
  n_models = length(models)
  
  if (quant == "Vol") {
    col_truth = 3
    col_pred = 2
  } else if (quant == "VaR1") {
    col_truth = 4
    col_pred = 4
  } else if (quant == "VaR2") {
    col_truth = 5
    col_pred = 5
  }
  
  path = paste0("./Simulations/", process)
  truth_file = paste0(path, "/", process, "_truth.csv")
  truth = read.csv(truth_file)[-c(1:10), col_truth] # 1% VaR
  
  data = list()
  
  for (j in 1:n_models) {
    
    pred_file = paste0(path, "/", process, "_", models[j], "_pred.csv")
    pred = read.csv(pred_file)[-c(1:10), col_pred]   # 1% VaR
    
    data[[j]] = data.frame("truth"=truth, "pred"=pred)
  }
  
  res_plots = list()
  range_low = min(unlist(data))
  range_high = max(unlist(data))
  
  n_col = ceiling(sqrt(n_models))
  n_row = ceiling(n_models/n_col)
  
  # create individual plots
  for (j in 1:n_models) {
    
    if ((j-1) %% n_col == 0) {
      y_lab = "predicted"
    } else {
      y_lab = ""
    }
    
    if ((j+n_col) > n_models) {
      x_lab = "true"
    } else {
      x_lab = ""
    }
    
    res_plots[[j]] = ggplot(data[[j]], aes(x=truth, y=pred)) +
      geom_point() +
      geom_abline(intercept=0, slope=1, linetype=1, col="red", linewidth=0.8) +
      labs(title = models[j]) +
      xlab(x_lab) +
      ylab(y_lab) +
      scale_x_log10(limits=c(range_low, range_high)) +
      scale_y_log10(limits=c(range_low, range_high))
    
  }
  
  # arrange in grid
  res_grid = ggarrange(plotlist=res_plots, nrow=n_row, ncol=n_col)
  
  # add title & caption
  #res_grid = annotate_figure(res_grid, top = text_grob("Value-at-Risk: true vs. fitted", face = "bold", size = 16),
  #                           bottom = text_grob(paste("Process:", process)))
  
  # save
  #png("./Plots/ResidualPlots.png", res_grid)
  
  if (save == TRUE) {
    filename = paste0("./Plots/", process, "_", quant, "-ResidualsPlot.png")
    png(filename, width=n_col*2000, height=n_row*1500, res=300)
    plot(res_grid)
    dev.off()
  }
  
  
  print(res_grid)
}

## create residual plots
models = c("GARCH", "Mix-GARCH", "MS-GARCH", "Bay-GARCH", "BayMix-GARCH", "BayMS-GARCH")
models = "all"

res_plots(process="GJR(N)", models=models, quant="Vol", save=TRUE)
res_plots(process="GJR(t)", models=models, quant="Vol", save=TRUE)
res_plots(process="Exp(N)", models=models, quant="Vol", save=TRUE)
res_plots(process="Exp(t)", models=models, quant="Vol", save=TRUE)
res_plots(process="MS(N)",  models=models, quant="Vol", save=TRUE)
res_plots(process="MS(t)",  models=models, quant="Vol", save=TRUE)

res_plots(process="GJR(N)", models=models, quant="VaR1", save=TRUE)
res_plots(process="GJR(t)", models=models, quant="VaR1", save=TRUE)
res_plots(process="Exp(N)", models=models, quant="VaR1", save=TRUE)
res_plots(process="Exp(t)", models=models, quant="VaR1", save=TRUE)
res_plots(process="MS(N)",  models=models, quant="VaR1", save=TRUE)
res_plots(process="MS(t)",  models=models, quant="VaR1", save=TRUE)

res_plots(process="GJR(N)", models=models, quant="VaR2", save=TRUE)
res_plots(process="GJR(t)", models=models, quant="VaR2", save=TRUE)
res_plots(process="Exp(N)", models=models, quant="VaR2", save=TRUE)
res_plots(process="Exp(t)", models=models, quant="VaR2", save=TRUE)
res_plots(process="MS(N)",  models=models, quant="VaR2", save=TRUE)
res_plots(process="MS(t)",  models=models, quant="VaR2", save=TRUE)



## plot VaR of innovation distributions
VaR_plot = function(save=FALSE) {
  
  grid = seq(0, 0.1, length=1000)
  df = data.frame("alpha"=rep(grid, each=2),
                  "VaR"  =rep(NA, 2*length(grid)),
                  "Distr"=rep(c("Norm", "t (5)"), 1000))
  df$VaR[df$Distr=="Norm"]   = (-1)*qnorm(p=grid, mean=0, sd=1)
  df$VaR[df$Distr=="t (5)"] = (-1)*q_std_t(p=grid, mu=0, sd=1, df=5)
  
  plot = ggplot(data=df, aes(x=alpha, y=VaR, color=Distr)) +
    geom_line() +
    geom_vline(xintercept=0.001, lwd=0.25, linetype="dashed") +
    geom_vline(xintercept=0.01, lwd=0.25, linetype="dashed") +
    geom_vline(xintercept=0.05, lwd=0.25, linetype="dashed") +
    xlab(expression(alpha)) +
    ylab(expression(alpha*"-VaR")) +
    scale_x_reverse()
  
  print(plot)
  
  if (save == TRUE) {
    filename = paste0("./Plots/VaRPlot.png")
    png(filename, width=1500, height=1000, res=200)
    plot(plot)
    dev.off()
  }
}







## barplot of errors

#plot_errors = function(errors_Vol, errors_VaR, true_process, true_err) {

# reshape for ggplot
#  errors_Vol_gg = errors_Vol
#  errors_VaR_gg = errors_VaR

#  errors_Vol_gg["Model"] = model_names
#  errors_VaR_gg["Model"] = model_names

#  errors_Vol_gg = errors_Vol_gg %>%
#    pivot_longer(cols = c("Vol_RMSE", "Vol_MAD"),
#                 names_to = "Metric",
#                 values_to = "Error")
#  errors_Vol_gg = data.frame(errors_Vol_gg)
#  errors_Vol_gg$Metric = errors_Vol_gg$Metric %>% stringr::str_remove("Vol_")
#  errors_Vol_gg["SD"] = errors_Vol_gg$Vol_RMSE_sd
#  errors_Vol_gg$SD[seq_len(nrow(errors_Vol_gg)) %% 2 == 0] = errors_Vol_gg$Vol_MAD_sd[seq_len(nrow(errors_Vol_gg)) %% 2 == 0]
#  errors_Vol_gg = errors_Vol_gg[-c(1:2)]

#  errors_VaR_gg = errors_VaR_gg %>%
#    pivot_longer(cols = c("VaR_RMSE", "VaR_MAD"),
#                 names_to = "Metric",
#                 values_to = "Error")
#  errors_VaR_gg = data.frame(errors_VaR_gg)
#  errors_VaR_gg$Metric = errors_VaR_gg$Metric %>% stringr::str_remove("VaR_")
#  errors_VaR_gg["SD"] = errors_VaR_gg$VaR_RMSE_sd
#  errors_VaR_gg$SD[seq_len(nrow(errors_VaR_gg)) %% 2 == 0] = errors_VaR_gg$VaR_MAD_sd[seq_len(nrow(errors_VaR_gg)) %% 2 == 0]
#  errors_VaR_gg = errors_VaR_gg[-c(1:2)]

# plot result

#  vol_plot = ggplot(errors_Vol_gg, aes(y=fct_inorder(Model), x=Error, fill=Metric)) +
#    geom_barh(stat = "identity",
#              color = "black",
#              position = "dodge") +
#    geom_errorbar(aes(xmin=Error-SD, xmax=Error+SD), width=.2,
#                  position=position_dodge(.9)) +
#    labs(title = "Volatility Estimation Error",
#         subtitle = paste("Truth:", true_process, "w", true_err, "innovations"),
#         caption = paste0("#Sim = ", n_sim, "; n = ", len)) +
#    xlab("Estimation Error") +
#    ylab("") +
#    theme(legend.position="right") +
#    scale_fill_brewer(palette="Paired")

#  var_plot = ggplot(errors_VaR_gg, aes(y=fct_inorder(Model), x=Error, fill=Metric)) +
#    geom_barh(stat = "identity",
#              color = "black",
#              position = "dodge") +
#    geom_errorbar(aes(xmin=Error-SD, xmax=Error+SD), width=.2,
#                  position=position_dodge(.9)) +
#    labs(title = "VaR Estimation Error",
#         subtitle = paste("Truth:", true_process, "w", true_err, "innovations"),
#         caption = paste0("#Sim = ", n_sim, "; n = ", len)) +
#    xlab("Estimation Error") +
#    ylab("") +
#    theme(legend.position="right") +
#    scale_fill_brewer(palette="Paired")

#  print(vol_plot)
#  print(var_plot)
#}
