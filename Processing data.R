#####################
## Process data from HCC

library(MASS)
library(dplyr)
library(tidyr)
library(statmod)
library(Matrix)
library(matrixcalc)
library(expm)
library(splines)
library(splines2)
library(ggplot2)


setwd('C:/Users/yuezh/OneDrive - University of Nebraska Medical Center/Research/Joint_model/joint_likelihood')
mode = 3
source('Initial parameter.R')
source('Numeric Diff.R')

row_names <- c("Mean_coef", "Mean_Bias", "Std_coef", "Mean_std","CP","2.5%pct","97.5%pct")
col_names <- colnames(parameters)
num_para <- length(parameters)

result = matrix(0,nrow=length(row_names),ncol = num_para,
                dimnames = list(row_names,col_names))

csv_dir <- "C:/Users/yuezh/Research/Joint_model/Result"

all_csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)

# Separate files by type
coef_files <- all_csv_files[grepl("coef", all_csv_files)]
std_files <- all_csv_files[grepl("std", all_csv_files)]
CP_files <- all_csv_files[grepl("CP_i", all_csv_files)]
time_files <- all_csv_files[grepl("time", all_csv_files)]
# data_files <- all_csv_files[grepl("data", all_csv_files)]
causal_S_files <- all_csv_files[grepl("causal_S", all_csv_files)]
true_NE_S_files <-all_csv_files[grepl("true_NE_S", all_csv_files)]
sd_NE_S_files <-all_csv_files[grepl("sd_NE_S", all_csv_files)]
CP_NE_S_files <-all_csv_files[grepl("CP_NE_S", all_csv_files)]

causal_M_files <- all_csv_files[grepl("causal_M", all_csv_files)]
true_NE_M_files <-all_csv_files[grepl("true_NE_M", all_csv_files)]
sd_NE_M_files <-all_csv_files[grepl("sd_NE_M", all_csv_files)]
CP_NE_M_files <-all_csv_files[grepl("CP_NE_M", all_csv_files)]
# Combine 'mean' files
coef <- coef_files %>%
  lapply(read.csv) %>%    # Read each file
  bind_rows()             # Combine rows

trials = nrow(coef)
# Combine 'std' files
stddev <- std_files %>%
  lapply(read.csv) %>%
  bind_rows()
# Coverage Probability
CP <- CP_files %>%
  lapply(read.csv) %>%
  bind_rows()
# Running time
time <- time_files %>%
  lapply(read.csv) %>%
  bind_rows()
# Causal mediation
causal_M <- causal_M_files %>%
  lapply(read.csv) %>%
  bind_rows()

true_NE_M <- true_NE_M_files %>%
  lapply(read.csv) %>%
  bind_rows()

sd_NE_M <- sd_NE_M_files %>%
  lapply(read.csv) %>%
  bind_rows()

CP_NE_M <- CP_NE_M_files %>%
  lapply(read.csv) %>%
  bind_rows()

causal_S <- causal_S_files %>%
  lapply(read.csv) %>%
  bind_rows()

true_NE_S <- true_NE_S_files %>%
  lapply(read.csv) %>%
  bind_rows()

sd_NE_S <- sd_NE_S_files %>%
  lapply(read.csv) %>%
  bind_rows()

CP_NE_S <- CP_NE_S_files %>%
  lapply(read.csv) %>%
  bind_rows()
# Data
# full_data <- data_files %>%
#   lapply(read.csv)
# 
# full_data2 <- bind_rows(
#   lapply(seq_along(full_data), function(i) {
#     full_data[[i]] |>
#       mutate(dataset_id = i)
#   })
# )

# full_data3 <- full_data2 %>%
#   group_by(dataset_id, id) |>               # unique person = (dataset_id, id)
#   mutate(id = cur_group_id()) |>        # assign 1..N IDs
#   ungroup() |> 
#   select(-dataset_id)


#Draw piecewise/spline hazard plot
if (mode == 2){
  hazard_files <- all_csv_files[grepl("hazard", all_csv_files)]
  time_grid <- seq(0, J, 0.01)
  plot_matrix <- matrix(0, nrow = length(time_grid), ncol = trials)
  j = 1
  for (i in hazard_files){
    data.i = read.csv(i)
    # knots = data.i$knots
    hazards = data.i$hazard
    # cumulative_hazard <- cumsum(hazards * c(0,diff(knots)))
    time_expand <- sort(unique(c(data.i$knots, time_grid)))
    # hazard_expanded <- full_join(data.i, time_expanded)
    hazard_df <- data.frame(time = time_expand)
    hazard_df$hazard = NA
    hazard_df$hazard[hazard_df$time %in% data.i$knots] <- hazards
    hazard_df <- hazard_df %>%
      fill(everything(), .direction = "down")
    cumulative_hazard <- cumsum(c(0,diff(hazard_df$time))*hazard_df$hazard)
    hazard_df <- cbind(hazard_df, cumulative_hazard)
    hazard_plot <- hazard_df[hazard_df$time %in% time_grid,]$cumulative_hazard
    plot_matrix[, j] <- hazard_plot
    j <- j + 1
  }
  hazard_mean <- apply(plot_matrix, 1, mean)
  hazard_sd <- apply(plot_matrix, 1, sd)
  hazard_2.5 <- apply(plot_matrix, 1, function(x) quantile(x, probs=0.025))
  hazard_97.5 <- apply(plot_matrix, 1, function(x) quantile(x, probs=0.975))
  plot(time_grid, hazard_2.5, type = "l", col = 'red',ylab = 'hazard', xlab = "time", ylim = c(0, max(hazard_97.5)))
  lines(time_grid, hazard_97.5, lty = 1, col = 'blue')
  lines(time_grid, hazard_mean, lty = 1, col = 'green')
  lines(time_grid, (time_grid/par[1])^par[2], lty = 1)
  legend("bottomright", legend = c("2.5% percentile", "97.5% percentile", "Mean", "True"), col = c("red", "blue", "green", "black"), lwd = 2)
}
if (mode == 3){
  hazard_files <- all_csv_files[grepl("hazard", all_csv_files)]
  time_grid <- seq(0, J -0.8, 0.01)
  plot_matrix <- matrix(0, nrow = trials, ncol = length(time_grid))
  for (i in 1:length(hazard_files)){
    data.i = read.csv(hazard_files[i])
    knots = unlist(data.i %>% dplyr::select(contains('knots')))
    alpha = data.i %>% dplyr::select(contains('alpha'))
    # boundary = as.numeric(data.i %>% dplyr::select(contains('boundary')))
    boundary_knots = c(0, J-0.8)
    coef.i = read.csv(coef_files[i])
    theta.a.hat = as.matrix(coef.i %>% dplyr::select(contains("theta.a")))
    Sigma.a.hat = matrix(0, nrow=K,ncol=K)
    upper_vec_a = as.matrix(coef.i %>% dplyr::select(contains("Sigma.a")))
    Sigma.a.hat[upper.tri(Sigma.a.hat,diag = TRUE)] = upper_vec_a
    if (K>1){
      Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
    }
    # hazards = data.i$hazard
    # cumulative_hazard <- cumsum(hazards * c(0,diff(knots)))
    # time_expand <- sort(unique(c(data.i$knots, time_grid)))
    # # hazard_expanded <- full_join(data.i, time_expanded)
    # hazard_df <- data.frame(time = time_expand)
    # hazard_df$hazard = NA
    # hazard_df$hazard[hazard_df$time %in% data.i$knots] <- hazards
    # hazard_df <- hazard_df %>%
    #   fill(everything(), .direction = "down")
    cumulative_hazard <- spline_cumulative_hazard(time_grid, knots, t(alpha), boundary_knots, degree)
    # cumulative_hazard <- evaluate_Ispline(time_grid, knots, t(alpha), boundary_knots = c(0, boundary), degree)
    plot_matrix[i, ] <- t(cumulative_hazard) # exp(as.numeric(theta.a.hat%*%Sigma.a.hat%*%t(theta.a.hat)))
  }
  hazard_mean <- apply(plot_matrix, 2, mean)
  # hazard_sd <- apply(plot_matrix, 2, sd)
  hazard_2.5 <- apply(plot_matrix, 2, function(x) quantile(x, probs=0.025))
  hazard_97.5 <- apply(plot_matrix, 2, function(x) quantile(x, probs=0.975))
  # plot(time_grid, hazard_2.5, type = "l", col = 'red',ylab = 'hazard', xlab = "time", ylim = c(0, max(hazard_97.5)))
  # lines(time_grid, hazard_97.5, lty = 1, col = 'blue')
  # lines(time_grid, hazard_mean, lty = 1, col = 'green')
  # lines(time_grid, (time_grid/par[1])^par[2], lty = 1)
  # legend("bottomright", legend = c("2.5% percentile", "97.5% percentile", "Mean", "True"), col = c("red", "blue", "green", "black"), lwd = 2)
  # ggplot2
  true_hazard = (time_grid/par[1])^par[2] # as.numeric(exp(t(theta.a)%*%Sigma.a%*%theta.a))
  plot_df <- data.frame(
    time = time_grid,
    mean = hazard_mean,
    lower_ci = hazard_2.5,
    upper_ci = hazard_97.5,
    true = true_hazard
  )

  # --- NEW: Generate the plot using ggplot2 ---
  ggplot(plot_df, aes(x = time)) +
    # Add the 95% confidence interval as a shaded ribbon
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "skyblue", alpha = 0.5) +

    # Add the mean cumulative hazard line
    geom_line(aes(y = mean, color = "Mean"), linewidth = 1) +

    # Add the true cumulative hazard line
    geom_line(aes(y = true, color = "True"), linetype = "dashed", linewidth = 1) +

    # Customize labels, colors, and title
    labs(
      title = "Simulated vs. True Cumulative Hazard",
      x = "Time",
      y = "Cumulative Hazard",
      color = "Legend" # Title for the legend
    ) +
    scale_color_manual(
      name = NULL, # Hide the legend title if you want
      values = c("Mean" = "blue", "True" = "red"),
      labels = c("Mean Estimate", "True Value")
    ) +

    # Apply a clean theme
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "bottom",
      legend.text = element_text(size = 12)
    )
}
# Combine results to one matrix
result[1,]=apply(coef,2,mean,na.rm=TRUE)
result[2,]=apply(coef,2,mean,na.rm=TRUE) - as.matrix(parameters)
result[3,]=apply(coef,2,sd,na.rm=TRUE)
result[4,]=apply(stddev,2,mean,na.rm=TRUE)
result[5,]=apply(CP,2,mean,na.rm=TRUE)
CI = apply(coef,2, function(x) quantile(x,probs = c(0.025,0.975), na.rm=TRUE))

result[6,]=CI[1,]
result[7,]=CI[2,]

if (mode == 1){
  path_out = 'C:/Users/yuezh/OneDrive - University of Nebraska Medical Center/Research/Joint_model/Result/Constant/'
  filename = paste(path_out,"K = ",K," p = ",p," n = ",n.id," trials = ",trials," ",Sys.Date(),".csv",sep="")
}
if (mode == 2){
  path_out = 'C:/Users/yuezh/OneDrive - University of Nebraska Medical Center/Research/Joint_model/Result/Piecewise/'
  filename = paste(path_out,"K = ",K," p = ",p," n = ",n.id," trials = ",trials," knots = ",knots," ",Sys.Date(),".csv",sep="")
}
if (mode == 3){
  path_out = 'C:/Users/yuezh/OneDrive - University of Nebraska Medical Center/Research/Joint_model/Result/B-spline/'
  filename = paste(path_out,"K = ",K," p = ",p," n = ",n.id," trials = ",trials," B_spline,knots = ",num_knots," ",Sys.Date(),".csv",sep="")
}
write.csv(result,filename)

# run_time=mean(as.matrix(time),na.rm=TRUE)

######################################
# Causation
source('Causal Mediation.R')
true_NDE_S = mean(true_NE_S$NDE)
true_NIE_S = mean(true_NE_S$NIE)

bias_NDE_S = mean(causal_S[,1]-true_NDE_S)
bias_NIE_S = mean(causal_S[,2]-true_NIE_S)

sd_NDE_S = sd(causal_S[,1])
sd_NIE_S = sd(causal_S[,2])

asd_NDE_S = mean(sd_NE_S[,1])
asd_NIE_S = mean(sd_NE_S[,2])

CP_NDE_S = mean(CP_NE_S[,1])
CP_NIE_S = mean(CP_NE_S[,2])


true_NDE_M = apply(true_NE_M[, 1:K], 2, mean)
true_NIE_M = apply(true_NE_M[, (K+1):(2*K)], 2, mean)

bias_NDE_M = apply(causal_M[,1:K], 2, mean)-true_NDE_M
bias_NIE_M = apply(causal_M[,(K+1):(2*K)], 2, mean)-true_NIE_M

sd_NDE_M = apply(causal_M[,1:K], 2, sd)
sd_NIE_M = apply(causal_M[,(K+1):(2*K)], 2, sd)

asd_NDE_M = apply(sd_NE_M[,1:K], 2, mean)
asd_NIE_M = apply(sd_NE_M[,(K+1):(2*K)], 2, mean)

CP_NDE_M = apply(CP_NE_M[,1:K], 2, mean)
CP_NIE_M = apply(CP_NE_M[,(K+1):(2*K)], 2, mean)

row_names <- c("Est_NE", "True_NE", "Bias_NE", "Std_NE", "Mean_std_NE", "CP")
col_names <- c("NDE_S","NIE_S",paste0("NDE_M",1:K),paste0("NIE_M",1:K))

result_NE = matrix(0,nrow=length(row_names),ncol = length(col_names),
                dimnames = list(row_names,col_names))

result_NE[1,1] = mean(causal_S[,1])
result_NE[1,2] = mean(causal_S[,2])
result_NE[1,3:(2+K)] = t(apply(causal_M[,1:K], 2, mean))
result_NE[1,(3+K):(2+2*K)] = t(apply(causal_M[,(K+1):(2*K)], 2, mean))
result_NE[2,1] = true_NDE_S
result_NE[2,2] = true_NIE_S
result_NE[2,3:(2+K)] = t(true_NDE_M)
result_NE[2,(3+K):(2+2*K)] = t(true_NIE_M)
result_NE[3,1] = bias_NDE_S
result_NE[3,2] = bias_NIE_S
result_NE[3,3:(2+K)] = t(bias_NDE_M)
result_NE[3,(3+K):(2+2*K)] = t(bias_NIE_M)
result_NE[4,1] = sd_NDE_S
result_NE[4,2] = sd_NIE_S
result_NE[4,3:(2+K)] = t(sd_NDE_M)
result_NE[4,(3+K):(2+2*K)] = t(sd_NIE_M)
result_NE[5,1] = asd_NDE_S
result_NE[5,2] = asd_NIE_S
result_NE[5,3:(2+K)] = t(asd_NDE_M)
result_NE[5,(3+K):(2+2*K)] = t(asd_NIE_M)
result_NE[6,1] = CP_NDE_S
result_NE[6,2] = CP_NIE_S
result_NE[6,3:(2+K)] = t(CP_NDE_M)
result_NE[6,(3+K):(2+2*K)] = t(CP_NIE_M)

filename = paste(path_out,"Causal ", "K = ",K," p = ",p," n = ",n.id," trials = ",trials," B_spline,knots = ",num_knots," ",Sys.Date(),".csv",sep="")
write.csv(result_NE,filename)

file.remove(all_csv_files)


