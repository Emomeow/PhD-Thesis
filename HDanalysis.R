## Real data analysis
###############################
## data cleaning and processing
library(dplyr)
library(tidyr)
library(MASS)
library(statmod)
library(Matrix)
library(matrixcalc)
library(expm)
library(splines)
library(splines2)

set.seed(1234)
data<-read.csv("C:/Users/yuezh/Research/Joint_model/data/20160219_predicthd_studywide_best_scan_release.csv")

data <- data[data$status == 'case',]

data<- data[, c("subject","visit","visit_diagnosis", "age1", "age","educ_yrs","gender","sex","ethnicity","race","cap","dx17", "sydigtot","stroopwo")]

data <- data[order(data$subject, data$visit),]


summary(data)
# delete rows all NA
data <- data[!apply(is.na(data), 1, all), ]

summary(data)
# After check that educ_yrs are missing for whole individual
data <- data[!is.na(data$educ_yrs),]

summary(data)
# remove those who diagnosed at first observation
data <- data |>
  group_by(subject) |>
  filter(first(visit_diagnosis) != "Diagnosed") |>
  ungroup()

# fill missing value for sydigtot, stroopwo and verflcor
data_notdiagnosed = data[data$visit_diagnosis == "",]

data_notdiagnosed <- data_notdiagnosed |>
  arrange(subject, visit) |>
  group_by(subject) |>
  fill(c(sydigtot, stroopwo), .direction = "downup") |>
  ungroup()

summary(data_notdiagnosed)
data_notdiagnosed[is.na(data_notdiagnosed$stroopwo)|is.na(data_notdiagnosed$sydigtot), ]
# 84 934 2593 2600 are completely missing
# no missing biomarker value now in not diagnosed data

data_diagnosed = data[data$visit_diagnosis == "Diagnosed",]
data_diagnosed <- data_diagnosed |>
  arrange(subject, visit) |>
  group_by(subject) |>
  fill(c(sydigtot, stroopwo), .direction = "updown") |>
  ungroup()
# There are some individuals who has no observation for some biomarkers
# after diagnosed, we need to delete them
summary(data_diagnosed)
data_diagnosed[is.na(data_diagnosed$stroopwo)|is.na(data_diagnosed$sydigtot), ]
#900 1012 are completely missing
data = rbind(data_diagnosed, data_notdiagnosed)

data <- data |>
  arrange(subject, visit) |>
  group_by(subject) |>
  filter(!subject %in% c(84, 900, 934, 1012, 2593, 2600)) |>
  ungroup()

summary(data)

missing_dx17 <- data |>
  group_by(subject) |>
  filter(any(is.na(dx17)) & any(visit_diagnosis == "Diagnosed")) |>
  ungroup()

not_missing_dx17 <- data |>
  group_by(subject) |>
  filter(!(any(is.na(dx17)) & any(visit_diagnosis == "Diagnosed"))) |>
  ungroup()
# Some subjects have dx17 missing before diagnosed
missing_dx17 <- missing_dx17 |>
  arrange(subject, visit) |>
  group_by(subject) |>
  mutate(
    first_diagnosed_visit = min(visit[visit_diagnosis == "Diagnosed"], na.rm = TRUE)
  ) |>
  filter(
    !(is.na(dx17) & visit == max(visit[visit < first_diagnosed_visit]))
  ) |>
  ungroup()

missing_dx17 <- missing_dx17 |>
  arrange(subject, visit) |>
  group_by(subject) |>
  mutate(
    first_diagnosed_visit = min(visit[visit_diagnosis == "Diagnosed"], na.rm = TRUE)
  ) |>
  filter(
    !(is.na(dx17) & visit == max(visit[visit < first_diagnosed_visit]))
  ) |>
  ungroup()

missing_dx17 <- missing_dx17 |>
  arrange(subject, visit) |>
  group_by(subject) |>
  mutate(
    first_diagnosed_visit = min(visit[visit_diagnosis == "Diagnosed"], na.rm = TRUE)
  ) |>
  filter(
    !(is.na(dx17) & visit == max(visit[visit < first_diagnosed_visit]))
  ) |>
  ungroup()

missing_dx17 <- missing_dx17 |>
  dplyr::select(-first_diagnosed_visit)
data <- rbind(missing_dx17, not_missing_dx17)

# those who has diagnose is event, others are censor
# Define V, U, visit, visits, 
data <- data |>
  arrange(subject, visit) |>
  group_by(subject) |>
  mutate(status = as.numeric(any(visit_diagnosis == "Diagnosed"))) |>
  ungroup()

data <- data |>
  group_by(subject) |>
  mutate(visittime = age - age1) |>
  ungroup()

case_data <- data |>
  group_by(subject) |>
  filter(status == 1) |>
  mutate(preetime = max(visittime[which(visit_diagnosis == "")]),
         etime = min(visittime[which(visit_diagnosis == "Diagnosed")]),
         visit = row_number(),
         visits = n(),
         befdiag = sum(visittime <= preetime),
         afterdiag = sum(visittime > preetime)) |>
  ungroup()

censor_data <- data |>
  group_by(subject) |>
  filter(status == 0) |>
  mutate(preetime = max(visittime),
         etime = max(visittime),
         visit = row_number(),
         visits = n(),
         befdiag = sum(visittime <= preetime),
         afterdiag = sum(visittime > preetime)) |>
  ungroup()

data <- rbind(case_data, censor_data)
data <- data |>
  filter(visits > 1)

# Define baseline cap and id
data <- data |>
  arrange(subject, visit) |>
  group_by(subject) |>
  mutate(base_cap = first(cap)/100) |>
  ungroup()

data <- data |>
  mutate(id = as.integer(factor(subject))) |>
  dplyr::select(-subject) |>
  relocate(id, .before = 1)

data <- data |>
  mutate(gender = as.numeric(factor(gender))-1,
         race = as.numeric(race == 5))

filename = 'C:/Users/yuezh/Research/Joint_model/data/cleaned data.csv'
write.csv(data, filename)
# sex: f=0 m=1
# Complete the construction of data

summary(data)


#############################
data = read.csv('C:/Users/yuezh/Research/Joint_model/data/cleaned data.csv')
# Now give initial guess of parameters:
library(lme4)
library(JMbayes2)
library(survival)
# Fit 2 longitudinal models for 2 biomarkers
data_longitudinal <- data |>
  group_by(id) |>
  mutate(gammaterm = as.numeric(visittime >= etime)*(visittime - (etime+preetime)/2)) |>
  ungroup()

model1_long = lme(fixed = sydigtot ~ visittime + gammaterm + age1 + educ_yrs + gender + base_cap,
                          random = ~ 1|id,
                   data = data_longitudinal)

model2_long = lme(fixed = stroopwo ~ visittime + gammaterm + age1 + educ_yrs + gender + base_cap,
                  random = ~ 1|id,
                   data = data_longitudinal)

summary(model1_long)
summary(model2_long)

# Get residual and random effect
resid1 <- residuals(model1_long, type = "response")
ranef1 <-ranef(model1_long)

resid2 <- residuals(model2_long, type = "response")
ranef2 <-ranef(model2_long)

#estimate of sigma.e, sigma.a
sigma.e_1_1 = as.numeric(model1_long$sigma)^2
sigma.e_2_2 = as.numeric(model2_long$sigma)^2
sigma.e_1_2 = cov(resid1, resid2)
sigma.a_1_2 = as.numeric(cov(ranef1, ranef2))
sigma.a_1_1 = as.numeric(VarCorr(model1_long)[1,"Variance"])
sigma.a_2_2 = as.numeric(VarCorr(model2_long)[1,"Variance"])


library(tibble)
ranef_df1 <- ranef(model1_long) |>
  rownames_to_column("id") |>
  mutate(id = as.integer(id))
ranef_df2 <- ranef(model2_long) |>
  rownames_to_column("id") |>
  mutate(id = as.integer(id))

data_survival <- data |>
  group_by(id) |>
  slice(1) |>
  mutate(est_eventtime = (preetime+etime)/2) |>
  ungroup() |>
  left_join(ranef_df1, by = "id") |>
  left_join(ranef_df2, by = "id") |>
  rename(ranef1 = `(Intercept).x`,
         ranef2 = `(Intercept).y`)

# Get theta.a
model_ranef = coxph(Surv(est_eventtime, status)~age1 + educ_yrs + gender + base_cap + ranef1 + ranef2, data=data_survival, x=TRUE)
base_surv <- survfit(model_ranef)
plot(base_surv, xlab = "Time", ylab = "Baseline Survival Probability",
     main = "Baseline Survival Curve", col = "blue", lwd = 2)

base_haz = basehaz(model_ranef)
hazard_df <- base_haz |>
  mutate(
    hazard_rate = c(NA, diff(hazard) / diff(time))
  )
plot(
  hazard_df$time, hazard_df$hazard,
  type = "l", col = "blue", lwd = 2,
  xlab = "Time", ylab = "Cumulative Hazard",
  main = "Cumulative Hazard Function"
)
model_hazard = lm(log(hazard)~log(time), data = hazard_df)
# coef: 1.49

model_cox = coxph(Surv(est_eventtime, status)~age1 + educ_yrs + gender + base_cap, data=data_survival, x=TRUE)
summary(model_cox)

theta.a_1 = as.numeric(model_ranef$coefficients[5])
theta.a_2 = as.numeric(model_ranef$coefficients[6])
joint_model <- jm(
  model_cox,
  list(model1_long, model2_long),
  time_var = "visittime",         # longitudinal time variable
  id_var = "id",        # subject ID
  control = list(n_iter = 1000)  # control settings (optional)
)

summary(joint_model)

# Set initial parameters
p=4
K=2
theta.x = as.matrix(coef(joint_model)$gammas)
rownames(theta.x) = paste("theta.x_",1:length(theta.x),sep="")
theta.m = as.matrix(coef(joint_model)$association)
rownames(theta.m) = paste("theta.m_",1:length(theta.m),sep="")
theta.a = as.matrix(c(theta.a_1, theta.a_2))
rownames(theta.a) = paste("theta.a_",1:length(theta.a),sep="")
beta.0 = as.matrix(c(summary(joint_model)$Outcome1$Mean[1], summary(joint_model)$Outcome2$Mean[1]))
rownames(beta.0) = paste("beta.0_",1:K,sep = "")
beta.1 = as.matrix(c(summary(joint_model)$Outcome1$Mean[2], summary(joint_model)$Outcome2$Mean[2]))
rownames(beta.1) = paste("beta.1_",1:K,sep = "")
beta.2 = matrix(c(summary(joint_model)$Outcome1$Mean[4:7], summary(joint_model)$Outcome2$Mean[4:7]),nrow=p,ncol=K)
beta.2.vec = matrix(as.vector(t(beta.2)))
index_grid <- expand.grid(row = 1:K, col = 1:p)
rownames(beta.2.vec) = paste("beta.2",index_grid$row,index_grid$col,sep = "_")
gamma = matrix(c(summary(joint_model)$Outcome1$Mean[3], summary(joint_model)$Outcome2$Mean[3]))
rownames(gamma) = paste("gamma_",1:K,sep = "")
Sigma.e.vec = matrix(c(sigma.e_1_1, sigma.e_1_2,sigma.e_2_2))
pos_labels_e <- sapply(1:K, function(i) {
  sapply(i:K, function(j) paste0("Sigma.e_", i, "_", j))
})
pos_labels_e <- unlist(pos_labels_e)
rownames(Sigma.e.vec) <- pos_labels_e

Sigma.a.vec = matrix(c(sigma.a_1_1, sigma.a_1_2,sigma.a_2_2))
pos_labels_a <- sapply(1:K, function(i) {
  sapply(i:K, function(j) paste0("Sigma.a_", i, "_", j))
})
pos_labels_a <- unlist(pos_labels_a)
rownames(Sigma.a.vec) <- pos_labels_a
#Use constant hazard/piecewise constant hazard
mode = 3
if (mode == 1){
  lambda0 = matrix(mean(hazard_df$hazard_rate, na.rm = TRUE))
  rownames(lambda0) = paste("lambda0")
  parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),t(lambda0),
                        t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))
} else if (mode == 2){
  parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),
                          t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))
} else if (mode == 3){
  parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),
                          t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))
}
num_para = length(parameters)
col_names <- colnames(parameters)
# Change data column name
long.data <- data |>
  rename(xval_1 = age1,
         xval_2 = educ_yrs,
         xval_3 = gender,
         xval_4 = base_cap,
         ynew_1 = sydigtot,
         ynew_2 = stroopwo)

# Define Gaussian nodes
no.gauss = 20
out.gauss=gauss.quad(no.gauss, kind="hermit")
nodes=as.matrix(do.call(expand.grid,replicate(K,out.gauss$nodes,simplify = FALSE)))
if (K>1){
  weight=as.vector(do.call(outer,replicate(K,out.gauss$weights,simplify = FALSE)))
} else {
  weight=as.vector(out.gauss$weights)
}
###############################################
setwd("C:/Users/yuezh/OneDrive - University of Nebraska Medical Center/Research/Joint_model/joint_likelihood")
source('Numeric Diff.R')
source('NR.R')
# Begin algorithm
# using constant likelihood
if (mode == 1) {
  ini.parameters = parameters
} else if (mode == 2){
  knots = 3
  long.data.id = long.data[!duplicated(long.data$id),]
  time = c(long.data.id$preetime,long.data.id$etime)
  d = c(0, quantile(time, probs = seq(1, knots)/knots, na.rm = TRUE))
  index <- findInterval(d, hazard_df$time) + 1
  matched_times <- ifelse(index <= length(hazard_df$time), hazard_df$time[index], NA)
  cumulative_hazard = c(hazard_df$hazard[hazard_df$time %in% matched_times], tail(hazard_df$hazard, 1))
  ini.hazard = data.frame(t(diff(cumulative_hazard)/diff(d)))
  colnames(ini.hazard) = paste("hazard",1:knots,sep = "_")
  ini.parameters = cbind(parameters, ini.hazard)
} else if (mode == 3){
  num_knots = 6
  degree = 3
  leg <- gauss.quad(no.gauss, kind = "legendre")
  z <- leg$nodes    # length = 20
  w <- leg$weights  # length = 20
  long.data.id = long.data[!duplicated(long.data$id),]
  # end_times <- pmin(long.data.id$etime, long.data.id$ctime)
  
  # Combine all end times with only the start times that are different
  time <- c(long.data.id$etime, long.data.id$preetime[long.data.id$preetime != long.data.id$etime])
  knots = quantile(time, probs = seq(1, num_knots)/(num_knots+1), na.rm = TRUE)
  
  # set boundary knots
  boundary_knots = c(0, max(max(long.data$visittime, na.rm = TRUE), knots))
  surv_data <- long.data.id %>%
    mutate(time_cox = (preetime + etime)/2)
  
  surv_object <- Surv(time = surv_data$time_cox, event = surv_data$status)
  hazard_fit <- estimate_bspline_hazard(surv_object, knots, degree)
  alpha = data.frame(t(hazard_fit[["coefficients"]][["estimate"]]))
  colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
  ini.parameters = cbind(parameters, alpha)
}

hat.parameters = ini.parameters

hat.parameters = NR_spline(long.data, hat.parameters, knots)

# end.time = Sys.time()
################################
coef = hat.parameters
# run.time = end.time - start.time
as.numeric(run.time)
filename = "C:/Users/yuezh/Research/Joint_model/data/coefficients_spline_knots=6.csv"
write.csv(coef, filename)
# plot B-spline cumulative hazard vs true 
alpha.hat = hat.parameters %>% dplyr::select(contains('alpha'))
time_grid = seq(0, max(time), 0.01)
cumulative_hazard = spline_cumulative_hazard(time_grid, knots, t(alpha.hat), boundary_knots, degree)

plot(
  hazard_df$time, hazard_df$hazard,
  type = "l", col = "blue", lwd = 2,
  xlab = "Time", ylab = "Cumulative Hazard",
  main = "Cumulative Hazard Function",
  ylim = c(0, max(cumulative_hazard))
)
lines(time_grid, cumulative_hazard)

# cumulative_hazard <- cumsum(hazards * c(0,diff(knots)))
time_expand <- sort(unique(c(d[-length(d)], time_grid)))
# hazard_expanded <- full_join(data.i, time_expanded)
hazard_df <- data.frame(time = time_expand)
hazard_df$hazard = NA
hazard_df$hazard[hazard_df$time %in% as.vector(d[-length(d)])] <- hazards
hazard_df <- hazard_df %>%
  fill(everything(), .direction = "down")
cumulative_hazard <- cumsum(c(0,diff(hazard_df$time))*hazard_df$hazard)
hazard_df <- cbind(hazard_df, cumulative_hazard)
hazard_plot <- hazard_df[hazard_df$time %in% time_grid,]$cumulative_hazard
plot(time_grid, hazard_plot)

# 8 minutes
# sd using hessian
Hessian.MLE = diff.likelihood.vec(long.data, hat.parameters)$Hessian
stddev = data.frame(t(sqrt(-diag(Hessian.MLE))))
stddev = stddev[,1:num_para]
colnames(stddev) = col_names


##############################################
# read existing coefficients
hat.parameters = read.csv("C:/Users/yuezh/Research/Joint_model/data/coefficients_spline.csv")
coef = hat.parameters %>% dplyr::select(!contains('alpha'))
# using Bootstrap
B = 50
bootstrap = 1
bootstrap_sample <- function(data) {
  # sampled_ids <- sample(unique(data$id), replace = TRUE)  # Sample IDs with replacement
  # resampled_data <- do.call(rbind, lapply(sampled_ids, function(id) data[data$id == id,]))
  
  unique_ids <- unique(data$id)
  boot_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
  
  # Step 3: Pull all rows for each sampled ID and assign a new sequential ID
  bootstrap_data <- purrr::map_dfr(seq_along(boot_ids), function(i) {
    data |>
      filter(id == boot_ids[i]) |>
      mutate(id = i)  # new ID from 1 to n
  })
  return(bootstrap_data)
}

bootstrap_fun = function(data, B, mode, coef){
  coef_bootstrap = matrix(NA, nrow = B, ncol = ncol(coef))
  
  for (b in 1:B){
    bootstrap_data = bootstrap_sample(data)
    if (mode == 1){
      hat.parameters = coef
    }
    else if (mode == 2){
      long.data.id = bootstrap_data[!duplicated(bootstrap_data$id),]
      time = c(long.data.id$preetime,long.data.id$etime)
      d = c(0, quantile(time, probs = seq(1, knots)/knots, na.rm = TRUE))
      
      bootstrap_data_survival <- bootstrap_data |>
        group_by(id) |>
        slice(1) |>
        mutate(est_eventtime = (preetime+etime)/2) |>
        ungroup()
      
      # Get theta.a
      model_cox_bootstrap = coxph(Surv(est_eventtime, status)~xval1 + xval2 + xval3 + xval4, data=bootstrap_data_survival, x=TRUE)
      
      base_haz = basehaz(model_cox_bootstrap)
      hazard_df <- base_haz |>
        mutate(
          hazard_rate = c(NA, diff(hazard) / diff(time))
        )
      index <- findInterval(d, hazard_df$time) + 1
      matched_times <- ifelse(index <= length(hazard_df$time), hazard_df$time[index], NA)
      cumulative_hazard = c(hazard_df$hazard[hazard_df$time %in% matched_times], tail(hazard_df$hazard, 1))
      ini.hazard = data.frame(t(diff(cumulative_hazard)/diff(d)))
      colnames(ini.hazard) = paste("hazard",1:knots,sep = "_")
      ini.parameters = cbind(coef, ini.hazard)
      hat.parameters = ini.parameters
    }
    #Calculate MLE for bootstrap data
    error = 0
    diff.theta = 1000
    tol = 1e-3
    loop = 1
    maxloop = 200
    try({
      while (diff.theta>tol & loop<=maxloop) {
        old.parameters <- hat.parameters
        error=0
        diff <- diff.likelihood.vec(bootstrap_data, hat.parameters)
        score <- diff$Score
        lr <- 1
        Hessian <- diff$Hessian
        if (any(is.na(Hessian))||any(is.infinite(Hessian))){
          error <- 1
        }
        step <- as.vector(Hessian%*%score)
        if (mode == 2){
          likelihood.hat <- likelihood.piecewise(bootstrap_data, hat.parameters)$ll
        }
        else if (mode == 1){
          likelihood.hat <- likelihood.vec(bootstrap_data, hat.parameters)$ll
        }
        while (1){
          est.parameters <- hat.parameters - lr*t(step)
          if (mode == 2){
          hazard.est <- as.numeric(est.parameters %>% dplyr::select(contains("hazard")))
          sigma.e.est <- est.parameters %>% dplyr::select(contains(paste("Sigma.e",1:K,1:K,sep="_")))
          sigma.a.est <- est.parameters %>% dplyr::select(contains(paste("Sigma.a",1:K,1:K,sep="_")))
          if (any(hazard.est<=0)|any(sigma.e.est<=0)|any(sigma.a.est<=0)){
            lr <- lr/2
            next
          }
          likelihood.est <- likelihood.piecewise(bootstrap_data, est.parameters)$ll
          }
          else if (mode == 1){
            lambda0.est = as.numeric(est.parameters %>% dplyr::select(contains("lambda0")))
            sigma.e.est = est.parameters %>% dplyr::select(contains(paste("Sigma.e",1:K,1:K,sep="_")))
            sigma.a.est = est.parameters %>% dplyr::select(contains(paste("Sigma.a",1:K,1:K,sep="_")))
            if (lambda0.est<=0|any(sigma.e.est<=0)|any(sigma.a.est<=0)){
              lr = lr/2
              next
            }
            likelihood.est = likelihood.vec(bootstrap_data, est.parameters)$ll
          }
          if (lr*max(abs(step))<1e-8){
            temp.parameters <- hat.parameters
            hat.parameters <- est.parameters
            break
          }
          if (likelihood.est<likelihood.hat){
            lr <- lr/2
            next
          } else {
            temp.parameters <- hat.parameters
            hat.parameters <- est.parameters
            break
          }
        }
        if (error == 1){
          hat.parameters <- matrix(NA, nrow=1, ncol=length(temp.parameters))
          break
        }
        diff.theta <- max(abs(hat.parameters-old.parameters))
        loop <- loop + 1
      }
      
      hat.parameters <- as.matrix(hat.parameters)
      coef_bootstrap[b, ] <- hat.parameters[, 1:length(parameters)]
    })
  }
  
  stddev = t(as.matrix(apply(coef_bootstrap, 2, sd,na.rm=TRUE)))
  return(stddev)
}

if (error == 0){
  # hat.hazard = hat.parameters %>% dplyr::select(contains("hazard"))
  # piecewise_hazard = cbind(as.vector(d[-length(d)]),t(hat.hazard))
  # colnames(piecewise_hazard) = c("knots","hazard")
  # coef = hat.parameters %>% dplyr::select(!contains("hazard"))
  # # bias = coef-parameters
  # colnames(coef) = col_names
  # # colnames(bias) = col_names
  col_names <- colnames(parameters)
  
  ####for Bootstrap, I would suggest add try function so if only a few bootstrap sample fail, you still can obtain the results
  if (bootstrap == 1){
    stddev = bootstrap_fun(long.data, B, mode, coef)
  }
  if (bootstrap == 0){
    hessian = diff.likelihood.vec(long.data, hat.parameters)$Hessian
    stddev = data.frame(t(as.matrix(sqrt(-diag(hessian)))))[,1:length(parameters)]
  }
  # Hessian.MLE = diff.likelihood.vec(long.data, hat.parameters)$Hessian
  # inv.obs.hessian = solve(obs.hessian)
  # sandwich.hessian = inv.obs.hessian%*%solve(-Hessian.MLE)%*%inv.obs.hessian
  # stddev = data.frame(t(sqrt(-diag(Hessian.MLE))))
  # stddev = stddev[,-((length(stddev)-knots+1):length(stddev))]
  colnames(stddev) = col_names
  # CP = as.matrix((parameters >= coef-1.96*stddev)&(parameters <= coef+1.96*stddev))
}

CIL = coef - 1.96*stddev
CIU = coef + 1.96*stddev
CI = rbind(CIL, CIU)
rownames(CI) = c("2.5% pct", "97.5% pct")
W=as.matrix(abs(coef)/stddev)
p_value = 2 * (1 - pnorm(W))


csv_dir <- "C:/Users/yuezh/Research/Joint_model/Result"

all_csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)

coef_files <- all_csv_files[grepl("coef", all_csv_files)]
hazard_files <- all_csv_files[grepl("hazard", all_csv_files)]
censor_files <- all_csv_files[grepl("censor", all_csv_files)]
causal_files <- all_csv_files[grepl("causal", all_csv_files)]

coef_bootstrap <- coef_files %>%
  lapply(read.csv) %>%    # Read each file
  bind_rows()             # Combine rows

causal <- causal_files %>%
  lapply(read.csv) %>%
  bind_rows()

sd_HD = apply(coef_bootstrap, 2, sd)
mean_HD = apply(coef_bootstrap, 2, mean)
outliers = rep(0, nrow(coef_bootstrap))
for (i in 1:nrow(coef_bootstrap)){
  outliers[i] = any((coef_bootstrap[i, ] > mean_HD + 3*sd_HD)+(coef_bootstrap[i, ] < mean_HD - 3*sd_HD))
}

outliers_index = which(outliers == 1)
coef_bootstrap = coef_bootstrap[-outliers_index, ]
sd_HD = apply(coef_bootstrap, 2, sd)


sd_causal = apply(causal, 2, sd)
coef = read.csv("C:/Users/yuezh/Research/Joint_model/data/coefficients_spline.csv")
coef = coef[,1:num_para]
CIL = coef - 1.96*sd_HD
CIU = coef + 1.96*sd_HD
CI = rbind(CIL, CIU)
rownames(CI) = c("2.5% pct", "97.5% pct")
W=as.matrix(abs(coef)/sd_HD)
p_value = 2 * (1 - pnorm(W))

hazard_files = hazard_files[-outliers_index]
plot_matrix <- matrix(0, nrow = length(hazard_files), ncol = length(time_grid))
for (i in 1:length(hazard_files)){
  # time_grid <- seq(0, J -0.8, 0.01)
  data.i = read.csv(hazard_files[i])
  knots = unlist(data.i %>% dplyr::select(contains('knots')))
  alpha = data.i %>% dplyr::select(contains('alpha'))
  boundary = as.numeric(data.i %>% dplyr::select(contains('boundary')))
  # hazards = data.i$hazard
  # cumulative_hazard <- cumsum(hazards * c(0,diff(knots)))
  # time_expand <- sort(unique(c(data.i$knots, time_grid)))
  # # hazard_expanded <- full_join(data.i, time_expanded)
  # hazard_df <- data.frame(time = time_expand)
  # hazard_df$hazard = NA
  # hazard_df$hazard[hazard_df$time %in% data.i$knots] <- hazards
  # hazard_df <- hazard_df %>%
  #   fill(everything(), .direction = "down")
  cumulative_hazard <- spline_cumulative_hazard(time_grid, knots, t(alpha), boundary_knots = c(0, boundary), degree)
  plot_matrix[i, ] <- t(cumulative_hazard)
}
hazard_2.5 <- apply(plot_matrix, 2, function(x) quantile(x, probs=0.025))
hazard_97.5 <- apply(plot_matrix, 2, function(x) quantile(x, probs=0.975))
plot_df <- data.frame(
  time = time_grid,
  estimate = cumulative_hazard,
  lower_ci = hazard_2.5,
  upper_ci = hazard_97.5
)
ggplot(plot_df, aes(x = time)) +
  # Add the 95% confidence interval as a shaded ribbon
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "skyblue", alpha = 0.5) +
  
  # Add the mean cumulative hazard line
  geom_line(aes(y = estimate, color = "Estimate"), linewidth = 1) +
  
  # Add the true cumulative hazard line
  # geom_line(aes(y = true, color = "True"), linetype = "dashed", linewidth = 1) +
  
  # Customize labels, colors, and title
  labs(
    title = "Estimated Cumulative Hazard and CI",
    x = "Time",
    y = "Cumulative Hazard",
    color = "Legend" # Title for the legend
  ) +
  scale_color_manual(
    name = NULL, # Hide the legend title if you want
    values = c("Estimate" = "red"),
    labels = c("Estimate")
  ) +
  
  # Apply a clean theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )


source('Causal Mediation.R')
t = 8
index = 4
x1 = quantile(long.data.id$xval_4, 0.25)
x2 = quantile(long.data.id$xval_4, 0.75)

NE = NE_M2(long.data, t, x1, x2, hat.parameters, index)
NE = data.frame(NE)
CIL = NE - 1.96*sd_causal
CIU = NE + 1.96*sd_causal
CI = rbind(CIL, CIU)
rownames(CI) = c("2.5% pct", "97.5% pct")
W=as.matrix(abs(NE)/sd_causal)
p_value = 2 * (1 - pnorm(W))

file.remove(all_csv_files)
