#Cox proportional Hazard using Lungs dataset

#with package
library(survival)
data(lung)
head(lung)
coxmodel =coxph(Surv(time,status)~age+sex+ph.ecog,data=lung)
cox_model

#without package
cox_partial_likelihood = function(beta,data){
  time = data$time
  status = data$status
  covariates = as.matrix(data[,c("age","sex","ph.ecog")])
  linear_pred =covariates%*%beta
  loglikelihood = 0
  for (i in 1:length(time)){
    if(status[i]==1){
      risk_set = which(time>=time[i])
      risk_set_linear_pred = linear_pred[risk_set]
      loglikelihood = loglikelihood + linear_pred[i]-log(sum(exp(risk_set_linear_pred)))
    }
  }
  return(-loglikelihood)
}

initial_beta = rep(0,ncol(lung[,c("age","sex","ph.ecog")]))
cox_fit = optim(initial_beta, cox_partial_likelihood,data=lung, method="BFGS")
cox_fit$par


#practical 13
# Load necessary libraries
set.seed(0)  # For reproducibility
library(survival)

# Simulation parameters
n <- 1000
lambda1 <- 0.02
lambda2 <- 0.03

# Simulate failure times
T1 <- rexp(n, rate = lambda1)
T2 <- rexp(n, rate = lambda2)

# Simulate censoring times
C <- runif(n, min = 0, max = 100)

# Calculate observed times and event indicators
T <- pmin(T1, T2, C)
delta <- ifelse(T == T1, 1, ifelse(T == T2, 2, 0))

# Create a data frame for the results
data <- data.frame(T = T, delta = delta)

# Create a Surv object for survival analysis
surv_obj <- Surv(data$T, data$delta)

# Fit the survival curves for each cause
fit1 <- survfit(surv_obj ~ (data$delta == 1))
fit2 <- survfit(surv_obj ~ (data$delta == 2))

# Estimate cause-specific hazards
hazard1 <- -diff(c(1, fit1$surv)) / diff(c(0, fit1$time))
hazard2 <- -diff(c(1, fit2$surv)) / diff(c(0, fit2$time))

# Calculate the average hazard rates
estimated_lambda1 <- mean(hazard1, na.rm = TRUE)
estimated_lambda2 <- mean(hazard2, na.rm = TRUE)

# Print the estimated hazards
cat("Estimated λ1:", estimated_lambda1, "True λ1:", lambda1, "\n")
cat("Estimated λ2:", estimated_lambda2, "True λ2:", lambda2, "\n")
