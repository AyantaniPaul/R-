Practical 1
For the “Nile” dataset available in R, obtain the results given below by writing a suitable R program

a)	Autocorrelation up to order three
b)	Partial autocorrelation
c)	Also plot these using suitable diagram

CODE:
  library(datasets) 
data=Nile 
plot(data) 
require(graphics)
data1<-ts(data,start=c(1871,1), end=c(1970,12), frequency=12)
decomp<-decompose(data1) # to decompose the data into different components 
plot(decomp)

# (a) Autocorrelation up to order three 
Autocorr<-acf(data1,lag.max = 3,plot=FALSE)
Autocorr

# (b) Partial Autocorrelation
Partial_Autocorr<-acf(data1,lag.max = 3,type=c("partial"),plot=FALSE)
Partial_Autocorr

# (c) Plot the autocorrelation and partial autocorrelation functions 
par(mfrow = c(2, 1)) # Set layout for two plots in one column 
acf(data1, lag.max = 10, main = "Autocorrelation Function") 
pacf(data1, main = "Partial Autocorrelation Function")


MSMS – 302 Statistical Machine Learning
Practical 1
1.	For the following dataset, obtain kernel density estimate and Naïve density estimator. Also plot both the estimators.

5.65746599	5.38283914	2.79892121 2.85423660 2.95252721 5.42626667
7.66239113	-0.18001073	0.65083500 2.40276530 -0.09929884 6.32619215
5.03650752	2.07470777	1.78019174 6.12891558 4.05352439 2.02686971

3.50834853 -2.76449768	
4.98428763 3.01292677 2.82448038 3.98110437
5.09371862 5.97961648	4.56968496 -0.48814532 5.08736697 2.41757609


CODE:
  # Dataset
data <- c(5.65746599, 5.38283914, 2.79892121, 2.85423660, 2.95252721, 5.42626667,
            7.66239113,-0.18001073, 0.65083500, 2.40276530,-0.09929884, 6.32619215,
            5.03650752, 2.07470777, 1.78019174, 6.12891558, 4.05352439, 2.02686971,
            3.50834853,-2.76449768, 4.98428763, 3.01292677, 2.82448038, 3.98110437,
            5.09371862, 5.97961648, 4.56968496,-0.48814532, 5.08736697, 2.41757609)

# Kernel Density Estimate 
kde <- density(data)
kde

# Naïve Density Estimator 
naive_density <- function(x, data, h) { n <- length(data)
return(sum((abs(x- data) <= h) / (2 * h)) / n)
}

# Naïve Density Estimator for a range of x
h <- 1.0 # Bandwidth for Naïve Density Estimation
x_range <- seq(min(data)- 1, max(data) + 1, length.out = 100) 
naive_estimates <- sapply(x_range, naive_density, data = data, h = h)

# Plotting both
plot(kde, main = "Kernel Density Estimate vs Naïve Density Estimator",
     col = "blue", lwd = 2, ylim = c(0, max(kde$y, naive_estimates)), xlab = "X", ylab = "Density") 
lines(x_range, naive_estimates, col = "red", lwd = 2)
legend("topright", legend = c("KDE", "Naïve Density"), col = c("blue", "red"), lty = 1, lwd = 2)


Practical 2
2.	Hosmer & Lerneshow (1989) give a dataset (“birthwt” available in R MASS library) on 189 births at a US hospital, with the main interest being in low birth weight. The main variable of interest is low birth weight, a binary response variable low. You can use variable “low” as binary response variable and remining variables as regressor variable. Divide the whole dataset into training and test dataset as solve perform following task
a)	Learn logistic classification model from training dataset and predict response using test dataset predictors.
b)	Obtain specificity, sensitivity, positive predictive value, negative predictive value of the model using test data set.

#Code:
data("birthwt", package = "MASS") 
birthwt <- birthwt
birthwt$low <- as.numeric(birthwt$low)

# Split data into training and test sets
train_indices <- sample(1:nrow(birthwt), size = 0.7 * nrow(birthwt)) # 70% training 
train_data <- birthwt[train_indices, ]
test_data <- birthwt[-train_indices, ]

# Logistic function 
logistic <- function(x) { 1 / (1 + exp(-x))
}

# Fit logistic regression using Newton-Raphson 
fit_logistic <- function(X, y, max_iter = 100, tol = 1e-6) { 
  n <- nrow(X)
p <- ncol(X)

# Initialize beta coefficients
beta <- matrix(0, nrow = p, ncol = 1) 
epsilon <- 1e-8

for (i in 1:max_iter) {
  eta <- X %*% beta # Linear predictor
  mu <- logistic(eta) # Predicted probabilities 
  mu <- pmin(pmax(mu, epsilon), 1- epsilon)
  
  # Compute gradient and Hessian
  W <- diag(as.vector(mu * (1- mu)), n, n)
  z <- eta + (y- mu) / (mu * (1- mu))
  
  Hessian <- t(X) %*% W %*% X 
  gradient <- t(X) %*% (y- mu)
  beta_new <- beta + solve(Hessian + diag(epsilon, p)) %*% gradient # Regularized update
  
  if (any(is.na(beta_new))) 
    stop("Divergence detected in Newton-Raphson method.") 
  else if (max(abs(beta_new- beta)) < tol) break
  
  beta <- beta_new
}
return(beta)
}

# Addding intercept to training data
X_train <- cbind(1, as.matrix(train_data[,-1])) # Intercept and regressors 
y_train <- train_data$low

# Fitting logistic regression model 
beta <- fit_logistic(X_train, y_train)


# Predict on test data
X_test <- cbind(1, as.matrix(test_data[,-1])) 
y_test <- test_data$low
log_odds <- X_test %*% beta 
probs <- logistic(log_odds)
predictions <- ifelse(probs > 0.5, 1, 0)

# Calculating confusion matrix
confusion_matrix <- table(Predicted = predictions, Actual = y_test) 
cat("Confusion Matrix:\n")
print(confusion_matrix)

# Extracting metrics from confusion 
true_positive <- confusion_matrix["1", "1"] 
true_negative <- confusion_matrix["0", "0"] 
false_positive <- confusion_matrix["1", "0"] 
false_negative <- confusion_matrix["0", "1"]

sensitivity <- true_positive / (true_positive + false_negative) 
specificity <- true_negative / (true_negative + false_positive)
positive_predictive_value <- true_positive / (true_positive + false_positive) 
negative_predictive_value <- true_negative / (true_negative + false_negative) 
# Printing evolution metrics
cat("\nModel Performance Metrics:\n")
cat("Sensitivity (True Positive Rate):", round(sensitivity, 4), "\n")
cat("Specificity (True Negative Rate):", round(specificity, 4), "\n")
cat("Positive Predictive Value (PPV):", round(positive_predictive_value, 4), "\n") 
cat("Negative Predictive Value (NPV):", round(negative_predictive_value, 4), "\n")




Practical 3
Generate 100 observations from a gamma distribution with shape parameter = 2 and scale parameter = 5. Obtain kernel density estimates using the Gaussian kernel and plot these estimates. [bandwidth may be 0.1, 0.2, 0.5]
Code:
  #Generating 100 observations from a Gamma distribution 
shape <- 2
scale <- 5
data <- rgamma(100, shape = shape, scale = scale) 
bandwidths <- c(0.1, 0.2, 0.5)

plot(NULL, xlim = range(data), ylim = c(0, 0.3), xlab = "Value", ylab = "Density",
     main = "Kernel Density Estimates with Gaussian Kernel")

# Colors for different bandwidths 
colors <- c("red", "blue", "green")

# Plotting kernel density estimates for each bandwidth 
for (i in 1:length(bandwidths)) {
bw <- bandwidths[i]
density_estimate <- density(data, kernel = "gaussian", bw = bw) 
lines(density_estimate, col = colors[i], lwd = 2)

}

# Adding a legend to the plot
legend("topright", legend = paste("Bandwidth =", bandwidths), col = colors, lwd = 2, cex = 1)






MSMS – 303 Multivariate Analysis
Practical 1
Find MLE of Σ , 𝜇 𝑎𝑛𝑑 𝜌 for the data given in table and also find the result given below.

Head Length,
First Son (𝑥1)	Head Breadth,
First Son (𝑥2)	Head Length,
Second Son (𝑥3)	Head Breadth,
Second Son (𝑥4)
191	155	179	145
195	149	201	152
181	148	185	149
183	153	188	149
176	144	171	142
208	157	192	152
189	150	190	149
197	159	189	152
188	152	197	159
192	150	187	151
179	158	186	148
183	147	174	147
174	150	185	152
190	159	195	157
188	151	187	158
163	137	161	130
195	155	183	158
186	153	173	148
181	145	182	146
175	140	165	137
192	154	185	152
174	143	178	147
176	139	176	143
197	167	200	158
190	163	187	150

A)	Find the estimates of parameters of conditional distribution of (𝑥3, 𝑥4) given (𝑥1, 𝑥2) i.e. find 𝑆21𝑆11–1 𝑎𝑛𝑑 𝑆22.1 = 𝑆22 − 𝑆21𝑆11–1 𝑆12
B)	Find the partial correlation 𝑟34.12
C)	Use Fisher’s Z to find a confidence interval for 𝜌34.12 with confidence 0.95
D)	Find the sample multiple correlation coefficients between 𝑥3 and (𝑥1, 𝑥2) and between 𝑥4 and (𝑥1, 𝑥2)
E)	Test the hypothesis that 𝑥3 is independent of (𝑥1, 𝑥2) and 𝑥4 is independent of (𝑥1, 𝑥2)

CODE:
  # Data Entry
  head_length_first_son <- c(191, 195, 181, 176, 208, 189, 188, 192, 179, 183, 190, 188, 163, 186, 181,
                             192)
head_breadth_first_son <- c(155, 149, 148, 144, 157, 150, 152, 152, 158, 147, 159, 151, 137, 153,
                            140, 154)
head_length_second_son <- c(179, 201, 185, 171, 192, 190, 197, 186, 187, 174, 195, 187, 161, 173,
                            182, 185)
head_breadth_second_son <- c(145, 152, 149, 142, 152, 149, 159, 151, 148, 147, 157, 158, 130, 148,
                             146, 152)

# Creating Matrix 
data <- cbind(
x1 = head_length_first_son, x2 = head_breadth_first_son,
x3 = head_length_second_son, x4 = head_breadth_second_son
)

# 1. Mean Vector (MLE of μ) n <- nrow(data)
mean_vector <- apply(data, 2, function(col) sum(col) / n) # Manual mean calculation 
cat("Mean Vector (MLE of μ):\n")
print(mean_vector)

# 2. Covariance Matrix (MLE of Σ) 
cov_matrix <- matrix(0, ncol = 4, nrow = 4) 
for (i in 1:4) {
for (j in 1:4) {
  cov_matrix[i, j] <- sum((data[, i] - mean_vector[i]) * (data[, j] - mean_vector[j])) / (n - 1)
}
}
cat("Covariance Matrix (MLE of Σ):\n") 
print(cov_matrix)

# 3. Correlation Matrix (ρ)
cor_matrix <- matrix(0, ncol = 4, nrow = 4) 
for (i in 1:4) {
  for (j in 1:4) {
    cor_matrix[i, j] <- cov_matrix[i, j] / sqrt(cov_matrix[i, i] * cov_matrix[j, j])
  }
}
cat("Correlation Matrix (ρ):\n") 
print(cor_matrix)

# 4. Conditional Distribution Parameters
S11 <- cov_matrix[1:2, 1:2]	# Sub-matrix for (x1, x2)

S12 <- cov_matrix[1:2, 3:4]	# Sub-matrix between (x1, x2) and (x3, x4) 
S21 <- t(S12)	# Transpose of S12
S22 <- cov_matrix[3:4, 3:4]	# Sub-matrix for (x3, x4)
# Solving for S11 Inverse using Manual Inversion (2x2 matrix) 
inv_S11 <- matrix(0, 2, 2)
det_S11 <- S11[1, 1] * S11[2, 2] - S11[1, 2] * S11[2, 1] 
inv_S11[1, 1] <- S11[2, 2] / det_S11
inv_S11[2, 2] <- S11[1, 1] / det_S11 
inv_S11[1, 2] <- -S11[1, 2] / det_S11 
inv_S11[2, 1] <- -S11[2, 1] / det_S11

# Conditional Covariance: S22.1 = S22 - S21 * inv(S11) * S12 
S22_1 <- S22 - S21 %*% inv_S11 %*% S12
cat("Conditional Covariance Matrix (S22.1):\n") 
print(S22_1)

# 5. Partial Correlation Between x3 and x4 Given x1, x2
inv_cov <- solve(cov_matrix) # Manually inverted covariance matrix using built-in solve() 
r34.12 <- -inv_cov[3, 4] / sqrt(inv_cov[3, 3] * inv_cov[4, 4])
cat("Partial Correlation r34.12:\n", r34.12, "\n")

# 6. Fisher's Z for Confidence Interval
z_value <- 0.5 * log((1 + r34.12) / (1 - r34.12)) # Fisher's Z-transform 
se <- 1 / sqrt(n - 4)
lower <- tanh(z_value - 1.96 * se) 
upper <- tanh(z_value + 1.96 * se)
cat("95% Confidence Interval for r34.12:\n")
cat("Lower Bound:", lower, "\nUpper Bound:", upper, "\n")

# 7. Sample Multiple Correlation Coefficient # Manual Calculation for x3 ~ (x1, x2)
X <- as.matrix(cbind(1, data[, 1:2])) # Adding intercept term 
Y_x3 <- data[, 3]
beta_x3 <- solve(t(X) %*% X) %*% t(X) %*% Y_x3 # Regression Coefficients 
Y_hat_x3 <- X %*% beta_x3
RSS_x3 <- sum((Y_x3 - Y_hat_x3)^2) 
TSS_x3 <- sum((Y_x3 - mean(Y_x3))^2) 
R_x3 <- sqrt(1 - RSS_x3 / TSS_x3)
cat("Multiple Correlation Coefficient (x3 ~ x1, x2):\n", R_x3, "\n")

#Calculation for x4 ~ (x1, x2) 
Y_x4 <- data[, 4]
beta_x4 <- solve(t(X) %*% X) %*% t(X) %*% Y_x4 
Y_hat_x4 <- X %*% beta_x4
RSS_x4 <- sum((Y_x4 - Y_hat_x4)^2) 
TSS_x4 <- sum((Y_x4 - mean(Y_x4))^2) 
R_x4 <- sqrt(1 - RSS_x4 / TSS_x4)
cat("Multiple Correlation Coefficient (x4 ~ x1, x2):\n", R_x4, "\n")



Biostatistics
1.	Imagine that the incidence of gun violence is compared in two cities, one with relaxed gun laws (A), the other with strict gun laws (B). In the city with relaxed gun laws, there were 50 shootings in a population of 100,000 and in the other city, 10 shootings in a population of 100,000.
a)	What is the relative risk of gun violence in the city with relaxed gun laws (A)?
b)	What is the relative risk of gun violence in the city with strict gun laws (B)?
c)	What questions need to be asked before concluding that there is an association between shootings and gun laws?
  
Code:
  # (a) Relative risk for city A (relaxed gun laws) 
  shootings_A <- 50
  population_A <- 100000
  risk_A <- shootings_A / population_A
  
  shootings_B <- 10
  population_B <- 100000
  risk_B <- shootings_B / population_B
  
  
  # Relative risk of gun violence in city A 
  relative_risk_A <- risk_A / risk_B
  cat("Relative Risk of gun violence in city A:", relative_risk_A, "\n")
  
  
  # (b) Relative risk for city B (strict gun laws) 
  relative_risk_B <- risk_B / risk_A
  cat("Relative Risk of gun violence in city B:", relative_risk_B, "\n")
  cat("This means the risk of gun violence in city B is 20% of the risk in city A, or 5 times lower.:\n")
  
  
  
  # (c) Additional questions to consider:
  cat("Questions to consider before concluding association:\n")
  cat("1. Are the populations comparable in demographics and economic factors?\n")
  cat("2. Could other factors (like law enforcement or economic conditions) influence gun violence rates?\n")
  cat("3. Are the data collection periods the same for both cities?\n") 
  cat("4. Is the gun law the only difference between the cities?\n")
  cat("5. How is gun violence measured (e.g., homicides, accidents, suicides)?\n")
  

2.	A study looking at breast cancer in women compared cases with non-cases, and found that 75/100 cases did not use calcium supplements compared with 25/100 of the non-cases.
a)	Develop a table to display the data.
b)	Calculate the odds of exposure in cases and non-cases.
c)	Calculate the odds ratio using the cross-product ratio
d)	How does the difference between the two prevalence of breast cancer (75% vs 25%) compare to the odds ratio?
 
  
   CODE:
cases_no_calcium <- 75
cases_use_calcium <- 25
non_cases_no_calcium <- 25
non_cases_use_calcium <- 75


odds_cases <- cases_no_calcium / cases_use_calcium
# (b) Calculate odds of exposure in cases and non-cases odds_cases <- cases_no_calcium / cases_use_calcium
odds_non_cases <- non_cases_no_calcium / non_cases_use_calcium
  
cat("Odds of exposure in cases:", odds_cases, "\n") 
cat("Odds of exposure in non-cases:", odds_non_cases, "\n")
# (c) Calculate the odds ratio using cross-product
odds_ratio <- (cases_no_calcium * non_cases_use_calcium) / (cases_use_calcium * non_cases_no_calcium)
cat("Odds Ratio:", odds_ratio, "\n")
# (d) Compare prevalence ratio and odds ratio
  
prevalence_ratio <- (cases_no_calcium / 100) / (non_cases_no_calcium / 100) 
cat("Prevalence Ratio:", prevalence_ratio, "\n")
cat("The odds ratio magnifies the difference compared to the prevalence ratio.\n")
  
  #Key Difference:
cat("The prevalence ratio (3) compares relative exposure percentages,
while the odds ratio (9) compares the likelihood of exposure between groups.
The odds ratio magnifies differences, especially for common events, making it a stronger measure of association.")
  


3.	Let us consider the relationship between smoking and lung cancer. Suppose exposure to cigarette smoke increases the incidence of lung cancer by 20% (i.e. the relative risk is 1.2). Lung cancer has a base line incidence of 3% per year (in the non- exposed group). Suppose as well that baseline incidence in obese individuals is 1/3 less (i.e. 1%/yr.), and the relative risk associated with the exposure is also 1.2. You follow up 1000 non-obese and 1000 obese subjects with the exposure, and an equivalent number without the exposure. The study lasts 25 years. Work with 25-year cumulative incidence and a denominator of1000.

a)	Create a table to show the data for obese and non-obese subjects.
b)	Calculate the odds ratio of disease in the exposed group in relation to those who are not exposed.
c)	Compare the odds ratio with the relative risk of 1.2.

#CODE:
  # Given data are: 
non_obese_non_exposed <- 0.03
non_obese_exposed <- non_obese_non_exposed * 1.2 
obese_non_exposed <- non_obese_non_exposed * (2 / 3) 
obese_exposed <- obese_non_exposed * 1.2

# Create a data frame to display incidence rates 
incidence_table <- data.frame(

Group = c("Non-Obese, Non-Exposed", "Non-Obese, Exposed", "Obese, Non-Exposed", "Obese, Exposed"),
Incidence_Rate = c(non_obese_non_exposed, non_obese_exposed, obese_non_exposed, obese_exposed)
)
print(incidence_table)

# Calculating Odds and Odds Ratio for Non-Obese group as an example 
odds_non_exposed <- non_obese_non_exposed / (1- non_obese_non_exposed) 
odds_exposed <- non_obese_exposed / (1- non_obese_exposed)
odds_ratio <- odds_exposed / odds_non_exposed

# Displaying Odds Ratio and Relative Risk 
cat("Odds Ratio (Non-Obese):", odds_ratio, "\n") 
cat("Relative Risk (given): 1.2\n")



4.	Use the following table to calculate the attributable risk associated with taking a supplement containing folate during pregnancy:
  
  
  #CODE:
no_folate_neural <- 631		 # Neural tube defects (No folate) 
folate_neural <- 24		# Neural tube defects (Folate) 
no_folate_premature <- 727 # Premature births (No folate) 
folate_premature <- 563	# Premature births (Folate)

# Attributable Risk (AR) calculations 
AR_neural <- no_folate_neural- folate_neural
AR_premature <- no_folate_premature- folate_premature

cat("Attributable Risk for Neural Tube Defects:", AR_neural, "per 100,000\n") 
cat("Attributable Risk for Premature Births:", AR_premature, "per 100,000\n")





LIFETIME DATA ANALYSIS
1.	Generate 100 observations from Weibull distribution with shape parameter 3 and scale parameter
10. Hence obtain the ML estimation of its parameters. Also draw the two-dimensional likelihood plot of Weibull model for the given dataset. Finally obtain the ML estimate of Mean failure time and compare it with sample mean.

#CODE:
  
shape_param <- 3	# Shape 
scale_param <- 10 # Scale
n_samples <- 100	# Number of observations

# Weibull-distributed data
data <- scale_param * (-log(runif(n_samples)))^(1 / shape_param)

# Negative log-likelihood function 
neg_log_likelihood <- function(alpha, beta) {
if (alpha <= 0 || beta <= 0) return(Inf) # Ensure parameters are positive n <- length(data)
log_likelihood <-
  n * log(alpha) - n * log(beta) + (alpha - 1) * sum(log(data)) - sum((data / beta)^alpha)
return(-log_likelihood) # Negative log-likelihood
}

# Generate 2D grid for shape (alpha) and scale (beta) 
alpha_vals <- seq(2, 4, length.out = 50) # Shape parameter grid 
beta_vals <- seq(8, 12, length.out = 50) # Scale parameter grid
likelihood_matrix <- matrix(0, nrow = length(alpha_vals), ncol = length(beta_vals))

# Computing log-likelihood for each combination 
for (i in 1:length(alpha_vals)) {
for (j in 1:length(beta_vals)) {
  
  likelihood_matrix[i, j] <- -neg_log_likelihood(alpha_vals[i], beta_vals[j])
}}
# 2D likelihood plot 
filled.contour(
x = alpha_vals, y = beta_vals, z = likelihood_matrix, xlab = "Shape (alpha)", ylab = "Scale (beta)",
main = "2D Likelihood Plot of Weibull Model"
)
# MLE Estimation
result <- optim(par = c(2, 5), fn = function(params) neg_log_likelihood(params[1], params[2]), method = "L-BFGS-B", lower = c(1e-5, 1e-5))
mle_shape <- result$par[1] 
mle_scale <- result$par[2]

# Mean failure time
mean_failure_time <- mle_scale * gamma(1 + 1 / mle_shape) # Using Weibull formula 
sample_mean <- mean(data) # Sample mean
#results
cat("MLE Shape (alpha):", mle_shape, "\n") 
cat("MLE Scale (beta):", mle_scale, "\n")
cat("MLE Mean Failure Time:", mean_failure_time, "\n") 
cat("Sample Mean Failure Time:", sample_mean, "\n")






2.	fifty leukaemia patients were subjected to a test and the test is terminated when 35 patients were failed. Their lifetimes (in weeks) are given below:
  22.3 26.8 30.3 31.9 32.1 33.3 33.7 33.9 34.7 36.1 36.4 36.5 36.6, 37.1 37.6 38.2 38.5 38.7 38.7 38.9
38.9 39.1 41.1 41.1 41.4 42.4 43.6 43.8 44.0 45.3 45.8 50.4 51.3 51.4 51.5

Assume lifetimes follow lognormal distribution and estimate the two parameters of the distribution. Also estimate mean time to failure and median time to failure. Draw survival and hazard curve.

CODE:
  # Given data
  lifetimes <- c(22.3, 26.8, 30.3, 31.9, 32.1, 33.3, 33.7, 33.9, 34.7, 36.1, 36.4, 36.5,
                 36.6, 37.1, 37.6, 38.2, 38.5, 38.7, 38.7, 38.9, 38.9, 39.1, 41.1, 41.1,
                 41.4, 42.4, 43.6, 43.8, 44.0, 45.3, 45.8, 50.4, 51.3, 51.4, 51.5)
# Log-transform the data 
log_lifetimes <- log(lifetimes)

# Estimate parameters (mu and sigma) using MLE 
mu_hat <- mean(log_lifetimes)
sigma_hat <- sd(log_lifetimes)

# Calculate Mean and Median Time to Failure 
mean_time_to_failure <- exp(mu_hat + (sigma_hat^2 / 2)) 
median_time_to_failure <- exp(mu_hat)

# Survival Function 
survival_function <- function(t) {
1- pnorm((log(t)- mu_hat) / sigma_hat)
}
# Hazard Function 
hazard_function <- function(t) {
dnorm((log(t)- mu_hat) / sigma_hat) / (t * survival_function(t))
}
# Generate a sequence of times for plotting
time_seq <- seq(min(lifetimes), max(lifetimes), length.out = 100) # Survival and Hazard values
survival_vals <- survival_function(time_seq) 
hazard_vals <- hazard_function(time_seq)

# Plot Survival Curve
plot(time_seq, survival_vals, type = "l", col = "blue", lwd = 2, xlab = "Time (weeks)", ylab = "Survival Probability",
     main = "Survival Curve") # Plot Hazard Curve
plot(time_seq, hazard_vals, type = "l", col = "red", lwd = 2, xlab = "Time (weeks)", ylab = "Hazard Rate",
     main = "Hazard Curve") # results
cat("Estimated Parameters:\n")
cat("Mu (mean of log-lifetimes):", mu_hat, "\n") 
cat("Sigma (std. dev. of log-lifetimes):", sigma_hat, "\n\n")
cat("Mean Time to Failure (MTTF):", mean_time_to_failure, "weeks\n") 
cat("Median Time to Failure:", median_time_to_failure, "weeks\n")





3.	The recorded death times of 15 patients were 7.35, 8.69, 8.80, 9.63, 9.63, 9.89, 9.98, 10.24, 10.36,
10.37, 10.48, 11.33, 11.39, 12.02 and 13.12 days, 10 patients whose are alive were removed from the test at 20 days. Suppose recorded time follows Weibull distribution, then
a)	Find maximum likelihood estimates of parameter.
b)	Using estimates of part 1 draw survival and hazard rate curve.
c)	Comment on behaviour of hazard rate.

CODE:
  death_times <- c(7.35, 8.69, 8.80, 9.63, 9.63, 9.89, 9.98, 10.24, 10.36, 10.37,
                   10.48, 11.33, 11.39, 12.02, 13.12) # Complete data
censored_times <- rep(20, 10) # Right-censored data # Combine data
all_times <- c(death_times, censored_times)
status <- c(rep(1, length(death_times)), rep(0, length(censored_times))) # 1=death, 0=censored

# Negative log-likelihood function 
neg_log_likelihood <- function(params) { 
alpha <- params[1] # Shape parameter 
beta <- params[2] # Scale parameter
if (alpha <= 0 || beta <= 0) return(Inf)
log_f <- log(alpha / beta) + (alpha- 1) * log(death_times / beta)- (death_times / beta)^alpha 
log_S <--(censored_times / beta)^alpha
total_log_likelihood <- sum(log_f) + sum(log_S)
return(-total_log_likelihood) # Negative log-likelihood for minimization
}
# Initial guesses for alpha and beta 
initial_guess <- c(1, 15)

# MLE using optim
result <- optim(par = initial_guess, fn = neg_log_likelihood, method = "L-BFGS-B", lower = c(1e-5, 1e-5))
mle_alpha <- result$par[1] 
mle_beta <- result$par[2]

# Survival function 
survival_function <- function(t) { exp(-(t / mle_beta)^mle_alpha)
}
# Hazard function 
hazard_function <- function(t) {
(mle_alpha / mle_beta) * (t / mle_beta)^(mle_alpha- 1)
}
# Generate time sequence for plotting 
time_seq <- seq(0, 25, length.out = 100) # Survival and hazard values
survival_vals <- survival_function(time_seq) 
hazard_vals <- hazard_function(time_seq) # Plotting Survival Curve
plot(time_seq, survival_vals, type = "l", col = "blue", lwd = 2,
     xlab = "Time (days)", ylab = "Survival Probability", main = "Survival Curve")

# Plotting Hazard Rate Curve
plot(time_seq, hazard_vals, type = "l", col = "red", lwd = 2,
     xlab = "Time (days)", ylab = "Hazard Rate", main = "Hazard Curve")

# Results
cat("MLE Estimates:\n")
cat("Shape (alpha):", mle_alpha, "\n")
cat("Scale (beta):", mle_beta, "\n")





