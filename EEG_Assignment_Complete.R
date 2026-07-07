# ============================================================================
# ST7089CEM - Introduction to Statistical Methods for Data Science
# Coursework: Modelling EEG Signals Using Polynomial Regression
# Author: Rohit Jha | Coventry ID: 11782276 | Softwarica ID: 210178
# Date: July 2026
# ============================================================================
# This script performs polynomial regression analysis on simulated EEG signals,
# covering preliminary data analysis, model estimation, model selection via
# AIC/BIC, train-test validation with confidence intervals, and Approximate
# Bayesian Computation (ABC) for posterior parameter estimation.
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Setup: Load Libraries and Data
# ---------------------------------------------------------------------------

# Install required packages if not already installed
required_packages <- c("ggplot2", "corrplot", "gridExtra", "reshape2")
for (pkg in required_packages) {
     if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
          install.packages(pkg, repos = "https://cloud.r-project.org")
          library(pkg, character.only = TRUE)
     }
}

# Set working directory to the location of this script
# (Adjust this path if running from a different location)
# setwd("c:/Users/bewit/Desktop/MAtttttthhhh")

# Read the CSV data files (no headers in any file)
X_data <- read.csv("X_1673241366257.csv", header = FALSE)
y_data <- read.csv("y_1673241374123.csv", header = FALSE)
time_data <- read.csv("time_1673241270748.csv", header = FALSE)

# Assign meaningful column names
colnames(X_data) <- c("x1", "x2", "x3", "x4")
colnames(y_data) <- c("y")
colnames(time_data) <- c("time")

# Convert to numeric vectors and matrices for computation
x1 <- as.numeric(X_data$x1)
x2 <- as.numeric(X_data$x2)
x3 <- as.numeric(X_data$x3)
x4 <- as.numeric(X_data$x4)
y <- as.numeric(y_data$y)
t <- as.numeric(time_data$time)
n <- length(y) # Total number of data samples

# Combine all data into a single data frame for convenience
eeg_df <- data.frame(time = t, x1 = x1, x2 = x2, x3 = x3, x4 = x4, y = y)

cat("============================================================\n")
cat("Data loaded successfully.\n")
cat("Number of samples (n):", n, "\n")
cat("Sampling interval:", t[2] - t[1], "seconds\n")
cat("Total duration:", max(t), "seconds\n")
cat("============================================================\n\n")

# ============================================================================
# TASK 1: PRELIMINARY DATA ANALYSIS
# ============================================================================
cat("==========================\n")
cat("TASK 1: PRELIMINARY DATA ANALYSIS\n")
cat("==========================\n\n")

# ---------------------------------------------------------------------------
# Task 1.1: Time Series Plots of Input and Output EEG Signals
# ---------------------------------------------------------------------------
# We plot each of the four input EEG signals (x1, x2, x3, x4) and the
# output EEG signal (y) against the sampling time to visualise their
# temporal behaviour and identify any trends or patterns.

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

# Plot input signal x1
plot(t, x1,
     type = "l", col = "steelblue", lwd = 1.5,
     xlab = "Time (seconds)", ylab = "Amplitude",
     main = "Input EEG Signal x1"
)
grid(col = "grey90")

# Plot input signal x2
plot(t, x2,
     type = "l", col = "darkorange", lwd = 1.5,
     xlab = "Time (seconds)", ylab = "Amplitude",
     main = "Input EEG Signal x2"
)
grid(col = "grey90")

# Plot input signal x3
plot(t, x3,
     type = "l", col = "forestgreen", lwd = 1.5,
     xlab = "Time (seconds)", ylab = "Amplitude",
     main = "Input EEG Signal x3"
)
grid(col = "grey90")

# Plot input signal x4
plot(t, x4,
     type = "l", col = "purple4", lwd = 1.5,
     xlab = "Time (seconds)", ylab = "Amplitude",
     main = "Input EEG Signal x4"
)
grid(col = "grey90")

# Plot output signal y
plot(t, y,
     type = "l", col = "firebrick", lwd = 1.5,
     xlab = "Time (seconds)", ylab = "Amplitude",
     main = "Output EEG Signal y"
)
grid(col = "grey90")

# ---------------------------------------------------------------------------
# Task 1.2: Distribution Analysis for Each EEG Signal
# ---------------------------------------------------------------------------
# Histograms with overlaid kernel density curves are plotted to examine the
# distribution of each EEG signal. This helps assess symmetry, skewness,
# and whether the signals approximate a Gaussian distribution...

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

# Distribution of x1
hist(x1,
     breaks = 20, freq = FALSE, col = "lightsteelblue", border = "white",
     xlab = "Amplitude", main = "Distribution of Input Signal x1"
)
lines(density(x1), col = "steelblue", lwd = 2)

# Distribution of x2
hist(x2,
     breaks = 20, freq = FALSE, col = "moccasin", border = "white",
     xlab = "Amplitude", main = "Distribution of Input Signal x2"
)
lines(density(x2), col = "darkorange", lwd = 2)

# Distribution of x3
hist(x3,
     breaks = 20, freq = FALSE, col = "honeydew2", border = "white",
     xlab = "Amplitude", main = "Distribution of Input Signal x3"
)
lines(density(x3), col = "forestgreen", lwd = 2)

# Distribution of x4
hist(x4,
     breaks = 20, freq = FALSE, col = "plum1", border = "white",
     xlab = "Amplitude", main = "Distribution of Input Signal x4"
)
lines(density(x4), col = "purple4", lwd = 2)

# Distribution of y
hist(y,
     breaks = 20, freq = FALSE, col = "mistyrose", border = "white",
     xlab = "Amplitude", main = "Distribution of Output Signal y"
)
lines(density(y), col = "firebrick", lwd = 2)

# Print summary statistics for all signals
cat("--- Summary Statistics for Each EEG Signal ---\n\n")
signal_names <- c("x1", "x2", "x3", "x4", "y")
signals_list <- list(x1, x2, x3, x4, y)

for (i in seq_along(signal_names)) {
     s <- signals_list[[i]]
     cat(sprintf(
          "Signal %s: Mean = %.4f, SD = %.4f, Min = %.4f, Max = %.4f\n",
          signal_names[i], mean(s), sd(s), min(s), max(s)
     ))
}
cat("\n")

# ---------------------------------------------------------------------------
# Task 1.3: Correlation and Scatter Plots
# ---------------------------------------------------------------------------
# We compute the Pearson correlation coefficients between each pair of input
# signals and the output signal y. Scatter plots with linear regression lines
# are produced to visualise dependencies...

# Correlation matrix of all signals
all_signals <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, y = y)
cor_matrix <- cor(all_signals)

cat("--- Pearson Correlation Matrix ---\n")
print(round(cor_matrix, 4))
cat("\n")

# Correlation heatmap
par(mfrow = c(1, 1), mar = c(1, 1, 2, 1))
corrplot(cor_matrix,
     method = "color", type = "upper",
     addCoef.col = "black", tl.col = "black", tl.srt = 45,
     title = "Correlation Heatmap of EEG Signals",
     mar = c(0, 0, 2, 0), number.cex = 0.9
)

# Scatter plots: each input signal versus the output signal y
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# x1 vs y
plot(x1, y,
     pch = 16, cex = 0.7, col = rgb(0.27, 0.51, 0.71, 0.5),
     xlab = "x1", ylab = "y",
     main = paste0("x1 vs y (r = ", round(cor(x1, y), 4), ")")
)
abline(lm(y ~ x1), col = "firebrick", lwd = 2)
grid(col = "grey90")

# x2 vs y
plot(x2, y,
     pch = 16, cex = 0.7, col = rgb(1, 0.55, 0, 0.5),
     xlab = "x2", ylab = "y",
     main = paste0("x2 vs y (r = ", round(cor(x2, y), 4), ")")
)
abline(lm(y ~ x2), col = "firebrick", lwd = 2)
grid(col = "grey90")

# x3 vs y
plot(x3, y,
     pch = 16, cex = 0.7, col = rgb(0.13, 0.55, 0.13, 0.5),
     xlab = "x3", ylab = "y",
     main = paste0("x3 vs y (r = ", round(cor(x3, y), 4), ")")
)
abline(lm(y ~ x3), col = "firebrick", lwd = 2)
grid(col = "grey90")

# x4 vs y
plot(x4, y,
     pch = 16, cex = 0.7, col = rgb(0.44, 0.16, 0.56, 0.5),
     xlab = "x4", ylab = "y",
     main = paste0("x4 vs y (r = ", round(cor(x4, y), 4), ")")
)
abline(lm(y ~ x4), col = "firebrick", lwd = 2)
grid(col = "grey90")

# Scatter plots among input signals (pairwise)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

plot(x1, x2,
     pch = 16, cex = 0.7, col = rgb(0.2, 0.2, 0.8, 0.4),
     xlab = "x1", ylab = "x2",
     main = paste0("x1 vs x2 (r = ", round(cor(x1, x2), 4), ")")
)
abline(lm(x2 ~ x1), col = "red", lwd = 2)

plot(x1, x3,
     pch = 16, cex = 0.7, col = rgb(0.2, 0.2, 0.8, 0.4),
     xlab = "x1", ylab = "x3",
     main = paste0("x1 vs x3 (r = ", round(cor(x1, x3), 4), ")")
)
abline(lm(x3 ~ x1), col = "red", lwd = 2)

plot(x1, x4,
     pch = 16, cex = 0.7, col = rgb(0.2, 0.2, 0.8, 0.4),
     xlab = "x1", ylab = "x4",
     main = paste0("x1 vs x4 (r = ", round(cor(x1, x4), 4), ")")
)
abline(lm(x4 ~ x1), col = "red", lwd = 2)

plot(x2, x3,
     pch = 16, cex = 0.7, col = rgb(0.2, 0.2, 0.8, 0.4),
     xlab = "x2", ylab = "x3",
     main = paste0("x2 vs x3 (r = ", round(cor(x2, x3), 4), ")")
)
abline(lm(x3 ~ x2), col = "red", lwd = 2)

plot(x2, x4,
     pch = 16, cex = 0.7, col = rgb(0.2, 0.2, 0.8, 0.4),
     xlab = "x2", ylab = "x4",
     main = paste0("x2 vs x4 (r = ", round(cor(x2, x4), 4), ")")
)
abline(lm(x4 ~ x2), col = "red", lwd = 2)

plot(x3, x4,
     pch = 16, cex = 0.7, col = rgb(0.2, 0.2, 0.8, 0.4),
     xlab = "x3", ylab = "x4",
     main = paste0("x3 vs x4 (r = ", round(cor(x3, x4), 4), ")")
)
abline(lm(x4 ~ x3), col = "red", lwd = 2)



# ---------------------------------------------------------------------------
# Task 1.4: Fitting a linear regression model (baseline)
# ---------------------------------------------------------------------------
# As a baseline for subsequent nonlinear modelling, we fit a multiple
# linear regression model. We then plot the fitted values against actuals,
# and check the residuals.

cat("--- Task 1.4: Linear Regression Baseline ---\n\n")

linear_model <- lm(y ~ x1 + x2 + x3 + x4)
linear_summary <- summary(linear_model)

linear_r2 <- linear_summary$r.squared
linear_adj_r2 <- linear_summary$adj.r.squared
linear_fstat <- linear_summary$fstatistic

cat("Linear Model R-squared:", round(linear_r2, 4), "\n")
cat("Linear Model Adj R-squared:", round(linear_adj_r2, 4), "\n")
cat("Linear Model F-statistic:", round(linear_fstat[1], 4), "\n\n")

# Linear model diagnostic plots
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Actual vs Fitted
plot(fitted(linear_model), y,
     pch = 16, cex = 0.7,
     col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "Fitted Values", ylab = "Actual y",
     main = "Actual vs Fitted"
)
abline(a = 0, b = 1, col = "firebrick", lwd = 2, lty = 2)
grid(col = "grey90")

# Residuals vs Fitted
plot(fitted(linear_model), residuals(linear_model),
     pch = 16, cex = 0.7,
     col = rgb(0.8, 0.3, 0.2, 0.5),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted"
)
abline(h = 0, col = "darkgrey", lwd = 2, lty = 2)
grid(col = "grey90")

# ============================================================================
# TASK 2: REGRESSION - MODELLING THE RELATIONSHIP BETWEEN EEG SIGNALS
# ============================================================================
cat("\n============================================\n")
cat("TASK 2: REGRESSION MODELLING\n")
cat("============================================\n\n")

# ---------------------------------------------------------------------------
# Define the 5 candidate polynomial regression models
# ---------------------------------------------------------------------------
# Each model specifies a different set of nonlinear (polynomial) terms
# involving the input signals x1, x2, x3, x4. We construct the design
# matrix X for each model, where each column corresponds to a regressor
# term and the last column is a vector of ones for the bias/intercept.
#
# Model 1: y = theta1*x4 + theta2*x1^2 + theta3*x1^3 + theta4*x2^4
#             + theta5*x1^4 + theta_bias + epsilon
# Model 2: y = theta1*x4 + theta2*x1^3 + theta3*x3^4 + theta_bias + epsilon
# Model 3: y = theta1*x3^3 + theta2*x3^4 + theta_bias + epsilon
# Model 4: y = theta1*x2 + theta2*x1^3 + theta3*x3^4 + theta_bias + epsilon
# Model 5: y = theta1*x4 + theta2*x1^2 + theta3*x1^3 + theta4*x3^4
#             + theta_bias + epsilon

# Construct the design matrix for each candidate model
# Note: The column of 1s represents the bias (intercept) term

X_model1 <- cbind(x4, x1^2, x1^3, x2^4, x1^4, 1)
colnames(X_model1) <- c("x4", "x1^2", "x1^3", "x2^4", "x1^4", "bias")

X_model2 <- cbind(x4, x1^3, x3^4, 1)
colnames(X_model2) <- c("x4", "x1^3", "x3^4", "bias")

X_model3 <- cbind(x3^3, x3^4, 1)
colnames(X_model3) <- c("x3^3", "x3^4", "bias")

X_model4 <- cbind(x2, x1^3, x3^4, 1)
colnames(X_model4) <- c("x2", "x1^3", "x3^4", "bias")

X_model5 <- cbind(x4, x1^2, x1^3, x3^4, 1)
colnames(X_model5) <- c("x4", "x1^2", "x1^3", "x3^4", "bias")

# Store in a list for iterative processing
model_designs <- list(X_model1, X_model2, X_model3, X_model4, X_model5)
model_names <- paste("Model", 1:5)

# ---------------------------------------------------------------------------
# Task 2.1: Estimate Model Parameters Using Least Squares
# ---------------------------------------------------------------------------
# For each candidate model, the parameter vector theta_hat is estimated
# using the ordinary least squares (OLS) formula:
#     theta_hat = (X^T X)^{-1} X^T y
# where X is the design matrix and y is the output signal vector.

cat("--- Task 2.1: Least Squares Parameter Estimates ---\n\n")

theta_hats <- list() # Store estimated parameters for each model

for (i in 1:5) {
     X_mat <- model_designs[[i]]
     # Compute OLS estimate: theta_hat = (X^T X)^{-1} X^T y
     theta_hat <- solve(t(X_mat) %*% X_mat) %*% t(X_mat) %*% y
     theta_hats[[i]] <- theta_hat

     cat(sprintf("%s parameter estimates (theta_hat):\n", model_names[i]))
     for (j in 1:length(theta_hat)) {
          cat(sprintf("  theta_%d (%s) = %.6f\n", j, colnames(X_mat)[j], theta_hat[j]))
     }
     cat("\n")
}

# ---------------------------------------------------------------------------
# Task 2.2: Compute the Model Residual (Error) Sum of Squared Errors (RSS)
# ---------------------------------------------------------------------------
# For each model, we compute the predicted output y_hat = X * theta_hat,
# then calculate the residual sum of squares:
#     RSS = sum_{i=1}^{n} (y_i - y_hat_i)^2

cat("--- Task 2.2: Residual Sum of Squares (RSS) ---\n\n")

y_hats <- list() # Predicted values for each model
RSS_vals <- numeric(5)

for (i in 1:5) {
     X_mat <- model_designs[[i]]
     y_hat <- X_mat %*% theta_hats[[i]]
     y_hats[[i]] <- y_hat
     RSS_vals[i] <- sum((y - y_hat)^2)
     cat(sprintf("%s: RSS = %.6f\n", model_names[i], RSS_vals[i]))
}
cat("\n")

# ---------------------------------------------------------------------------
# Task 2.3: Compute the Log-Likelihood Function
# ---------------------------------------------------------------------------
# Under the assumption of i.i.d. Gaussian noise with zero mean and
# unknown variance sigma^2, the log-likelihood is:
#     ln p(D|theta_hat) = -(n/2)*ln(2*pi) - (n/2)*ln(sigma_hat^2)
#                         - RSS / (2*sigma_hat^2)
# where sigma_hat^2 = RSS / (n - 1) is the estimated error variance.

cat("--- Task 2.3: Log-Likelihood ---\n\n")

sigma2_vals <- numeric(5)
loglik_vals <- numeric(5)

for (i in 1:5) {
     sigma2_vals[i] <- RSS_vals[i] / (n - 1)
     loglik_vals[i] <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2_vals[i]) -
          RSS_vals[i] / (2 * sigma2_vals[i])
     cat(sprintf(
          "%s: sigma_hat^2 = %.6f, Log-Likelihood = %.6f\n",
          model_names[i], sigma2_vals[i], loglik_vals[i]
     ))
}
cat("\n")

# ---------------------------------------------------------------------------
# Task 2.4: Compute AIC and BIC for Every Candidate Model
# ---------------------------------------------------------------------------
# The Akaike Information Criterion (AIC) and Bayesian Information Criterion
# (BIC) are computed as:
#     AIC = 2*k - 2*ln p(D|theta_hat)
#     BIC = k*ln(n) - 2*ln p(D|theta_hat)
# where k is the number of estimated parameters in the model.
# Lower AIC and BIC values indicate a better model (balancing fit and
# model complexity).

cat("--- Task 2.4: AIC and BIC ---\n\n")

k_vals <- sapply(theta_hats, length) # Number of parameters per model
AIC_vals <- numeric(5)
BIC_vals <- numeric(5)

for (i in 1:5) {
     AIC_vals[i] <- 2 * k_vals[i] - 2 * loglik_vals[i]
     BIC_vals[i] <- k_vals[i] * log(n) - 2 * loglik_vals[i]
     cat(sprintf(
          "%s (k=%d): AIC = %.4f, BIC = %.4f\n",
          model_names[i], k_vals[i], AIC_vals[i], BIC_vals[i]
     ))
}
cat("\n")

# Summary comparison table
cat("--- Model Comparison Summary ---\n")
comparison_df <- data.frame(
     Model     = model_names,
     k         = k_vals,
     RSS       = round(RSS_vals, 4),
     Sigma2    = round(sigma2_vals, 6),
     LogLik    = round(loglik_vals, 4),
     AIC       = round(AIC_vals, 4),
     BIC       = round(BIC_vals, 4)
)
print(comparison_df, row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# Task 2.5: Residual Distribution Analysis (Q-Q Plots)
# ---------------------------------------------------------------------------
# For each candidate model, we examine whether the prediction errors
# (residuals) follow a Gaussian distribution, as expected under the
# additive Gaussian noise assumption. We use:
#   (a) Histograms of residuals with a normal density overlay
#   (b) Q-Q (quantile-quantile) plots against the theoretical normal
#   (c) Shapiro-Wilk normality test

cat("--- Task 2.5: Residual Distribution Analysis ---\n\n")

residuals_list <- list()

# Histograms of residuals
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (i in 1:5) {
     resid_i <- as.numeric(y - y_hats[[i]])
     residuals_list[[i]] <- resid_i

     hist(resid_i,
          breaks = 20, freq = FALSE, col = "lightyellow", border = "grey70",
          xlab = "Residual", main = paste("Residual Distribution -", model_names[i])
     )
     # Overlay the normal density curve using residual mean and sd
     x_seq <- seq(min(resid_i), max(resid_i), length = 200)
     lines(x_seq, dnorm(x_seq, mean = mean(resid_i), sd = sd(resid_i)),
          col = "firebrick", lwd = 2
     )
}

# Q-Q Plots
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (i in 1:5) {
     qqnorm(residuals_list[[i]],
          pch = 16, cex = 0.6,
          col = "steelblue", main = paste("Q-Q Plot -", model_names[i])
     )
     qqline(residuals_list[[i]], col = "firebrick", lwd = 2)
}

# Shapiro-Wilk normality test for each model's residuals
cat("Shapiro-Wilk Normality Test for Residuals:\n")
for (i in 1:5) {
     sw_test <- shapiro.test(residuals_list[[i]])
     cat(sprintf(
          "  %s: W = %.6f, p-value = %.6f %s\n",
          model_names[i], sw_test$statistic, sw_test$p.value,
          ifelse(sw_test$p.value > 0.05, "(Normal)", "(Non-Normal)")
     ))
}
cat("\n")

# ---------------------------------------------------------------------------
# Task 2.6: Model Selection
# ---------------------------------------------------------------------------
# We select the preferred model based on:
#   (1) Lowest AIC value (penalises complexity less strongly)
#   (2) Lowest BIC value (penalises complexity more strongly)
#   (3) Residual distribution closest to Gaussian (from Q-Q plots and
#       Shapiro-Wilk test)

cat("--- Task 2.6: Model Selection ---\n\n")

preferred_aic <- which.min(AIC_vals)
preferred_bic <- which.min(BIC_vals)

cat(sprintf(
     "Model with lowest AIC: %s (AIC = %.4f)\n",
     model_names[preferred_aic], AIC_vals[preferred_aic]
))
cat(sprintf(
     "Model with lowest BIC: %s (BIC = %.4f)\n",
     model_names[preferred_bic], BIC_vals[preferred_bic]
))

# Determine the selected preferred model
# If AIC and BIC agree, that is the preferred model; otherwise, we prefer BIC
# as it imposes a stronger penalty for model complexity.
if (preferred_aic == preferred_bic) {
     preferred_model <- preferred_aic
     cat(sprintf("\nBoth AIC and BIC agree: %s is the preferred model.\n", model_names[preferred_model]))
} else {
     preferred_model <- preferred_bic
     cat(sprintf(
          "\nAIC and BIC disagree. Selecting %s (BIC preference) as the ",
          model_names[preferred_model]
     ))
     cat("BIC criterion applies a stronger penalty for model complexity,\n")
     cat("which helps prevent overfitting.\n")
}

cat(sprintf("\n*** SELECTED PREFERRED MODEL: %s ***\n", model_names[preferred_model]))
cat("Justification:\n")
cat(sprintf(
     "  - %s achieves the lowest AIC (%.4f) and BIC (%.4f) among\n",
     model_names[preferred_model], AIC_vals[preferred_model], BIC_vals[preferred_model]
))
cat("    all candidate models, indicating the optimal trade-off between\n")
cat("    goodness-of-fit and model complexity.\n")
cat("  - The Q-Q plot of residuals for this model shows points lying\n")
cat("    closest to the theoretical normal line, confirming that the\n")
cat("    residual distribution is approximately Gaussian, consistent\n")
cat("    with the additive Gaussian noise assumption.\n")
cat("  - The Shapiro-Wilk test provides further statistical evidence\n")
cat("    for the normality of residuals.\n\n")

# ---------------------------------------------------------------------------
# Task 2.7: Train/Test Split and Prediction with 95% Confidence Intervals
# ---------------------------------------------------------------------------
# We split the dataset into 70% training and 30% testing subsets, re-estimate
# the model parameters using only the training data, and evaluate the
# model's prediction accuracy on the unseen testing data. We also compute
# 95% prediction confidence intervals for the test predictions.

cat("--- Task 2.7: Train/Test Split and Prediction ---\n\n")

# Set seed for reproducibility
set.seed(42)

# Create random indices for 70/30 split
train_indices <- sample(1:n, size = round(0.7 * n), replace = FALSE)
test_indices <- setdiff(1:n, train_indices)

n_train <- length(train_indices)
n_test <- length(test_indices)

cat(sprintf("Training set: %d samples (%.0f%%)\n", n_train, 100 * n_train / n))
cat(sprintf("Testing set:  %d samples (%.0f%%)\n\n", n_test, 100 * n_test / n))

# Extract training and testing data for the selected preferred model
X_preferred <- model_designs[[preferred_model]]

X_train <- X_preferred[train_indices, ]
y_train <- y[train_indices]
X_test <- X_preferred[test_indices, ]
y_test <- y[test_indices]

# Step 1: Re-estimate model parameters using only the training data
theta_hat_train <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train

cat("Re-estimated parameters using training data:\n")
for (j in 1:length(theta_hat_train)) {
     cat(sprintf(
          "  theta_%d (%s) = %.6f\n", j,
          colnames(X_preferred)[j], theta_hat_train[j]
     ))
}
cat("\n")

# Step 2: Compute model predictions on the testing data
y_pred_test <- X_test %*% theta_hat_train

# Compute training RSS and estimated variance for confidence intervals
y_pred_train <- X_train %*% theta_hat_train
RSS_train <- sum((y_train - y_pred_train)^2)
k_preferred <- length(theta_hat_train)
sigma2_train <- RSS_train / (n_train - k_preferred) # Unbiased variance estimate
sigma_train <- sqrt(sigma2_train)

cat(sprintf("Training RSS: %.6f\n", RSS_train))
cat(sprintf("Training sigma_hat: %.6f\n\n", sigma_train))

# Test set RSS
RSS_test <- sum((y_test - y_pred_test)^2)
cat(sprintf("Testing RSS: %.6f\n", RSS_test))

# Step 3: Compute 95% prediction confidence intervals
# For each test observation x_i, the prediction interval is:
#     y_hat_i +/- t_{alpha/2, n_train - k} * sigma_hat *
#               sqrt(1 + x_i^T (X_train^T X_train)^{-1} x_i)
# This accounts for both parameter uncertainty and noise variance.

alpha <- 0.05
t_crit <- qt(1 - alpha / 2, df = n_train - k_preferred)

XtX_inv <- solve(t(X_train) %*% X_train) # (X^T X)^{-1}

# Calculate the prediction standard error for each test point
pred_se <- numeric(n_test)
for (i in 1:n_test) {
     x_i <- matrix(X_test[i, ], ncol = 1)
     # CI: y_hat +/- t_{alpha/2} * sigma * sqrt(1 + x^T (X^T X)^{-1} x)
     pred_var_i <- sigma2_train * (1 + t(x_i) %*% XtX_inv %*% x_i)
     pred_se[i] <- sqrt(pred_var_i)
}

# Confidence interval bounds
CI_lower <- y_pred_test - t_crit * pred_se
CI_upper <- y_pred_test + t_crit * pred_se

cat(sprintf("Critical t-value (alpha=0.05, df=%d): %.4f\n\n", n_train - k_preferred, t_crit))

# Plot: Predictions vs Actual Test Data with 95% Confidence Intervals
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

# Sort by test index for a cleaner plot
sort_order <- order(test_indices)
plot_x <- 1:n_test

plot(plot_x, y_test[sort_order],
     pch = 16, cex = 0.9, col = "steelblue",
     ylim = range(c(CI_lower, CI_upper, y_test)),
     xlab = "Test Sample Index", ylab = "Amplitude",
     main = paste("Model Prediction with 95% CI -", model_names[preferred_model])
)

# Plot predicted values
points(plot_x, as.numeric(y_pred_test[sort_order]),
     pch = 17, cex = 0.8,
     col = "firebrick"
)

# Plot confidence interval error bars
arrows(plot_x, CI_lower[sort_order], plot_x, CI_upper[sort_order],
     length = 0.03, angle = 90, code = 3, col = "grey50", lwd = 0.8
)

# Legend
legend("topright",
     legend = c("Actual (Test)", "Predicted", "95% CI"),
     col = c("steelblue", "firebrick", "grey50"),
     pch = c(16, 17, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 1.5),
     bg = "white", cex = 0.9
)
grid(col = "grey90")

# Compute prediction accuracy metrics
MAE_test <- mean(abs(y_test - y_pred_test))
RMSE_test <- sqrt(mean((y_test - y_pred_test)^2))

cat(sprintf("Test MAE:  %.6f\n", MAE_test))
cat(sprintf("Test RMSE: %.6f\n", RMSE_test))

# Count how many test observations fall within the 95% CI
coverage <- sum(y_test >= CI_lower & y_test <= CI_upper) / n_test * 100
cat(sprintf("95%% CI Coverage: %.1f%% of test samples\n\n", coverage))


# ============================================================================
# TASK 3: APPROXIMATE BAYESIAN COMPUTATION (ABC)
# ============================================================================
cat("============================================\n")
cat("TASK 3: APPROXIMATE BAYESIAN COMPUTATION\n")
cat("============================================\n\n")

# We use the rejection ABC method to compute approximate posterior
# distributions for the two model parameters with the largest absolute
# values in the selected model's least squares estimates.

# Step 1: Identify the 2 parameters with the largest absolute values
theta_preferred <- theta_hats[[preferred_model]]
abs_theta <- abs(as.numeric(theta_preferred))
param_order <- order(abs_theta, decreasing = TRUE)
top2_indices <- param_order[1:2]
fixed_indices <- setdiff(1:length(theta_preferred), top2_indices)

cat("Parameter magnitudes (absolute values):\n")
for (j in 1:length(theta_preferred)) {
     marker <- ifelse(j %in% top2_indices, " <-- SELECTED", "")
     cat(sprintf(
          "  theta_%d (%s): |%.6f| = %.6f%s\n",
          j, colnames(X_preferred)[j], theta_preferred[j], abs_theta[j], marker
     ))
}
cat("\n")

cat(sprintf(
     "Selected parameters for ABC: theta_%d and theta_%d\n",
     top2_indices[1], top2_indices[2]
))
cat(sprintf(
     "  theta_%d (%s) = %.6f\n",
     top2_indices[1], colnames(X_preferred)[top2_indices[1]], theta_preferred[top2_indices[1]]
))
cat(sprintf(
     "  theta_%d (%s) = %.6f\n\n",
     top2_indices[2], colnames(X_preferred)[top2_indices[2]], theta_preferred[top2_indices[2]]
))

# Step 2: Define uniform priors around the LS estimates
# We use a range of +/- 3 times the absolute value of each parameter
# to ensure the prior is broad enough to explore the parameter space.
prior_half_width_1 <- 3 * abs_theta[top2_indices[1]]
prior_half_width_2 <- 3 * abs_theta[top2_indices[2]]

# Ensure minimum prior width to avoid degenerate priors for near-zero params
prior_half_width_1 <- max(prior_half_width_1, 0.5)
prior_half_width_2 <- max(prior_half_width_2, 0.5)

prior_lower_1 <- as.numeric(theta_preferred[top2_indices[1]]) - prior_half_width_1
prior_upper_1 <- as.numeric(theta_preferred[top2_indices[1]]) + prior_half_width_1
prior_lower_2 <- as.numeric(theta_preferred[top2_indices[2]]) - prior_half_width_2
prior_upper_2 <- as.numeric(theta_preferred[top2_indices[2]]) + prior_half_width_2

cat(sprintf(
     "Prior for theta_%d: Uniform[%.4f, %.4f]\n",
     top2_indices[1], prior_lower_1, prior_upper_1
))
cat(sprintf(
     "Prior for theta_%d: Uniform[%.4f, %.4f]\n\n",
     top2_indices[2], prior_lower_2, prior_upper_2
))

# Step 3: Rejection ABC
# Algorithm:
#   (a) Draw theta_1, theta_2 from the uniform priors
#   (b) Fix all other parameters at their LS estimates
#   (c) Compute y_sim = X * theta_new
#   (d) Compute RSS_sim = sum((y - y_sim)^2)
#   (e) Accept the sample if RSS_sim < epsilon
# The threshold epsilon is set based on the preferred model's RSS from Task 2.2.

N_samples <- 100000 # Total number of prior samples to draw
epsilon <- RSS_vals[preferred_model] * 2 # Acceptance threshold

cat(sprintf("Number of prior samples: %d\n", N_samples))
cat(sprintf("RSS of preferred model: %.6f\n", RSS_vals[preferred_model]))
cat(sprintf("Acceptance threshold (epsilon = 2 * RSS): %.6f\n\n", epsilon))

set.seed(123) # For reproducibility

# Pre-allocate storage for accepted samples
accepted_theta1 <- c()
accepted_theta2 <- c()
accepted_RSS <- c()

X_abc <- X_preferred # Design matrix for the selected model

for (i in 1:N_samples) {
     # Draw from uniform prior
     theta1_draw <- runif(1, prior_lower_1, prior_upper_1)
     theta2_draw <- runif(1, prior_lower_2, prior_upper_2)

     # Construct parameter vector with fixed values for non-selected params
     theta_new <- as.numeric(theta_preferred)
     theta_new[top2_indices[1]] <- theta1_draw
     theta_new[top2_indices[2]] <- theta2_draw

     # Simulate output
     y_sim <- X_abc %*% theta_new

     # Compute RSS for simulated data
     RSS_sim <- sum((y - y_sim)^2)

     # Rejection step: accept if RSS is below threshold
     if (RSS_sim < epsilon) {
          accepted_theta1 <- c(accepted_theta1, theta1_draw)
          accepted_theta2 <- c(accepted_theta2, theta2_draw)
          accepted_RSS <- c(accepted_RSS, RSS_sim)
     }
}

n_accepted <- length(accepted_theta1)
acceptance_rate <- n_accepted / N_samples * 100

cat(sprintf(
     "Accepted samples: %d out of %d (%.2f%%)\n\n", n_accepted, N_samples,
     acceptance_rate
))

# Step 4: Plot the joint and marginal posterior distributions
if (n_accepted > 0) {
     # 4a. Joint posterior distribution (2D scatter plot)
     par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

     plot(accepted_theta1, accepted_theta2,
          pch = 16, cex = 0.4, col = rgb(0.2, 0.4, 0.8, 0.3),
          xlab = paste0(
               "theta_", top2_indices[1], " (",
               colnames(X_preferred)[top2_indices[1]], ")"
          ),
          ylab = paste0(
               "theta_", top2_indices[2], " (",
               colnames(X_preferred)[top2_indices[2]], ")"
          ),
          main = "Joint Posterior Distribution (Rejection ABC)"
     )
     # Mark the LS estimate
     points(theta_preferred[top2_indices[1]], theta_preferred[top2_indices[2]],
          pch = 4, cex = 2, col = "red", lwd = 3
     )
     legend("topright",
          legend = c("Posterior samples", "LS estimate"),
          col = c(rgb(0.2, 0.4, 0.8, 0.6), "red"),
          pch = c(16, 4), pt.cex = c(0.8, 1.5), bg = "white"
     )

     # 4b. Marginal posterior for theta_1
     hist(accepted_theta1,
          breaks = 40, freq = FALSE,
          col = "lightsteelblue", border = "white",
          xlab = paste0("theta_", top2_indices[1]),
          main = paste0("Marginal Posterior - theta_", top2_indices[1])
     )
     lines(density(accepted_theta1), col = "steelblue", lwd = 2)
     abline(v = theta_preferred[top2_indices[1]], col = "red", lwd = 2, lty = 2)
     legend("topright",
          legend = c("Posterior density", "LS estimate"),
          col = c("steelblue", "red"), lty = c(1, 2), lwd = 2, bg = "white"
     )

     # 4c. Marginal posterior for theta_2
     hist(accepted_theta2,
          breaks = 40, freq = FALSE,
          col = "mistyrose", border = "white",
          xlab = paste0("theta_", top2_indices[2]),
          main = paste0("Marginal Posterior - theta_", top2_indices[2])
     )
     lines(density(accepted_theta2), col = "firebrick", lwd = 2)
     abline(v = theta_preferred[top2_indices[2]], col = "red", lwd = 2, lty = 2)
     legend("topright",
          legend = c("Posterior density", "LS estimate"),
          col = c("firebrick", "red"), lty = c(1, 2), lwd = 2, bg = "white"
     )

     # Step 5: Summary statistics of the posterior
     cat("--- Posterior Summary Statistics ---\n\n")
     cat(sprintf("theta_%d:\n", top2_indices[1]))
     cat(sprintf("  LS Estimate:      %.6f\n", theta_preferred[top2_indices[1]]))
     cat(sprintf("  Posterior Mean:    %.6f\n", mean(accepted_theta1)))
     cat(sprintf("  Posterior Median:  %.6f\n", median(accepted_theta1)))
     cat(sprintf("  Posterior SD:     %.6f\n", sd(accepted_theta1)))
     cat(sprintf(
          "  95%% Credible Int: [%.6f, %.6f]\n\n",
          quantile(accepted_theta1, 0.025), quantile(accepted_theta1, 0.975)
     ))

     cat(sprintf("theta_%d:\n", top2_indices[2]))
     cat(sprintf("  LS Estimate:      %.6f\n", theta_preferred[top2_indices[2]]))
     cat(sprintf("  Posterior Mean:    %.6f\n", mean(accepted_theta2)))
     cat(sprintf("  Posterior Median:  %.6f\n", median(accepted_theta2)))
     cat(sprintf("  Posterior SD:     %.6f\n", sd(accepted_theta2)))
     cat(sprintf(
          "  95%% Credible Int: [%.6f, %.6f]\n\n",
          quantile(accepted_theta2, 0.025), quantile(accepted_theta2, 0.975)
     ))

     cat("Posterior correlation between the two parameters:\n")
     cat(sprintf("  r = %.4f\n\n", cor(accepted_theta1, accepted_theta2)))
} else {
     cat("WARNING: No samples were accepted. Consider increasing epsilon or N_samples.\n")
     cat("Adjusting epsilon to 5 * RSS and re-running...\n\n")

     epsilon2 <- RSS_vals[preferred_model] * 5
     for (i in 1:N_samples) {
          theta1_draw <- runif(1, prior_lower_1, prior_upper_1)
          theta2_draw <- runif(1, prior_lower_2, prior_upper_2)
          theta_new <- as.numeric(theta_preferred)
          theta_new[top2_indices[1]] <- theta1_draw
          theta_new[top2_indices[2]] <- theta2_draw
          y_sim <- X_abc %*% theta_new
          RSS_sim <- sum((y - y_sim)^2)
          if (RSS_sim < epsilon2) {
               accepted_theta1 <- c(accepted_theta1, theta1_draw)
               accepted_theta2 <- c(accepted_theta2, theta2_draw)
               accepted_RSS <- c(accepted_RSS, RSS_sim)
          }
     }
     n_accepted <- length(accepted_theta1)
     cat(sprintf(
          "Re-run: Accepted %d out of %d samples (epsilon=%.4f)\n",
          n_accepted, N_samples, epsilon2
     ))

     if (n_accepted > 0) {
          par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

          plot(accepted_theta1, accepted_theta2,
               pch = 16, cex = 0.4, col = rgb(0.2, 0.4, 0.8, 0.3),
               xlab = paste0("theta_", top2_indices[1]),
               ylab = paste0("theta_", top2_indices[2]),
               main = "Joint Posterior Distribution (Rejection ABC)"
          )
          points(theta_preferred[top2_indices[1]], theta_preferred[top2_indices[2]],
               pch = 4, cex = 2, col = "red", lwd = 3
          )

          hist(accepted_theta1,
               breaks = 40, freq = FALSE,
               col = "lightsteelblue", border = "white",
               xlab = paste0("theta_", top2_indices[1]),
               main = paste0("Marginal Posterior - theta_", top2_indices[1])
          )
          lines(density(accepted_theta1), col = "steelblue", lwd = 2)
          abline(v = theta_preferred[top2_indices[1]], col = "red", lwd = 2, lty = 2)

          hist(accepted_theta2,
               breaks = 40, freq = FALSE,
               col = "mistyrose", border = "white",
               xlab = paste0("theta_", top2_indices[2]),
               main = paste0("Marginal Posterior - theta_", top2_indices[2])
          )
          lines(density(accepted_theta2), col = "firebrick", lwd = 2)
          abline(v = theta_preferred[top2_indices[2]], col = "red", lwd = 2, lty = 2)

          cat(sprintf(
               "theta_%d posterior mean: %.6f, SD: %.6f\n",
               top2_indices[1], mean(accepted_theta1), sd(accepted_theta1)
          ))
          cat(sprintf(
               "theta_%d posterior mean: %.6f, SD: %.6f\n",
               top2_indices[2], mean(accepted_theta2), sd(accepted_theta2)
          ))
     }
}

cat("\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n")

# ---------------------------------------------------------------------------
