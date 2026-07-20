# ============================================================================
# ST7089CEM - Introduction to Statistical Methods for Data Science
# Coursework: Modelling EEG Signals Using Polynomial Regression
# Author: Rohit Jha | Coventry ID: 11782276 | Softwarica ID: 210178
# Module Leader: Shrawan Thakur
# Date: July 2026
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Setup: Load Libraries and Data
# ---------------------------------------------------------------------------

# Clear any broken or hanging graphics devices
# graphics.off()

# Install required packages if not already installed
required_packages <- c("ggplot2", "corrplot", "gridExtra", "reshape2")
for (pkg in required_packages) {
     if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
          install.packages(pkg, repos = "https://cloud.r-project.org")
          library(pkg, character.only = TRUE)
     }
}

# Set working directory to the location of this script
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
# Sample data of X inputs and time
# ---------------------------------------------------------------------------
cat("--- Sample data of X inputs and time ---\n\n")

# Display sample data (first 6 rows)
sample_df <- head(eeg_df, 6)
print(sample_df)
cat("\n")

# ---------------------------------------------------------------------------
# Task 1.1: Time Series Plots of Input and Output EEG Signals
# ---------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(reshape2)

# Define a custom attractive theme for all plots
custom_theme <- theme_minimal(base_size = 12) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold", color = "#2C3E50", size = 14),
          axis.title = element_text(face = "bold", color = "#34495E"),
          axis.text = element_text(color = "#7F8C8D"),
          panel.grid.major = element_line(color = "#D5DBDB", linetype = "dotted", linewidth = 0.8),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(15, 15, 15, 15),
          legend.text = element_text(size = 12, color = "#34495E")
     )

# Function to generate individual attractive plots
create_plot <- function(data, y_var, color, title) {
     ggplot(data, aes(x = time, y = .data[[y_var]])) +
          geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
          geom_line(color = color, linewidth = 0.8) +
          labs(title = title, x = "Time (seconds)", y = "Amplitude") +
          custom_theme
}

# 1. Overlaid Plot (Position 1)
df_melt <- melt(eeg_df[, c("time", "x1", "x2", "x3", "x4")], id.vars = "time", variable.name = "Signal", value.name = "Value")
p1 <- ggplot(df_melt, aes(x = time, y = Value, color = Signal)) +
     geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
     geom_line(linewidth = 0.8) +
     labs(title = "Time series comparison of X with time", x = "Time (seconds)", y = "Amplitude") +
     scale_color_manual(values = c("x1" = "steelblue", "x2" = "darkorange", "x3" = "forestgreen", "x4" = "purple4"), labels = c("X1", "X2", "X3", "X4")) + # nolint # nolint line_length_linter.
     custom_theme +
     theme(legend.position = "top", legend.title = element_blank())

# 2-6. Individual Plots
p2 <- create_plot(eeg_df, "x1", "steelblue", "Input EEG Signal x1")
p3 <- create_plot(eeg_df, "x2", "darkorange", "Input EEG Signal x2")
p4 <- create_plot(eeg_df, "x3", "forestgreen", "Input EEG Signal x3")
p5 <- create_plot(eeg_df, "x4", "purple4", "Input EEG Signal x4")
eeg_df$y_col <- y
p6 <- create_plot(eeg_df, "y_col", "firebrick", "Output EEG Signal y")

# Save the individual X inputs (x1 to x4) to a 2x2 grid
grid.arrange(p2, p3, p4, p5, ncol = 2)

# Display them sequentially in the R viewer
print(p1)
grid.arrange(p2, p3, p4, p5, ncol = 2)
print(p6)

# ---------------------------------------------------------------------------
# Task 1.2: Distribution Analysis for Each EEG Signal
# ---------------------------------------------------------------------------

# Create overlaid distribution plot for all X inputs
x_df <- data.frame(X1 = x1, X2 = x2, X3 = x3, X4 = x4)
x_melt <- melt(x_df, variable.name = "Inputs", value.name = "Signal")

p_dist_overlaid <- ggplot(x_melt, aes(x = Signal, fill = Inputs, color = Inputs)) +
     geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.4, bins = 20) +
     geom_density(alpha = 0.1, size = 1) +
     geom_rug(sides = "b", color = "black", alpha = 0.2) +
     scale_fill_brewer(palette = "Pastel1") +
     scale_color_brewer(palette = "Pastel1") +
     labs(title = "Distribution of X inputs", x = "X Signal", y = "Density") +
     theme_minimal(base_size = 14) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save the overlaid distribution plot
print(p_dist_overlaid)

# Individual Distributions using ggplot2
p_dist_x1 <- ggplot(data.frame(Signal = x1), aes(x = Signal)) +
     geom_histogram(aes(y = ..density..), fill = "lightpink", color = "white", bins = 20, alpha = 0.7) +
     geom_density(color = "indianred", size = 1) +
     geom_rug(sides = "b", color = "black", alpha = 0.2) +
     labs(title = "Distribution of X1", x = "Amplitude", y = "Density") +
     theme_minimal(base_size = 12) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_dist_x2 <- ggplot(data.frame(Signal = x2), aes(x = Signal)) +
     geom_histogram(aes(y = ..density..), fill = "lightblue", color = "white", bins = 20, alpha = 0.7) +
     geom_density(color = "steelblue", size = 1) +
     geom_rug(sides = "b", color = "black", alpha = 0.2) +
     labs(title = "Distribution of X2", x = "Amplitude", y = "Density") +
     theme_minimal(base_size = 12) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_dist_x3 <- ggplot(data.frame(Signal = x3), aes(x = Signal)) +
     geom_histogram(aes(y = ..density..), fill = "lightgreen", color = "white", bins = 20, alpha = 0.7) +
     geom_density(color = "forestgreen", size = 1) +
     geom_rug(sides = "b", color = "black", alpha = 0.2) +
     labs(title = "Distribution of X3", x = "Amplitude", y = "Density") +
     theme_minimal(base_size = 12) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_dist_x4 <- ggplot(data.frame(Signal = x4), aes(x = Signal)) +
     geom_histogram(aes(y = ..density..), fill = "thistle", color = "white", bins = 20, alpha = 0.7) +
     geom_density(color = "purple4", size = 1) +
     geom_rug(sides = "b", color = "black", alpha = 0.2) +
     labs(title = "Distribution of X4", x = "Amplitude", y = "Density") +
     theme_minimal(base_size = 12) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save the 2x2 grid of individual X distributions
grid.arrange(p_dist_x1, p_dist_x2, p_dist_x3, p_dist_x4, ncol = 2)

# Distribution of y

y_df <- data.frame(Signal = y)
p_dist_y <- ggplot(y_df, aes(x = Signal)) +
     geom_histogram(aes(y = ..density..), fill = "lightcoral", color = "white", bins = 15, alpha = 0.4) +
     geom_density(color = "lightcoral", size = 1) +
     geom_rug(sides = "b", color = "black", alpha = 0.2) +
     labs(title = "Distribution of Y", x = "Y Signal", y = "Density") +
     theme_minimal(base_size = 14) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_dist_y)

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

df_scatter <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, y = y)

create_scatter <- function(data, x_var, color_val, title) {
     r_val <- round(cor(data[[x_var]], data$y), 4)
     ggplot(data, aes(x = .data[[x_var]], y = y)) +
          geom_point(alpha = 0.4, color = color_val, size = 2) +
          geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 0.8, se = FALSE) +
          geom_smooth(method = "loess", color = color_val, linewidth = 1.2, se = FALSE) +
          labs(title = paste0(title, " vs Y (r = ", r_val, ")"), x = paste("Input", title), y = "Output Y") +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_s1 <- create_scatter(df_scatter, "x1", "#e41a1c", "X1")
p_s2 <- create_scatter(df_scatter, "x2", "#377eb8", "X2")
p_s3 <- create_scatter(df_scatter, "x3", "#4daf4a", "X3")
p_s4 <- create_scatter(df_scatter, "x4", "#984ea3", "X4")

grid.arrange(p_s1, p_s2, p_s3, p_s4, ncol = 2)

# Scatter plots among input signals (pairwise)
create_scatter_pair <- function(data, x_var, y_var, color_val) {
     r_val <- round(cor(data[[x_var]], data[[y_var]]), 4)
     ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
          geom_point(alpha = 0.3, color = color_val, size = 1.5) +
          geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 0.8, se = FALSE) +
          geom_smooth(method = "loess", color = color_val, linewidth = 1.2, se = FALSE) +
          labs(
               title = paste0(toupper(x_var), " vs ", toupper(y_var), " (r = ", r_val, ")"),
               x = toupper(x_var), y = toupper(y_var)
          ) +
          theme_minimal(base_size = 11) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_p1 <- create_scatter_pair(df_scatter, "x1", "x2", "steelblue")
p_p2 <- create_scatter_pair(df_scatter, "x1", "x3", "steelblue")
p_p3 <- create_scatter_pair(df_scatter, "x1", "x4", "steelblue")
p_p4 <- create_scatter_pair(df_scatter, "x2", "x3", "steelblue")
p_p5 <- create_scatter_pair(df_scatter, "x2", "x4", "steelblue")
p_p6 <- create_scatter_pair(df_scatter, "x3", "x4", "steelblue")

grid.arrange(p_p1, p_p2, p_p3, p_p4, p_p5, p_p6, ncol = 2)

# ---------------------------------------------------------------------------
# Task 1.4: Fitting a linear regression model (baseline)
# ---------------------------------------------------------------------------

cat("--- Task 1.4: Linear Regression Baseline ---\n\n")

linear_model <- lm(y ~ x1 + x2 + x3 + x4)
linear_summary <- summary(linear_model)

linear_r2 <- linear_summary$r.squared
linear_adj_r2 <- linear_summary$adj.r.squared
linear_fstat <- linear_summary$fstatistic

cat("Linear Model R-squared:", round(linear_r2, 4), "\n")
cat("Linear Model Adj R-squared:", round(linear_adj_r2, 4), "\n")
cat("Linear Model F-statistic:", round(linear_fstat[1], 4), "\n\n")

# Linear model diagnostic plots (Redesigned with ggplot2)
df_lm <- data.frame(
     Fitted = fitted(linear_model),
     Actual = y,
     Residuals = residuals(linear_model)
)

p_lm1 <- ggplot(df_lm, aes(x = Fitted, y = Actual)) +
     geom_point(alpha = 0.5, color = "steelblue", size = 2) +
     geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 0.8) +
     geom_smooth(method = "loess", color = "firebrick", linewidth = 1.2, se = FALSE) +
     labs(title = "Actual vs Fitted", x = "Fitted Values", y = "Actual Y") +
     theme_minimal(base_size = 12) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_lm2 <- ggplot(df_lm, aes(x = Fitted, y = Residuals)) +
     geom_point(alpha = 0.5, color = "indianred", size = 2) +
     geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
     geom_smooth(method = "loess", color = "firebrick", linewidth = 1.2, se = FALSE) +
     labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
     theme_minimal(base_size = 12) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))

grid.arrange(p_lm1, p_lm2, ncol = 2)

# ============================================================================
# ============================================================================
cat("\n============================================\n")
cat("TASK 2: REGRESSION MODELLING\n")
cat("============================================\n\n")

# ---------------------------------------------------------------------------
# Define the 5 candidate polynomial regression models
# ---------------------------------------------------------------------------

# Construct the design matrix for each candidate model

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
#     theta_hat = (X^T X)^{-1} X^T y

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
#                         - RSS / (2*sigma_hat^2)

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
#     AIC = 2*k - 2*ln p(D|theta_hat)
#     BIC = k*ln(n) - 2*ln p(D|theta_hat)

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

cat("--- Task 2.5: Residual Distribution Analysis ---\n\n")

residuals_list <- list()

library(grid)

colors_fill <- c("#e0e2ff", "lightgreen", "lightcoral", "thistle", "moccasin")
colors_line <- c("deepskyblue4", "forestgreen", "firebrick", "purple4", "darkorange3")

# Histograms of residuals (Task 1.2 Redesign)
hist_plots <- list()
for (i in 1:5) {
     resid_i <- as.numeric(y - y_hats[[i]])
     residuals_list[[i]] <- resid_i

     df_resid <- data.frame(Residual = resid_i)

     p_hist <- ggplot(df_resid, aes(x = Residual)) +
          geom_histogram(aes(y = after_stat(density)), fill = colors_fill[i], color = "white", bins = 20, alpha = 0.7) +
          stat_function(fun = dnorm, args = list(mean = mean(resid_i), sd = sd(resid_i)), color = colors_line[i], linewidth = 1) +
          geom_rug(sides = "b", color = "black", alpha = 0.2) +
          labs(title = paste("Model", i), x = "Residual", y = "Density") +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))

     hist_plots[[i]] <- p_hist
}

# Save histogram grid as PNG
p_hist_grid <- grid.arrange(grobs = hist_plots, ncol = 2, top = textGrob("Histograms of Residuals", gp = gpar(fontsize = 16, fontface = "bold")))

# Q-Q Plots (Task 1.2 Redesign)
qq_plots <- list()
for (i in 1:5) {
     df_resid <- data.frame(Residual = residuals_list[[i]])

     p_qq <- ggplot(df_resid, aes(sample = Residual)) +
          stat_qq(color = colors_line[i], alpha = 0.7, size = 1.5) +
          stat_qq_line(color = "black", linewidth = 1, linetype = "dashed") +
          labs(title = paste("Model", i), x = "Theoretical Quantiles", y = "Sample Quantiles") +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))

     qq_plots[[i]] <- p_qq
}

# Save Q-Q grid as PNG
p_qq_grid <- grid.arrange(grobs = qq_plots, ncol = 2, top = textGrob("Q-Q Plots of Residuals", gp = gpar(fontsize = 16, fontface = "bold")))

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
if (preferred_aic == preferred_bic) {
     preferred_model <- preferred_aic
     cat(sprintf("\nBoth AIC and BIC agree: %s is the preferred model.\n", model_names[preferred_model]))
} else {
     preferred_model <- preferred_bic
     cat(sprintf(
          "\nAIC and BIC disagree. Selecting %s (BIC preference) as the ",
          model_names[preferred_model]
     ))
}

cat(sprintf("\n*** SELECTED PREFERRED MODEL: %s ***\n\n", model_names[preferred_model]))
# ---------------------------------------------------------------------------
# Task 2.7: Train/Test Split and Prediction with 95% Confidence Intervals
# ---------------------------------------------------------------------------

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
#     y_hat_i +/- t_{alpha/2, n_train - k} * sigma_hat *

alpha <- 0.05
t_crit <- qt(1 - alpha / 2, df = n_train - k_preferred)

XtX_inv <- solve(t(X_train) %*% X_train) # (X^T X)^{-1}

pred_se <- numeric(n_test)
for (i in 1:n_test) {
     x_i <- matrix(X_test[i, ], ncol = 1)
     pred_var_i <- sigma2_train * (1 + t(x_i) %*% XtX_inv %*% x_i)
     pred_se[i] <- sqrt(pred_var_i)
}

# Confidence interval bounds
CI_lower <- y_pred_test - t_crit * pred_se
CI_upper <- y_pred_test + t_crit * pred_se

cat(sprintf("Critical t-value (alpha=0.05, df=%d): %.4f\n\n", n_train - k_preferred, t_crit))

sort_order <- order(test_indices)

df_pred <- data.frame(
     Index = 1:n_test,
     Actual = y_test[sort_order],
     Predicted = as.numeric(y_pred_test[sort_order]),
     CI_Lower = CI_lower[sort_order],
     CI_Upper = CI_upper[sort_order]
)

p_pred <- ggplot(df_pred, aes(x = Index)) +
     geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, color = "95% CI"), width = 0.3, linewidth = 0.6, alpha = 0.8) +
     geom_line(aes(y = Predicted, color = "Predicted"), linewidth = 1) +
     geom_point(aes(y = Actual, color = "Actual (Test)", shape = "Actual (Test)"), size = 2.5, alpha = 0.9) +
     geom_point(aes(y = Predicted, color = "Predicted", shape = "Predicted"), size = 2.5, stroke = 1.2) +
     scale_color_manual(name = "", values = c("Actual (Test)" = "#2C3E50", "Predicted" = "#E74C3C", "95% CI" = "#BDC3C7")) +
     scale_shape_manual(name = "", values = c("Actual (Test)" = 16, "Predicted" = 4, "95% CI" = 32)) +
     labs(
          title = paste("Model Prediction with 95% CI -", model_names[preferred_model]),
          subtitle = "Actual Test Data vs Predictions",
          x = "Test Sample Index",
          y = "Amplitude"
     ) +
     theme_minimal(base_size = 14) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold", color = "#2C3E50", size = 16),
          plot.subtitle = element_text(hjust = 0.5, color = "#7F8C8D", size = 12),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12, color = "#34495E"),
          axis.title = element_text(face = "bold", color = "#34495E"),
          axis.text = element_text(color = "#7F8C8D"),
          panel.grid.major = element_line(color = "#D5DBDB", linetype = "dotted", linewidth = 0.8),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
     )

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
prior_half_width_1 <- 3 * abs_theta[top2_indices[1]]
prior_half_width_2 <- 3 * abs_theta[top2_indices[2]]

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
     df_abc <- data.frame(
          theta1 = accepted_theta1,
          theta2 = accepted_theta2
     )

     param1_name <- paste0("theta_", top2_indices[1])
     param2_name <- paste0("theta_", top2_indices[2])

     p_joint <- ggplot(df_abc, aes(x = theta1, y = theta2)) +
          geom_point(color = "steelblue", alpha = 0.3, size = 1.5) +
          geom_point(aes(x = theta_preferred[top2_indices[1]], y = theta_preferred[top2_indices[2]]),
               color = "red", shape = 4, size = 4, stroke = 2
          ) +
          labs(title = "Joint Posterior Distribution (Rejection ABC)", x = param1_name, y = param2_name) +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))

     # 4b. Marginal posterior for theta_1
     p_marg1 <- ggplot(df_abc, aes(x = theta1)) +
          geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "lightsteelblue", color = "white") +
          geom_density(color = "steelblue", linewidth = 1) +
          geom_vline(xintercept = theta_preferred[top2_indices[1]], color = "red", linetype = "dashed", linewidth = 1) +
          labs(title = paste("Marginal Posterior -", param1_name), x = param1_name, y = "Density") +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))

     # 4c. Marginal posterior for theta_2
     p_marg2 <- ggplot(df_abc, aes(x = theta2)) +
          geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "mistyrose", color = "white") +
          geom_density(color = "firebrick", linewidth = 1) +
          geom_vline(xintercept = theta_preferred[top2_indices[2]], color = "red", linetype = "dashed", linewidth = 1) +
          labs(title = paste("Marginal Posterior -", param2_name), x = param2_name, y = "Density") +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))

     # Combine using gridExtra
     p_abc_grid <- grid.arrange(p_joint, p_marg1, p_marg2,
          layout_matrix = rbind(c(1, 1), c(2, 3)),
          top = grid::textGrob("ABC Joint and Marginal Posteriors", gp = grid::gpar(fontsize = 16, fontface = "bold"))
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
          df_abc <- data.frame(
               theta1 = accepted_theta1,
               theta2 = accepted_theta2
          )

          param1_name <- paste0("theta_", top2_indices[1])
          param2_name <- paste0("theta_", top2_indices[2])

          mean_theta1 <- mean(df_abc$theta1)
          mean_theta2 <- mean(df_abc$theta2)
          y_range <- max(df_abc$theta2) - min(df_abc$theta2)

          p_joint <- ggplot(df_abc, aes(x = theta1, y = theta2)) +
               # 1. Filled density contours (continuous blue gradient)
               stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.9) +
               scale_fill_distiller(palette = "Blues", direction = 1) +
               # 2. Red contour lines
               geom_density_2d(colour = "red", linewidth = 0.5) +
               # 3. Scatter points (semi-transparent blue)
               geom_point(color = "steelblue", alpha = 0.4, size = 1) +
               # 4. Dashed lines intersecting at the mean
               geom_vline(xintercept = mean_theta1, linetype = "dashed", color = "gray50") +
               geom_hline(yintercept = mean_theta2, linetype = "dashed", color = "gray50") +
               # 5. Posterior Mean point (Yellow circle with black border)
               geom_point(aes(x = mean_theta1, y = mean_theta2),
                    fill = "gold", color = "black", shape = 21, size = 4, stroke = 1.2
               ) +
               # 6. Text label for Posterior Mean
               annotate("text",
                    x = mean_theta1, y = mean_theta2 + y_range * 0.04,
                    label = "Posterior Mean", fontface = "bold"
               ) +
               labs(
                    title = "Joint Posterior Distribution",
                    subtitle = "Accepted samples obtained using the Rejection ABC algorithm",
                    x = param1_name,
                    y = param2_name,
                    fill = "Posterior\nDensity"
               ) +
               theme_classic(base_size = 15) +
               theme(
                    plot.title.position = "plot",
                    plot.title = element_text(
                         face = "bold",
                         size = 18,
                         hjust = 0.5
                    ),
                    plot.subtitle = element_text(
                         size = 14,
                         hjust = 0.5
                    ),
                    axis.title = element_text(
                         face = "bold"
                    ),
                    axis.text = element_text(
                         colour = "black"
                    ),
                    axis.line = element_line(
                         colour = "black"
                    ),
                    legend.position = "right",
                    legend.title = element_text(
                         face = "bold"
                    )
               )

          p_marg1 <- ggplot(df_abc, aes(x = theta1)) +
               geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "lightsteelblue", color = "white") +
               geom_density(color = "steelblue", linewidth = 1) +
               geom_vline(xintercept = theta_preferred[top2_indices[1]], color = "red", linetype = "dashed", linewidth = 1) +
               labs(title = paste("Marginal Posterior -", param1_name), x = param1_name, y = "Density") +
               theme_minimal(base_size = 15) +
               theme(
                    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                    axis.title = element_text(face = "bold"),
                    axis.text = element_text(colour = "black")
               )

          p_marg2 <- ggplot(df_abc, aes(x = theta2)) +
               geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "mistyrose", color = "white") +
               geom_density(color = "firebrick", linewidth = 1) +
               geom_vline(xintercept = theta_preferred[top2_indices[2]], color = "red", linetype = "dashed", linewidth = 1) +
               labs(title = paste("Marginal Posterior -", param2_name), x = param2_name, y = "Density") +
               theme_minimal(base_size = 15) +
               theme(
                    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                    axis.title = element_text(face = "bold"),
                    axis.text = element_text(colour = "black")
               )

          p_abc_grid <- grid.arrange(p_joint, p_marg1, p_marg2,
               layout_matrix = rbind(c(1, 1), c(2, 3)),
               top = grid::textGrob("ABC Joint and Marginal Posteriors", gp = grid::gpar(fontsize = 20, fontface = "bold"))
          )

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
