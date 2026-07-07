# EEG Assignment Complete

This project contains one main R script, [EEG_Assignment_Complete.R](EEG_Assignment_Complete.R), that analyses simulated EEG data using polynomial regression.

The explanation below is written for beginners and follows the script step by step. It focuses on what the code is doing, in the same order as the script, using simple language.

## 1. What this project does

The script:

1. Loads the EEG dataset from the CSV files in `Data_Set/`.
2. Plots the signals so you can see how they behave over time.
3. Checks the distribution and correlation of the signals.
4. Fits a simple linear regression model as a starting point.
5. Tries 5 polynomial regression models.
6. Compares those models using RSS, AIC, and BIC.
7. Chooses the best model.
8. Splits the data into training and testing sets.
9. Measures prediction accuracy and builds 95% confidence intervals.
10. Runs Approximate Bayesian Computation (ABC) for the most important parameters.

## 2. Files in the project

- [EEG_Assignment_Complete.R](EEG_Assignment_Complete.R) - the main analysis script.
- `Data_Set/time_1673241270748.csv` - time values for the samples.
- `Data_Set/X_1673241366257.csv` - input signals `x1`, `x2`, `x3`, `x4`.
- `Data_Set/y_1673241374123.csv` - output signal `y`.
- `Output_Images/` - place where you can save plots.
- `Table_Screenshot/` - place where you can save tables.

## 3. How to run it

1. Open the project folder in RStudio or another R environment.
2. Make sure the working directory is the project root folder.
3. Install the required packages if they are missing: `ggplot2`, `corrplot`, `gridExtra`, and `reshape2`.
4. Run the script from top to bottom.
5. Read the console output and inspect the plots that appear.

If the CSV files are not found, check that the script is being run from the same folder that contains the `Data_Set/` directory.

## 4. Step-by-step explanation of the script

### Step 0. Setup

The script first loads the required R packages. If a package is missing, the script installs it and then loads it.

What the code is doing:

- `required_packages <- c(...)` creates a list of package names.
- `for (pkg in required_packages)` goes through that list one package at a time.
- `require(...)` checks whether the package is already available.
- `install.packages(...)` installs the package only if it is missing.
- `library(...)` loads the package so its functions can be used later.

After that, it reads the three CSV files and stores them in memory:

- `X_data` for the input signals.
- `y_data` for the output signal.
- `time_data` for the time values.

The script then renames the columns to make them easier to understand and converts the data into numeric form.

What the code is doing:

- `read.csv(...)` loads the raw files.
- `colnames(...) <- ...` gives the columns meaningful names.
- `as.numeric(...)` converts each column into numeric values.
- `n <- length(y)` counts how many samples are in the data.
- `data.frame(...)` combines everything into one clean table called `eeg_df`.

### Step 1. Preliminary data analysis

This part helps you understand the data before modelling.

#### 1.1 Plot the signals over time

The script creates line plots for `x1`, `x2`, `x3`, `x4`, and `y` against time.

What the code is doing:

- `par(mfrow = c(3, 2), ...)` arranges the plotting area into multiple panels.
- `plot(t, x1, type = "l", ...)` draws a line plot for signal `x1`.
- The same idea is repeated for `x2`, `x3`, `x4`, and `y`.
- `grid(...)` adds a light grid to make the graph easier to read.

Why this matters:

- You can see whether the signals are stable or changing.
- You can spot unusual spikes or patterns.

#### 1.2 Check distributions

The script draws histograms with density curves for each signal.

What the code is doing:

- `hist(...)` shows how often values appear in each signal.
- `freq = FALSE` makes the histogram show density instead of raw counts.
- `density(...)` estimates a smooth curve for the distribution.
- `lines(density(...))` draws that smooth curve on top of the histogram.

Why this matters:

- You can see whether the data looks roughly normal.
- You can compare the spread of the input and output signals.

It also prints summary statistics:

- Mean
- Standard deviation
- Minimum value
- Maximum value

#### 1.3 Check correlation

The script calculates the Pearson correlation between all signals and shows:

- A correlation matrix
- A heatmap
- Scatter plots between `y` and each input signal
- Scatter plots between input signals

Why this matters:

- Correlation shows whether two variables move together.
- The scatter plots help you see whether a linear or curved relationship may exist.

#### 1.4 Fit a linear regression baseline

The script fits a basic multiple linear regression model:

`y ~ x1 + x2 + x3 + x4`

What the code is doing:

- `lm(...)` fits a standard linear regression model.
- `summary(...)` gives useful model results such as R-squared.
- `fitted(...)` gets the predicted values from the model.
- `residuals(...)` gets the differences between actual and predicted values.
- The plots compare the fitted values with the real values and check whether residuals look random.

This is only a starting point.

Why this matters:

- It gives a baseline result.
- Later polynomial models can be compared against it.
- The script also shows fitted vs actual values and residuals vs fitted values.

## 5. Task 2: Polynomial regression modelling

This is the main modelling part of the assignment.

### 2.0 Build 5 candidate models

The script creates 5 different design matrices.

What the code is doing:

- `cbind(...)` joins several predictor columns together.
- `x1^2`, `x1^3`, `x2^4`, and similar terms create polynomial features.
- The final column of `1` adds the intercept term.
- `colnames(...) <- ...` labels each column so the output is easier to read.
- `model_designs <- list(...)` stores all 5 models together for easy looping.

Each model uses different combinations of polynomial terms such as:

- `x1^2`
- `x1^3`
- `x2^4`
- `x3^4`
- `x4`

Why this matters:

- EEG relationships are often not purely linear.
- Polynomial terms can capture curved patterns.

Each model also includes a bias term, which acts like an intercept.

### 2.1 Estimate the parameters

For each model, the script uses ordinary least squares to estimate the coefficients.

Simple idea:

- The model tries to choose coefficient values that make the prediction errors as small as possible.

What the code is doing:

- `solve(t(X) %*% X) %*% t(X) %*% y` is the ordinary least squares formula.
- `theta_hat` stores the estimated coefficients for one model.
- The loop repeats this for all 5 candidate models.

The script prints every estimated parameter.

### 2.2 Compute RSS

RSS means Residual Sum of Squares.

Formula idea:

`RSS = sum((actual - predicted)^2)`

What the code is doing:

- `y_hat <- X_mat %*% theta_hats[[i]]` computes the fitted values.
- `sum((y - y_hat)^2)` adds up the squared errors.
- Smaller RSS means the model fits better.

Why this matters:

- Lower RSS means the model fits the data better.
- It is one of the first measures used to compare models.

### 2.3 Compute log-likelihood

The script assumes the errors are Gaussian and calculates the log-likelihood for each model.

What the code is doing:

- `sigma2_vals[i] <- RSS_vals[i] / (n - 1)` estimates the error variance.
- `loglik_vals[i] <- ...` computes the log-likelihood formula.
- The loop stores the result for every candidate model.

Why this matters:

- Log-likelihood is used later in AIC and BIC.
- It helps measure how well the model explains the data.

### 2.4 Compute AIC and BIC

The script calculates:

- AIC: Akaike Information Criterion
- BIC: Bayesian Information Criterion

Why this matters:

- Both measures reward good fit.
- Both also penalise overly complex models.
- Lower values are better.

The script then prints a comparison table for all 5 models.

What the code is doing:

- `k_vals <- sapply(theta_hats, length)` counts the number of parameters in each model.
- `AIC_vals[i] <- 2 * k_vals[i] - 2 * loglik_vals[i]` calculates AIC.
- `BIC_vals[i] <- k_vals[i] * log(n) - 2 * loglik_vals[i]` calculates BIC.
- `comparison_df <- data.frame(...)` creates a simple table for comparison.

### 2.5 Check residual normality

The script studies the residuals using:

- Histograms
- Q-Q plots
- Shapiro-Wilk test

Why this matters:

- Good regression models often have residuals that look approximately normal.
- If residuals are badly non-normal, the model may be missing something.

### 2.6 Select the preferred model

The script chooses the best model using AIC and BIC.

Main idea:

- If AIC and BIC agree, that model is selected.
- If they differ, the script prefers the BIC choice because BIC penalises complexity more strongly.

The script also uses the residual plots and Shapiro-Wilk test as supporting evidence.

What the code is doing:

- `which.min(AIC_vals)` finds the model with the smallest AIC.
- `which.min(BIC_vals)` finds the model with the smallest BIC.
- `if (preferred_aic == preferred_bic)` checks whether both criteria agree.
- If they disagree, the script chooses the BIC winner because it is stricter.

### 2.7 Train/test split and prediction intervals

The data is split into:

- 70% training data
- 30% testing data

The chosen model is re-fit using only the training set.

What the code is doing:

- `set.seed(42)` makes the random split repeatable.
- `sample(...)` selects training indices.
- `setdiff(...)` finds the remaining test indices.
- `theta_hat_train <- ...` re-fits the model on the training data only.
- `X_test %*% theta_hat_train` produces predictions for unseen samples.

Then the script:

- Predicts the test values
- Computes test RSS
- Calculates MAE and RMSE
- Draws 95% prediction intervals

Why this matters:

- Training accuracy alone is not enough.
- Testing shows how well the model works on unseen data.
- Confidence intervals show how uncertain the predictions are.

## 6. Task 3: Approximate Bayesian Computation (ABC)

ABC is used here to estimate approximate posterior distributions for the two most important parameters in the selected model.

### 3.1 Choose the two strongest parameters

The script finds the two coefficients with the largest absolute values.

What the code is doing:

- `abs(as.numeric(theta_preferred))` gets the size of each coefficient.
- `order(..., decreasing = TRUE)` sorts them from largest to smallest.
- The top two parameters are selected for the ABC step.

Why this matters:

- These are treated as the most influential parameters in the model.

### 3.2 Set uniform priors

The script creates broad uniform prior ranges around the least squares estimates.

What the code is doing:

- `prior_half_width_1 <- 3 * abs_theta[...]` makes the prior wide enough to explore values around the estimate.
- `max(..., 0.5)` prevents the prior from becoming too narrow.
- `prior_lower_1`, `prior_upper_1`, `prior_lower_2`, and `prior_upper_2` define the allowed ranges.

Why this matters:

- A prior gives a plausible starting range for the parameter values.

### 3.3 Run rejection ABC

The script repeatedly:

1. Randomly draws values for the two selected parameters.
2. Keeps the other parameters fixed.
3. Simulates the output.
4. Computes the simulated RSS.
5. Accepts the draw if the RSS is small enough.

Why this matters:

- Accepted draws act like samples from an approximate posterior distribution.

What the code is doing:

- `runif(...)` draws random parameter values from the prior ranges.
- The selected coefficients are replaced in `theta_new`.
- `y_sim <- X_abc %*% theta_new` creates a simulated output.
- `RSS_sim <- sum((y - y_sim)^2)` measures how close the simulation is to the real data.
- `if (RSS_sim < epsilon)` keeps only good samples.

### 3.4 Plot posterior results

If samples are accepted, the script shows:

- A joint posterior scatter plot
- Marginal posterior histograms
- Posterior summary statistics

These results help you understand how the selected parameters may vary around their least squares estimates.

What the code is doing:

- The scatter plot shows how the two accepted parameters move together.
- The histograms show the distribution of each parameter separately.
- `mean(...)`, `median(...)`, `sd(...)`, and `quantile(...)` summarize the posterior samples.
- The red line marks the original least squares estimate for comparison.

## 7. Beginner summary of the whole workflow

In simple words, the script does this:

1. Load the EEG data.
2. Look at the data with plots and summary statistics.
3. Try several polynomial regression models.
4. Compare the models using fit quality and complexity penalties.
5. Choose the best model.
6. Test the chosen model on unseen data.
7. Use ABC to study the uncertainty of the most important parameters.

If you want to read the script like a beginner, follow it in this order:

1. Load packages and data.
2. Plot and inspect the signals.
3. Check correlation and distributions.
4. Fit the simple linear model.
5. Build and compare the 5 polynomial models.
6. Choose the best model.
7. Test the model on unseen data.
8. Run ABC on the two most important parameters.

## 8. Helpful notes

- If you want the plots to be saved automatically, you can add `png()` or `pdf()` commands before the plotting sections.
- If the working directory causes file-loading problems, set it to the project root before running the script.
- The script may take time during the ABC step because it uses 100,000 samples.

## 9. Output folders

You can use these folders to organise results:

- `Output_Images/` for figures and plots
- `Table_Screenshot/` for tables or printed summaries
- `Console_Screenshoot/` for console captures

## 10. Final note

This project is mostly about understanding the full statistical workflow, not just getting a final number. The important idea is to compare models carefully, check assumptions, and then validate the chosen model on test data.