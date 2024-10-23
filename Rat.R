# This code applies Bayesian Stability Selection to the rat microarray dataset
# The isolation forest algorithm is used to identify outliers
# The LASSO algorithm is used to identify significant predictors
# For reproducing the results reported in the manuscript,
# you need to set b = 1000. It takes a few minutes to run.
# For faster results, you may set b = 10.
#################  Install and Load Necessary Libraries #################
# Don't display warnings
options(warn=-1)
# Install necessary libraries if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
library(BiocManager) # For installing GEOquery
if(!require("GEOquery")){
  BiocManager::install("GEOquery")
}
if(!require("glmnet")){install.packages("glmnet")}
if(!require("isotree")){install.packages("isotree")}
# Load necessary libraries
library(GEOquery) # For accessing Gene Expression Omnibus (GEO) data
library(glmnet)  # For LASSO
library(isotree)  # For Isolation Forest


# Session Info
#R version 4.3.2 (2023-10-31)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS 15.0.1

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] isotree_0.6.1-1     glmnet_4.1-8        Matrix_1.6-1.1      GEOquery_2.70.0     Biobase_2.62.0     
#[6] BiocGenerics_0.48.1 BiocManager_1.30.25

#loaded via a namespace (and not attached):
#[1] jsonlite_1.8.9    limma_3.58.1      dplyr_1.1.4       compiler_4.3.2    tidyselect_1.2.1 
#[6] Rcpp_1.0.13       xml2_1.3.6        parallel_4.3.2    tidyr_1.3.1       splines_4.3.2    
#[11] statmod_1.5.0     lattice_0.22-6    readr_2.1.5       R6_2.5.1          generics_0.1.3   
#[16] shape_1.4.6.1     iterators_1.0.14  tibble_3.2.1      pillar_1.9.0      tzdb_0.4.0       
#[21] rlang_1.1.4       utf8_1.2.4        cli_3.6.3         magrittr_2.0.3    foreach_1.5.2    
#[26] grid_4.3.2        rstudioapi_0.16.0 hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5      
#[31] glue_1.8.0        data.table_1.16.0 codetools_0.2-20  survival_3.7-0    fansi_1.0.6      
#[36] purrr_1.0.2       tools_4.3.2       pkgconfig_2.0.3 

#################  Preprocessing the Data #################
# Load the GSE5680 dataset (takes a few seconds)
gset <- getGEO("GSE5680", GSEMatrix =TRUE, getGPL=FALSE) # takes a few seconds
# Select the first platform (GPL1355) if multiple platforms are present
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
# Extract the expression data
gset <- gset[[idx]]
# Extract the expression matrix
ex <- exprs(gset)
# Select TRIM32 gene expression values (probe 1389163_at) as the response variable
response_var <- ex["1389163_at", ]
# Remove the TRIM32 gene expression values from the expression matrix
ex <- ex[rownames(ex) != "1389163_at", ]
# Transpose the expression matrix for easier manipulation
ex <- data.frame(t(ex))
# remove unnecessary objects
rm(gset, idx);gc()

# Remove probes with max expression less than 25th percentile of all samples
max_expr <- as.numeric(apply(ex, 2, max))
threshold <- quantile(max_expr, 0.25)
filtered_ex <- ex[, max_expr > threshold]

# Remove probes with a range of expression less than 2
range_expr <- as.numeric(apply(filtered_ex, 2, function(x) diff(range(x))))
filtered_ex <- filtered_ex[, range_expr >= 2]
# Combine the response variable and filtered expression matrix
data <- cbind(Y = response_var, filtered_ex)
# Remove unnecessary objects
rm(ex, filtered_ex, max_expr, range_expr, response_var, threshold);gc()

set.seed(26)  # Set seed for reproducibility

# Fit the isolation forest model
iso_forest <- isolation.forest(data)

# Predict anomaly scores
anomaly_scores <- predict(iso_forest, data)

# Define a threshold for outliers (e.g., top 5% most anomalous points)
threshold_iso <- quantile(anomaly_scores, 0.95)

# Identify outliers
outliers_iso <- anomaly_scores > threshold_iso
# Get the indices of outliers
#which(outliers_iso)
# Remove outliers from the data
data <- data[!outliers_iso, ]
# Remove unnecessary objects
rm(iso_forest, anomaly_scores, outliers_iso, threshold_iso);gc()

#################  Bayesian Stability Selection #################

# Function to update Beta distribution parameters
# Here, success is defined as the predictor being selected and,
# failure is defined as the predictor not being selected
update_beta <- function(i , successes, failures, prior_alphas, prior_betas) {
  alpha_posterior <- prior_alphas[i] + successes
  beta_posterior <- prior_betas[i] + failures
  return(list(alpha_posterior = alpha_posterior, beta_posterior = beta_posterior))
}
# Function to calculate mean and 95% credible interval of Beta distribution
beta_stats <- function(alpha, beta) {
  # Mean of the Beta distribution
  mean <- alpha / (alpha + beta)
  
  # Lower and Upper bounds of 95% credible interval
  lower <- qbeta(0.025, alpha, beta)
  upper <- qbeta(0.975, alpha, beta)
  # Return results as a list
  return(list(mean = mean, lower = lower, upper = upper))
}

p <- ncol(data) -1 # Number of predictors
b <- 10 # Number of subsamples

x <- data[,-1]; y <- data[,1] # Separate predictors and response
x <- scale(x) # Standardise the predictors
# Cross-validated LASSO
cv_lasso <- cv.glmnet(as.matrix(x), y, nfolds = 10, alpha = 1)
# Get the optimal lambda value
lambda <- cv_lasso$lambda.1se

# Initialize selection matrix S
S <- matrix(data = 0, nrow = b, ncol = p)

# Stability Selection
for (i in 1:b) {
  # Sub-sampling
  model_data <- data[sample(1:nrow(data), nrow(data)/2, replace = FALSE), ]
  # Prepare the response variable
  x <- model_data[,-1]
  y <- model_data[,1]
  x <- scale(x) # Standardise the predictors
  lasso_model <- glmnet(x, y, alpha = 1, lambda = lambda)
  # Determine significant predictors (excluding intercept)
  significant_predictors <- ifelse(coef(lasso_model) != 0, 1, 0)[-1]
  # Store significant predictors in selection matrix S
  S[i, ] <- significant_predictors
}

# Store selection frequency predictors
Selection_Frequency <- colMeans(S)
# Set name for selection frequency elements
names(Selection_Frequency) <- colnames(x)
# Display predictors with selection frequency greater than 0.3
#Selection_Frequency[Selection_Frequency > 0.2]

# Define prior parameters for Beta distributions
prior_alphas <- rep(1, p)
prior_betas <- rep(1, p)

# Set name for the selection frequency columns
colnames(S) <- colnames(x)

# Initialize lists to store posterior parameters for each predictor
alpha_posterior_list <- list()
beta_posterior_list <- list()
#A vector to store the posterior means
m <- c()
#Two vectors to store the lower and upper bounds of the 95% credible intervals
L <- c(); U <- c()

# Perform Bayesian updating for each predictor
for (j in 1:ncol(S)) {
  # Count number of successes (1's) and failures (0's) for predictor j
  successes <- sum(S[, j])
  failures <- nrow(S) - successes
  
  # Update Beta distribution parameters
  posterior_params <- update_beta(j , successes, failures, prior_alphas, prior_betas)
  
  # Store posterior parameters for predictor j
  alpha_posterior_list[[j]] <- posterior_params$alpha_posterior
  beta_posterior_list[[j]] <- posterior_params$beta_posterior
  stats <- beta_stats(alpha_posterior_list[[j]], beta_posterior_list[[j]])
  m <- c(m, stats$mean)
  L <- c(L, stats$lower)
  U <- c(U, stats$upper)
}

# Create a data frame to display the results
output <- data.frame(matrix(nrow = p, ncol = 6))
# Set row names corresponding to predictor names
rownames(output) <- colnames(S)
# Set column names
colnames(output) <- c('Selection Frequency', 'Posterior Mean', 
                      'Posterior alpha', 'Posterior beta', 
                      'Lower Bound', 'Upper Bound')
output$`Selection Frequency` <- Selection_Frequency
output$`Posterior Mean` <- m
output$`Posterior alpha` <- unlist(alpha_posterior_list)
output$`Posterior beta` <- unlist(beta_posterior_list)
output$`Lower Bound` <- L
output$`Upper Bound` <- U
# Round the output to 3 decimal places
output <- round(output, 3)
# Order the output by posterior mean in descending order
output <- output[order(output$`Posterior Mean`, decreasing = TRUE), ]
# Display the top 10 predictors with the highest posterior mean
head(output, 10)


