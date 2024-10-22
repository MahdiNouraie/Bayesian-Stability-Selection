# This code generates synthetic data and 
#performs Bayesian Stability Selection using Elastic Net regression.
# n_dataset denotes the number of datasets generated.
# For getting the results similar to those reported in the manuscript,
# you need to set n_datasets = 100. It takes a few minutes to run.
# For faster results, you may set n_datasets = 10.

#################  Install and Load Necessary Libraries #################
# Install necessary libraries if not already installed
if(!require("MASS")){install.packages("MASS")}
if(!require("glmnet")){install.packages("glmnet")}

# Load necessary libraries
library(MASS)  # For multivariate Normal sampling
library(glmnet)  # For Elastic Net regression

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
#[1] glmnet_4.1-8   Matrix_1.6-1.1 MASS_7.3-60   

#loaded via a namespace (and not attached):
#[1] compiler_4.3.2    tools_4.3.2       survival_3.7-0    rstudioapi_0.16.0 Rcpp_1.0.13      
#[6] splines_4.3.2     codetools_0.2-20  grid_4.3.2        iterators_1.0.14  foreach_1.5.2    
#[11] shape_1.4.6.1     lattice_0.22-6   

################# Set Parameters #################
p <- 500 # Number of predictors
n <- 50 # Number of observations
b <- 100 # Number of subsamples
n_datasets <- 10 # Number of datasets

# Initialize matrix for selection frequencies over datasets
Selection_Frequency <- matrix(0, nrow = n_datasets, ncol = p)

# Create the covariance matrix Σ
Sigma <- diag(p)  # Identity matrix of size p x p
# Specify the non-zero off-diagonal elements
Sigma[1, 2] <- Sigma[2, 1] <- 0.8
Sigma[3, 4] <- Sigma[4, 3] <- 0.8
Sigma[3, 5] <- Sigma[5, 3] <- 0.8
Sigma[4, 5] <- Sigma[5, 4] <- 0.8

# Define the coefficient vector
beta <- c(rep(0.9, 2), rep(0.7, 3), 1.5, rep(0, p - 6))  

################# Data Generation and Bayesian Stability Selection #################
# Generate synthetic datasets and loop through to calculate selection frequencies
for (i in 1:n_datasets) {
  set.seed(26 + i)  # Set a different seed for each dataset
  # Sample predictors X from multivariate normal distribution N(0, Σ)
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # Generate the response Y using the first 5 predictors and adding Gaussian noise
  epsilon <- rnorm(n, mean = 0, sd = 2)  # Gaussian noise with variance 4 (sd = 2)
  Y <- X %*% beta + epsilon
  
  # Combine predictors and response into a single data frame
  data <- data.frame(X, Y)
  
  set.seed(26)  # Set seed for reproducibility
  # Perform cross-validated Elastic Net
  cv_elastic <- cv.glmnet(X, Y, nfolds = 10, alpha = 0.2)
  
  # Get the optimal lambda value
  lambda <- cv_elastic$lambda.1se
  
  # Initialize matrix S for stability selection for current dataset
  S <- matrix(data = 0, nrow = b, ncol = p)  
  for (j in 1:b) {  # Number of subsamples
    # Sub-sampling
    model_data <- data[sample(1:nrow(data), nrow(data) / 2, replace = FALSE), ]
    
    # Prepare the predictor matrix x and response vector y for the subsample
    x_sub <- as.matrix(model_data[, 1:ncol(model_data) - 1])
    y_sub <- model_data$Y
    
    # Fit the final LASSO model with the optimal lambda value
    elastic_model <- glmnet(x_sub, y_sub, alpha = 0.2, lambda = lambda)
    
    # Determine significant predictors (excluding intercept)
    significant_predictors <- ifelse(coef(elastic_model) != 0, 1, 0)[-1]
    
    # Store significant predictors in matrix S
    S[j, ] <- significant_predictors
  }
  
  # Store Inclusion Probability for the dataset in the corresponding row
  Selection_Frequency[i, ] <- colMeans(S)
}

# The number of times the first 6 predictors were selected
n_j <- round(colMeans(Selection_Frequency)[1:6] * 100)
# Using Priors for the first 6 predictors
# Assuming the answer to the first question is 50% for the first 6 predictors and,
# Assuming the perceived importance of the first 6 predictors is 70%
alpha_j <- rep(70, 6)
beta_j <- rep(30, 6)

# The average selection frequency of the first 6 predictors over all datasets
cat("Selection Frequency: ", round(colMeans(Selection_Frequency)[1:6], 3),'\n')
# Posterior Inclusion Probability
cat("Inclusion Probability: ", round((n_j + alpha_j) / (b + alpha_j + beta_j), 3))
