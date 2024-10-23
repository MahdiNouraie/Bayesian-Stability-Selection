# This code applies Bayesian Stability Selection to the Riboflavin data set
# The RLARS algorithm is used to identify significant predictors
# For reproducing the results reported in Table 2 of the manuscript,
# you need to set b = 1000. It takes a few minutes to run.
# For faster results, you may set b = 10.
#################  Install and Load Necessary Libraries #################
# Don't display warnings
options(warn=-1)
# Install necessary libraries if not already installed
if(!require("hdi")){install.packages("hdi")}
if(!require("robustHD")){install.packages("robustHD")}

# Load necessary libraries
library(hdi)  # For getting riboflavin dataset
library(robustHD)  # For RLARS

# Session Info
# Session Info
#R version 4.3.2 (2023-10-31)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS 15.0.1

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#attached base packages:
#[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] robustHD_0.8.1      robustbase_0.99-4-1 perry_0.3.1         ggplot2_3.5.1       hdi_0.1-9          
#[6] scalreg_1.0.1       lars_1.3           

#loaded via a namespace (and not attached):
#[1] Matrix_1.6-1.1    glmnet_4.1-8      gtable_0.3.5      dplyr_1.1.4       compiler_4.3.2   
#[6] tidyselect_1.2.1  Rcpp_1.0.13       splines_4.3.2     scales_1.3.0      lattice_0.22-6   
#[11] R6_2.5.1          linprog_0.9-4     generics_0.1.3    shape_1.4.6.1     iterators_1.0.14 
#[16] MASS_7.3-60       tibble_3.2.1      munsell_0.5.1     pillar_1.9.0      rlang_1.1.4      
#[21] utf8_1.2.4        cli_3.6.3         withr_3.0.1       magrittr_2.0.3    foreach_1.5.2    
#[26] grid_4.3.2        rstudioapi_0.16.0 lifecycle_1.0.4   DEoptimR_1.1-3    vctrs_0.6.5      
#[31] lpSolve_5.6.21    glue_1.8.0        codetools_0.2-20  survival_3.7-0    fansi_1.0.6      
#[36] colorspace_2.1-1  tools_4.3.2       pkgconfig_2.0.3  

#################  Bayesian Stability Selection #################
# Load the Riboflavin dataset
data(riboflavin)
# Convert the dataset to a data frame
riboflavin <- as.data.frame(cbind(Y=riboflavin$y - 1, X=riboflavin$x))

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

# Set seed for reproducibility
set.seed(26)

p <- ncol(riboflavin) -1 # Number of predictors
b <- 1000 # Number of subsamples

# Initialize selection matrix S
S <- matrix(data = 0, nrow = b, ncol = p)

# Stability Selection
for (i in 1:b) {
  # Sub-sampling
  model_data <- riboflavin[sample(1:nrow(riboflavin), nrow(riboflavin)/2, replace = FALSE), ]
  # Prepare the response variable
  x <- model_data[,-1]
  y <- model_data[,1]
  # Standardise the predictors
  x <- scale(x)
  # Fit the RLARS model
  rlars_model <- rlars(x, y)
  # Get the active predictors
  active <- rlars_model$active
  # Create a binary vector of significant predictors
  significant_predictors <- rep(0, p)
  # Set the significant predictors to 1
  significant_predictors[active] <- 1
  # Store significant predictors in selection matrix S
  S[i, ] <- significant_predictors
}

# Selection frequency of each predictor
Selection_Frequency <- colMeans(S)
# set names on the selection frequency elements
names(Selection_Frequency) <- colnames(x)
# Display predictors with selection frequency greater than 0.5
#Selection_Frequency[Selection_Frequency > 0.5]

# set names on the selection frequency columns
colnames(S) <- colnames(x)

# 41 genes identified by Arashi et al (2021)
genes <- c('ARGF_at', 'DNAJ_at', 'GAPB_at', 'XHLB_at', 'YACN_at',
           'YBFI_at', 'LYSC_at', 'PKSA_at', 'PRIA_at', 'YCDH_at',
           'YCGO_at', 'YCKE_at', 'SPOIIAA_at', 'SPOVAA_at', 'THIK_at',
           'YCLB_at', 'YCLF_at', 'YDDH_at', 'YDDK_at', 'YEBC_at', 
           'YFHE_r_at', 'YLXW_at', 'YMFE_at', 'YOAB_at', 'YFII_at',
           'YFIO_at', 'YFIR_at', 'YPGA_at', 'YQJT_at', 'YQJU_at',
           'YHDS_r_at', 'YKBA_at', 'YKVJ_at', 'YRVJ_at', 'YTGB_at',
           'YUID_at', 'YURQ_at', 'YXLD_at', 'YXLE_at', 'YYBG_at', 'YYDA_at')
# Get the column indices of the 41 genes
column_index <- which(colnames(x) %in% genes)
# determining the prior parameters for the Beta distributions
# We assume that the for the 41 previously identified genes,
# the answer to the first question is 50% and their perceived importance is 70%
# For the remaining predictors, we use non-informative priors
prior_alphas <- rep(1, p)
prior_betas <- rep(1, p)
prior_alphas[column_index] <- c(rep(700, length(genes)))
prior_betas[column_index] <- c(rep(300, length(genes)))

# Initialize lists to store posterior parameters for each predictor
alpha_posterior_list <- list()
beta_posterior_list <- list()
# A vector to store the posterior mean of each predictor
m <- c()
# Two vectors to store the lower and upper bounds of the 95% credible interval for each predictor
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
  # Store the mean, lower and upper bounds of the posterior distribution for predictor j
  m <- c(m, stats$mean)
  L <- c(L, stats$lower)
  U <- c(U, stats$upper)
}
# create a data frame to store the output
output <- data.frame(matrix(nrow = p, ncol = 6))
# Set names on rows corresponding to predictors
rownames(output) <- colnames(S)
# Set names on columns
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
