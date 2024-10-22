# Bayesian Stability Selection Paper Repository

This repository contains the accompanying code for the paper titled *"Bayesian Stability Selection and Inference on Inclusion Probabilities"*. All scripts are written in R. 
**We recomment to start by running the `Synthetic_data.R` file**.

### Overview of Files:

**1. `Synthetic_data.R`**  
This script generates synthetic data sets, applies stability selection using LASSO, and tracks the selection frequencies of variables across data sets. At the end, the Bayesian method is used to infer inclusion probabilities for variables.

**2. `Riboflavin.R`**  
This script applies stability selection using RLARS on the riboflavin data set. It incorporates prior information from a related study and employs the Bayesian tool to estimate inclusion probabilities for the variables.

**3. `Rat.R`**  
This script performs stability selection using LASSO on the rat microarray data, following appropriate preprocessing steps. The Bayesian approach is then applied to estimate inclusion probabilities, providing insights into the significance of the variables.

---







