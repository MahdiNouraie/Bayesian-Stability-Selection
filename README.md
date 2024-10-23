# Bayesian Stability Selection Paper Repository

This repository contains the accompanying code for the paper titled *"Bayesian Stability Selection and Inference on Inclusion Probabilities"*. All scripts are written in R. 

**Disclaimer: For reproducing the results reported in the paper, we can only guarantee consistency if the code is run on posit Cloud. The results may differ if executed on other systems due to variations in system settings** ([Posit link](https://posit.cloud/content/9064090)).



**The sequence of results presented in the paper follows the order: first, the `Synthetic_data.R` file, followed by the `Riboflavin.R` file, and finally with the `Rat.R` file.**.

### Overview of Files:

**1. `Synthetic_data.R`**  
This script generates synthetic the synthetic data sets in the paper, applies stability selection using elastic net, and tracks the selection frequencies of variables across data sets. At the end, the Bayesian stability selection is used to infer inclusion probabilities for variables.

**2. `Riboflavin.R`**  
This script applies stability selection using RLARS on the riboflavin data set. It incorporates prior information from a Arashi rt al (2021) and employs the Bayesian stability selection to infer inclusion probabilities for the genes.

**3. `Rat.R`**  
This script performs stability selection using LASSO on the rat microarray data, following appropriate preprocessing steps. The Bayesian approach is then applied to estimate inclusion probabilities, providing insights into the significance of the variables.

---







