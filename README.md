# Bayesian Stability Selection Repository

This repository accompanies the paper *"Bayesian Stability Selection and Inference on Inclusion Probabilities"*([preprint version](https://arxiv.org/pdf/2410.21914)), providing the R scripts necessary for reproducing the analysis presented. 

**Note:** Reproducibility of results is assured when the code is executed on Posit Cloud. System differences may result in variation, particularly with dependent package installations. We have observed that certain packages may not install on all systems. [Project on Posit Cloud](https://posit.cloud/content/9074125)

---

### File Descriptions:

**1. `Synthetic_data.R`**  
This script generates synthetic data sets and applies stability selection using the elastic net. It monitors variable selection frequencies across the generated data sets. Bayesian stability selection is then implemented to infer the inclusion probabilities of variables.

**2. `Riboflavin_1.R`**  
This script applies Bayesian stability selection using RLARS and non-informative priors on the riboflavin data set. It reproduces the results presented in Table 1 of the paper.

**3. `Riboflavin_2.R`**  
In this script, Bayesian stability selection is performed using RLARS with informative priors on the riboflavin data set. It reproduces the results presented in Table 2 of the paper.

**4. `Rat_1.R`**  
This script uses Bayesian stability selection with RLARS and non-informative priors on the rat microarray data. It reproduces the results presented in Table 3 of the paper.

**5. `Rat_2.R`**  
This script begins by removing outliers from the rat microarray data using the isolation tree algorithm. Bayesian stability selection is then applied with LASSO and non-informative priors to the cleaned data. It reproduces the results presented in Table 4 of the paper.

---








