# Bayesian Stability Selection

This repository accompanies the paper *"Bayesian Stability Selection and Inference on Selection Probabilities"*, providing the R scripts necessary for reproducing the analysis presented. [Preprint Version](https://arxiv.org/pdf/2410.21914)

**Note:** Reproducibility of results is assured when the code is executed on Posit Cloud. System differences may result in variation, particularly with dependent package installations. We have observed that certain packages may not install on all systems. [Project on Posit Cloud](https://posit.cloud/content/11667810)

---

### File Descriptions:

**1. `Figure_1_2_3.R`**  
This script genrates plots presented in the methodology section of the paper.

**2. `Scenario_1.R`**  
This script genrates results and plots presented for the first synthetic scenario presented in the paper.

**3. `Scenario_2.R`**  
This script genrates results and plots presented for the second synthetic scenario presented in the paper.

**4. `Riboflavin_1.R`**  
This script applies Bayesian stability selection using RLARS and uniform priors on the riboflavin data set. It reproduces the results presented in Table 1 of the paper.

**5. `Riboflavin_2.R`**  
In this script, Bayesian stability selection is performed using RLARS with informative priors on the riboflavin data set. It reproduces the results presented in Table 2 of the paper.

**6. `Rat_1.R`**  
This script uses Bayesian stability selection with RLARS and uniform priors on the rat microarray data. It reproduces the results presented in Table 3 of the paper.

**7. `Rat_2.R`**  
This script begins by removing outliers from the rat microarray data using the isolation tree algorithm. Bayesian stability selection is then applied with LASSO and uniform priors to the cleaned data. It reproduces the results presented in Table 4 of the paper.

---








