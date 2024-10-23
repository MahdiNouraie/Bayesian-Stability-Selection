# Bayesian Stability Selection Repository

This repository contains the accompanying code for the paper titled *"Bayesian Stability Selection and Inference on Inclusion Probabilities"*. All scripts are written in R. 

**Disclaimer: Reproducibility of the results presented in the paper can only be ensured if the code is executed on Posit Cloud. Results may vary if run on different systems due to variations in system configurations. In particular, we have observed that certain dependent packages may fail to install on some machines.** ([Posit link](https://posit.cloud/content/9064090))



**The order of results presented in the paper proceeds as follows: first, the `Synthetic_data.R` file, followed by the `Riboflavin.R` file, and concluding with the `Rat.R` file.**

### Overview of Files:

**1. `Synthetic_data.R`**  
This script generates the synthetic data sets used in the paper, applies stability selection via elastic net, and tracks variable selection frequencies across data sets. Bayesian stability selection is employed at the end to infer the inclusion probabilities of variables.

**2. `Riboflavin.R`**  
This script performs stability selection using RLARS on the riboflavin data set. It incorporates prior information from [Arashi et al. (2021)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0245376) and uses Bayesian stability selection to infer inclusion probabilities for the genes. The script reproduces the results reported in Table 2 of the paper.

**3. `Rat.R`**  
This script begins by removing outlier samples from the rat microarray data before applying stability selection using LASSO on the cleaned data set. Bayesian stability selection is then used to infer inclusion probabilities for the probes. The script reproduces the results reported in Table 4 of the paper.


---







