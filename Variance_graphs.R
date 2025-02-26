# This script generates variance graphs presented in the 
# methodology section of the paper.
rm(list=ls()); gc() # Clean up
options(warn=-1) # Don't display warnings
if(!require("ggplot2")){install.packages("ggplot2")} 
if(!require("latex2exp")){install.packages("latex2exp")}
if(!require("reshape2")){install.packages("reshape2")}
if(!require("plotly")){install.packages("plotly")}

library(ggplot2)
library(latex2exp)
library(reshape2)
library(plotly)

#sessionInfo()
#R version 4.4.2 (2024-10-31)
#Platform: x86_64-pc-linux-gnu
#Running under: Ubuntu 20.04.6 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0

#locale:
#  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C          
#[3] LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8    
#[5] LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#[7] LC_PAPER=C.UTF-8       LC_NAME=C             
#[9] LC_ADDRESS=C           LC_TELEPHONE=C        
#[11] LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

#time zone: UTC
#tzcode source: system (glibc)

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods  
#[7] base     

#other attached packages:
#  [1] plotly_4.10.4   reshape2_1.4.4  latex2exp_0.9.6
#[4] ggplot2_3.5.1  

#loaded via a namespace (and not attached):
#  [1] gtable_0.3.6      jsonlite_1.9.0    dplyr_1.1.4      
#[4] compiler_4.4.2    tidyselect_1.2.1  Rcpp_1.0.14      
#[7] stringr_1.5.1     tidyr_1.3.1       scales_1.3.0     
#[10] fastmap_1.2.0     R6_2.6.1          plyr_1.8.9       
#[13] generics_0.1.3    htmlwidgets_1.6.4 tibble_3.2.1     
#[16] munsell_0.5.1     pillar_1.10.1     rlang_1.1.5      
#[19] stringi_1.8.4     lazyeval_0.2.2    viridisLite_0.4.2
#[22] cli_3.6.4         withr_3.0.2       magrittr_2.0.3   
#[25] digest_0.6.37     grid_4.4.2        lifecycle_1.0.4  
#[28] vctrs_0.6.5       glue_1.8.0        data.table_1.17.0
#[31] colorspace_2.1-1  purrr_1.0.4       httr_1.4.7       
#[34] tools_4.4.2       pkgconfig_2.0.3   htmltools_0.5.8.1

# Define the function to return both variance values
# V1 is variance of posterior selection probability for non-informative prior
# V2 is variance of posterior selection probability for informative prior
variance_values <- function(B, n, a) {
  V1 <- ((1 + n) * (1 + B - n)) / ((2 + B)^2 * (3 + B))
  V2 <- ((a + n) * (2*B - a - n)) / (4 * B^2 * (2 * B + 1))
  return(list(V1 = V1, V2 = V2))
}

# Figure 1a
# Fix B and n
B <- 100
n <- 30

# Create a sequence of a values from 0 to 1000
a_values <- seq(0, 100, by = 1)

# Compute V1 and V2 for each a value
variance_results <- lapply(a_values, function(a) variance_values(B, n, a))

# Convert the list to a data frame
df <- data.frame(
  a = a_values,
  Non_Informative = sapply(variance_results, function(x) x$V1),
  Informative = sapply(variance_results, function(x) x$V2)
)
# Convert data to long format for ggplot
df_long <- melt(df, id.vars = "a", variable.name = "Mode", value.name = "Value")

# Plot density curves with thicker lines
ggplot(df_long, aes(x = a, y = Value, colour = Mode)) +
  geom_density(stat = "identity", size = 1.2) +  # Increase line thickness
  labs(
    title = "Variance of Posterior Selection Probability",
    x = expression(alpha[j]),
    y = "Variance",
    colour = "Mode"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),  # Center title, increase font size
    axis.title = element_text(size = 28),  # Increase font size for axis labels
    axis.text = element_text(size = 28),  # Increase font size for axis text
    legend.title = element_text(size = 26),  # Increase legend title font size
    legend.text = element_text(size = 24)  # Increase legend text font size
  )

# Figure 1b
# Fix B and n
B <- 100
n <- 50

# Create a sequence of a values from 0 to 1000
a_values <- seq(0, 100, by = 1)

# Compute V1 and V2 for each a value
variance_results <- lapply(a_values, function(a) variance_values(B, n, a))

# Convert the list to a data frame
df <- data.frame(
  a = a_values,
  Non_Informative = sapply(variance_results, function(x) x$V1),
  Informative = sapply(variance_results, function(x) x$V2)
)
# Convert data to long format for ggplot
df_long <- melt(df, id.vars = "a", variable.name = "Mode", value.name = "Value")

# Plot density curves with thicker lines
ggplot(df_long, aes(x = a, y = Value, colour = Mode)) +
  geom_density(stat = "identity", size = 1.2) +  # Increase line thickness
  labs(
    title = "Variance of Posterior Selection Probability",
    x = expression(alpha[j]),
    y = "Variance",
    colour = "Mode"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),  # Center title, increase font size
    axis.title = element_text(size = 28),  # Increase font size for axis labels
    axis.text = element_text(size = 28),  # Increase font size for axis text
    legend.title = element_text(size = 26),  # Increase legend title font size
    legend.text = element_text(size = 24)  # Increase legend text font size
  )

# Figure 1c
# Fix B and n
B <- 100
n <- 70

# Create a sequence of a values from 0 to 1000
a_values <- seq(0, 100, by = 1)

# Compute V1 and V2 for each a value
variance_results <- lapply(a_values, function(a) variance_values(B, n, a))

# Convert the list to a data frame
df <- data.frame(
  a = a_values,
  Non_Informative = sapply(variance_results, function(x) x$V1),
  Informative = sapply(variance_results, function(x) x$V2)
)
# Convert data to long format for ggplot
df_long <- melt(df, id.vars = "a", variable.name = "Mode", value.name = "Value")

# Plot density curves with thicker lines
ggplot(df_long, aes(x = a, y = Value, colour = Mode)) +
  geom_density(stat = "identity", size = 1.2) +  # Increase line thickness
  labs(
    title = "Variance of Posterior Selection Probability",
    x = expression(alpha[j]),
    y = "Variance",
    colour = "Mode"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),  # Center title, increase font size
    axis.title = element_text(size = 28),  # Increase font size for axis labels
    axis.text = element_text(size = 28),  # Increase font size for axis text
    legend.title = element_text(size = 26),  # Increase legend title font size
    legend.text = element_text(size = 24)  # Increase legend text font size
  )

# Figure 2a
# Fix B and n
B <- 100
n <- 10

# Create a sequence of a values from 0 to 1000
a_values <- seq(0, 100, by = 1)

# Compute V1 and V2 for each a value
variance_results <- lapply(a_values, function(a) variance_values(B, n, a))

# Convert the list to a data frame
df <- data.frame(
  a = a_values,
  Non_Informative = sapply(variance_results, function(x) x$V1),
  Informative = sapply(variance_results, function(x) x$V2)
)
# Convert data to long format for ggplot
df_long <- melt(df, id.vars = "a", variable.name = "Mode", value.name = "Value")

# Plot density curves with thicker lines
ggplot(df_long, aes(x = a, y = Value, colour = Mode)) +
  geom_density(stat = "identity", size = 1.2) +  # Increase line thickness
  labs(
    title = "Variance of Posterior Selection Probability",
    x = expression(alpha[j]),
    y = "Variance",
    colour = "Mode"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),  # Center title, increase font size
    axis.title = element_text(size = 28),  # Increase font size for axis labels
    axis.text = element_text(size = 28),  # Increase font size for axis text
    legend.title = element_text(size = 26),  # Increase legend title font size
    legend.text = element_text(size = 24)  # Increase legend text font size
  )


# Figure 2b
# Fix B and n
B <- 100
n <- 90

# Create a sequence of a values from 0 to 1000
a_values <- seq(0, 100, by = 1)

# Compute V1 and V2 for each a value
variance_results <- lapply(a_values, function(a) variance_values(B, n, a))

# Convert the list to a data frame
df <- data.frame(
  a = a_values,
  Non_Informative = sapply(variance_results, function(x) x$V1),
  Informative = sapply(variance_results, function(x) x$V2)
)
# Convert data to long format for ggplot
df_long <- melt(df, id.vars = "a", variable.name = "Mode", value.name = "Value")

# Plot density curves with thicker lines
ggplot(df_long, aes(x = a, y = Value, colour = Mode)) +
  geom_density(stat = "identity", size = 1.2) +  # Increase line thickness
  labs(
    title = "Variance of Posterior Selection Probability",
    x = expression(alpha[j]),
    y = "Variance",
    colour = "Mode"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),  # Center title, increase font size
    axis.title = element_text(size = 28),  # Increase font size for axis labels
    axis.text = element_text(size = 28),  # Increase font size for axis text
    legend.title = element_text(size = 26),  # Increase legend title font size
    legend.text = element_text(size = 24)  # Increase legend text font size
  )


# Figure 3a , 3b
# Define the function to return both V1 and V2

# Fix B
B <- 100

# Create a grid of a and n values from 200 to 800
a_values <- seq(0, 100, by = 1)
n_values <- seq(0, 100, by = 1)

# Expand the grid
grid <- expand.grid(a = a_values, n = n_values)

# Calculate V1 and V2 for each combination of a and n
variance_results <- mapply(variance_values, grid$a, grid$n, MoreArgs = list(B = B), SIMPLIFY = FALSE)

# Extract V1, V2, and compute the difference
grid$V1 <- sapply(variance_results, function(x) x$V1)
grid$V2 <- sapply(variance_results, function(x) x$V2)

# Reshape data into matrices for surface plotting
matrix_V1 <- matrix(grid$V1, nrow = length(a_values), ncol = length(n_values), byrow = TRUE)
matrix_V2 <- matrix(grid$V2, nrow = length(a_values), ncol = length(n_values), byrow = TRUE)

# Define axis labels
x_label <- list(title = "n")  # Now n is on the x-axis
y_label <- list(title = "alpha")  # alpha is on the y-axis

# Create 3D surface plots
plot1 <- plot_ly(x = n_values, y = a_values, z = matrix_V1, type = 'surface') %>%
  layout(title = "Surface Plot of V1",
         scene = list(xaxis = x_label,
                      yaxis = y_label,
                      zaxis = list(title = "Variance")))

plot2 <- plot_ly(x = n_values, y = a_values, z = matrix_V2, type = 'surface') %>%
  layout(title = "Surface Plot of V2",
         scene = list(xaxis = x_label,
                      yaxis = y_label,
                      zaxis = list(title = "Variance")))

# Display plots
# Figure 3a
plot1
# Figure 3b
plot2
rm(list=ls()); gc() # Clean up




