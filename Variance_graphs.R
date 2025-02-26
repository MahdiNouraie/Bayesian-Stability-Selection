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




