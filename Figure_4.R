# This scripts reproduces Figure 4 of the paper
# For Figure 4(a) set n = 20 in line 15
# For Figure 4(b) set n = 50
# For Figure 4(c) set n = 80

rm(list=ls()); gc() # Clean up
options(warn=-1) # Don't display warnings
if(!require("plotly")){install.packages("plotly")}

library(plotly)

# Parameters
B <- 100
threshold <- 0.6
n <- 80

# Axes
alpha_vals <- seq(0, B, length.out = 100)
gamma_vals <- seq(0, B, length.out = 100)

# Full grid
ALPHA_GRID <- matrix(rep(alpha_vals, each = length(gamma_vals)), nrow = length(gamma_vals))
GAMMA_GRID <- matrix(rep(gamma_vals, times = length(alpha_vals)), nrow = length(gamma_vals))

# Fraction matrix: only for gamma >= alpha
FRAC <- (ALPHA_GRID + n) / (GAMMA_GRID + B)
FRAC[ALPHA_GRID > GAMMA_GRID] <- NA  # Fraction undefined below diagonal

# Threshold plane 
THRESHOLD <- matrix(threshold, nrow = nrow(FRAC), ncol = ncol(FRAC))

# Fraction surface
p <- plot_ly() |>
  add_surface(x = ~ALPHA_GRID, y = ~GAMMA_GRID, z = ~FRAC,
              colorscale = "Viridis",
              showscale = TRUE,
              name = "",            
              colorbar = list(title = ""))

# Red threshold plane
p <- p |>
  add_surface(x = ~ALPHA_GRID, y = ~GAMMA_GRID, z = ~THRESHOLD,
              surfacecolor = matrix(1, nrow = nrow(THRESHOLD), ncol = ncol(THRESHOLD)),
              colorscale = list(c(0,1), c("red","red")),
              showscale = FALSE,
              opacity = 0.5,
              showlegend = FALSE)

# Layout
p <- p |> layout(
  scene = list(
    xaxis = list(title = "alpha"),
    yaxis = list(title = "gamma"),
    zaxis = list(title = "Posterior Mean")
  ),
  showlegend = FALSE
)

p





