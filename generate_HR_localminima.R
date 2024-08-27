# This is the final version of generating HR (July 17th 2024)

generate_HR <- function(L1, L2) {
  rounds <- nrow(L1)
  hr_mv_der <- numeric(rounds)
  hr_mv <- numeric(rounds)
  nh_var <- numeric(rounds)
  hr_lpm <- numeric(rounds)
  x1_lpm <- numeric(rounds)
  x2_lpm <- numeric(rounds)
  mean_value <- numeric(rounds)
  mean_value2 <- numeric(rounds)
  nh_var1 <- numeric(rounds)
  nh_var2 <- numeric(rounds)
  skewness_L1 <- numeric(rounds)
  skewness_L2 <- numeric(rounds)
  
  for (i in 1:rounds) {
    p <- L1[i, 3]
    x_h <- L1[i, 2]
    x_l <- L1[i, 4]
    y_h <- L2[i, 4]
    y_l <- L2[i, 2]
    
    # Calculate mean values
    mean_value[i] <- (1 - p) * x_h + p * x_l
    mean_value2[i] <- (1 - p) * y_l + p * y_h
    
    # Calculate variances 
    nh_var1[i] <- (1 - p) * (x_h - mean_value[i])^2 + p * (x_l - mean_value[i])^2
    nh_var2[i] <- p * (y_h - mean_value2[i])^2 + (1 - p) * (y_l - mean_value2[i])^2
    
    # Calculate skewness for L1 and L2
    skewness_L1[i] <- ((1 - p) * (x_h - mean_value[i])^3 + p * (x_l - mean_value[i])^3) / sqrt(nh_var1[i])^3
    skewness_L2[i] <- (p * (y_h - mean_value2[i])^3 + (1 - p) * (y_l - mean_value2[i])^3) / sqrt(nh_var2[i])^3
    
    
    result <- optim(par = 0, fn = varianceObjective, 
                    p = p, x_h = x_h, x_l = x_l, y_h = y_h, y_l = y_l,
                    method = "BFGS")
    hr_mv[i] <- result$par
    
    expected_value <- (1 - p) * (x_h + hr_mv[i] * y_l) + p * (x_l + hr_mv[i] * y_h)
    nh_var[i] <- (1 - p) * ((x_h + hr_mv[i] * y_l) - expected_value)^2 + p * ((x_l + hr_mv[i] * y_h) - expected_value)^2
    
    hr_mv_der[i] <- sqrt(nh_var1[i]) / sqrt(nh_var2[i])
    
   
    # LPM2 minimization
    # Create a sequence of h values
    h_values <- seq(0, 20, by = 0.05)  
    lpm_values <- sapply(h_values, function(h) lpmObjective(h, p, x_h, x_l, y_h, y_l))
    
    # Find the minimum LPM2 value
    min_lpm <- min(lpm_values)
    
    if (min_lpm < 1e-10) {  # Consider LPM2 is zero
      # Find all h values that give LPM2 = 0
      zero_lpm_indices <- which(lpm_values < 1e-10)
      
      # Choose the smallest h among these
      hr_lpm[i] <- h_values[min(zero_lpm_indices)]
    } else {
      # If min_lpm is not 0, solve the minimization problem
      min_lpm_index <- which.min(lpm_values)
      hr_lpm[i] <- h_values[min_lpm_index]
    }
    
    x1_lpm[i] <- (1 - p) * max(0 - x_h, 0)^2 + p * max(0 - x_l, 0)^2
    x2_lpm[i] <- p * max(0 - y_h, 0)^2 + (1 - p) * max(0 - y_l, 0)^2
  }
  
  # Return the calculated hedge ratios
  return(list(hr_mv_der = hr_mv_der,
              hr_mv = hr_mv,
              hr_lpm = hr_lpm,
              nh_var = nh_var,
              nh_var1 = nh_var1,
              nh_var2 = nh_var2,
              x1_lpm = x1_lpm,
              x2_lpm = x2_lpm,
              mean_value = mean_value,
              mean_value2 = mean_value2,
              skewness_L1 = skewness_L1,
              skewness_L2 = skewness_L2))
}

varianceObjective <- function(h, p, x_h, x_l, y_h, y_l) {
  expected_value <- (1 - p) * (x_h + h * y_l) + p * (x_l + h * y_h)
  variance <- (1 - p) * ((x_h + h * y_l) - expected_value)^2 + p * ((x_l + h * y_h) - expected_value)^2
  return(variance)
}

lpmObjective <- function(h, p, x_h, x_l, y_h, y_l) {
  xbar <- 0
  lpm <- (1 - p) * max(xbar - (x_h + h * y_l), 0)^2 + p * max(xbar - (x_l + h * y_h), 0)^2
  return(lpm)
}

scenarios <- list(
  S1 = list(
    p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    x1_l_values = c(-6.0, -4.0, -3.05, -2.449, -2.0, -1.633, -1.309, -1.0, -0.667),
    x1_h_values = c(0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.1, 4.0, 6.0),
    x2_h_values = c(0.39, 0.7, 1.06, 1.49, 2.0, 2.59, 3.3, 4.3, 6.28),
    x2_l_values = c(-6.28, -4.3, -3.3, -2.59, -2.0, -1.49, -1.06, -0.7, -0.39)
  ),
  S2 = list(
    p_values = rep(0.5, 9),
    x1_l_values = -8:0,
    x1_h_values = 1:9,
    x2_h_values = 1:9,
    x2_l_values = -8:0
  ),
  S3 = list(
    p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    x1_l_values = c(-6.0, -4.0, -3.05, -2.449, -2.0, -1.633, -1.309, -1.0, -0.667),
    x1_h_values = c(0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.1, 4.0, 6.0),
    x2_h_values = c(6.4, 4, 3.1, 2.4, 2.0, 1.6, 1.3, 1.0, 0.7),
    x2_l_values = c(-0.667, -1.0, -1.309, -1.633, -2.0, -2.449, -3.055, -4.0, -6.0)
  ),
  
  S4 = list(
    p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    x1_l_values =c(-6.0, -4.0, -3.05, -2.449, -2.0, -1.633, -1.309, -1.0, -0.667),
    x1_h_values = c(0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.1, 4.0, 6.0),
    x2_h_values = rep(0.4, 9),
    x2_l_values = c(-6.3, -4.6, -4.0, -3.7, -3.6, -3.7, -4.0, -4.6, -6.3)
  ),
  S5 = list(
    p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    x1_l_values = c(-6.0, -4.0, -3.05, -2.449, -2.0, -1.633, -1.309, -1.0, -0.667),
    x1_h_values = c(0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.1, 4.0, 6.0),
    x2_h_values = c(6.3, 4.6, 4.0, 3.7, 3.6, 3.7, 4.0, 4.6, 6.3),
    x2_l_values = rep(-0.4, 9)
  )
)

#   S4 = list(
#     p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#     x1_l_values = rep(-0.4, 9),
#     x1_h_values = c(6.3, 4.6, 4.0, 3.7, 3.6, 3.7, 4.0, 4.6, 6.3),
#     x2_h_values = c(6.3, 4.6, 4.0, 3.7, 3.6, 3.7, 4.0, 4.6, 6.3),
#     x2_l_values = rep(-0.4, 9)
#   ),
#   S5 = list(
#     p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#     x1_l_values = c(-6.0, -4.0, -3.05, -2.449, -2.0, -1.633, -1.309, -1.0, -0.667),
#     x1_h_values = c(0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.1, 4.0, 6.0),
#     x2_h_values = c(6.0, 4.0, 3.06, 2.45, 2.0, 1.63, 1.31, 1.0, 0.67),
#     x2_l_values = c(-0.667, -1.0, -1.309, -1.633, -2.0, -2.449, -3.055, -4.0, -6.0)
#   ),
#     S8 = list(
#     p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#     x1_l_values = c(-6.0, -4.0, -3.055, -2.449, -2.0, -1.633, -1.309, -1.0, -0.667),
#     x1_h_values = c(0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.1, 4.0, 6.0),
#     x2_h_values = c(10.00, 7.40, 6.40, 5.95, 5.83, 5.95, 6.40, 7.40, 10.00),
#     x2_l_values = rep(-0.4, 9)
#   ),
#   S9 = list(
#     p_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#     x1_l_values = c(-6.0, -4.0, -3.05, -2.449, -2.0, -1.633, -1.309, -1.0, -0.667),
#     x1_h_values = c(0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.1, 4.0, 6.0),
#     x2_h_values = rep(0.4, 9),
#     x2_l_values = c(-10.00, -7.40, -6.40, -5.95, -5.83, -5.95, -6.40, -7.40, -10.00)
#   )
# )

# Function to create L1 and L2 for a given scenario
create_L1_L2 <- function(scenario) {
  L1 <- cbind(1 - scenario$p_values, scenario$x1_h_values, scenario$p_values, scenario$x1_l_values)
  L2 <- cbind(1-scenario$p_values, scenario$x2_l_values, scenario$p_values, scenario$x2_h_values)
  list(L1 = L1, L2 = L2)
}

# Create L1 and L2 for all scenarios
all_L1_L2 <- lapply(scenarios, create_L1_L2)


generate_HR(all_L1_L2$S1$L1,all_L1_L2$S1$L2)
generate_HR(all_L1_L2$S2$L1,all_L1_L2$S2$L2)
generate_HR(all_L1_L2$S3$L1,all_L1_L2$S3$L2)
generate_HR(all_L1_L2$S4$L1,all_L1_L2$S4$L2)
generate_HR(all_L1_L2$S5$L1,all_L1_L2$S5$L2)
#generate_HR(all_L1_L2$S6$L1,all_L1_L2$S6$L2)
#generate_HR(all_L1_L2$S7$L1,all_L1_L2$S7$L2)
#generate_HR(all_L1_L2$S8$L1,all_L1_L2$S8$L2)
#generate_HR(all_L1_L2$S9$L1,all_L1_L2$S9$L2)


