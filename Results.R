# Function to collect scenario results
collect_scenario_results <- function(L1, L2, scenario_num) {
  hr <- generate_HR(L1, L2)
  rounds <- nrow(L1)
  
  df <- data.frame(
    '1-p' = 1 - L1[, 3],
    L1_high = L1[,2],
    p = L1[, 3],
    L1_low = L1[, 4],
    L1_mean = hr$mean_value,
    L1_variance = hr$nh_var1,
    L1_skewness = hr$skewness_L1,
    p = L1[, 3],
    L2_high = L2[, 4],
    '1-p' = 1 - L1[, 3],
    L2_low = L2[, 2],
    L2_mean = hr$mean_value2,
    L2_variance = hr$nh_var2,
    L2_skewness = hr$skewness_L2,
    hr_mv = hr$hr_mv,
    hr_lpm = hr$hr_lpm
  )
  
  return(df)
}

# List of scenarios
scenarios <- list(
  list(L1 = all_L1_L2$S1$L1, L1 = all_L1_L2$S1$L1),
  list(L1 = all_L1_L2$S2$L1, L2 = all_L1_L2$S2$L2),
  list(L1 = all_L1_L2$S3$L1, L2 = all_L1_L2$S3$L2),
  list(L1 = all_L1_L2$S4$L1, L2 = all_L1_L2$S4$L2),
  list(L1 = all_L1_L2$S5$L1, L2 = all_L1_L2$S5$L2),
  list(L1 = all_L1_L2$S6$L1, L2 = all_L1_L2$S6$L2),
  list(L1 = all_L1_L2$S7$L1, L2 = all_L1_L2$S7$L2),
  list(L1 = all_L1_L2$S8$L1, L2 = all_L1_L2$S8$L2),
  list(L1 = all_L1_L2$S9$L1, L2 = all_L1_L2$S9$L2)
)

# Create a new workbook
wb <- createWorkbook()

# Loop through scenarios and add each to the workbook
for (i in 1:length(scenarios)) {
  scenario_results <- collect_scenario_results(scenarios[[i]]$L1, scenarios[[i]]$L2, i)
  addWorksheet(wb, paste("Scenario", i))
  writeData(wb, sheet = paste("Scenario", i), scenario_results)
}

# Save the workbook to an Excel file
saveWorkbook(wb, "HR_Scenarios3.xlsx", overwrite = TRUE)