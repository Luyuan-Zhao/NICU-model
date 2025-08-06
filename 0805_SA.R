# Codes for One-at-a-Time (OAT) Sensitivity Analysis
# luyuan Zhao - 3020120Z
# 6th August
# This script performs a OAT sensitivity analysis using bootstrapped Spearman correlations 
# between model input parameters and key outcomes, and visualizes the results as tornado plots.
#------------------------------------------------

# Load necessary libraries for data manipulation and plotting
library(tidyverse)
library(ggplot2)

# Define the file paths for sensitivity analysis results for each parameter
files <- list(
  p_empiric_abx = "sensitivity_p_empiric_abx.csv",
  abx_course = "sensitivity_abx_course.csv",
  p_decontam_dev = "sensitivity_p_decontam_dev.csv",
  p_decontam_hcw = "sensitivity_p_decontam_hcw.csv",
  p_clean_env = "sensitivity_p_clean_env.csv"
)

# Define model outcomes of interest and the name of the parameter column in the CSV files
outcomes <- c("cumulative_I_R", "peak_infection", "selection_events", "transmission_R")
param_column <- "param" 

# Function to compute Spearman's correlation and estimate 95% confidence intervals via bootstrapping
calculate_ci <- function(data, param_name, outcome, n_samples = 2000, sample_size = 1000) {
  # Repeat sampling process n_samples times (default: 2000) for CI estimation
  spearman_list <- replicate(n_samples, {
    sample_df <- data %>% sample_n(min(sample_size, nrow(data)))
    
    x <- as.numeric(sample_df[[param_name]])
    y <- as.numeric(sample_df[[outcome]])
    
    # Remove NA values to ensure complete cases for correlation
    valid_idx <- complete.cases(x, y)
    x <- x[valid_idx]
    y <- y[valid_idx]
    
    if (length(x) >= 10 && length(x) == length(y)) {
      suppressWarnings(cor(x, y, method = "spearman"))
    } else {
      NA
    }
  })
  
  # Return median and 95% CI (2.5% and 97.5% quantiles) of Spearman rho values
  ci <- quantile(spearman_list, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  tibble(
    parameter = param_name,
    outcome = outcome,
    rho_median = ci[2],
    rho_lower = ci[1],
    rho_upper = ci[3]
  )
}

# Store the results of correlation analysis for each file and each outcome
all_results <- list()
# Loop through each CSV file and compute CI for all outcomes
for (param in names(files)) {
  df <- read_csv(files[[param]])
  
  for (outcome in outcomes) {
    result <- calculate_ci(df, param_column, outcome)
    result$parameter <- param  
    all_results[[length(all_results) + 1]] <- result
  }
}

# Combine all results into a single data frame
summary_df <- bind_rows(all_results)

# Print the summary results to console
print(summary_df)

# Prepare and visualize tornado plots of correlation results

plot_df <- summary_df %>%
  mutate(parameter = fct_reorder(parameter, abs(rho_median))) 

parameter_labels <- c(
  p_empiric_abx   = "Empiric Antibiotic Use",
  abx_course      = "Antibiotic Course Duration",
  p_decontam_dev  = "Device Cleaning",
  p_decontam_hcw  = "Hand Hygiene Compliance",
  p_clean_env     = "Environmental Cleaning"
)

outcome_labels <- c(
  cumulative_I_R    = "Cumulative Resistance Infections",
  peak_infection    = "Peak Infections",
  selection_events  = "Selection Events",
  transmission_R    = "Transmission Resistance Events"
)

plot_df <- plot_df %>%
  mutate(
    parameter = recode(as.character(parameter), !!!parameter_labels),
    outcome   = recode(as.character(outcome), !!!outcome_labels),
    direction = factor(
      ifelse(rho_median >= 0, "Positive Correlation", "Negative Correlation"),
      levels = c("Positive Correlation", "Negative Correlation")
    )
  ) %>%
  drop_na(parameter, outcome, rho_median)

ggplot(plot_df, aes(x = parameter, y = rho_median, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = rho_lower, ymax = rho_upper), width = 0.2, color = "black") +
  facet_wrap(~ outcome, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Positive Correlation" = "#3498db",  # 蓝色
      "Negative Correlation" = "#e74c3c"   # 红色
    ),
    name = "Correlation Direction"
  ) +
  labs(
    x = "Model Parameter",
    y = "Spearman’s ρ (95% CI)",
    title = "Tornado Plot of Sensitivity Analysis"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    strip.text  = element_text(size = 14, face = "bold"),
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )
# Prepare plot data: relabel variables and classify direction of correlation
export_df <- plot_df %>%
  select(parameter, outcome, rho_median, rho_lower, rho_upper, direction) %>%
  arrange(outcome, desc(abs(rho_median)))  

write_csv(export_df, "sensitivity_analysis_summary.csv")






