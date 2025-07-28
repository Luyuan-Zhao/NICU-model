# Results for Sensitivity Analysis
# Luyuan Zhao, 14 July
# This is the R scripts for statistical analysis and visualisation in the Results part.
library(tidyverse)

# Load all the SA results here
files_info <- list(
  p_empiric_abx = "sensitivity_p_empiric_abx.csv",
  abx_course = "sensitivity_abx_course.csv",
  p_decontam_hcw = "sensitivity_p_decontam_hcw.csv",
  p_decontam_dev = "sensitivity_p_decontam_dev.csv",
  p_clean_env = "sensitivity_p_clean_env.csv"
)

# Initialize a blank dataframe
all_data <- tibble()

# Loop
for (param_name in names(files_info)) {
  df <- read_csv(files_info[[param_name]])
  df <- df %>%
    rename(param_value = param) %>%
    mutate(param_name = param_name)
  all_data <- bind_rows(all_data, df)
}

# Reorganize the table
long_df <- all_data %>%
  pivot_longer(cols = c("peak_infection", "selection_events", "transmission_R", "cumulative_I_R"),
               names_to = "metric", values_to = "value")

# Plotting for every metrics
ggplot(long_df, aes(x = param_value, y = value)) +
  geom_point(alpha = 0.4, color = "darkblue") +
  stat_summary(fun = mean, geom = "line", color = "black", size = 1) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.02) +
  facet_grid(metric ~ param_name, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(
    title = "One-at-a-Time Sensitivity Analysis",
    x = "Parameter Value",
    y = "Model Output"
  )

# Spearman analysis
cor_summary <- long_df %>%
  group_by(param_name, metric) %>%
  summarise(
    spearman_rho = cor(param_value, value, method = "spearman"),
    p_value = cor.test(param_value, value, method = "spearman")$p.value,
    .groups = "drop"
  )

# Print results
print(cor_summary)

# Store Statistical results
write_csv(cor_summary, "OAT_spearman_summary.csv")
