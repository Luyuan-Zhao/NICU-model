
# Load necessary packages
library(ggplot2)
library(dplyr)
library(FSA)
library(patchwork)
library(readr)

# Read my data
df <- read_csv("D:/Github Files/NICU-model/all_scenarios_combined.csv")

# Order scenarios
df$scenario <- factor(df$scenario, levels = c("Baseline", LETTERS[1:8]))

# Select indicators we need 
metrics <- c("cumulative_I_R", "selection_events", "transmission_R", "peak_infection")

# Initialize the plotting lists
plot_list <- list()

# Loop 
for (metric in metrics) {
  cat("\n==============================\n")
  cat(paste("Analyzing:", metric, "\n"))
  
  # 1. Normality test
  cat("\n[1] Shapiro-Wilk Normality test:\n")
  shapiro_results <- by(df[[metric]], df$scenario, shapiro.test)
  print(shapiro_results)
  
  # 2. Calculate the mean value ± 95% CI
  summary_df <- df %>%
    group_by(scenario) %>%
    summarise(
      mean_val = mean(get(metric), na.rm = TRUE),
      sd_val = sd(get(metric), na.rm = TRUE),
      n = n(),
      ci_lower = mean_val - 1.96 * sd_val / sqrt(n),
      ci_upper = mean_val + 1.96 * sd_val / sqrt(n),
      .groups = "drop"
    )
  
  # 3. Plotting
  p <- ggplot(summary_df, aes(x = scenario, y = mean_val, fill = scenario)) +
    geom_col(alpha = 0.8) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "black") +
    labs(title = paste("Mean ± 95% CI of", metric, "across Scenarios"),
         x = "Scenario", y = metric) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Save plots in the list
  plot_list[[metric]] <- p
  
  # 4. Kruskal-Wallis test
  cat("\n[2] Kruskal-Wallis test:\n")
  kw <- kruskal.test(df[[metric]] ~ df$scenario)
  print(kw)
  
  # 5. Dunn tests（Bonferroni validation）
  cat("\n[3] Dunn test (Bonferroni correction):\n")
  dunn_res <- dunnTest(df[[metric]] ~ df$scenario, method = "bonferroni")
  print(dunn_res)
}

# Conjunct and gather the figures in a （2x2） matrix
combined_plot <- (plot_list[["cumulative_I_R"]] | plot_list[["selection_events"]]) /
  (plot_list[["transmission_R"]] | plot_list[["peak_infection"]])

# Visualisation
print(combined_plot)


write.csv(summary_df, "ScenarioSummaryTable.csv", row.names = FALSE)

