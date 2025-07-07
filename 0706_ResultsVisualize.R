
# Use 3 types of Visulisations here: boxplot+jitterplot, bar plot+ CI and heatmap

# 1. Boxplots with jitter
# loading Necessary Pkgs
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)

# Read the dataset
df <- read_csv("all_scenarios_full.csv")

# Reorganize the table
df_long <- df %>%
  pivot_longer(cols = c(peak_infection, selection_events, transmission_R, cumulative_I_R),
               names_to = "metric", values_to = "value") %>%
  mutate(
    strategy = factor(strategy, levels = c("None", "A", "B", "C", "D", "E", "F", "G", "H")),
    baseline = factor(baseline, levels = c("Low", "Mid", "High")),
    metric = recode(metric,
                    peak_infection = "Peak infections",
                    selection_events = "Number of Selection Events",
                    transmission_R = "Number of Resistance Transmission Events",
                    cumulative_I_R = "Number of Cummulative Resistant Infections")
  )

# Color choosing
paired_palette <- brewer.pal(n = 12, name = "Paired")
baseline_colors <- setNames(paired_palette[c(1, 3, 5)], c("Low", "Mid", "High"))  

# Output path
output_dir <- "figures/boxplots"
dir.create(output_dir, showWarnings = FALSE)

# Loop for every output metric
for (m in unique(df_long$metric)) {
  df_plot <- df_long %>% filter(metric == m)
  
  p <- ggplot(df_plot, aes(x = strategy, y = value, fill = baseline, color = baseline)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                shape = 21, size = 1.5, stroke = 0.2, alpha = 0.4) +
    scale_fill_manual(values = baseline_colors) +
    scale_color_manual(values = baseline_colors) +
    labs(
      title = paste(m, "by Strategy and Baseline"),
      x = "Intervention Strategy",
      y = m,
      fill = "Baseline",
      color = "Baseline"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  print(p)
}






# 2. Barplots with 95% CI
# Load Pkgs
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)

# Read the dataset
df <- read_csv("all_scenarios_full.csv")

# 
df_long <- df %>%
  pivot_longer(cols = c(peak_infection, selection_events, transmission_R, cumulative_I_R),
               names_to = "metric", values_to = "value") %>%
  mutate(
    strategy = factor(strategy, levels = c("None", "A", "B", "C", "D", "E", "F", "G", "H")),
    baseline = factor(baseline, levels = c("Low", "Mid", "High")),
    metric = recode(metric,
                    peak_infection = "Peak infections",
                    selection_events = "Number of Selection Events",
                    transmission_R = "Number of Resistance Transmission Events",
                    cumulative_I_R = "Number of Cummulative Resistant Infections")
  )

# Colours
paired_palette <- brewer.pal(n = 12, name = "Paired")
baseline_colors <- setNames(paired_palette[c(1, 3, 5)], c("Low", "Mid", "High"))

# Create Bar plots with Mean + 95% CI
for (m in unique(df_long$metric)) {
  df_plot <- df_long %>% filter(metric == m)
  
  # Calculate Mean Value and 95% CI
  df_summary <- df_plot %>%
    group_by(strategy, baseline) %>%
    summarise(
      mean_val = mean(value),
      sd_val = sd(value),
      n = n(),
      se = sd_val / sqrt(n),
      ci = qt(0.975, df = n - 1) * se,
      .groups = "drop"
    )
  
  # Plotting 
  p <- ggplot(df_summary, aes(x = strategy, y = mean_val, fill = baseline)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.65, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_val - ci, ymax = mean_val + ci),
                  position = position_dodge(0.8), width = 0.2) +
    scale_fill_manual(values = baseline_colors) +
    labs(
      title = paste(m, "(mean ± 95% CI) by Strategy and Baseline"),
      x = "Intervention Strategy",
      y = m,
      fill = "Baseline"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  print(p)
}





# 3.1  Heatmap for simple comparison
# Pkgs
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)
library(scales)

# Read the data
df <- read_csv("all_scenarios_full.csv")

# 
df_long <- df %>%
  pivot_longer(cols = c(peak_infection, selection_events, transmission_R, cumulative_I_R),
               names_to = "metric", values_to = "value") %>%
  mutate(
    strategy = factor(strategy, levels = c("None", "A", "B", "C", "D", "E", "F", "G", "H")),
    baseline = factor(baseline, levels = c("Low", "Mid", "High")),
    metric = recode(metric,
                    peak_infection = "Peak infections",
                    selection_events = "Number of Selection Events",
                    transmission_R = "Number of Resistance Transmission Events",
                    cumulative_I_R = "Number of Cummulative Resistant Infections")
  )

# Colours
heat_palette <- scale_fill_distiller(palette = "YlGnBu", direction = 1)

# Aggregate mean values
df_heat <- df_long %>%
  group_by(strategy, baseline, metric) %>%
  summarise(mean_value = mean(value), .groups = "drop")

# Plotting for every metrics
for (m in unique(df_heat$metric)) {
  df_plot <- df_heat %>% filter(metric == m)
  
  p <- ggplot(df_plot, aes(x = strategy, y = baseline, fill = mean_value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    heat_palette +
    geom_text(aes(label = round(mean_value, 1)), size = 4) +
    labs(
      title = paste(m, "— Mean Values Across Baseline × Strategy"),
      x = "Strategy",
      y = "Baseline",
      fill = m
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  print(p)
}




# 3.2  Processed Heatmap for ro evaluate the effectiveness of strategies
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)  
library(purrr)

# Main plotting function
plot_combined_improvements <- function(file_path = "all_scenarios_full.csv") {
  df <- read_csv(file_path)
  
  outcome_vars <- c(  "Peak infections (I_S + I_R)" = "peak_infection",
                      "Resistance selection events (C_S → C_R)" = "selection_events",
                      "Transmission-induced resistance (S → C_R)" = "transmission_R",
                      "Cumulative resistant infections (I_R)" = "cumulative_I_R")
  plot_list <- list()
  
  for (label in names(outcome_vars)) {
    var <- outcome_vars[[label]]
    
    df_summary <- df %>%
      group_by(baseline, strategy) %>%
      summarise(mean_val = mean(.data[[var]]), .groups = "drop")
    
    baseline_vals <- df_summary %>%
      filter(strategy == "None") %>%
      select(baseline, mean_val) %>%
      rename(baseline_val = mean_val)
    
    df_relative <- df_summary %>%
      left_join(baseline_vals, by = "baseline") %>%
      mutate(improvement = 100 * (baseline_val - mean_val) / baseline_val) %>%
      filter(strategy != "None")
    
    p <- ggplot(df_relative, aes(x = strategy, y = baseline, fill = improvement)) +
      geom_tile(color = "white") +
      geom_text(aes(label = paste0(round(improvement), "%")), color = "black", size = 3.5) +
      scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850", midpoint = 0,
                           name = "Improvement\n(%)") +
      labs(title = label, x = "Strategy", y = "Baseline") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 10)
      )
    
    plot_list[[var]] <- p
  }
  
  # Matrix (2X2)
  combined_plot <- (plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]])
  return(combined_plot)
}

# Call for Results
final_plot <- plot_combined_improvements("all_scenarios_full.csv")
print(final_plot)  
# ggsave("combined_outcomes_heatmap.png", final_plot, width = 12, height = 8)






