# Results 
# Luyuan Zhao, 14 July
# This is the R scripts for statistical analysis and visualisation in the Results part.

# ---------------1. Bar pots---------------
# Load necessary packages here
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

# Prepare names for all strategies (A-H, 8 in total)
strategy_names <- c(
  "None" = "None",
  "A" = "A",
  "B" = "B",
  "C" = "C",
  "D" = "D",
  "E" = "E",
  "F" = "F",
  "G" = "G",
  "H" = "H"
)

# Upload the dataset and reshape it for better fitting
df <- read_csv("all_scenarios_full.csv")

df_long <- df %>%
  pivot_longer(cols = c(peak_infection, selection_events,
                        transmission_R, cumulative_I_R),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(
    strategy_label = factor(strategy_names[as.character(strategy)],
                            levels = strategy_names),
    baseline = factor(baseline, levels = c("Low", "Mid", "High")),
    metric = recode(metric,
                    peak_infection   = "Peak infections",
                    selection_events = "Number of Selection Events",
                    transmission_R   = "Number of Resistance Transmission Events",
                    cumulative_I_R   = "Number of Cumulative Resistant Infections")
  )

# Color palette 
pal            <- brewer.pal(12, "Paired")
baseline_cols  <- setNames(pal[c(1, 3, 5)], c("Low", "Mid", "High"))

#  Helper function to build a single subplot
build_barplot <- function(metric_name, show_x_axis = TRUE) {
  
  df_sum <- df_long %>%
    filter(metric == metric_name) %>%
    group_by(strategy_label, baseline) %>%
    summarise(
      mean_val = mean(value),
      sd_val   = sd(value),
      n        = n(),
      se       = sd_val / sqrt(n),
      ci       = qt(0.975, df = n - 1) * se,
      .groups  = "drop"
    )
  
  p <- ggplot(df_sum,
              aes(x = strategy_label, y = mean_val, fill = baseline)) +
    geom_bar(stat = "identity",
             position = position_dodge(0.8),
             width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(ymin = mean_val - ci, ymax = mean_val + ci),
                  position = position_dodge(0.8),
                  width = 0.35, size = 0.6, colour = "black") +
    scale_fill_manual(values = baseline_cols) +
    labs(
      title = paste0(metric_name, " \n(Mean ± 95% CI)"),
      x = if (show_x_axis) "Intervention Strategy" else NULL,
      y = metric_name,
      fill = "Baseline"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title    = element_text(face = "bold", size = 16),
      axis.text.x   = element_text(face = "bold", angle = 0, hjust = 1, size = 14),
      axis.title.y  = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  if (!show_x_axis) {
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x  = element_blank())
  }
  return(p)
}

# Build four plots (top row without x-axis and keep the bottom one)
metric_order <- c("Peak infections",
                  "Number of Selection Events",
                  "Number of Resistance Transmission Events",
                  "Number of Cumulative Resistant Infections")

p1 <- build_barplot(metric_order[1], show_x_axis = FALSE)  # top-left
p2 <- build_barplot(metric_order[2], show_x_axis = FALSE)  # top-right
p3 <- build_barplot(metric_order[3], show_x_axis = TRUE)   # bottom-left
p4 <- build_barplot(metric_order[4], show_x_axis = TRUE)   # bottom-right

# Stitch 2×2 panel & collect legend
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")  # apply ONLY legend rule globally

# Export (PNG + PDF) 
dir.create("figures", showWarnings = FALSE)

ggsave("figures/panel_sharedX_bottomOnly.png",
       combined_plot, width = 35, height = 28, units = "cm", dpi = 600)

ggsave("figures/panel_sharedX_bottomOnly.pdf",
       combined_plot, width = 35, height = 28, units = "cm")




# ---------------2. Processed Heatmap for improvement evaluate the effectiveness of strategies---------------
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
        # label x- and y-axis
        axis.title.x = element_text(size = 14, face = "bold"),  
        axis.title.y = element_text(size = 14, face = "bold"),
        
        # tick labels
        axis.text.x  = element_text(size = 12, face = "bold"), 
        axis.text.y  = element_text(size = 12, face = "bold"),
      
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 14)
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
ggsave("combined_outcomes_heatmap.png", final_plot, width = 12, height = 8)





# ---------------3. Synergy Assessment of NICU Intervention Strategies---------------
# 3.1 Datapreparation
library(dplyr)
library(readr)

# read the dataset
df <- read_csv("all_scenarios_full.csv")

# create space to store the data
metrics <- c("cumulative_I_R", "peak_infection", "selection_events", "transmission_R")

improvements <- df %>%
  group_by(baseline, strategy) %>%
  summarise(across(all_of(metrics), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
  group_by(baseline, metric) %>%
  mutate(
    baseline_value = value[strategy == "None"],
    improvement_pct = round(100 * (baseline_value - value) / baseline_value, 1)
  ) %>%
  filter(strategy != "None") %>%
  select(baseline, strategy, metric, improvement_pct)

# save as .csv
write_csv(improvements, "intervention_improvements.csv")


# 3.2 Visualisaion
# load necessary packages
library(tidyverse)
library(patchwork)

# define what bundles I want to study
combo_list <- list(
  "F vs A+B"      = list(target = "F", parts = c("A", "B")),
  "G vs C+D+E"    = list(target = "G", parts = c("C", "D", "E")),
  "H vs A+B+C+D+E"= list(target = "H", parts = c("A", "B", "C", "D", "E"))
)

# read the dataset ("improvements")
df <- read_csv("intervention_improvements.csv")

# extract relevant data for both independent and combination bundles
sum_parts_df <- purrr::imap_dfr(combo_list, ~{
  comp <- .y
  parts <- .x$parts
  df %>%
    filter(strategy %in% parts) %>%
    group_by(baseline, metric) %>%
    summarise(improvement_pct = sum(improvement_pct), .groups = "drop") %>%
    mutate(strategy = paste(parts, collapse = "+"),
           group = "sum_of_parts",
           comparison = comp)
})

combined_df <- purrr::imap_dfr(combo_list, ~{
  comp <- .y
  tgt <- .x$target
  df %>%
    filter(strategy == tgt) %>%
    mutate(group = "combined",
           comparison = comp)
})

plot_df <- bind_rows(sum_parts_df, combined_df)

# justify whether it's synergy, additive or antagonism
synergy_df <- plot_df %>%
  select(baseline, metric, comparison, group, improvement_pct) %>%
  pivot_wider(names_from = group, values_from = improvement_pct) %>%
  mutate(diff = combined - sum_of_parts,
         class = case_when(
           abs(diff) < 1e-6 ~ "additive",
           diff > 0 ~ "synergy",
           diff < 0 ~ "antagonism"
         ))

# to beautify a bit
fill_cols <- c("sum_of_parts" = "#E3822F",
               "combined"     = "#495A84")

# setting for labels
metric_labels <- c(
  cumulative_I_R   = "Cumulative\nresistant infections",
  peak_infection   = "Peak infection\ncount",
  selection_events = "Resistance selection\nevents",
  transmission_R   = "Transmission-induced\nresistance"
)

# plot for each key metric
make_metric_plot <- function(metric_code) {
  pretty_name <- metric_labels[[metric_code]]
  
  ggplot(plot_df %>% filter(metric == metric_code),
         aes(x = baseline, y = improvement_pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.65), width = 0.55) +
    geom_text(data = synergy_df %>% filter(metric == metric_code),
              aes(x = baseline, y = combined,
                  label = paste0("Δ = ", round(diff, 1), "%\n", class)),
              vjust = -0.6, size = 4.2, fontface = "bold",
              inherit.aes = FALSE)+
    facet_wrap(~ comparison, ncol = 3) +
    scale_y_continuous(name = pretty_name,
                       labels = scales::label_percent(scale = 1),
                       expand = expansion(mult = c(0.05, 0.30))) +
    scale_fill_manual(values = fill_cols,
                      breaks = names(fill_cols),
                      labels = c("Sum of parts", "Combined")) +
    labs(x = NULL, fill = NULL) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "top",
      strip.text      = element_text(size = 15, face = "bold"),         
      axis.title.y    = element_text(size = 15, face = "bold"),         
      axis.text       = element_text(size = 13, face = "bold"),         
      plot.margin     = margin(10, 10, 10, 10)
    )
}

# get the four individual plots
p1 <- make_metric_plot("cumulative_I_R")
p2 <- make_metric_plot("peak_infection")
p3 <- make_metric_plot("selection_events")
p4 <- make_metric_plot("transmission_R")

# shut down the legend and captions using patchwork
final_fig <- ((p1 + p2 + plot_layout(guides = "collect")) /
                (p3 + p4 + plot_layout(guides = "collect"))) &
  theme(legend.position = "top") &
  plot_annotation(
    title = "Synergy assessment of NICU intervention strategies",
    caption = "Δ = Combined − Sum of parts; class: synergy (>0), additive (=0), antagonism (<0)"
  )

# download the plots
ggsave("synergy_barplot_split_clean.png",
       plot = final_fig,
       width = 18, height = 12, dpi = 600)





