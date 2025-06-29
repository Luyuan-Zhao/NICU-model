# Package loading

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# Files reading

files <- list(
  "p_empiric_abx" = "C:/Users/Willi/Desktop/Independent Project/0629 Methods+SimuResults/sensitivity_p_empiric_abx.csv",
  "abx_course" = "C:/Users/Willi/Desktop/Independent Project/0629 Methods+SimuResults/sensitivity_abx_course.csv",
  "p_decontam_hcw" = "C:/Users/Willi/Desktop/Independent Project/0629 Methods+SimuResults/sensitivity_p_decontam_hcw.csv",
  "p_decontam_dev" = "C:/Users/Willi/Desktop/Independent Project/0629 Methods+SimuResults/sensitivity_p_decontam_dev.csv",
  "p_clean_env" = "C:/Users/Willi/Desktop/Independent Project/0629 Methods+SimuResults/sensitivity_p_clean_env.csv"
)

all_data <- bind_rows(lapply(names(files), function(pname) {
  df <- read_csv(files[[pname]])
  df$parameter <- pname
  return(df)
}))

# Output selection

outputs <- c("cumulative_I_R", "selection_events", "transmission_R", "peak_infection")

# Main loop for Sensitivity analysis

for (p in unique(all_data$parameter)) {
  df_sub <- all_data %>% filter(parameter == p)
  plot_list <- list()
  
  cat("\n=========================================\n")
  cat(paste("Parameter analysis：", p, "\n"))
  
  for (metric in outputs) {
    # Spearman tests
    spearman_res <- cor.test(df_sub[[metric]], as.numeric(df_sub$param), method = "spearman")
    rho <- round(spearman_res$estimate, 3)
    pval <- signif(spearman_res$p.value, 3)
    
    cat(paste0("Spearman: ", metric, " ~ ", p,
               " | rho = ", rho, ", p = ", pval, "\n"))
    
    # Mean ± SD
    df_summary <- df_sub %>%
      group_by(param) %>%
      summarise(mean_val = mean(.data[[metric]], na.rm = TRUE),
                sd_val = sd(.data[[metric]], na.rm = TRUE),
                .groups = "drop")
    
    # plotting
    p_plot <- ggplot(df_summary, aes(x = as.factor(param), y = mean_val, group = 1)) +
      geom_line(color = "#1b9e77", linewidth = 1.2) +
      geom_point(color = "#d95f02", size = 3) +
      geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                    width = 0.2, color = "gray30") +
      labs(title = paste0(metric, "\nρ = ", rho, ", p = ", pval),
           x = paste0(p, " (tested values)"), y = "Mean ± SD") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_list[[metric]] <- p_plot
  }
  
  # 2x2 Matrix for visualisation
  print((plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]]))
}