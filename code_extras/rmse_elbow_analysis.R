# rmse analysis to find optimal sample size - designing shearing experiments

# this script estimates the most efficient sample size for our wet lab work
# by calculating RMSE for log-log models across increasing sample sizes.
# 
# - fits log-log linear models to subsets of each dataset
# - finds the smallest n where improvement steps drop below 10%
# - generates a ggplot summary with cutoff shown at n=optimal_n_all

# bpfeiffer@covaris.com, 8/4/25

library(ggplot2)
library(dplyr)
library(readr)

# config 
base_path <- "C:/Users/bpfeiffer/Desktop/programming_projects/Protocol_Builder_RShiny/"
# kris' R230 data on different volumes in the same consumable
datasets <- c("elbow_plot_10", "elbow_plot_20", "elbow_plot_40", "elbow_plot_80")

# function to calculate rmse for different subset sizes
calculate_rmse_curve <- function(file_path, dataset_name) {
  # read and prepare data
  df <- read_csv(file_path, col_names = c("x", "y"), skip = 1)
  df <- df[order(-df$x), ]
  
  # inner function to fit model and calc rmse for a given subset size
  fit_subset_model <- function(df, size) {
    subset_df <- df[1:size, ]
    model <- lm(log(y) ~ log(x), data = subset_df)
    preds <- predict(model, newdata = df) 
    rmse <- sqrt(mean((log(df$y) - preds)^2))  # rmse in log space
    return(data.frame(n = size, rmse = rmse, dataset = dataset_name))
  }
  
  # calculate rmse for different subset sizes
  results <- lapply(3:nrow(df), function(n) fit_subset_model(df, n))
  return(bind_rows(results))
}

# process all datasets
all_results <- data.frame()
for (dataset in datasets) {
  file_path <- paste0(base_path, dataset, ".csv")
  results <- calculate_rmse_curve(file_path, dataset)
  all_results <- rbind(all_results, results)
}

# calculate percent improvement for each step
# we're using a bootstrap elbow plot type system, as our curves are too smooth to
# use the standard methods to find an elbow
improvement_analysis <- all_results %>%
  group_by(dataset) %>%
  arrange(n) %>%
  mutate(pct_improvement = (lag(rmse) - rmse) / lag(rmse) * 100) %>%
  filter(!is.na(pct_improvement))

# find where each condition reaches <10% improvement
individual_optimal <- improvement_analysis %>%
  group_by(dataset) %>%
  filter(pct_improvement < 10) %>%
  slice(1) %>%  # take the first one that meets criteria
  select(dataset, optimal_n = n)

# take the max to be conservative
optimal_n_all <- max(individual_optimal$optimal_n)

# make the plot
# using manual colors because the default ones are kind of ugly
plot <- ggplot(all_results, aes(x = n, y = rmse, color = dataset)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, alpha = 0.8) +  # alpha makes overlapping points clearer
  
  # custom colors - used the covaris colors
  scale_color_manual(
    values = c("elbow_plot_10" = "#cc1543", "elbow_plot_20" = "#00aae7", 
               "elbow_plot_40" = "#00827a", "elbow_plot_80" = "#800f69"),
    labels = c("Condition 10", "Condition 20", "Condition 40", "Condition 80")
  ) +
  
  # add padding to y axis so points don't touch edges
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  
  # labels - trying to be descriptive but not too wordy
  labs(
    title = "Model Prediction Error vs. Sample Size",
    subtitle = "Finding the optimal number of data points for reliable predictions",
    x = "Number of Data Points Used",
    y = "Root Mean Square Error (RMSE)",
    color = "Experimental\nCondition" 
  ) +
  
  # theme tweaks
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    axis.title = element_text(face = "bold")
  )

# add line showing cutoff
plot <- plot + 
  geom_vline(xintercept = optimal_n_all, linetype = "dashed", 
             color = "black", linewidth = 1, alpha = 0.7) +
  annotate("text", x = optimal_n_all + 0.5, y = max(all_results$rmse) * 0.85, 
           label = paste0("Optimal Range: n = ", optimal_n_all, "\n(<10% improvement)"), 
           color = "black", size = 3.5, hjust = 0, fontface = "bold")

print(plot)

# summary output for console
cat("\n=== rmse analysis summary ===\n")
cat("Recommended sample size for reliable predictions: n =", optimal_n_all, "\n")
cat("(based on <10% rmse improvement threshold)\n\n")

cat("individual dataset analysis:\n")
individual_detail <- improvement_analysis %>%
  group_by(dataset) %>%
  filter(pct_improvement < 10) %>%
  slice(1) %>%
  select(dataset, n, rmse, pct_improvement)
print(individual_detail)
