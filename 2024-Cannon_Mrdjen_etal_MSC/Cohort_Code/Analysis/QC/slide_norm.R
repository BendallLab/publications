install.packages("tidyverse", "dplyr")
library(tidyverse)
library(dplyr)

# Reading the CSV file
# original, before ctrlNorm
data <- read_csv("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/QC_attempts/qc_tma_metrics_original/combined_csv_files/combined_nonzero_mean_stats.csv")
# total intensity trlNorm corrected
data <- read_csv("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/qc_tma_metrics/combined_csv_files/combined_total_intensity_stats.csv")
# ctrlNorm corrected
data <- read_csv("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/qc_tma_metrics_ctrlNorm/combined_data_nzm.csv")
# tma7 corrected
data <- read_csv("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/qc_tma_metrics_norm_TMA7/combined_data_nzm.csv")
# slide median corrected
data <- read_csv("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/qc_tma_metrics_hc_median_norm/combined_data_nzm.csv")
# slide median 1M corrected
data <- read_csv("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/qc_tma_metrics_hc_median_norm_1M/combined_data_nzm.csv")


################## First looks ################## 
# Processing the data
processed_data <- data %>%
  filter(TMA != 0 & Region == "HC-CA1") %>%
  #filter(TMA != 0 & Region != "AD-CA1" & Region != "HC-CA1") %>%
  group_by(channel) %>%
  mutate(z_score = scale(`Non-zero mean intensity`)) %>%
  ungroup() %>%
  mutate(TMA_Region_Core = paste(TMA, Region, Core, sep = "_")) %>%
  group_by(TMA_Region_Core, channel) %>%
  summarize(median_z_score = median(z_score)) %>%
  #summarize(avg_z_score = mean(z_score, na.rm = TRUE)) %>%
  #mutate(exp_avg_z_score = 2**avg_z_score) %>%
  ungroup()

# Plotting the heatmap
library(ggplot2)
ggplot(processed_data, aes(x = TMA_Region_Core, y = channel, fill = median_z_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Median Z-score")

# Plotting the boxplot
library(ggplot2)
ggplot(processed_data, aes(x = TMA_Region_Core, y = z_score)) +
  geom_boxplot() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  ylim(-3,3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Median Z-score")

# Save the plot if needed
plot_save_path <-"/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/qc_tma_metrics/combined_csv_files/"
ggsave("non-zero_mean_heatmap_plot_Cores_AVG.png", width = 10, height = 10, dpi=300, bg="white")
# Save the processed data csv
write_csv(processed_data, "data/control_percents.csv")




################## percent norm filtering ################## 

# Simplified processing for normalization.
processed_data <- data %>%
  filter(TMA != 0 & Region == "HC-CA1") %>%
  group_by(TMA) %>%
  summarize(median_nzm_expr = median(`Non-zero mean intensity`, na.rm = TRUE)) %>%
  # mutate(z_score = scale(median_nzm_expr)) %>%
  # mutate(percentile_z = pnorm(z_score)) %>%
  # Use mutate to create the new column Y with 1 - the percent difference
  mutate(normalize_factor = 1 - (median_nzm_expr - min(median_nzm_expr, na.rm = TRUE)) / min(median_nzm_expr, na.rm = TRUE)) %>%
  ungroup()

processed_data$TMA = as.factor(processed_data$TMA)
# Save the processed data csv
write_csv(processed_data, "data/control_percents.csv")



################## Slide nzmexpr median norm filtering ################## 

# Processing the data for Slide Mean norm
processed_data <- data %>%
  #filter(TMA != 0) %>%
  filter(TMA != 0 & Region == "HC-CA1") %>%
  group_by(TMA) %>%
  mutate(nzm_median = 1000**median(`Non-zero mean intensity`)) %>%
  mutate(new_nzm = `Non-zero mean intensity`/nzm_median) %>%
  ungroup() %>%
  group_by(TMA) %>%
  summarize(x = median(nzm_median))
  # below was used for other iterations.
  group_by(TMA, channel) %>%
  mutate(z_score = scale(new_nzm)) %>%
  #mutate(adj_z_score = 2**scale(nzm_median)) %>%
  ungroup() %>%
  #
  mutate(TMA_Region_Core = paste(TMA, Region, Core, sep = "_")) %>%
  group_by(channel, TMA) %>%
  summarize(new_median_z_score = median(z_score)) %>%
  #summarize(avg_z_score = mean(z_score, na.rm = TRUE)) %>%
  #mutate(exp_avg_z_score = 2**avg_z_score) %>%
  ungroup()

# checking slide mean norm
processed_data2 <- data %>%
  filter(TMA != 0 & Region != "HC-CA1" & Region != "AD-CA1") %>%
  left_join(processed_data, by = "TMA") %>%
  mutate(new_nzm = `Non-zero mean intensity`/ x) %>%
  ungroup() %>%
  group_by(TMA, channel) %>%
  mutate(z_score = scale(new_nzm)) %>%
  #mutate(adj_z_score = 2**scale(nzm_median)) %>%
  ungroup() %>%
  #
  mutate(TMA_Region_Core = paste(TMA, Region, Core, sep = "_")) %>%
  group_by(channel, TMA_Region_Core) %>%
  summarize(new_median_z_score = median(z_score)) %>%
  #summarize(avg_z_score = mean(z_score, na.rm = TRUE)) %>%
  #mutate(exp_avg_z_score = 2**avg_z_score) %>%
  ungroup()
  


# Plotting the heatmap
library(ggplot2)
ggplot(processed_data2, aes(x = TMA_Region_Core, y = channel, fill = new_median_z_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "New Median Z-score")

# Plotting the boxplot
library(ggplot2)
ggplot(processed_data2, aes(x = TMA_Region_Core, y = z_score)) +
  geom_boxplot() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Median Z-score")

processed_data$TMA = as.factor(processed_data$TMA)
# Save the processed data csv
write_csv(processed_data, "data/control_medians_1M.csv")

