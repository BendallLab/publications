data_path = "/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data"
name_acta_brain_paper_1_cell_table = "ADD_Partial_HIRes_Hippocampus_cell_table.csv"
abp1_cell_table = data.table::fread(file=paste(data_path, name_acta_brain_paper_1_cell_table, sep = "/"))
# get rid of zeros
abp1_cell_table <- abp1_cell_table %>% mutate(Iba1_normalized = Iba1 + 0.0001)

## Violin plot for establishing mglia / neuron Iba1 baselinbe expression values - Reviewer plot
ggplot(abp1_cell_table) +
  aes(x = Meta, y = Iba1_normalized, fill = Meta) +
  geom_violin(adjust = 1L, scale = "area") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "white", position = position_dodge(width = 0.9)) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 5), size = 20), 
               vjust = -7.5, hjust = -0.2, color = "black", position = position_dodge(width = 0.9)) +
  scale_fill_hue(direction = 1) +
  scale_y_continuous(trans = "log10") +
  theme_minimal()

## for cohort density
density_mglia = mglia_to_save %>% select(c("label", "patient_label", "Disease_Status", "fov_row", "fov_col", "Region")) %>%
  mutate(point_id = paste0(fov_row, "_", fov_col))
write_csv(density_mglia, paste0(data_path, "/density_mglia_large_cohort.csv"))
