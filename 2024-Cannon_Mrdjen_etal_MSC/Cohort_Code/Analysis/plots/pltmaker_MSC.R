library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

data_path = "/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data"
trajectory_values = readRDS(paste(data_path, "myFits_mean.rds", sep = "/"))
name_H_AD = "myABCs_eml_region_pval (1).csv"
name_H_only = "myABCs_eml_region_pval-CA1_vs_Caudate.csv"
volcano_sig_H_AD = data.table::fread(file=paste(data_path, name_H_AD, sep = "/"))
volcano_sig_H_only = data.table::fread(paste(data_path, name_H_only, sep = "/"))

# reorder early/mid/late
volcano_sig_H_only$pseudotime_bins_eml = factor(volcano_sig_H_only$pseudotime_bins_eml, levels=c("early", "middle", "late"))
volcano_sig_H_AD$pseudotime_bins_eml = factor(volcano_sig_H_AD$pseudotime_bins_eml, levels=c("early", "middle", "late"))

# rename markers
old_names = c("ApoE..pan.", "Beta.amyloid.1.40", "Beta.amyloid.1.42", "Beta.amyloid.17.24", "Biotin", "CleavedCaspase3", "MLKL", "PHF.tau", "PHF.tau.AT8", "UbiquitinK48", "X8OHGuano")
new_names = c("ApoE", "Abeta40", "Abeta42", "Abeta17-24", "TREM2", "Cleaved Caspase 3", "pMLKL", "PHF tau", "PHF tau AT8", "PolyubiquitinK48", "8OHGuanosine")

# Changing values using mutate and recode
trajectory_values = trajectory_values %>%
  mutate(feature = recode(feature, !!!setNames(new_names, old_names)))

volcano_sig_H_only = volcano_sig_H_only %>%
  mutate(feature = recode(feature, !!!setNames(new_names, old_names)))

volcano_sig_H_AD = volcano_sig_H_AD %>%
  mutate(feature = recode(feature, !!!setNames(new_names, old_names)))

# Colors
Region_colors = c(CA1 = "#F8766C" , Caudate = "#00BFC4")
Disease_Status_colors = c(Healthy = "#A6E002", AD = "#172869")

# sig markers to include in trajectory plots
traj_markers_H_only_immune = c("CD68", "CD44", "P2RY12", "HLADR", "ApoE", "Iba1", "CD45", "TREM2", "CD33", "TMEM119")
traj_markers_H_only_other = c("PolyubiquitinK63", "PolyubiquitinK48", "Abeta42", "PINK1", "ACC", "Fe", "Cleaved Caspase 3", 
                              "pMLKL", "Parvalbumin", "PHF tau")
traj_markers_H_AD_CA1_immune = c("Iba1", "CD45", "CD44", "ApoE", "HLADR", "TMEM119", "CD33", "P2RY12", "TREM2")
traj_markers_H_AD_CA1_other = c("Cleaved Caspase 3", "PSD95", "Parvalbumin", "ACC", "Abeta40", "PHF tau")
traj_markers_H_AD_Caudate_immune = c("CD45", "ApoE")
traj_markers_H_AD_Caudate_other = c("Fe", "Parvalbumin", "ACC", "Abeta40", "Abeta42", "8OHGuanosine", "PINK1")


p1.1 = ggplot(trajectory_values %>% filter(Disease_Status == "Healthy" & feature %in% traj_markers_H_only_immune), aes(x=x, y=y, group=Region, color=Region)) + 
  geom_line(linewidth = 2, aes(color = Region)) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', linewidth = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', linewidth = 1, color = 'darkgray') +
  #facet_grid(rows=vars(feature), scales="free_y") +
  facet_wrap(~feature, scales="free_y") +
  scale_color_manual(values = Region_colors) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = Disease_Status_colors["Healthy"], 
                                    fill = NA, 
                                    size = 1))

p1.2 = ggplot(trajectory_values %>% filter(Disease_Status == "Healthy" & feature %in% traj_markers_H_only_other), aes(x=x, y=y, group=Region, color=Region)) + 
  geom_line(linewidth = 2, aes(color = Region)) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', linewidth = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', linewidth = 1, color = 'darkgray') +
  #facet_grid(rows=vars(feature), scales="free_y") +
  facet_wrap(~feature, scales="free_y") +
  scale_color_manual(values = Region_colors) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = Disease_Status_colors["Healthy"], 
                                    fill = NA, 
                                    size = 1))

p1.3 = ggplot(trajectory_values %>% filter(Disease_Status == "Healthy" & (feature %in% traj_markers_H_only_immune | feature %in% traj_markers_H_only_other)), aes(x=x, y=y, group=Region, color=Region)) + 
  geom_line(linewidth = 2, aes(color = Region)) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', linewidth = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', linewidth = 1, color = 'darkgray') +
  #facet_grid(rows=vars(feature), scales="free_y") +
  facet_wrap(~feature, scales="free_y") +
  scale_color_manual(values = Region_colors) +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(color = Disease_Status_colors["Healthy"], 
                                    fill = NA, 
                                    size = 1))

p2.1 = ggplot(trajectory_values %>% filter(Region == "CA1" & feature %in% c(traj_markers_H_AD_CA1_immune, traj_markers_H_AD_Caudate_immune)),
              aes(x=x, y=y, group=Disease_Status, color=Disease_Status)) + 
  geom_line(linewidth = 2) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', size = 1, color = 'darkgray') +
  facet_grid(rows=vars(feature), cols=vars(Region), scales="free_y") +
  scale_color_manual(values = Disease_Status_colors) +
  theme_minimal(base_size = 12)  +
  theme(panel.border = element_rect(color = Region_colors["CA1"], 
                                    fill = NA, 
                                    size = 1))

p2.2 = ggplot(trajectory_values %>% filter(Region == "CA1" & feature %in% c(traj_markers_H_AD_CA1_other, traj_markers_H_AD_Caudate_other)), 
              aes(x=x, y=y, group=Disease_Status, color=Disease_Status)) + 
  geom_line(linewidth = 2) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', size = 1, color = 'darkgray') +
  facet_grid(rows=vars(feature), cols=vars(Region), scales="free_y") +
  scale_color_manual(values = Disease_Status_colors) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = Region_colors["CA1"], 
                                    fill = NA, 
                                    size = 1))

p2AND3 = ggplot(trajectory_values %>% filter(feature %in% c(traj_markers_H_AD_CA1_immune, traj_markers_H_AD_Caudate_immune, traj_markers_H_AD_CA1_other, traj_markers_H_AD_Caudate_other)), 
              aes(x=x, y=y, group=Disease_Status, color=Disease_Status)) + 
  geom_line(linewidth = 2) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', size = 1, color = 'darkgray') +
  #facet_grid(rows=vars(feature), scales="free_y") +
  facet_wrap(vars(feature, Region), scales="free_y") +
  scale_color_manual(values = Disease_Status_colors) +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

p3.1 = ggplot(trajectory_values %>% filter(Region == "Caudate" & feature %in% c(traj_markers_H_AD_CA1_immune, traj_markers_H_AD_Caudate_immune)),
                                           aes(x=x, y=y, group=Disease_Status, color=Disease_Status)) + 
  geom_line(linewidth = 2) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', size = 1, color = 'darkgray') +
  facet_grid(rows=vars(feature), cols=vars(Region), scales="free_y") +
  scale_color_manual(values = Disease_Status_colors) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(Region_colors["Caudate"], 
                                    fill = NA, 
                                    size = 1))

p3.2 = ggplot(trajectory_values %>% filter(Region == "Caudate" & feature %in% c(traj_markers_H_AD_CA1_other, traj_markers_H_AD_Caudate_other)),
              aes(x=x, y=y, group=Disease_Status, color=Disease_Status)) + 
  geom_line(linewidth = 2) +
  geom_ribbon(aes(ymin=y-std, ymax=y+std), alpha=0, linetype=2, linewidth = 1) +
  geom_vline(xintercept = 0.33333333, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0.66666666, linetype ='dashed', size = 1, color = 'darkgray') +
  facet_grid(rows=vars(feature), cols=vars(Region), scales="free_y") +
  scale_color_manual(values = Disease_Status_colors) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = Region_colors["Caudate"], 
                                    fill = NA, 
                                    size = 1))

# volcano plots below

p4 = ggplot(volcano_sig_H_only, aes(x = aBC_Direction, y = -1*log10(padj_threshold)))+
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = 'transparent', linewidth = 1)) +
  geom_hline(yintercept = 1, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0, size = 1, linetype = 'dashed', color = 'darkgray') +
  geom_point(aes(color = aBC_Direction_cat), size = 4)   +
  xlab('area between curve: CA1 vs Caudate') +
  scale_y_continuous(limits = c(0,3.5)) +
  ggrepel::geom_text_repel(data = volcano_sig_H_only %>% 
                           filter(-1*log10(padj_threshold) >= 1), 
                           aes(label = feature, color = aBC_Direction_cat), 
                           size = 5, force = 20, max.overlaps=20
                           ) +
  scale_color_manual(values = Region_colors) +
  facet_wrap(~pseudotime_bins_eml)

p5 = ggplot(volcano_sig_H_AD %>% filter(Region=="CA1"), aes(x = aBC_Direction, y = -1*log10(padj_threshold)))+
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = 'transparent', size = 1)) +
  geom_hline(yintercept = 1, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0, size = 1, linetype = 'dashed', color = 'darkgray') +
  geom_point(aes(color = aBC_Direction_cat), size = 4)   +
  xlab('area between curve: Healthy vs AD - CA1') +
  scale_y_continuous(limits = c(0,3.5)) +
  ggrepel::geom_text_repel(data = volcano_sig_H_AD %>% filter(Region=="CA1") %>% 
                           filter(-1*log10(padj_threshold) >= 1), 
                           aes(label = feature, color = aBC_Direction_cat), 
                           size = 5, force = 20, max.overlaps=20
  ) +
  scale_color_manual(values = Disease_Status_colors) +
  facet_wrap(~pseudotime_bins_eml)

p6 = ggplot(volcano_sig_H_AD %>% filter(Region=="Caudate"), aes(x = aBC_Direction, y = -1*log10(padj_threshold))) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = 'transparent', size =1)) +
  geom_hline(yintercept = 1, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0 ,size = 1, linetype = 'dashed', color = 'darkgray') +
  geom_point(aes(color = aBC_Direction_cat), size = 4)   +
  xlab('area between curve: Healthy vs AD - Caudate') +
  scale_y_continuous(limits = c(0,3.5)) +
  ggrepel::geom_text_repel(data = volcano_sig_H_AD %>% filter(Region=="Caudate") %>% 
                           filter(-1*log10(padj_threshold) >= 1), 
                           aes(label = feature, color = aBC_Direction_cat),
                           size = 5, force = 20, max.overlaps=20
  ) +
  scale_color_manual(values = Disease_Status_colors) +
  facet_wrap(~pseudotime_bins_eml)


#save_path
save_path = "/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/Papers/mglia_revisions_2024/Plots_Images/personal_MSC_plots"
# p1_1
ggsave(paste(save_path,"p1_1_traj_HOnly_Imm.pdf", sep = "/"), plot = p1.1, width=6, height=25, device="pdf", dpi=300)
# p1_2
ggsave(paste(save_path,"p1_2_traj_HOnly_Other.pdf", sep = "/"), plot = p1.2, width=6, height=25, device="pdf", dpi=300)
# p1_3
ggsave(paste(save_path,"p1_3_traj_HOnly_ImmOther.pdf", sep = "/"), plot = p1.3, width=20, height=10, device="pdf", dpi=300)
# p2_1
ggsave(paste(save_path,"p2_1_traj_HAD_CA1_Imm.pdf", sep = "/"), plot = p2.1, width=6, height=25, device="pdf", dpi=300)
# p2_2
ggsave(paste(save_path,"p2_2_traj_HAD_CA1_Other.pdf", sep = "/"), plot = p2.2, width=6, height=25, device="pdf", dpi=300)
# p2_3
ggsave(paste(save_path,"p2AND3_traj_HAD_CA1Caudate_ImmOther.pdf", sep = "/"), plot = p2AND3, width=35, height=20, device="pdf", dpi=300)
# p3_1
ggsave(paste(save_path,"p3_1_traj_HAD_Caudate_Other.pdf", sep = "/"), plot = p3.1, width=6, height=25, device="pdf", dpi=300)
# p3_2
ggsave(paste(save_path,"p3_2_traj_HAD_Caudate_Other.pdf", sep = "/"), plot = p3.2, width=6, height=25, device="pdf", dpi=300)
# p4
ggsave(paste(save_path,"p4_volc_HOnly_split.pdf", sep = "/"), plot = p4, width=15, height=5, device="pdf", dpi=300)
# p5
ggsave(paste(save_path,"p5_volc_HAD_CA1_split.pdf", sep = "/"), plot = p5, width=15, height=5, device="pdf", dpi=300)
# p6
ggsave(paste(save_path,"p6_volc_HAD_Caudate_split.pdf", sep = "/"), plot = p6, width=15, height=5, device="pdf", dpi=300)
  

#####
# other - for reviewer comments (R1)
merged_temp = merge(mglia_to_save, og_data[c("pseudotime_bins_eml", "label")], by="label")
merged_temp$pseudotime_bins_eml = factor(merged_temp$pseudotime_bins_eml, levels=c("early", "middle", "late"))
# reorder early/mid/late
merged_temp$pseudotime_bins_eml = factor(merged_temp$pseudotime_bins_eml, levels=c("low", "medium", "high"))

# Calculate the 99th percentile of HLADR, Iba1
HLADR_percentile_99 <- quantile(merged_temp$HLADR, 0.99, na.rm = TRUE)
Iba1_percentile_99 <- quantile(merged_temp$Iba1, 0.99, na.rm = TRUE)

# Normalize HLADR, Iba1 by the 99th percentile
merged_temp <- merged_temp %>%
  mutate(HLADR_normalized = HLADR / HLADR_percentile_99)
merged_temp <- merged_temp %>%
  mutate(Iba1_normalized = Iba1 / Iba1_percentile_99)

p_reviewer_6 = ggplot(merged_temp) +
       aes(x = pseudotime_bins_eml, y = HLADR_normalized, fill = Disease_Status) +
       geom_violin(scale = "area") +
       stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "white", position = position_dodge(width = 0.9)) +
       stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 5), size = 20), 
                                       vjust = -7.5, hjust = -0.2, color = "black", position = position_dodge(width = 0.9)) +
       scale_fill_manual(values = Disease_Status_colors) +
       scale_y_log10() +  # Apply log10 transformation to the y-axis
       scale_x_discrete(labels = c("early" = "low", "middle" = "medium", "late" = "high")) +  # Rename x-axis labels
       theme_minimal(base_size = 20)

p_reviewer_6 = ggplot(merged_temp) +
  aes(x = pseudotime_bins_eml, y = Iba1_normalized, fill = Disease_Status) +
  geom_violin(scale = "area") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "white", position = position_dodge(width = 0.9)) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 3), size = 20), 
               vjust = -5.5, hjust = 0, color = "black", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = Disease_Status_colors) +
  scale_y_log10() +  # Apply log10 transformation to the y-axis
  scale_x_discrete(labels = c("early" = "low", "middle" = "medium", "late" = "high")) +  # Rename x-axis labels
  geom_hline(yintercept = 0.2731, linetype = "dotted", color = "black", linewidth = 1) +  # Horizontal line for Microglia
  geom_hline(yintercept = 0.003, linetype = "dotted", color = "black", linewidth = 1) +   # Horizontal line for Neurons
  #annotate("text", x = Inf, y = 0.2731, label = "Microglia: average Iba1 expr", hjust = 1.1, vjust = -0.5, color = "gray") +  # Label for Microglia line
  #annotate("text", x = Inf, y = 0.003, label = "Neurons: average Iba1 expr", hjust = 1.1, vjust = -0.5, color = "gray") +    # Label for Neurons line
  facet_wrap(vars(Region)) +
  theme_minimal(base_size = 20)

ggsave(paste(save_path,"p_reviewer_6.pdf", sep = "/"), plot = p_reviewer_6, width=15, height=5, device="pdf", dpi=300)



######
# Calculate means and standard errors for each group
summary_stats <- merged_temp %>%
  group_by(pseudotime_bins_eml, Disease_Status, Region) %>%
  summarise(
    mean_HLADR = mean(Iba1_normalized, na.rm = TRUE),
    sd_HLADR = sd(HLADR_normalized, na.rm = TRUE),
    se_HLADR = sd(HLADR_normalized, na.rm = TRUE) / sqrt(n()),
    mean_Iba1 = mean(Iba1_normalized, na.rm = TRUE),
    sd_Iba1 = sd(Iba1_normalized, na.rm = TRUE),
    se_Iba1 = sd(Iba1_normalized, na.rm = TRUE) / sqrt(n())
  )

# Now plot using the normalized HLADR values
ggplot(summary_stats) +
  aes(x = pseudotime_bins_eml, y = mean_HLADR, fill = Disease_Status) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_HLADR - se_HLADR, ymax = mean_HLADR + se_HLADR), 
                width = 0.3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = Disease_Status_colors) +
  
  scale_x_discrete(labels = c("early" = "low", "middle" = "medium", "late" = "high")) +  # Rename x-axis labels
  theme_minimal() +
  facet_wrap(vars(Region))

# Now plot using the normalized Iba1 values
ggplot(summary_stats) +
  aes(x = pseudotime_bins_eml, y = mean_Iba1, fill = Disease_Status) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_Iba1 - se_Iba1, ymax = mean_Iba1 + se_Iba1), 
                width = 0.3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = Disease_Status_colors) +
  
  scale_x_discrete(labels = c("early" = "low", "middle" = "medium", "late" = "high")) +  # Rename x-axis labels
  theme_minimal() +
  facet_wrap(vars(Region))



