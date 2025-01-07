og_data = read.csv("/Users/bryjc/Downloads/microglia_data_for_bryan_pseudotime_bins_eml.csv")

cols_to_keep = c("label","fov","mask_type","Diagnosis","Pathology","Disease_Status","TMA","Region","Core","fov_row","fov_col",
                 "patient_label","Thal_Phase","Braak_Stage","Path_other","MMSE_Total","Age_at_Death","pca1","pca2","trnsfm",
                 "quantile_value","pst_path1","PSEUDOTIME_TO_DBPN","x","y","y_normalized","pseudotime_bin_01_min","pseudotime_bin_01_max",
                 "cumulative.sum","pseudotime_density_min","pseudotime_density_max","pseudotime_bins_density_based","PSEUDOTIME_NORMALIZED",
                 "pseudotime_bins_eml"
)

mean_neighborhoods = read.csv("/Users/bryjc/Downloads/CellNeighborhood.csv")
mean_neighborhoods = rename(mean_neighborhoods, label = CellID)

updated_mglia_mean_neighbor = inner_join(
  og_data %>% select(cols_to_keep), mean_neighborhoods, by = "label")

# write csv
write.csv(updated_mglia_mean_neighbor, file="/Users/bryjc/Downloads/mglia_with_neighbors.csv")

###
devtools::install_github("mamouzgar/trajiq",dependencies = T)

#3.
library(tidyverse)
library(trajiq)

#4. If using single-cell data (not pseudobulk), Optionally downsample your data using the built-in function in trajic to help things run faster. This function downsamples by group, and if there are not n # of cells in a group, it'll take the maximum available cells. Here we downsample to a maximum of 1000 cells per group. Check the groups are correct in the group_by part of this code. You might want to add the timepoint  column or any other group of interest too!


## load example data
glm_input=trajiq::glm_input
## real data
glm_input=updated_mglia_mean_neighbor

# reorder early/mid/late
glm_input$pseudotime_bins_eml = recode(glm_input$pseudotime_bins_eml, "early" = "low", "middle" = "medium", "late" = "high")
glm_input$pseudotime_bins_eml = factor(glm_input$pseudotime_bins_eml, levels=c("low", "medium", "high"))

## features of interest
my_features = c("PHF.tau.AT8","PHF.tau","Beta.amyloid.1.40","Beta.amyloid.1.42","Beta.amyloid.17.24","Clusterin", "CD31_CD105", "Fe", "GAD65", "GFAP", "Glutamine.Synthetase",
"MAG_MCNPase", "MBP", "NEFL_MAP2_NEFH", "Parvalbumin", "PSD95",
"SMA", "Synaptophysin", "VGAT", "VGlut1")
covariates_in_model = c('patient_label') ## covariates of interest
contrast_variables = c('Disease_Status') ## variable you want to compare

###############################################################################################
## (A) example for running analysis on glm_input features that look like celltype abundance fractions
###############################################################################################
## note: be sure to make covariates as factor variables when it's not a continuous variable. Eg, mouse ids as a covariate for matched analysis
glm_input = glm_input%>%
  mutate(patient_label = factor(patient_label)) 
dif_res_CA1 = trajiq::differential_analysis_program(glm_input %>% filter(Region == 'CA1'), outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm', SPLIT_BY_NAMES= 'pseudotime_bins_eml')
head(dif_res_CA1)
dif_res_Caudate =trajiq::differential_analysis_program(glm_input %>% filter(Region == 'Caudate'), outcome_features = my_features, contrast_variables = contrast_variables, covariates_in_model = covariates_in_model,intercept = TRUE,contrast_method = 'emm', SPLIT_BY_NAMES= 'pseudotime_bins_eml')
head(dif_res_Caudate)

###################################################

# Pull in t tests
t_test_results = read.csv('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/renamed_t_test_results.csv')

# Pull in separately calculated area between curves (matching previous models from Meelad)
abcurves = read.csv('/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/summarized_differences_with_direction.csv')

# rename abcurve 
abcurves = abcurves %>% dplyr::rename(c(pseudotime_bins_eml = "Group", Feature = "Biomarker"))

# Attach curve differences to model significance
combined = t_test_results %>% left_join(abcurves, by=c("Feature","pseudotime_bins_eml","Region"))
combined$pseudotime_bins_eml = factor(combined$pseudotime_bins_eml, levels=c("low", "medium", "high"))

# rename markers
old_names = c("ApoE..pan.", "Beta.amyloid.1.40", "Beta.amyloid.1.42", "Beta.amyloid.17.24", "Biotin", "CleavedCaspase3", "MLKL", "PHF.tau", "PHF.tau.AT8", "UbiquitinK48", "X8OHGuano")
new_names = c("ApoE", "Abeta40", "Abeta42", "Abeta17-24", "TREM2", "Cleaved Caspase 3", "pMLKL", "PHF tau", "PHF tau AT8", "PolyubiquitinK48", "8OHGuanosine")

# Changing values using mutate and recode
combined = combined  %>%
  mutate(Feature = recode(Feature, !!!setNames(new_names, old_names)))

# Assign directionality
combined = combined %>% mutate(Direction_cat = as.factor(if_else(Direction > 0, "AD", "Healthy")))

# Add percentile cat
combined = combined %>% mutate(Effect_cat = ifelse(abs(combined$Total_Difference) >= quantile(abs(combined$Total_Difference), 0.10), 'Y', 'N'))

# Add colors
Disease_Status_colors = c(Healthy = "#A6E002", AD = "#172869")


##
p1 = ggplot(combined %>% filter(Region == "CA1"), aes(x = Total_Difference, y = transf_p_adj)) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = 'transparent', size =1)) +
  geom_hline(yintercept = 1, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 12 ,size = 0.5, linetype = 'dashed', color = 'darkgray') +
  geom_vline(xintercept = -12 ,size = 0.5, linetype = 'dashed', color = 'darkgray') +
  geom_vline(xintercept = 0 ,size = 1, linetype = 'dashed', color = 'darkgray') +
  geom_point(aes(color = Direction_cat), size = 4)   +
  xlab('Effect estimate: Healthy vs AD - CA1') +
  scale_y_continuous(trans = "log10") +
  #scale_y_continuous(limits = c(0,3.5)) +
  ggrepel::geom_text_repel(data = combined %>%
                             filter(Region == "CA1" & transf_p_adj >= 1 & Effect_cat == "Y"), 
                           aes(label = Feature, color = Direction_cat),
                           size = 5, force = 20, max.overlaps=20
  ) +
  scale_color_manual(values = Disease_Status_colors) +
  facet_wrap(~pseudotime_bins_eml)

p2 = ggplot(combined %>% filter(Region == "Caudate"), aes(x = Total_Difference, y = transf_p_adj)) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = 'transparent', size =1)) +
  geom_hline(yintercept = 1, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 12 ,size = 0.5, linetype = 'dashed', color = 'darkgray') +
  geom_vline(xintercept = -12 ,size = 0.5, linetype = 'dashed', color = 'darkgray') +
  geom_vline(xintercept = 0 ,size = 1, linetype = 'dashed', color = 'darkgray') +
  geom_point(aes(color = Direction_cat), size = 4)   +
  xlab('Effect estimate: Healthy vs AD - Caudate') +
  scale_y_continuous(trans = "log10") +
  #scale_y_continuous(limits = c(0,3.5)) +
  ggrepel::geom_text_repel(data = combined %>%
                             filter(Region == "Caudate" & transf_p_adj >= 1 & Effect_cat == "Y"), 
                           aes(label = Feature, color = Direction_cat),
                           size = 5, force = 20, max.overlaps=20
  ) +
  scale_color_manual(values = Disease_Status_colors) +
  facet_wrap(~pseudotime_bins_eml)

save_path = "/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/Papers/mglia_revisions_2024/Plots_Images/personal_MSC_plots"
# p1
ggsave(paste(save_path,"p1_CA1_meanNeighborhood_volcano.pdf", sep = "/"), plot = p1, width=15, height=5, device="pdf", dpi=300)
# p2
ggsave(paste(save_path,"p1_Caudate_meanNighborhood_volcano.pdf", sep = "/"), plot = p2, width=15, height=5, device="pdf", dpi=300)

#
#
# find microglia for representative images
# normalize channels
channels_etc = c("X8OHGuano", "ACC", "ApoE..pan.", "Au", "Beta.amyloid.1.40", "Beta.amyloid.1.42", "Beta.amyloid.17.24", "Biotin", "C1q", "CD31_CD105", "CD33", 
                 "CD44", "CD45", "CD47", "CD68", "CleavedCaspase3", "Clusterin", "Fe" , "GAD65", "GFAP" , "Glutamine.Synthetase", "HH3_dsDNA",  "HLADR", "Iba1", 
                 "LAIR1", "MAG_MCNPase", "MBP", "MFN2", "MLKL", "NEFL_MAP2_NEFH", "P2RY12", "PHF.tau.AT8", "PHF.tau", "PINK1", "PSD95", "Parvalbumin", 
                 "PolyubiquitinK63", "SMA", "Synaptophysin", "TMEM119", "UbiquitinK48", "VGAT", "VGlut1", "chan_39", "chan_48"
)
normalize_99 <- function(x) {
  x / quantile(x, 0.99, na.rm = TRUE)
}

#### Finding mglia df ####
mglia_neighborhood_to_find = as.data.frame(updated_mglia_mean_neighbor)
mglia_neighborhood_to_find[channels_etc] = apply(mglia_neighborhood_to_find[channels_etc], 2, normalize_99)

# find them
hi_AD = mglia_neighborhood_to_find %>% filter(Disease_Status == "AD" & Region == "CA1" & GFAP > 0.7 & Parvalbumin > 0.7) %>%
  select(c(GFAP, Parvalbumin, PSD95, Glutamine.Synthetase, VGlut1, patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))
hi_Healthy = mglia_neighborhood_to_find %>% filter(Disease_Status == "AD" & Region == "CA1" & PSD95 > 0.7 & Glutamine.Synthetase > 0.7 & VGlut1 > 0.7) %>%
  select(c(GFAP, Parvalbumin, PSD95, Glutamine.Synthetase, VGlut1, patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))

################################################################################































# using tajiq
abcurves_CA1 = abcurves %>% filter(Region == "CA1")
abcurves_Caudate = abcurves %>% filter(Region == "Caudate")

# rename abcurve colnames
abcurves_CA1 = abcurves_CA1 %>% dplyr::rename(c(pseudotime_bins_eml = "Group", feature = "Biomarker"))
abcurves_Caudate = abcurves_Caudate %>% dplyr::rename(c(pseudotime_bins_eml = "Group", feature = "Biomarker"))

abcurves_CA1 = abcurves_CA1 %>% dplyr::rename(c(feature = "Feature"))
abcurves_Caudate = abcurves_Caudate %>% dplyr::rename(c(feature = "Feature"))

# Attach curve differences to model significance
dif_res_CA1 = dif_res_CA1 %>% left_join(abcurves_CA1, by=c("feature","pseudotime_bins_eml"))
dif_res_CA1$pseudotime_bins_eml = factor(dif_res_CA1$pseudotime_bins_eml, levels=c("low", "medium", "high"))
dif_res_Caudate = dif_res_Caudate %>% left_join(abcurves_Caudate, by=c("feature","pseudotime_bins_eml"))
dif_res_Caudate$pseudotime_bins_eml = factor(dif_res_Caudate$pseudotime_bins_eml, levels=c("low", "medium", "high"))

# rename markers
old_names = c("ApoE..pan.", "Beta.amyloid.1.40", "Beta.amyloid.1.42", "Beta.amyloid.17.24", "Biotin", "CleavedCaspase3", "MLKL", "PHF.tau", "PHF.tau.AT8", "UbiquitinK48", "X8OHGuano")
new_names = c("ApoE", "Abeta40", "Abeta42", "Abeta17-24", "TREM2", "Cleaved Caspase 3", "pMLKL", "PHF tau", "PHF tau AT8", "PolyubiquitinK48", "8OHGuanosine")

# Changing values using mutate and recode
dif_res_CA1  = dif_res_CA1  %>%
  mutate(feature = recode(feature, !!!setNames(new_names, old_names)))

# Assign directionality
dif_res_CA1 = dif_res_CA1 %>% mutate(Total_Difference_cat = as.factor(if_else(Total_Difference > 0, "AD", "Healthy")))
dif_res_Caudate = dif_res_Caudate %>% mutate(Total_Difference_cat = as.factor(if_else(Total_Difference > 0, "AD", "Healthy")))

# Add colors
Disease_Status_colors = c(Healthy = "#A6E002", AD = "#172869")


##
p1 = ggplot(dif_res_CA1, aes(x = Total_Difference, y = minus_log10padj)) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = 'transparent', size =1)) +
  geom_hline(yintercept = 1, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0 ,size = 1, linetype = 'dashed', color = 'darkgray') +
  geom_point(aes(color = Total_Difference_cat), size = 4)   +
  xlab('Effect estimate: Healthy vs AD - CA1') +
  scale_y_continuous(trans = "log10") +
  #scale_y_continuous(limits = c(0,3.5)) +
  ggrepel::geom_text_repel(data = dif_res_CA1 %>%
                             filter(minus_log10padj >= 1), 
                           aes(label = feature, color = Total_Difference_cat),
                           size = 5, force = 20, max.overlaps=20
  ) +
  scale_color_manual(values = Disease_Status_colors) +
  facet_wrap(~pseudotime_bins_eml)

p2 = ggplot(dif_res_Caudate, aes(x = estimate, y = minus_log10padj)) +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = 'transparent', size =1)) +
  geom_hline(yintercept = 1, linetype ='dashed', size = 1, color = 'darkgray') +
  geom_vline(xintercept = 0 ,size = 1, linetype = 'dashed', color = 'darkgray') +
  geom_point(aes(color = estimate_cat), size = 4)   +
  xlab('Effect estimate: Healthy vs AD - Caudate') +
  scale_y_continuous(limits = c(0,3.5)) +
  ggrepel::geom_text_repel(data = dif_res_Caudate %>%
                             filter(minus_log10padj >= 1), 
                           aes(label = feature, color = estimate_cat),
                           size = 5, force = 20, max.overlaps=20
  ) +
  scale_color_manual(values = Disease_Status_colors) +
  facet_wrap(~pseudotime_bins_eml)


