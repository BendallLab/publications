
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
BiocManager::install("esquisse")
library("esquisse")
library("sva")
library("ggplot2")

#### PCA plots ####
# Assuming df is your dataframe and X is the factor for splitting/coloring
# Remove the factor variable before performing PCA
df_numeric <- mglia_DM_annot %>% select(all_of(unlist(channels_etc)))
df_scaled <- scale(df_numeric)  # Scale numeric data

# Perform PCA - scale. = TRUE standardizes the data
pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)

library(ggplot2)

# Add the PCA scores to the original dataframe for plotting
mglia_DM_annot_pca <- mglia_DM_annot
mglia_DM_annot_pca$pca1 <- pca_result$x[,1]
mglia_DM_annot_pca$pca2 <- pca_result$x[,2]
mglia_DM_annot_pca$TMA <- as.factor(mglia_DM_annot_pca$TMA)



#### Size filtering ####
# Filter df based on markers
#TMEM119 < 0.0015, CD31_CD105 < 0.0015, C1q < 0.0015, HLA-DR < 0.0015, ApoE < 0.02, Biotin (TREM2) < 0.02 (all other MSC markers below 0.025)
#size filtering keeps nucleated cells and anucleated, larger projections.
mglia_DM_annot_pca_sized <- mglia_DM_annot_pca %>%
  filter(!if_any(any_of(channels_MSC), ~ . > 0.025)) %>%
  filter(TMEM119 < 0.0015 & CD31_CD105 < 0.0015 & C1q < 0.0015 & HLADR < 0.0015) %>%
  filter((HH3_dsDNA > 0.01 & cell_size < 4000) | (HH3_dsDNA < 0.01 & major_minor_axis_ratio >= 2))

# Assuming 'df' is your dataframe with columns 'A' and 'B'
ggplot(data = mglia_DM_annot, aes(x = CD45, y = TMEM119)) +
  geom_point() +  # This adds the points to the plot
  theme_minimal() +  # Optional: Adds a minimal theme for aesthetics
  labs(title = "Biaxial Plot of A vs. B", x = "Column A", y = "Column B")  # Adding labels

# dna plotting
ggplot(mglia_DM_annot) +
  aes(x = HH3_dsDNA) +
  geom_histogram(bins = 100L, fill = "#112446") +
  scale_x_continuous(trans = "log10") +
  theme_minimal()

#### Plotting PCAs ####
library(ggplot2)
library(dplyr)

# Filter df based on variables
mglia_DM_annot_pca_toplot <- mglia_DM_annot_pca_sized %>%
  #filter(TMA != 0 & Core == "HC")
  filter(TMA != 0 & Core != "HC" & Core != "AD")# & patient_label != "7-3")

# Define a vector of 12 colors
my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
               "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB")

ggplot(data = mglia_DM_annot_pca_toplot, aes(x = pca1, y = pca2)) +
  geom_point(data = mglia_DM_annot_pca_toplot %>% sample_n(5000) %>% select(-c("patient_label")), color = "grey", size = 0.5) +  # Plot all points in grey
  geom_point(aes(x = pca1, y = pca2, color = patient_label), size = 0.5) +  # Overlay points with color based on group
  #scale_color_manual(values = my_colors) +  # Use the custom colors [c(1,2,3)]
  facet_wrap(~patient_label) +  # Create a facet for each category of X
  theme_minimal() +
  #xlim(-50,100) +
  #ylim(-50,50) +
  labs(title = "PCA Plot Colored by pt_label", x = "Principal Component 1", y = "Principal Component 2") +
  theme(legend.position = "right")  # Adjust legend position


# Extracting loadings
loadings <- pca_result$rotation

# For PC1
loadings_PC1 <- sort(loadings[,1], decreasing = TRUE)

# For PC2
loadings_PC2 <- sort(loadings[,2], decreasing = TRUE)

# You can print or analyze the loadings_PC1 and loadings_PC2 to see the most contributing variables
print(loadings_PC1)
print(loadings_PC2)




#### Final filtering ####
mglia_to_save <- mglia_DM_annot_pca_sized %>%
  filter(TMA != 0 & Core != "HC" & Core != "AD") %>%
  filter(Disease_Status != "Resilient")

#### ComBat adjustment ####
  pheno = mglia_DM_annot %>% select(!all_of(unlist(channels_etc)))
  edata = t(as.data.frame(mglia_DM_annot %>% select(all_of(unlist(channels_etc))) %>% scale()))
  batch = pheno$patient_label
  modcombat = model.matrix(~1, data=pheno)

  
  # parametric adjustment
  combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
  
  # non-parametric adjustment, mean-only version
  combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

  
  # combat parametric returned to dataset:
  edata_post_combat = t(combat_edata1)
  # Perform PCA - scale. = TRUE standardizes the data
  pca_postcombat_result <- prcomp(edata_post_combat, center = TRUE, scale. = TRUE)
  
  
  # Add the PCA scores to the original dataframe for plotting
  mglia_DM_postcombat_pca <- cbind(pheno, edata_post_combat)
  mglia_DM_postcombat_pca$pca1 <- pca_postcombat_result$x[,1]
  mglia_DM_postcombat_pca$pca2 <- pca_postcombat_result$x[,2]
  mglia_DM_postcombat_pca$TMA <- as.factor(mglia_DM_postcombat_pca$TMA)
  
  
#### Processed data table ####
# export csv
write.csv(mglia_to_save, "/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/cell_table/microglia-DM_fixedV3_wAnucleated_table_size_transformed.csv")
# export healthy only csv
mglia_to_save_HONLY <- filter(mglia_to_save, Disease_Status == "Healthy")  
write.csv(mglia_to_save_HONLY, "/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/cell_table/microglia-DM_fixedV3_wAnucleated_table_size_transformed_HONLY.csv")
