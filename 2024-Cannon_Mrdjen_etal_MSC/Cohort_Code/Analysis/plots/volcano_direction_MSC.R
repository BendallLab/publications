# Establish MSC binning
trajectory_values = trajectory_values %>% 
  mutate(pseudotime_bins_eml = case_when((x >= 0 & x < 0.33) ~ "early",
                                         (x >= 0.33 & x <= 0.66) ~ "middle",
                                         (x > 0.66 & x <= 1) ~ "late"))
# Difference in Curve Calc
# group_by(Region, feature, pseudotime_bins_eml) -> subtract Disease_Status(AD) y - Disease_Status(H) y -> if positive, assign +1 (areaBetweenCurves), if negative, assign -1 (areaBetweenCurves)
trajectory_values = trajectory_values %>% group_by(x, Region, feature) %>% mutate(dValue_DS = y[Disease_Status == "AD"] - y)
trajectory_values = trajectory_values %>% group_by(x, Disease_Status, feature) %>% mutate(dValue_Reg = y[Region == "CA1"] - y)


# Caclulate mean curve differences
summed_traj_values_DiseaseStatus = trajectory_values %>% group_by(pseudotime_bins_eml, Region, feature) %>% summarise(avg_dValueDS = mean(dValue_DS)) 
summed_traj_values_Region = trajectory_values %>% group_by(pseudotime_bins_eml, Disease_Status, feature) %>% summarise(avg_dValueReg = mean(dValue_Reg))

# Join to sig df
volcano_sig_H_AD = volcano_sig_H_AD %>% left_join(summed_traj_values_DiseaseStatus, by = c("pseudotime_bins_eml", "Region", "feature"))
volcano_sig_H_only = volcano_sig_H_only %>% left_join(summed_traj_values_Region %>% filter(Disease_Status == "Healthy"), by = c("pseudotime_bins_eml", "feature"))

# Assign directionality to areaBetweenCurves
volcano_sig_H_AD = volcano_sig_H_AD %>% mutate(aBC_Direction = if_else(avg_dValueDS > 0, areaBetweenCurves, -1*areaBetweenCurves))
volcano_sig_H_only = volcano_sig_H_only %>% mutate(aBC_Direction = if_else(avg_dValueReg > 0, areaBetweenCurves, -1*areaBetweenCurves))

# Add a categorical element
volcano_sig_H_AD = volcano_sig_H_AD %>% mutate(aBC_Direction_cat = as.factor(if_else(aBC_Direction > 0, "AD", "Healthy")))
volcano_sig_H_only = volcano_sig_H_only %>% mutate(aBC_Direction_cat = as.factor(if_else(aBC_Direction > 0, "CA1", "Caudate")))

# Add percentile cat
volcano_sig_H_AD = volcano_sig_H_AD %>% mutate(Effect_cat = ifelse(abs(volcano_sig_H_AD$areaBetweenCurves) >= quantile(abs(volcano_sig_H_AD$areaBetweenCurves), 0.75), 'Y', 'N'))
volcano_sig_H_only = volcano_sig_H_only %>% mutate(Effect_cat = ifelse(abs(volcano_sig_H_only$areaBetweenCurves) >= quantile(abs(volcano_sig_H_only$areaBetweenCurves), 0.75), 'Y', 'N'))
