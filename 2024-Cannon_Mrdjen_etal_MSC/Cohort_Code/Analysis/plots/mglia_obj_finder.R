#### Variables ####
channels_etc = c("8OHGuano", "ACC", "ApoE (pan)", "Au", "Beta amyloid 1 40", "Beta amyloid 1 42", "Beta amyloid 17-24", "Biotin", "C1q", "CD31_CD105", "CD33", 
                 "CD44", "CD45", "CD47", "CD68", "CleavedCaspase3", "Clusterin", "Fe" , "GAD65", "GFAP" , "Glutamine Synthetase", "HH3_dsDNA",  "HLADR", "Iba1", 
                 "LAIR1", "MAG_MCNPase", "MBP", "MFN2", "MLKL", "NEFL_MAP2_NEFH", "Noodle", "P2RY12", "PHF-tau-AT8", "PHF-tau", "PINK1", "PSD95", "Parvalbumin", 
                 "PolyubiquitinK63", "SMA", "Synaptophysin", "TMEM119", "UbiquitinK48", "VGAT", "VGlut1", "chan_39", "chan_48"
)
normalize_99 <- function(x) {
  x / quantile(x, 0.99, na.rm = TRUE)
}

#### Finding mglia df ####
mglia_to_find = as.data.frame(mglia_to_save)
mglia_to_find[channels_etc] = apply(mglia_to_find[channels_etc], 2, normalize_99)

#### Filtering ####
# Hi MSC AD UP Program
hi_MSC_AD_CA1_mglia_UP = mglia_to_find %>% filter(Disease_Status == "AD" & Region == "CA1" & Iba1 > 0.5 & CD33 > 0.4 & HLADR < 0.2 & `ApoE (pan)` < 0.2 & P2RY12 < 0.2 & CD44 < 0.8) %>%
  select(c(Iba1, HLADR, `ApoE (pan)`, P2RY12, CD44, CD33, patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))

# Hi MSC Healthy UP Program
hi_MSC_H_CA1_mglia_UP = mglia_to_find %>% filter(Disease_Status == "Healthy" & Region == "CA1" & Iba1 > 0.5 & CD33 < 0.3 & HLADR > 0.3 & `ApoE (pan)` > 0.3 & P2RY12 > 0.3 & CD44 < 0.6) %>%
  select(c(Iba1, HLADR, `ApoE (pan)`, P2RY12, CD44, CD33,  patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))
#
#
#
# Hi MSC AD Down Program - Basic
hi_MSC_AD_CA1_mglia_UP_Base = mglia_to_find %>% filter(Disease_Status == "AD" & Region == "CA1" & Iba1 > 0.5 & (CD33 > 0.7 | CD44 > 0.7)) %>%
  select(c(Iba1, HLADR, `ApoE (pan)`, P2RY12, CD44, CD33, patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))

# Hi MSC Healthy UP Program - Basic
hi_MSC_H_CA1_mglia_UP_Base = mglia_to_find %>% filter(Disease_Status == "Healthy" & Region == "CA1" & Iba1 > 0.5 & (HLADR > 0.9 & `ApoE (pan)` > 0.9 & P2RY12 > 0.9)) %>%
  select(c(Iba1, HLADR, `ApoE (pan)`, P2RY12, CD44, CD33, patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))

#
table(hi_MSC_AD_CA1_mglia_UP_Base$patient_label, hi_MSC_AD_CA1_mglia_UP_Base$row_col)

#
#
#
# Hi Healthy only - CA1 UP Program
hi_MSC_Region_CA1_mglia_UP = mglia_to_find %>% filter(Disease_Status == "Healthy" & Region == "CA1" & Iba1 > 0.5 & HLADR > 0.7 & P2RY12 > 0.7 & Biotin < 0.6 & CD45 < 0.6 & CD68 < 0.6) %>%
  select(c(Iba1, HLADR, P2RY12, `ApoE (pan)`, CD44, Biotin, CD45, CD68, patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))

# Hi Healthy only - Caudate UP Program
hi_MSC_Region_Caudate_mglia_UP = mglia_to_find %>% filter(Disease_Status == "Healthy" & Region == "Caudate" & Iba1 > 0.5 & HLADR < 0.2 & P2RY12 < 0.2 & Biotin > 0.4 & CD45 > 0.4 & CD68 > 0.4) %>%
  select(c(Iba1, HLADR, P2RY12, `ApoE (pan)`, CD44, Biotin, CD45, CD68, patient_label, fov_row, fov_col)) %>% mutate(row_col = paste0(fov_row,"_",fov_col))

#
#
#

# Hi Healthy only - CA1 PSD95 Program
hi_MSC_excit_CA1_mglia_UP = mglia_to_find %>% filter(Disease_Status == "Healthy" & Region == "CA1" & Iba1 > 0.5 & PSD95 > 0.7 & Parvalbumin < 0.3 & `Beta amyloid 1 40` < 0.1 & `PHF-tau` < 0.1) %>%
  select(c(Iba1, PSD95, Parvalbumin, HLADR, P2RY12, `ApoE (pan)`, CD44, Biotin, CD45, CD68, CD33, `Beta amyloid 1 40`, `PHF-tau`, patient_label, fov_row, fov_col)) %>% 
  mutate(row_col = paste0(fov_row,"_",fov_col))

# Hi AD only - CA1 Parv Program
hi_MSC_inhib_CA1_mglia_UP = mglia_to_find %>% filter(Disease_Status == "AD" & Region == "CA1" & Iba1 > 0.5 & PSD95 < 0.3 & Parvalbumin > 0.7 & `Beta amyloid 1 40` > 0.1 & `PHF-tau` > 0.1) %>%
  select(c(Iba1, PSD95, Parvalbumin, HLADR, P2RY12, `ApoE (pan)`, CD44, Biotin, CD45, CD68, CD33, `Beta amyloid 1 40`, `PHF-tau`, patient_label, fov_row, fov_col)) %>% 
  mutate(row_col = paste0(fov_row,"_",fov_col))


