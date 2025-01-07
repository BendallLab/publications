install.packages("data.table", "dplyr")
library(data.table)
library(dplyr)
# import cell table
fname = "/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/cell_table/microglia-DM_table_size_normalized.csv"
mglia_DM_RAW <- data.table::fread(fname)
# import and add cohort data
cohort_annot <- data.table::fread("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/docs/cohort_info/final_merged_cohort_data.csv")
mglia_DM_annot <- mglia_DM_RAW %>% dplyr::left_join(cohort_annot, by = "fov")

# channels
channels_etc = c("8OHGuano", "ACC", "ApoE (pan)", "Au", "Beta amyloid 1 40", "Beta amyloid 1 42", "Beta amyloid 17-24", "Biotin", "C1q", "CD31_CD105", "CD33", 
                    "CD44", "CD45", "CD47", "CD68", "CleavedCaspase3", "Clusterin", "Fe" , "GAD65", "GFAP" , "Glutamine Synthetase", "HH3_dsDNA",  "HLADR", "Iba1", 
                    "LAIR1", "MAG_MCNPase", "MBP", "MFN2", "MLKL", "NEFL_MAP2_NEFH", "Noodle", "P2RY12", "PHF-tau-AT8", "PHF-tau", "PINK1", "PSD95", "Parvalbumin", 
                    "PolyubiquitinK63", "SMA", "Synaptophysin", "TMEM119", "UbiquitinK48", "VGAT", "VGlut1", "chan_39", "chan_48"
)

# channels MSC
channels_MSC = c("ApoE (pan)", "Biotin", "C1q", "CD44", "CD45", "CD68", "HLADR", "Iba1", "P2RY12", "TMEM119", "UbiquitinK48")

# antigen targets you want to include in (1) for the statistical analysis
channels_to_check = c("ACC", "ApoE (pan)", "Biotin", "C1q", "CD33", 
                      "CD44", "CD45", "CD68", "CleavedCaspase3", "Fe" , "HLADR", "Iba1", 
                      "LAIR1", "MLKL", "P2RY12", "PINK1", "PSD95", "Parvalbumin", 
                      "PolyubiquitinK63", "TMEM119", "UbiquitinK48", 
                      "8OHGuano", "Beta amyloid 1 40", "Beta amyloid 1 42", "Beta amyloid 17-24", "PHF-tau-AT8", "PHF-tau"
                      )
# morphology targets 
morphos_to_check = c('area', 'eccentricity', 'major_axis_length', 'minor_axis_length', 'perimeter', 
                     'major_minor_axis_ratio', 'perim_square_over_area', 
                     'major_axis_equiv_diam_ratio', 'num_concavities'
                     )

# Processed data table
# export csv
write.csv(mglia_DM_annot, "/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/cell_table/microglia-DM_annot_table_size_transformed.csv")
# Need to handoff caudate, healthy and ad samples, and healthy controls
# export csv with just healthy
