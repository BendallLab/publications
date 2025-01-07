# channels plus elements -> Gold, Potassium(39), Titanium(48), Iron, Noodle 
channels_plus_elements = list("8OHGuano", "ACC", "ApoE (pan)", "Au", "Beta amyloid 1 40", "Beta amyloid 1 42", "Beta amyloid 17-24", "Biotin", "C1q", "CD31_CD105", "CD33", "CD44", 
                              "CD45", "CD47", "CD68", "CleavedCaspase3", "Clusterin", "Fe" , "GAD65", "GFAP" , "Glutamine Synthetase", "HH3_dsDNA",  "HLADR", "Iba1", "LAIR1", 
                              "MAG_MCNPase", "MBP", "MFN2", "MLKL", "NEFL_MAP2_NEFH", "Noodle", "P2RY12", "PHF-tau-AT8", "PHF-tau", "PINK1", "PSD95", "Parvalbumin", "PolyubiquitinK63", 
                              "SMA", "Synaptophysin", "TMEM119", "UbiquitinK48", "VGAT", "VGlut1", "chan_39", "chan_48"
                              )
# panel channels alone
channels = list("8OHGuano", "ACC", "ApoE (pan)", "Beta amyloid 1 40", "Beta amyloid 1 42", "Beta amyloid 17-24", "Biotin", "C1q", "CD31_CD105", "CD33", "CD44", 
                "CD45", "CD47", "CD68", "CleavedCaspase3", "Clusterin", "GAD65", "GFAP" , "Glutamine Synthetase", "HH3_dsDNA",  "HLADR", "Iba1", "LAIR1", 
                "MAG_MCNPase", "MBP", "MFN2", "MLKL", "NEFL_MAP2_NEFH", "P2RY12", "PHF-tau-AT8", "PHF-tau", "PINK1", "PSD95", "Parvalbumin", "PolyubiquitinK63", 
                "SMA", "Synaptophysin", "TMEM119", "UbiquitinK48", "VGAT", "VGlut1"
                )

# annotation columns
annot_cols = list("fov", "Diagnosis", "Pathology", "Disease_Status", "TMA", "Region", "Core", "fov_row", "fov_col")

# common markers between DM and BJC panels
common_ch = list("8OHGuano", "ApoE (pan)", "Beta amyloid 1 40", "Beta amyloid 1 42", "Beta amyloid 17-24", "Biotin", "CD31_CD105", "CD44", 
                 "CD45", "CD68", "GFAP" , "HH3_dsDNA",  "HLADR", "Iba1", "MAG_MCNPase", "MFN2", "NEFL_MAP2_NEFH", "P2RY12", "PHF-tau", 
                 "TMEM119", "VGAT", "VGlut1"
                 )
