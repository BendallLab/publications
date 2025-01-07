import pandas as pd
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multitest
import numpy as np

# Load the uploaded CSV file
file_path = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/loess_values.csv'
df = pd.read_csv(file_path)

# Define the groups
bins = [0, 0.33, 0.66, 1]
group_labels = ['low', 'medium', 'high']
df['x_group'] = pd.cut(df['x'], bins=bins, labels=group_labels, include_lowest=True)

# Initialize a list to store the results
results = []

# Iterate over each region and biomarker
for (region, biomarker), group_data in df.groupby(['Region', 'Biomarker']):
    for group in group_labels:
        group_df = group_data[group_data['x_group'] == group]

        # Separate data into AD and Healthy
        ad_data = group_df[group_df['Status'] == 'AD']['y']
        healthy_data = group_df[group_df['Status'] == 'Healthy']['y']

        # Perform t-test if there are enough data points in both groups
        if len(ad_data) > 1 and len(healthy_data) > 1:
            t_stat, p_val = ttest_ind(ad_data, healthy_data)
            results.append({
                'Region': region,
                'Biomarker': biomarker,
                'Group': group,
                't_stat': t_stat,
                'p_value': p_val
            })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Apply Benjamini-Hochberg correction
results_df['p_adj'] = multitest.multipletests(results_df['p_value'], method='fdr_bh')[1]

# Calculate transformed p-values
results_df['transf_p_adj'] = -1 * np.log10(results_df['p_adj'])

# Rename the columns
results_df = results_df.rename(columns={
    'Group': 'pseudotime_bins_eml',
    'Biomarker': 'Feature'
})

# Save the updated DataFrame to a new CSV file
output_csv_path = '/Users/bryjc/Library/CloudStorage/GoogleDrive-bryjc@stanford.edu/My Drive/Bangelo_Lab/AD_Resilience_Project/inital_DM_analysis/data/renamed_t_test_results.csv'
results_df.to_csv(output_csv_path, index=False)

